//! Non-conforming P1 red refinement seam types (RFC-0001).
use crate::base::Point2;
use crate::constraints::AffineConstraints;
use crate::grid::{BoundaryEdgeGeometry, BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::integrators::{diffusion_reaction_p1, DiffusionTensor};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::{error, Error, Result};
use std::collections::{BTreeMap, BTreeSet};

fn edge_key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Walk sub-edge keys `(a, m1)`, `(m1, b)` of an interface edge for level
/// differences greater than one (RFC-0001 §3): emit one `HangingNode` per
/// split sub-edge key found in `edge_splits`. Iterative with a visited set so
/// deep chains cannot exhaust the stack and malformed split maps cannot loop.
fn collect_subedges(
    edge_splits: &BTreeMap<(VertexId, VertexId), VertexId>,
    hanging: &mut BTreeMap<VertexId, HangingNode>,
    a: VertexId,
    b: VertexId,
) -> Result<()> {
    let mut pending = vec![(a, b)];
    let mut seen = BTreeSet::new();
    while let Some((x, y)) = pending.pop() {
        let (lo, hi) = edge_key(x, y);
        if !seen.insert((lo, hi)) {
            continue;
        }
        if let Some(&mid) = edge_splits.get(&(lo, hi)) {
            if mid == lo || mid == hi {
                continue;
            }
            let node = HangingNode {
                vertex: mid,
                masters: [lo, hi],
            };
            if let Some(existing) = hanging.get(&mid) {
                if *existing != node {
                    return Err(error("inconsistent hanging node masters"));
                }
            } else {
                hanging.insert(mid, node);
            }
            pending.push((lo, mid));
            pending.push((mid, hi));
        }
    }
    Ok(())
}

/// Discover the hanging nodes of an active mesh from `edge_splits` alone
/// (RFC-0001 §3). An edge is an interface root iff its canonical key is both
/// incident to an active cell and present in `edge_splits`; no cell level or
/// `parent_cells` chain is consulted. Shared by `refine` and `coarsen` so a
/// pre-existing interface reappears when both neighbours become coarse again.
fn discover_hanging(
    edge_splits: &BTreeMap<(VertexId, VertexId), VertexId>,
    incidence: &BTreeMap<(VertexId, VertexId), Vec<usize>>,
) -> Result<Vec<HangingNode>> {
    let mut hanging = BTreeMap::new();
    let roots: BTreeSet<(VertexId, VertexId)> = incidence
        .keys()
        .copied()
        .filter(|edge| edge_splits.contains_key(edge))
        .collect();
    for (a, b) in roots {
        collect_subedges(edge_splits, &mut hanging, a, b)?;
    }
    Ok(hanging.into_values().collect())
}

/// A midpoint DoF on an unsplit parent edge, constrained to the arithmetic
/// mean of its endpoints for continuous P1 fields. All IDs belong to the
/// refined object, even though original vertex indices are retained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HangingNode {
    pub vertex: VertexId,
    pub masters: [VertexId; 2],
}

/// The active non-conforming mesh seam type (RFC-0001 §1): the flattened
/// active (leaf) cell list plus the cross-level provenance needed for
/// multi-level refinement. Typed IDs in `cells` and `parent_cells` and the
/// `edge_splits` keys and the `midpoint_parents` keys and values belong to
/// this mesh instance.
///
/// # Invariants
///
/// [`validate`](Self::validate) checks the structural invariants: original
/// vertices precede midpoints, `parent_cells` is parallel to `cells`, every
/// typed vertex id is in range, vertex coordinates are finite, edge endpoints
/// are distinct, and each `midpoint_parents` entry agrees with `edge_splits`.
/// Geometric conformity, nonoverlap and absence of hanging interfaces remain
/// caller obligations (as for [`Mesh2`]); use
/// [`try_new`](crate::grid::NonconformingRefinementHistory::try_new) for
/// untrusted input.
#[derive(Debug, Clone)]
pub struct NonconformingMesh {
    /// All vertices: original plus midpoints (compacted).
    pub vertices: Vec<Point2>,
    /// Active (leaf) cells only, canonically oriented and possibly
    /// nonconforming.
    pub cells: Vec<[VertexId; 3]>,
    /// Per active cell, the index of the source cell in the mesh that produced
    /// this level (opaque, level-scoped provenance for the frozen
    /// `NonconformingRefinement::parent_cell` accessor). These ids are not
    /// self-parent pointers and not a cell level; they must not be used for
    /// level comparisons or to index another level's cell list (RFC-0001 §1).
    pub parent_cells: Vec<CellId>,
    /// Canonical (sorted) split-edge key to its midpoint vertex: the single
    /// source of truth for "is this edge split, by which vertex" (RFC-0001
    /// §2).
    pub edge_splits: BTreeMap<(VertexId, VertexId), VertexId>,
    /// One entry per midpoint (slave) vertex, with its two masters.
    pub hanging: Vec<HangingNode>,
    /// Exterior edge vertex pairs only. Interface edges are never boundary
    /// edges.
    pub boundary_edges: Vec<[VertexId; 2]>,
    /// Unique exterior vertices; hanging nodes on internal edges are
    /// excluded.
    pub boundary_vertices: Vec<VertexId>,
    /// Count of level-0 vertices.
    pub original_vertices: usize,
    /// Per midpoint vertex, its two master vertices; inverse of
    /// `edge_splits` by master pair.
    pub midpoint_parents: BTreeMap<VertexId, (VertexId, VertexId)>,
}

impl NonconformingMesh {
    /// Check the structural invariants documented on [`NonconformingMesh`].
    ///
    /// This is a boundary check for meshes supplied by callers: it rejects out
    /// of range typed ids, nonfinite coordinates, degenerate edge keys and
    /// `midpoint_parents` entries that disagree with `edge_splits`. It does
    /// **not** verify geometric conformity or nonoverlap; those remain caller
    /// obligations, as for [`Mesh2`].
    pub fn validate(&self) -> Result<()> {
        let n = self.vertices.len();
        if self.original_vertices > n {
            return Err(error("original vertex count exceeds vertex list"));
        }
        if self.parent_cells.len() != self.cells.len() {
            return Err(error("parent_cells is not parallel to cells"));
        }
        if self
            .vertices
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite())
        {
            return Err(error("nonfinite vertex coordinate"));
        }
        let in_range = |v: VertexId| v.0 < n;
        for cell in &self.cells {
            if cell.iter().any(|&v| !in_range(v)) {
                return Err(error("cell vertex id out of range"));
            }
        }
        for (&(a, b), &mid) in &self.edge_splits {
            if a == b || !in_range(a) || !in_range(b) || !in_range(mid) {
                return Err(error("invalid edge split endpoints or midpoint"));
            }
        }
        for (&mid, &(a, b)) in &self.midpoint_parents {
            if a == b || !in_range(mid) || !in_range(a) || !in_range(b) {
                return Err(error("invalid midpoint parent ids"));
            }
            if self.edge_splits.get(&edge_key(a, b)) != Some(&mid) {
                return Err(error("midpoint parents disagree with edge splits"));
            }
        }
        for h in &self.hanging {
            if h.masters[0] == h.masters[1]
                || !in_range(h.vertex)
                || !in_range(h.masters[0])
                || !in_range(h.masters[1])
            {
                return Err(error("invalid hanging node ids"));
            }
        }
        for edge in &self.boundary_edges {
            if edge[0] == edge[1] || !in_range(edge[0]) || !in_range(edge[1]) {
                return Err(error("invalid boundary edge id"));
            }
        }
        if self.boundary_vertices.iter().any(|&v| !in_range(v)) {
            return Err(error("boundary vertex id out of range"));
        }
        Ok(())
    }
}

/// A single nonconforming refinement of a conforming triangular mesh: the
/// [`NonconformingMesh`] of the refined level plus transition provenance
/// (RFC-0001 §9).
///
/// Marked triangles are split into four; unmarked triangles are not split.
/// Midpoints are shared among marked neighbors. Hanging nodes are identified
/// *from the parent edge connectivity and the marking*, not by geometric
/// searches; boundary edges are derived from the parent boundary, not from
/// incidence in the nonconforming cells. Thus false internal boundary edges
/// are never exposed through this type. Parent input must be geometrically
/// conforming, nonoverlapping and without hanging nodes (not fully checked by
/// `Mesh2`). No conversion to `Mesh2` is provided: its boundary detection and
/// conforming-mesh algorithms would otherwise misinterpret the interface.
///
/// The refined-level mesh is exposed by [`mesh`](Self::mesh). Further levels
/// are produced by [`refine`](Self::refine); one level can be undone for a
/// selection of sibling groups by [`coarsen`](Self::coarsen). Multi-level
/// constraints compose through the iterative `AffineConstraints::close`.
/// Markers and solution vectors remain caller-managed.
#[derive(Debug, Clone)]
pub struct NonconformingRefinement {
    mesh: NonconformingMesh,
    children: Vec<Vec<CellId>>,
    /// Level-scoped provenance for coarsening: the midpoints newly created
    /// while refining each sibling group. Keys index this transition's
    /// `children`; a midpoint reused via `edge_splits` is not recorded here.
    introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>>,
    oriented_boundary_edges: Vec<[VertexId; 2]>,
    parent_boundary_edges: Vec<BoundaryEdgeId>,
    boundary_vertices: Vec<VertexId>,
    original_vertices: usize,
}
impl NonconformingRefinement {
    /// Refine exactly the marked parent cells. Duplicate and invalid IDs, and
    /// numerically invalid children, are rejected. An empty selection is valid.
    pub fn new(parent: &Mesh2, marked: &[CellId]) -> Result<Self> {
        let mut marks = BTreeSet::new();
        let mut split = BTreeSet::new();
        for &id in marked {
            if !marks.insert(id) {
                return Err(error("duplicate marked cell"));
            }
            let cell = parent
                .cells()
                .get(id.0)
                .ok_or_else(|| error("marked cell index out of bounds"))?;
            for i in 0..3 {
                split.insert(edge_key(cell[i], cell[(i + 1) % 3]));
            }
        }
        parent
            .cells()
            .len()
            .checked_mul(4)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = parent.vertices().to_vec();
        let mut edge_splits = BTreeMap::new();
        let prior = edge_splits.clone();
        let mut midpoint_parents = BTreeMap::new();
        for &(a, b) in &split {
            let p = vertices[a.0];
            let q = vertices[b.0];
            let midpoint = Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y)?;
            let id = VertexId(vertices.len());
            vertices.push(midpoint);
            midpoint_parents.insert(id, (a, b));
            edge_splits.insert((a, b), id);
        }
        // Record provenance after creation: every marked group owns the
        // midpoints of its three edges, including edges shared with a group.
        let mut introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>> = BTreeMap::new();
        for (i, &[a, b, c]) in parent.cells().iter().enumerate() {
            if !marks.contains(&CellId(i)) {
                continue;
            }
            let entry = introduced_midpoints.entry(CellId(i)).or_default();
            for (x, y) in [(a, b), (b, c), (c, a)] {
                let key = edge_key(x, y);
                if !prior.contains_key(&key) {
                    entry.insert(edge_splits[&key]);
                }
            }
        }
        let mut cells = Vec::new();
        let mut parent_cells = Vec::new();
        let mut children = Vec::with_capacity(parent.cells().len());
        let mut hanging = BTreeMap::new();
        for (i, &[a, b, c]) in parent.cells().iter().enumerate() {
            let mut ids = Vec::new();
            if marks.contains(&CellId(i)) {
                let ab = edge_splits[&edge_key(a, b)];
                let bc = edge_splits[&edge_key(b, c)];
                let ca = edge_splits[&edge_key(c, a)];
                for cell in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
                    ids.push(CellId(cells.len()));
                    cells.push(cell);
                    parent_cells.push(CellId(i));
                }
            } else {
                ids.push(CellId(cells.len()));
                cells.push([a, b, c]);
                parent_cells.push(CellId(i));
                for (x, y) in [(a, b), (b, c), (c, a)] {
                    let (v, w) = edge_key(x, y);
                    if let Some(&mid) = edge_splits.get(&(v, w)) {
                        hanging.insert(
                            (v, w),
                            HangingNode {
                                vertex: mid,
                                masters: [v, w],
                            },
                        );
                    }
                }
            }
            children.push(ids);
        }
        // Validate each triangle and conforming edge incidence. Mesh2 cannot
        // detect hanging interfaces; its inferred boundary is NOT used here.
        let checked = Mesh2::new(vertices, cells)?;
        let vertices = checked.vertices().to_vec();
        let cells = checked.cells().to_vec();
        let mut boundary = BTreeMap::new();
        for i in 0..parent.boundary_edges().len() {
            let id = BoundaryEdgeId(i);
            let edge = parent.boundary_edge_geometry(id)?.vertices;
            let [a, b] = edge;
            if let Some(&mid) = edge_splits.get(&edge_key(a, b)) {
                boundary.insert(edge_key(a, mid), ([a, mid], id));
                boundary.insert(edge_key(mid, b), ([mid, b], id));
            } else {
                boundary.insert(edge_key(a, b), ([a, b], id));
            }
        }
        let mut boundary_edges = Vec::new();
        let mut oriented_boundary_edges = Vec::new();
        let mut parent_boundary_edges = Vec::new();
        let mut boundary_vertices = BTreeSet::new();
        for ((a, b), (oriented, parent_id)) in boundary {
            boundary_edges.push([a, b]);
            oriented_boundary_edges.push(oriented);
            parent_boundary_edges.push(parent_id);
            boundary_vertices.insert(a);
            boundary_vertices.insert(b);
        }
        let boundary_vertices: Vec<VertexId> = boundary_vertices.into_iter().collect();
        let original_vertices = parent.vertices().len();
        let mesh = NonconformingMesh {
            vertices,
            cells,
            parent_cells,
            edge_splits,
            hanging: hanging.into_values().collect(),
            boundary_edges,
            boundary_vertices: boundary_vertices.clone(),
            original_vertices,
            midpoint_parents,
        };
        Ok(Self {
            mesh,
            children,
            introduced_midpoints,
            oriented_boundary_edges,
            parent_boundary_edges,
            boundary_vertices,
            original_vertices,
        })
    }

    /// The refined-level mesh: active cells, provenance, and the cross-level
    /// edge-split map (RFC-0001 §1).
    pub fn mesh(&self) -> &NonconformingMesh {
        &self.mesh
    }

    /// Refine further from this level, returning the next transition
    /// (RFC-0001 §1, §2, §3, §11). Empty markings are rejected; duplicate and
    /// out-of-bounds cell ids are rejected as in [`new`](Self::new).
    ///
    /// Splitting an edge always consults `edge_splits` first, so a midpoint
    /// created at an earlier level keeps its index. An active edge is a
    /// coarse/fine interface iff its canonical key is in `edge_splits` and is
    /// still incident to an active cell; no cell level is derived from
    /// `parent_cells`. Its midpoint (and every split sub-edge midpoint, for
    /// level differences greater than one) becomes a hanging node. Existing
    /// interfaces are re-tested and dropped once both sides have been refined.
    /// Unbalanced level differences are allowed; no neighbours are expanded,
    /// and continuity is enforced by constraints.
    pub fn refine(&self, marked: &[CellId]) -> Result<Self> {
        if marked.is_empty() {
            return Err(error("empty refinement marking"));
        }
        let mut marks = BTreeSet::new();
        let mut split = BTreeSet::new();
        for &id in marked {
            if !marks.insert(id) {
                return Err(error("duplicate marked cell"));
            }
            let cell = self
                .mesh
                .cells
                .get(id.0)
                .ok_or_else(|| error("marked cell index out of bounds"))?;
            for i in 0..3 {
                split.insert(edge_key(cell[i], cell[(i + 1) % 3]));
            }
        }

        let cell_capacity = self
            .mesh
            .cells
            .len()
            .checked_mul(4)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = self.mesh.vertices.clone();
        let mut edge_splits = self.mesh.edge_splits.clone();
        let prior = edge_splits.clone();
        let mut midpoint_parents = self.mesh.midpoint_parents.clone();
        for (a, b) in split {
            if edge_splits.contains_key(&(a, b)) {
                continue;
            }
            let p = vertices[a.0];
            let q = vertices[b.0];
            let point = Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y)?;
            let id = VertexId(vertices.len());
            vertices.push(point);
            midpoint_parents.insert(id, (a, b));
            edge_splits.insert((a, b), id);
        }
        // Only keys absent from the prior split map are newly introduced by
        // this group. A reused interface midpoint belongs to no group here.
        let mut introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>> = BTreeMap::new();
        for i in marks.iter().map(|id| id.0) {
            let [a, b, c] = self.mesh.cells[i];
            let entry = introduced_midpoints.entry(CellId(i)).or_default();
            for (x, y) in [(a, b), (b, c), (c, a)] {
                let key = edge_key(x, y);
                if !prior.contains_key(&key) {
                    entry.insert(edge_splits[&key]);
                }
            }
        }

        let mut cells = Vec::with_capacity(cell_capacity);
        let mut parent_cells = Vec::with_capacity(cell_capacity);
        let mut children = Vec::with_capacity(self.mesh.cells.len());
        for (i, &[a, b, c]) in self.mesh.cells.iter().enumerate() {
            let mut ids = Vec::new();
            if marks.contains(&CellId(i)) {
                let ab = edge_splits[&edge_key(a, b)];
                let bc = edge_splits[&edge_key(b, c)];
                let ca = edge_splits[&edge_key(c, a)];
                for cell in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
                    ids.push(CellId(cells.len()));
                    cells.push(cell);
                    parent_cells.push(CellId(i));
                }
            } else {
                ids.push(CellId(cells.len()));
                cells.push([a, b, c]);
                parent_cells.push(CellId(i));
            }
            children.push(ids);
        }

        // Edge -> active cell incidence on the new level. An edge is an
        // interface root iff its key is split and still incident to an active
        // cell; no level is derived from `parent_cells` (RFC-0001 §3).
        let mut incidence = BTreeMap::<(VertexId, VertexId), Vec<usize>>::new();
        for (i, cell) in cells.iter().enumerate() {
            for j in 0..3 {
                incidence
                    .entry(edge_key(cell[j], cell[(j + 1) % 3]))
                    .or_default()
                    .push(i);
            }
        }
        let hanging = discover_hanging(&edge_splits, &incidence)?;

        // Boundary edges inherit their parent edge id and orientation; a
        // split exterior edge contributes its two half-edges.
        let mut boundary_edges = Vec::new();
        let mut oriented_boundary_edges = Vec::new();
        let mut parent_boundary_edges = Vec::new();
        let mut boundary_vertices = BTreeSet::new();
        for (i, &[a, b]) in self.oriented_boundary_edges.iter().enumerate() {
            let parent_id = BoundaryEdgeId(i);
            let halves: Vec<[VertexId; 2]> = match edge_splits.get(&edge_key(a, b)) {
                Some(&mid) => vec![[a, mid], [mid, b]],
                None => vec![[a, b]],
            };
            for oriented in halves {
                let (lo, hi) = edge_key(oriented[0], oriented[1]);
                boundary_edges.push([lo, hi]);
                oriented_boundary_edges.push(oriented);
                parent_boundary_edges.push(parent_id);
                boundary_vertices.insert(oriented[0]);
                boundary_vertices.insert(oriented[1]);
            }
        }
        let boundary_vertices: Vec<VertexId> = boundary_vertices.into_iter().collect();

        let checked = Mesh2::new(vertices, cells)?;
        let mesh = NonconformingMesh {
            vertices: checked.vertices().to_vec(),
            cells: checked.cells().to_vec(),
            parent_cells,
            edge_splits,
            hanging,
            boundary_edges,
            boundary_vertices: boundary_vertices.clone(),
            original_vertices: self.mesh.original_vertices,
            midpoint_parents,
        };
        Ok(Self {
            mesh,
            children,
            introduced_midpoints,
            oriented_boundary_edges,
            parent_boundary_edges,
            boundary_vertices,
            original_vertices: self.mesh.original_vertices,
        })
    }

    /// Restore the parents of fully marked sibling groups, one level per call
    /// (RFC-0001 §8). The input is unchanged; a new transition wrapping a new
    /// [`NonconformingMesh`] is returned.
    ///
    /// A parent is restored only when all four of its children are marked,
    /// none of those children has descendants, no surviving active cell uses
    /// a midpoint newly introduced by that group (reused interface midpoints
    /// are never removed), and no input hanging node references any removed
    /// midpoint as a slave or a master. The result's hanging nodes are
    /// recomputed from the surviving split edges. Otherwise, an error is
    /// returned. Empty markings are rejected (§11); duplicate and
    /// out-of-bounds cell IDs are rejected as in [`refine`](Self::refine).
    ///
    /// # ID invalidation
    ///
    /// Survivors keep their original vertex order, removed midpoints are
    /// excised and indices are compacted. **All** `VertexId`, `CellId`, `DofId`
    /// and `BoundaryEdgeId` values from the input are invalidated by this
    /// structural change (ADR-0003 §4); no automatic remapping is provided.
    /// Use only the IDs exposed by the returned transition. Surviving
    /// `edge_splits` entries are renumbered with the same vertex map; entries
    /// whose midpoint or either endpoint was removed are dropped, cascading
    /// to sub-edges that referenced a removed vertex.
    pub fn coarsen(&self, marked: &[CellId]) -> Result<Self> {
        if marked.is_empty() {
            return Err(error("empty coarsening marking"));
        }
        let mut marks = BTreeSet::new();
        for &id in marked {
            if !marks.insert(id) {
                return Err(error("duplicate marked cell"));
            }
            if id.0 >= self.mesh.cells.len() {
                return Err(error("marked cell index out of bounds"));
            }
        }

        // Each parent group owns a disjoint set of active children. Groups of
        // one child are carried cells and cannot be coarsened.
        let mut group_of = vec![None; self.mesh.cells.len()];
        let mut covered: usize = 0;
        for (p, kids) in self.children.iter().enumerate() {
            covered = covered
                .checked_add(kids.len())
                .ok_or_else(|| error("sibling coverage overflow"))?;
            if kids.len() == 1 {
                continue;
            }
            if kids.len() != 4 {
                return Err(error("unexpected sibling group size"));
            }
            for &kid in kids {
                if kid.0 >= self.mesh.cells.len() {
                    return Err(error("stale child cell"));
                }
                if group_of[kid.0].is_some() {
                    return Err(error("cell in multiple sibling groups"));
                }
                group_of[kid.0] = Some(p);
            }
        }
        if covered != self.mesh.cells.len() {
            return Err(error("incomplete sibling coverage"));
        }
        let mut restore = BTreeSet::new();
        for &id in &marks {
            match group_of[id.0] {
                Some(p) => {
                    restore.insert(p);
                }
                None => return Err(error("marked cell has no coarsenable sibling group")),
            }
        }

        // Condition 1 (all children marked). The candidate removal set comes
        // from the transition's provenance, never from master membership: a
        // midpoint reused from an earlier level is not introduced by this
        // group and must survive (RFC-0001 §8).
        let mut parent_cell = BTreeMap::new();
        let mut removed = BTreeSet::new();
        for &p in &restore {
            let kids = &self.children[p];
            let mut verts = BTreeSet::new();
            for &kid in kids {
                if !marks.contains(&kid) {
                    return Err(error("partial sibling group"));
                }
                // Condition 2 is vacuous in the active-only model: a cell with
                // descendants is not a member of `mesh.cells`. It is verified
                // here by requiring every child to be an active cell owned by
                // this group (`group_of` above), rather than a stale or
                // level-scoped id borrowed from another transition.
                let cell = self
                    .mesh
                    .cells
                    .get(kid.0)
                    .ok_or_else(|| error("marked cell has descendants"))?;
                for &v in cell {
                    verts.insert(v);
                }
            }
            // Red refinement emits children in the order
            // [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]; the first
            // three corners recover the parent's original vertex order.
            let a = self.mesh.cells[kids[0].0][0];
            let b = self.mesh.cells[kids[1].0][1];
            let c = self.mesh.cells[kids[2].0][2];
            if !verts.contains(&a) || !verts.contains(&b) || !verts.contains(&c) {
                return Err(error("inconsistent sibling vertices"));
            }
            removed.extend(
                self.introduced_midpoints
                    .get(&CellId(p))
                    .cloned()
                    .unwrap_or_default(),
            );
            parent_cell.insert(p, [a, b, c]);
        }

        // Cells of the result, still in the old vertex index space: restored
        // parents replace their four children at the group's first position.
        let mut old_cells = Vec::with_capacity(self.mesh.cells.len());
        for (p, kids) in self.children.iter().enumerate() {
            if restore.contains(&p) {
                old_cells.push(parent_cell[&p]);
            } else {
                for &kid in kids {
                    old_cells.push(
                        *self
                            .mesh
                            .cells
                            .get(kid.0)
                            .ok_or_else(|| error("stale child cell"))?,
                    );
                }
            }
        }
        // `Mesh2` forbids unused vertices: a midpoint can become orphaned when
        // every cell that used it is restored to its parent (cascade).
        let mut used = BTreeSet::new();
        for cell in &old_cells {
            for &v in cell {
                used.insert(v);
            }
        }
        for i in 0..self.mesh.vertices.len() {
            let v = VertexId(i);
            if used.contains(&v) {
                continue;
            }
            if i < self.mesh.original_vertices || !self.mesh.midpoint_parents.contains_key(&v) {
                return Err(error("unused original vertex after coarsening"));
            }
            removed.insert(v);
        }

        // Conditions 3 and 4, against the final removal set.
        for (i, cell) in self.mesh.cells.iter().enumerate() {
            if marks.contains(&CellId(i)) {
                continue;
            }
            for &v in cell {
                if removed.contains(&v) {
                    return Err(error("surviving cell references a removed midpoint"));
                }
            }
        }
        for h in &self.mesh.hanging {
            if removed.contains(&h.vertex)
                || removed.contains(&h.masters[0])
                || removed.contains(&h.masters[1])
            {
                return Err(error("hanging node references a removed midpoint"));
            }
        }

        // Compact the survivors, preserving their original order.
        let n_vertices = self
            .mesh
            .vertices
            .len()
            .checked_sub(removed.len())
            .ok_or_else(|| error("coarsening size underflow"))?;
        let mut remap = vec![None; self.mesh.vertices.len()];
        let mut vertices = Vec::with_capacity(n_vertices);
        for (i, &point) in self.mesh.vertices.iter().enumerate() {
            if !removed.contains(&VertexId(i)) {
                remap[i] = Some(VertexId(vertices.len()));
                vertices.push(point);
            }
        }
        for slot in remap.iter().take(self.mesh.original_vertices) {
            if slot.is_none() {
                return Err(error("original vertex removed by coarsening"));
            }
        }
        let renumber = |v: VertexId| -> Result<VertexId> {
            remap
                .get(v.0)
                .and_then(|slot| *slot)
                .ok_or_else(|| error("removed vertex in coarsened cell"))
        };

        // Rebuild cells and the parent/child index space. `parent_cells` and
        // `children` share the parent-level group index `p`: every cell
        // produced from group `p` records `CellId(p)`, and `children[p]` holds
        // all of those cells in the new numbering (RFC-0001 §1, §9).
        let mut cells = Vec::with_capacity(old_cells.len());
        let mut parent_cells = Vec::with_capacity(old_cells.len());
        let mut children = vec![Vec::new(); self.children.len()];
        for (p, kids) in self.children.iter().enumerate() {
            let start = cells.len();
            if restore.contains(&p) {
                let [a, b, c] = parent_cell[&p];
                cells.push([renumber(a)?, renumber(b)?, renumber(c)?]);
            } else {
                for &kid in kids {
                    let cell = self
                        .mesh
                        .cells
                        .get(kid.0)
                        .ok_or_else(|| error("stale child cell"))?;
                    cells.push([renumber(cell[0])?, renumber(cell[1])?, renumber(cell[2])?]);
                }
            }
            let new_ids: Vec<CellId> = (start..cells.len()).map(CellId).collect();
            parent_cells.resize(cells.len(), CellId(p));
            children[p] = new_ids;
        }

        // Drop removed midpoints and cascade to any split key that contains a
        // removed endpoint; renumber the survivors with the same vertex map.
        let mut edge_splits = BTreeMap::new();
        for (&(a, b), &mid) in &self.mesh.edge_splits {
            if removed.contains(&a) || removed.contains(&b) || removed.contains(&mid) {
                continue;
            }
            edge_splits.insert(edge_key(renumber(a)?, renumber(b)?), renumber(mid)?);
        }
        let mut midpoint_parents = BTreeMap::new();
        for (&mid, &(a, b)) in &self.mesh.midpoint_parents {
            if removed.contains(&mid) || removed.contains(&a) || removed.contains(&b) {
                continue;
            }
            midpoint_parents.insert(renumber(mid)?, (renumber(a)?, renumber(b)?));
        }
        // Carry midpoint provenance through the vertex compaction: restored
        // groups lose all of their entries, surviving groups keep theirs.
        let mut introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>> = BTreeMap::new();
        for (&p, mids) in &self.introduced_midpoints {
            let kept = mids
                .iter()
                .filter(|m| !removed.contains(m))
                .map(|&m| renumber(m))
                .collect::<Result<BTreeSet<_>>>()?;
            if !kept.is_empty() {
                introduced_midpoints.insert(p, kept);
            }
        }

        // Recompute interfaces on the coarsened cells rather than carrying the
        // input's hanging list: an interface inactive while both neighbours
        // were fine reappears when one is restored (RFC-0001 §3, §8).
        let mut incidence = BTreeMap::<(VertexId, VertexId), Vec<usize>>::new();
        for (i, cell) in cells.iter().enumerate() {
            for j in 0..3 {
                incidence
                    .entry(edge_key(cell[j], cell[(j + 1) % 3]))
                    .or_default()
                    .push(i);
            }
        }
        let hanging = discover_hanging(&edge_splits, &incidence)?;

        // Collapse pairs of boundary half-edges back into their parent edge
        // when the shared midpoint was excised; keep orientation and the
        // previous-level boundary edge id.
        if self.oriented_boundary_edges.len() != self.parent_boundary_edges.len() {
            return Err(error("inconsistent boundary provenance"));
        }
        let mut boundary = BTreeMap::<(VertexId, VertexId), ([VertexId; 2], BoundaryEdgeId)>::new();
        for (i, &[x, y]) in self.oriented_boundary_edges.iter().enumerate() {
            let parent_id = self.parent_boundary_edges[i];
            let (ox, oy) = if removed.contains(&x) && removed.contains(&y) {
                continue;
            } else if removed.contains(&x) {
                let (a, b) = *self
                    .mesh
                    .midpoint_parents
                    .get(&x)
                    .ok_or_else(|| error("removed boundary vertex is not a midpoint"))?;
                let other = if a == y {
                    b
                } else if b == y {
                    a
                } else {
                    return Err(error("inconsistent boundary midpoint"));
                };
                (other, y)
            } else if removed.contains(&y) {
                let (a, b) = *self
                    .mesh
                    .midpoint_parents
                    .get(&y)
                    .ok_or_else(|| error("removed boundary vertex is not a midpoint"))?;
                let other = if a == x {
                    b
                } else if b == x {
                    a
                } else {
                    return Err(error("inconsistent boundary midpoint"));
                };
                (x, other)
            } else {
                (x, y)
            };
            if removed.contains(&ox) || removed.contains(&oy) {
                return Err(error("nested removed boundary midpoint"));
            }
            let nx = renumber(ox)?;
            let ny = renumber(oy)?;
            boundary
                .entry(edge_key(nx, ny))
                .or_insert(([nx, ny], parent_id));
        }
        let mut boundary_edges = Vec::new();
        let mut oriented_boundary_edges = Vec::new();
        let mut parent_boundary_edges = Vec::new();
        let mut boundary_vertices = BTreeSet::new();
        for ((a, b), (oriented, parent_id)) in boundary {
            boundary_edges.push([a, b]);
            oriented_boundary_edges.push(oriented);
            parent_boundary_edges.push(parent_id);
            boundary_vertices.insert(oriented[0]);
            boundary_vertices.insert(oriented[1]);
        }
        let boundary_vertices: Vec<VertexId> = boundary_vertices.into_iter().collect();

        let checked = Mesh2::new(vertices, cells)?;
        let mesh = NonconformingMesh {
            vertices: checked.vertices().to_vec(),
            cells: checked.cells().to_vec(),
            parent_cells,
            edge_splits,
            hanging,
            boundary_edges: boundary_edges.clone(),
            boundary_vertices: boundary_vertices.clone(),
            original_vertices: self.mesh.original_vertices,
            midpoint_parents,
        };
        Ok(Self {
            mesh,
            children,
            introduced_midpoints,
            oriented_boundary_edges,
            parent_boundary_edges,
            boundary_vertices,
            original_vertices: self.mesh.original_vertices,
        })
    }

    /// Vertex coordinates in this refined object's DoF order.
    pub fn vertices(&self) -> &[Point2] {
        &self.mesh.vertices
    }
    /// Positively oriented triangle connectivity, potentially nonconforming.
    pub fn cells(&self) -> &[[VertexId; 3]] {
        &self.mesh.cells
    }
    /// Exterior edges only, with canonically ordered endpoints. Edge IDs are
    /// indices into this list and belong to this object.
    pub fn boundary_edges(&self) -> &[[VertexId; 2]] {
        &self.mesh.boundary_edges
    }
    /// Unique exterior vertices; hanging nodes on internal edges are excluded.
    pub fn boundary_vertices(&self) -> &[VertexId] {
        &self.boundary_vertices
    }
    /// Midpoints constrained across a coarse/fine interface.
    pub fn hanging_nodes(&self) -> &[HangingNode] {
        &self.mesh.hanging
    }
    /// The children of a parent cell; four if marked, one otherwise.
    pub fn children(&self, parent: CellId) -> Result<&[CellId]> {
        self.children
            .get(parent.0)
            .map(Vec::as_slice)
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Parent ID of a refined cell.
    pub fn parent_cell(&self, child: CellId) -> Result<CellId> {
        self.mesh
            .parent_cells
            .get(child.0)
            .copied()
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Parent boundary edge for a refined exterior edge.
    pub fn parent_boundary_edge(&self, child: BoundaryEdgeId) -> Result<BoundaryEdgeId> {
        self.parent_boundary_edges
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Cell coordinates, in local vertex order.
    pub fn cell_points(&self, id: CellId) -> Result<[Point2; 3]> {
        let cell = self
            .mesh
            .cells
            .get(id.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok(cell.map(|v| self.mesh.vertices[v.0]))
    }
    /// Oriented exterior edge, length and outward normal inherited from the
    /// parent boundary orientation. IDs and coordinates are those of this mesh.
    pub fn boundary_edge_geometry(&self, id: BoundaryEdgeId) -> Result<BoundaryEdgeGeometry> {
        let vertices = *self
            .oriented_boundary_edges
            .get(id.0)
            .ok_or_else(|| error("boundary edge index out of bounds"))?;
        let points = vertices.map(|v| self.mesh.vertices[v.0]);
        let dx = points[1].x - points[0].x;
        let dy = points[1].y - points[0].y;
        let length = dx.hypot(dy);
        if !length.is_finite() || length <= 0.0 {
            return Err(error("invalid boundary edge length"));
        }
        Ok(BoundaryEdgeGeometry {
            vertices,
            points,
            outward_normal: [dy / length, -dx / length],
            length,
        })
    }
    /// Build midpoint relations for all hanging nodes. Add Dirichlet or other
    /// compatible constraints before closing the returned builder. Conflicting
    /// rules for a slave are rejected by `AffineConstraints::add`.
    pub fn hanging_constraints(&self) -> Result<AffineConstraints> {
        let mut builder = AffineConstraints::new(self.mesh.vertices.len());
        for h in &self.mesh.hanging {
            builder.add(
                h.vertex.0,
                0.0,
                &[(h.masters[0].0, 0.5), (h.masters[1].0, 0.5)],
            )?;
        }
        Ok(builder)
    }
    /// Interpolate a finite scalar function on the actual exterior vertices,
    /// yielding `(vertex index, value)` suitable for explicit Dirichlet rules.
    pub fn interpolate_boundary<F: Fn(Point2) -> f64>(&self, f: F) -> Result<Vec<(usize, f64)>> {
        self.boundary_vertices
            .iter()
            .map(|v| {
                let value = f(self.mesh.vertices[v.0]);
                if !value.is_finite() {
                    Err(error("boundary function returned a nonfinite value"))
                } else {
                    Ok((v.0, value))
                }
            })
            .collect()
    }
    /// Prolong a scalar P1 parent field by midpoint interpolation. The result
    /// satisfies all hanging-node constraints in exact arithmetic. Masters may
    /// themselves be midpoints from an earlier level: midpoint ids are
    /// assigned after their masters, so accumulating into `result` resolves
    /// the interpolation chain in one pass (RFC-0001 §7).
    pub fn prolongate_p1(&self, values: &[f64]) -> Result<Vec<f64>> {
        if values.len() != self.original_vertices || values.iter().any(|v| !v.is_finite()) {
            return Err(error("invalid parent P1 field"));
        }
        let mut result = values.to_vec();
        for &(a, b) in self.mesh.midpoint_parents.values() {
            let left = *result
                .get(a.0)
                .ok_or_else(|| error("invalid midpoint parent"))?;
            let right = *result
                .get(b.0)
                .ok_or_else(|| error("invalid midpoint parent"))?;
            let v = 0.5 * left + 0.5 * right;
            if !v.is_finite() {
                return Err(error("nonfinite interpolated P1 value"));
            }
            result.push(v);
        }
        Ok(result)
    }
    /// Assemble the unconstrained P1 volume system for `-Δu=f` by summing
    /// local cell contributions. Apply [`hanging_constraints`](Self::hanging_constraints)
    /// plus boundary rules and condense before solving. The mesh must be
    /// nonoverlapping and derived from a conforming parent. This assembles no
    /// Neumann terms and does not verify coercivity or SPD. Source is sampled
    /// by the degree-five triangle quadrature used by `diffusion_reaction_p1`.
    pub fn assemble_poisson<F: Fn(Point2) -> f64>(
        &self,
        source: F,
    ) -> Result<(CsrMatrix, Vec<f64>)> {
        let n = self.mesh.vertices.len();
        let mut matrix = CooBuilder::new(n, n);
        let mut rhs = vec![0.0; n];
        for (i, cell) in self.mesh.cells.iter().enumerate() {
            let local = diffusion_reaction_p1(
                self.cell_points(CellId(i))?,
                |_| DiffusionTensor::identity(),
                |_| 0.0,
                &source,
            )
            .map_err(|e| Error(format!("cell {}: {}", i, e)))?;
            for a in 0..3 {
                rhs[cell[a].0] += local.rhs[a];
                if !rhs[cell[a].0].is_finite() {
                    return Err(error("nonfinite assembled right hand side"));
                }
                for b in 0..3 {
                    matrix.add(cell[a].0, cell[b].0, local.matrix[a][b])?;
                }
            }
        }
        Ok((matrix.into_csr()?, rhs))
    }
}

/// Root-level transition view used by [`NonconformingRefinementHistory::refine`]
/// while the history is still at level zero.
///
/// A bare [`NonconformingMesh`] carries the mesh-level fields but not the
/// level-scoped transition provenance that [`NonconformingRefinement::refine`]
/// reads: oriented exterior edges and sibling groups. At a history root those
/// are derived from the mesh itself — boundary orientation comes from the
/// incident active cell, boundary edges parent to themselves, and every active
/// cell is its own (uncoarsenable) sibling group. `coarsen` is never called on
/// this view: the history removes levels by popping, not by mesh coarsening.
fn root_transition(mesh: &NonconformingMesh) -> Result<NonconformingRefinement> {
    let mut oriented_boundary_edges = Vec::with_capacity(mesh.boundary_edges.len());
    for &[a, b] in &mesh.boundary_edges {
        let oriented = mesh
            .cells
            .iter()
            .find_map(|cell| {
                [(cell[0], cell[1]), (cell[1], cell[2]), (cell[2], cell[0])]
                    .iter()
                    .copied()
                    .find(|&(x, y)| edge_key(x, y) == edge_key(a, b))
            })
            .ok_or_else(|| error("boundary edge has no incident active cell"))?;
        oriented_boundary_edges.push([oriented.0, oriented.1]);
    }
    let children = (0..mesh.cells.len()).map(|i| vec![CellId(i)]).collect();
    let parent_boundary_edges = (0..mesh.boundary_edges.len()).map(BoundaryEdgeId).collect();
    Ok(NonconformingRefinement {
        mesh: mesh.clone(),
        children,
        introduced_midpoints: BTreeMap::new(),
        oriented_boundary_edges,
        parent_boundary_edges,
        boundary_vertices: mesh.boundary_vertices.clone(),
        original_vertices: mesh.original_vertices,
    })
}

/// A stack of nonconforming transitions retaining every mesh level (RFC-0001 §9).
///
/// Mirrors [`crate::grid::RefinementHistory`] for the nonconforming seam: each
/// level is stored as a full [`NonconformingRefinement`] transition, so
/// level-scoped provenance (sibling groups, boundary ancestry, midpoint
/// parents) stays with its level. [`mesh`](Self::mesh) and
/// [`mesh_at`](Self::mesh_at) expose only the [`NonconformingMesh`] view.
///
/// [`refine`](Self::refine) appends one transition, rejecting empty selections
/// so successive levels stay distinguishable.
/// [`coarsen_last`](Self::coarsen_last) removes the entire latest level; it is
/// not selective element coarsening and cannot remove the base level. IDs,
/// vertex fields and constraints belong to their level and must not be reused
/// at another level.
#[derive(Debug, Clone)]
pub struct NonconformingRefinementHistory {
    base: NonconformingMesh,
    transitions: Vec<NonconformingRefinement>,
}

impl NonconformingRefinementHistory {
    /// Start with a root-level (level zero) mesh. Geometric conformity and
    /// nonoverlap remain caller obligations, as for [`Mesh2`]. Structural
    /// invariants are checked in debug builds; use
    /// [`try_new`](Self::try_new) for untrusted input.
    pub fn new(mesh: NonconformingMesh) -> Self {
        debug_assert!(mesh.validate().is_ok(), "invalid NonconformingMesh");
        Self {
            base: mesh,
            transitions: Vec::new(),
        }
    }

    /// Fallible constructor for a root-level (level zero) mesh: checks the
    /// structural invariants documented on [`NonconformingMesh`] before
    /// storing it. Geometric conformity is still a caller obligation.
    pub fn try_new(mesh: NonconformingMesh) -> Result<Self> {
        mesh.validate()?;
        Ok(Self {
            base: mesh,
            transitions: Vec::new(),
        })
    }

    /// Index of the current level; the mesh passed to [`new`](Self::new) is
    /// level zero.
    pub fn level(&self) -> usize {
        self.transitions.len()
    }

    /// Mesh at the current level. IDs are meaningful only at this level.
    pub fn mesh(&self) -> &NonconformingMesh {
        self.transitions
            .last()
            .map_or(&self.base, NonconformingRefinement::mesh)
    }

    /// Mesh at a retained level, from zero through [`level`](Self::level).
    pub fn mesh_at(&self, level: usize) -> Result<&NonconformingMesh> {
        if level == 0 {
            return Ok(&self.base);
        }
        self.transitions
            .get(level - 1)
            .map(NonconformingRefinement::mesh)
            .ok_or_else(|| error("refinement level out of bounds"))
    }

    /// Refine the current level, appending the resulting transition. Errors —
    /// including an empty selection — leave the history unchanged.
    pub fn refine(&mut self, marked: &[CellId]) -> Result<()> {
        self.mesh().validate()?;
        if marked.is_empty() {
            return Err(error("cannot append an empty refinement"));
        }
        let next = match self.transitions.last() {
            Some(current) => current.refine(marked)?,
            None => root_transition(&self.base)?.refine(marked)?,
        };
        self.transitions.push(next);
        Ok(())
    }

    /// Remove the entire latest refinement level. Errors at level zero leave
    /// the history unchanged. Retained old fields are caller-managed.
    pub fn coarsen_last(&mut self) -> Result<()> {
        if self.transitions.pop().is_none() {
            return Err(error("cannot coarsen the base level"));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RFC-0001 §1: `mesh()` exposes the seam fields; the refinement
    /// accessors delegate to the mesh or the retained transition provenance.
    #[test]
    fn mesh_exposes_seam_fields_and_delegates() {
        let parent = Mesh2::unit_square(1).unwrap();
        let r = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
        let mesh = r.mesh();
        assert_eq!(mesh.vertices.as_slice(), r.vertices());
        assert_eq!(mesh.cells.as_slice(), r.cells());
        assert_eq!(mesh.hanging.as_slice(), r.hanging_nodes());
        assert_eq!(mesh.boundary_edges.as_slice(), r.boundary_edges());
        assert_eq!(mesh.boundary_vertices.as_slice(), r.boundary_vertices());
        assert_eq!(mesh.original_vertices, parent.vertices().len());
        assert_eq!(mesh.cells.len(), 5);
        // Children of the marked parent point at it.
        assert!(mesh.parent_cells[0..4].iter().all(|&c| c == CellId(0)));
        // The carried cell keeps the historical parent convention.
        assert_eq!(mesh.parent_cells[4], CellId(1));
        // Transition provenance is retained per RFC-0001 §9.
        assert_eq!(r.children(CellId(0)).unwrap().len(), 4);
        assert_eq!(r.children(CellId(1)).unwrap().len(), 1);
        assert_eq!(r.oriented_boundary_edges.len(), r.boundary_edges().len());
    }

    /// RFC-0001 §2: `edge_splits` is the canonical key -> midpoint map and
    /// `midpoint_parents` is its inverse; midpoints sit at averaged masters.
    #[test]
    fn edge_splits_and_midpoint_parents_agree() {
        let parent = Mesh2::unit_square(1).unwrap();
        let r = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
        let mesh = r.mesh();
        assert_eq!(mesh.edge_splits.len(), 3);
        for (&(a, b), &mid) in &mesh.edge_splits {
            assert_eq!(mesh.midpoint_parents[&mid], (a, b));
            let u = mesh.vertices[a.0];
            let v = mesh.vertices[b.0];
            let p = mesh.vertices[mid.0];
            assert!((p.x - 0.5 * (u.x + v.x)).abs() < 1e-15);
            assert!((p.y - 0.5 * (u.y + v.y)).abs() < 1e-15);
        }
        // The single hanging node is the midpoint of the split interface edge.
        let h = r.hanging_nodes()[0];
        assert_eq!(
            mesh.edge_splits[&edge_key(h.masters[0], h.masters[1])],
            h.vertex
        );
    }

    /// RFC-0001 §9: the transition retains parent boundary edges per refined
    /// exterior edge, and `parent_boundary_edge` stays on the transition.
    #[test]
    fn parent_boundary_edges_are_retained_on_transition() {
        let parent = Mesh2::unit_square(1).unwrap();
        let r = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
        assert_eq!(r.parent_boundary_edges.len(), r.boundary_edges().len());
        for (i, id) in r.parent_boundary_edges.iter().enumerate() {
            assert_eq!(r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap(), *id);
            assert!(*id < BoundaryEdgeId(parent.boundary_edges().len()));
        }
    }

    /// RFC-0001 §1/§3 validation seam: a malformed base mesh must be rejected
    /// at `History::refine` entry instead of panicking while creating midpoints
    /// (the private fields let this unit test bypass the `new` debug assert).
    #[test]
    fn history_refine_rejects_malformed_base_mesh() {
        let parent = Mesh2::unit_square(1).unwrap();
        let mut base = NonconformingRefinement::new(&parent, &[])
            .unwrap()
            .mesh()
            .clone();
        base.cells[0][0] = VertexId(usize::MAX);
        assert!(base.validate().is_err());
        let mut history = NonconformingRefinementHistory {
            base,
            transitions: Vec::new(),
        };

        assert!(history.refine(&[CellId(1)]).is_err());
        assert_eq!(history.level(), 0);
        assert!(history.mesh().cells[0].contains(&VertexId(usize::MAX)));
    }
}
