//! One-level P1 red refinement with explicitly represented hanging nodes.
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

/// A midpoint DoF on an unsplit parent edge, constrained to the arithmetic
/// mean of its endpoints for continuous P1 fields. All IDs belong to the
/// refined object, even though original vertex indices are retained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HangingNode {
    pub vertex: VertexId,
    pub masters: [VertexId; 2],
}

/// A single nonconforming refinement of a conforming triangular mesh.
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
/// This type handles a single level only. It does not support refinement of
/// an already nonconforming mesh, coarsening, hanging-node chains, general
/// curved geometry, Neumann conditions, spatial lookup or error estimation.
/// Markers and solution vectors remain caller-managed.
#[derive(Debug, Clone)]
pub struct NonconformingRefinement {
    vertices: Vec<Point2>,
    cells: Vec<[VertexId; 3]>,
    parent_cells: Vec<CellId>,
    children: Vec<Vec<CellId>>,
    boundary_edges: Vec<[VertexId; 2]>,
    oriented_boundary_edges: Vec<[VertexId; 2]>,
    boundary_parents: Vec<BoundaryEdgeId>,
    boundary_vertices: Vec<VertexId>,
    original_vertices: usize,
    midpoint_parents: Vec<[VertexId; 2]>,
    hanging: Vec<HangingNode>,
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
        let mut midpoints = BTreeMap::new();
        let mut midpoint_parents = Vec::new();
        for &(a, b) in &split {
            let p = vertices[a.0];
            let q = vertices[b.0];
            let midpoint = Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y)?;
            let id = VertexId(vertices.len());
            vertices.push(midpoint);
            midpoint_parents.push([a, b]);
            midpoints.insert((a, b), id);
        }
        let mut cells = Vec::new();
        let mut parent_cells = Vec::new();
        let mut children = Vec::with_capacity(parent.cells().len());
        let mut hanging = BTreeMap::new();
        for (i, &[a, b, c]) in parent.cells().iter().enumerate() {
            let mut ids = Vec::new();
            if marks.contains(&CellId(i)) {
                let ab = midpoints[&edge_key(a, b)];
                let bc = midpoints[&edge_key(b, c)];
                let ca = midpoints[&edge_key(c, a)];
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
                    if let Some(&mid) = midpoints.get(&(v, w)) {
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
            if let Some(&mid) = midpoints.get(&edge_key(a, b)) {
                boundary.insert(edge_key(a, mid), ([a, mid], id));
                boundary.insert(edge_key(mid, b), ([mid, b], id));
            } else {
                boundary.insert(edge_key(a, b), ([a, b], id));
            }
        }
        let mut boundary_edges = Vec::new();
        let mut oriented_boundary_edges = Vec::new();
        let mut boundary_parents = Vec::new();
        let mut boundary_vertices = BTreeSet::new();
        for ((a, b), (oriented, parent_id)) in boundary {
            boundary_edges.push([a, b]);
            oriented_boundary_edges.push(oriented);
            boundary_parents.push(parent_id);
            boundary_vertices.insert(a);
            boundary_vertices.insert(b);
        }
        Ok(Self {
            vertices,
            cells,
            parent_cells,
            children,
            boundary_edges,
            oriented_boundary_edges,
            boundary_parents,
            boundary_vertices: boundary_vertices.into_iter().collect(),
            original_vertices: parent.vertices().len(),
            midpoint_parents,
            hanging: hanging.into_values().collect(),
        })
    }
    /// Vertex coordinates in this refined object's DoF order.
    pub fn vertices(&self) -> &[Point2] {
        &self.vertices
    }
    /// Positively oriented triangle connectivity, potentially nonconforming.
    pub fn cells(&self) -> &[[VertexId; 3]] {
        &self.cells
    }
    /// Exterior edges only, with canonically ordered endpoints. Edge IDs are
    /// indices into this list and belong to this object.
    pub fn boundary_edges(&self) -> &[[VertexId; 2]] {
        &self.boundary_edges
    }
    /// Unique exterior vertices; hanging nodes on internal edges are excluded.
    pub fn boundary_vertices(&self) -> &[VertexId] {
        &self.boundary_vertices
    }
    /// Midpoints constrained across a coarse/fine interface.
    pub fn hanging_nodes(&self) -> &[HangingNode] {
        &self.hanging
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
        self.parent_cells
            .get(child.0)
            .copied()
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Parent boundary edge for a refined exterior edge.
    pub fn parent_boundary_edge(&self, child: BoundaryEdgeId) -> Result<BoundaryEdgeId> {
        self.boundary_parents
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Cell coordinates, in local vertex order.
    pub fn cell_points(&self, id: CellId) -> Result<[Point2; 3]> {
        let cell = self
            .cells
            .get(id.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok(cell.map(|v| self.vertices[v.0]))
    }
    /// Oriented exterior edge, length and outward normal inherited from the
    /// parent boundary orientation. IDs and coordinates are those of this mesh.
    pub fn boundary_edge_geometry(&self, id: BoundaryEdgeId) -> Result<BoundaryEdgeGeometry> {
        let vertices = *self
            .oriented_boundary_edges
            .get(id.0)
            .ok_or_else(|| error("boundary edge index out of bounds"))?;
        let points = vertices.map(|v| self.vertices[v.0]);
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
        let mut builder = AffineConstraints::new(self.vertices.len());
        for h in &self.hanging {
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
                let value = f(self.vertices[v.0]);
                if !value.is_finite() {
                    Err(error("boundary function returned a nonfinite value"))
                } else {
                    Ok((v.0, value))
                }
            })
            .collect()
    }
    /// Prolong a scalar P1 parent field by midpoint interpolation. The result
    /// satisfies all hanging-node constraints in exact arithmetic.
    pub fn prolongate_p1(&self, values: &[f64]) -> Result<Vec<f64>> {
        if values.len() != self.original_vertices || values.iter().any(|v| !v.is_finite()) {
            return Err(error("invalid parent P1 field"));
        }
        let mut result = values.to_vec();
        for &[a, b] in &self.midpoint_parents {
            let v = 0.5 * values[a.0] + 0.5 * values[b.0];
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
        let n = self.vertices.len();
        let mut matrix = CooBuilder::new(n, n);
        let mut rhs = vec![0.0; n];
        for (i, cell) in self.cells.iter().enumerate() {
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
