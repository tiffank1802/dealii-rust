//! One-level conforming red-green refinement for marked P1 triangles.
use crate::base::Point2;
use crate::grid::{BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::{error, Result};
use std::collections::{BTreeMap, BTreeSet};

fn key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// One conforming local refinement of a triangular mesh.
///
/// Marked triangles are split at all three edge midpoints (red); unmarked
/// neighbors are split along whichever shared edges were marked (green).
/// Thus all resulting edges match, with no hanging nodes or affine constraints.
/// Unmarked cells away from marked edges remain unchanged. Parent and child
/// cell and boundary IDs belong to different meshes, even for unchanged cells.
/// Input geometry must already be conforming and nonoverlapping. Repeated
/// refinement is possible but no persistent hierarchy, coarsening, shape
/// regularity, or automatic physical-tag propagation is provided.
#[derive(Debug, Clone)]
pub struct LocalRefinement {
    mesh: Mesh2,
    original_vertices: usize,
    midpoint_parents: Vec<[VertexId; 2]>,
    parent_cells: Vec<CellId>,
    children: Vec<Vec<CellId>>,
    boundary_parents: Vec<BoundaryEdgeId>,
}
impl LocalRefinement {
    /// Refine marked cells and make adjacent cells conforming. Repeated or
    /// out-of-range cell IDs are rejected; an empty selection is allowed.
    /// Numerical degeneracy in any generated cell is an error.
    pub fn new(parent: &Mesh2, marked: &[CellId]) -> Result<Self> {
        let mut marked_ids = BTreeSet::new();
        let mut split_edges = BTreeSet::new();
        for &id in marked {
            if !marked_ids.insert(id) {
                return Err(error("duplicate marked cell"));
            }
            let cell = parent
                .cells()
                .get(id.0)
                .ok_or_else(|| error("marked cell index out of bounds"))?;
            for j in 0..3 {
                split_edges.insert(key(cell[j], cell[(j + 1) % 3]));
            }
        }
        parent
            .cells()
            .len()
            .checked_mul(4)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = parent.vertices().to_vec();
        let mut midpoint_parents = Vec::new();
        let mut midpoints = BTreeMap::new();
        for &(a, b) in &split_edges {
            let p = vertices[a.0];
            let q = vertices[b.0];
            let point = Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y)?;
            let id = VertexId(vertices.len());
            vertices.push(point);
            midpoint_parents.push([a, b]);
            midpoints.insert((a, b), id);
        }
        let mut cells = Vec::new();
        let mut parent_cells = Vec::new();
        let mut children = Vec::with_capacity(parent.cells().len());
        for (i, cell) in parent.cells().iter().enumerate() {
            let [a, b, c] = *cell;
            let mids = [
                midpoints.get(&key(a, b)).copied(),
                midpoints.get(&key(b, c)).copied(),
                midpoints.get(&key(c, a)).copied(),
            ];
            let mut local = Vec::with_capacity(4);
            match mids {
                [None, None, None] => local.push([a, b, c]),
                [Some(m), None, None] => local.extend([[a, m, c], [m, b, c]]),
                [None, Some(m), None] => local.extend([[b, m, a], [m, c, a]]),
                [None, None, Some(m)] => local.extend([[c, m, b], [m, a, b]]),
                [Some(ab), Some(bc), None] => local.extend([[b, bc, ab], [a, ab, c], [ab, bc, c]]),
                [None, Some(bc), Some(ca)] => local.extend([[c, ca, bc], [b, bc, a], [bc, ca, a]]),
                [Some(ab), None, Some(ca)] => local.extend([[a, ab, ca], [c, ca, b], [ca, ab, b]]),
                [Some(ab), Some(bc), Some(ca)] => {
                    local.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]])
                }
            }
            let mut ids = Vec::with_capacity(local.len());
            for triangle in local {
                ids.push(CellId(cells.len()));
                cells.push(triangle);
                parent_cells.push(CellId(i));
            }
            children.push(ids);
        }
        let mesh = Mesh2::new(vertices, cells)?;
        let mut boundary_map = BTreeMap::new();
        for (i, &[a, b]) in parent.boundary_edges().iter().enumerate() {
            let parent_id = BoundaryEdgeId(i);
            if let Some(&m) = midpoints.get(&key(a, b)) {
                boundary_map.insert(key(a, m), parent_id);
                boundary_map.insert(key(m, b), parent_id);
            } else {
                boundary_map.insert(key(a, b), parent_id);
            }
        }
        let boundary_parents = mesh
            .boundary_edges()
            .iter()
            .map(|&[a, b]| {
                boundary_map
                    .get(&key(a, b))
                    .copied()
                    .ok_or_else(|| error("refinement produced an unexpected boundary edge"))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            mesh,
            original_vertices: parent.vertices().len(),
            midpoint_parents,
            parent_cells,
            children,
            boundary_parents,
        })
    }
    /// Refined mesh; all of its IDs belong to this mesh only.
    pub fn mesh(&self) -> &Mesh2 {
        &self.mesh
    }
    /// Parent cell ID of a refined cell, including unchanged cells.
    pub fn parent_cell(&self, child: CellId) -> Result<CellId> {
        self.parent_cells
            .get(child.0)
            .copied()
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Children of one parent in deterministic order; one if unchanged.
    pub fn children(&self, parent: CellId) -> Result<&[CellId]> {
        self.children
            .get(parent.0)
            .map(Vec::as_slice)
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Parent boundary edge ID, whether or not that edge was split.
    pub fn parent_boundary_edge(&self, child: BoundaryEdgeId) -> Result<BoundaryEdgeId> {
        self.boundary_parents
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Prolong a scalar P1 field, preserving its piecewise affine function in
    /// exact arithmetic. Parent values must be finite and in parent vertex order.
    pub fn prolongate_p1(&self, parent_values: &[f64]) -> Result<Vec<f64>> {
        if parent_values.len() != self.original_vertices
            || parent_values.iter().any(|v| !v.is_finite())
        {
            return Err(error("invalid parent P1 field"));
        }
        let mut values = parent_values.to_vec();
        for &[a, b] in &self.midpoint_parents {
            let value = 0.5 * parent_values[a.0] + 0.5 * parent_values[b.0];
            if !value.is_finite() {
                return Err(error("nonfinite interpolated P1 value"));
            }
            values.push(value);
        }
        Ok(values)
    }
}
