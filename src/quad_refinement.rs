//! Uniform refinement of isoparametric Q1 quadrilaterals.
use std::collections::BTreeMap;

use crate::base::Point2;
use crate::grid::{BoundaryEdgeId, CellId, VertexId};
use crate::quad::QuadMesh;
use crate::{error, Result};

#[derive(Debug, Clone, Copy)]
enum NewVertex {
    Edge([VertexId; 2]),
    Center([VertexId; 4]),
}

/// One global conforming refinement of an immutable quadrilateral mesh.
/// Each parent is subdivided into four Q1 cells by edge midpoints and its
/// reference-cell center. The parent mapping is restricted to each reference
/// quadrant, including on nonaffine cells. Mesh conformity and nonoverlap
/// are preconditions, just as for [`QuadMesh::new`].
///
/// Original vertices retain their indices, while cell and boundary edge IDs
/// are mesh-specific. Relations describe only this one refinement level;
/// no persistent hierarchy, local refinement, or tag propagation is provided.
#[derive(Debug, Clone)]
pub struct UniformRefinement {
    mesh: QuadMesh,
    original_vertices: usize,
    new_vertices: Vec<NewVertex>,
    boundary_parents: Vec<BoundaryEdgeId>,
}

fn edge_key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

impl UniformRefinement {
    /// Construct four children per cell, sharing midpoints across cells.
    /// Errors are returned if refined geometry is invalid or cannot be
    /// represented, including coordinate rounding that degenerates a child.
    pub fn new(parent: &QuadMesh) -> Result<Self> {
        let count = parent
            .cells()
            .len()
            .checked_mul(4)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = parent.vertices().to_vec();
        let mut new_vertices = Vec::new();
        let mut midpoints = BTreeMap::new();
        let mut cells = Vec::with_capacity(count);
        for &[a, b, c, d] in parent.cells() {
            let mut midpoint = |a: VertexId, b: VertexId| -> Result<VertexId> {
                let key = edge_key(a, b);
                if let Some(&id) = midpoints.get(&key) {
                    return Ok(id);
                }
                let p = vertices[a.0];
                let q = vertices[b.0];
                let point = Point2::new(0.5 * p.x + 0.5 * q.x, 0.5 * p.y + 0.5 * q.y)?;
                let id = VertexId(vertices.len());
                vertices.push(point);
                new_vertices.push(NewVertex::Edge([key.0, key.1]));
                midpoints.insert(key, id);
                Ok(id)
            };
            let ab = midpoint(a, b)?;
            let bc = midpoint(b, c)?;
            let cd = midpoint(c, d)?;
            let da = midpoint(d, a)?;
            let pts = [a, b, c, d].map(|v| parent.vertices()[v.0]);
            let center = Point2::new(
                pts.iter().map(|p| p.x * 0.25).sum(),
                pts.iter().map(|p| p.y * 0.25).sum(),
            )?;
            let m = VertexId(vertices.len());
            vertices.push(center);
            new_vertices.push(NewVertex::Center([a, b, c, d]));
            cells.extend_from_slice(&[
                [a, ab, m, da],
                [ab, b, bc, m],
                [m, bc, c, cd],
                [da, m, cd, d],
            ]);
        }
        let mesh = QuadMesh::new(vertices, cells)?;
        let mut parent_of_edge = BTreeMap::new();
        for (i, &[a, b]) in parent.boundary_edges().iter().enumerate() {
            let midpoint = *midpoints
                .get(&edge_key(a, b))
                .ok_or_else(|| error("missing boundary midpoint"))?;
            for (x, y) in [(a, midpoint), (midpoint, b)] {
                parent_of_edge.insert(edge_key(x, y), BoundaryEdgeId(i));
            }
        }
        let boundary_parents = mesh
            .boundary_edges()
            .iter()
            .map(|edge| {
                parent_of_edge
                    .get(&edge_key(edge[0], edge[1]))
                    .copied()
                    .ok_or_else(|| error("refinement produced an unexpected boundary edge"))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            mesh,
            original_vertices: parent.vertices().len(),
            new_vertices,
            boundary_parents,
        })
    }
    /// Refined mesh; its cell and boundary IDs are not parent IDs.
    pub fn mesh(&self) -> &QuadMesh {
        &self.mesh
    }
    /// Parent of a child cell, in the original mesh.
    pub fn parent_cell(&self, child: CellId) -> Result<CellId> {
        if child.0 >= self.mesh.cells().len() {
            return Err(error("cell index out of bounds"));
        }
        Ok(CellId(child.0 / 4))
    }
    /// Four child cells in reference order: bottom-left, bottom-right,
    /// top-right, top-left. IDs refer to the refined mesh.
    pub fn children(&self, parent: CellId) -> Result<[CellId; 4]> {
        if parent.0 >= self.mesh.cells().len() / 4 {
            return Err(error("cell index out of bounds"));
        }
        let i = parent.0 * 4;
        Ok([CellId(i), CellId(i + 1), CellId(i + 2), CellId(i + 3)])
    }
    /// Original boundary edge containing the given refined boundary edge.
    pub fn parent_boundary_edge(&self, child: BoundaryEdgeId) -> Result<BoundaryEdgeId> {
        self.boundary_parents
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Prolong a continuous scalar Q1 nodal field onto the refined mesh.
    /// Edge values are averages of two parent nodal values; center values are
    /// averages of four. The same Q1 field is represented on every child in
    /// exact arithmetic, including on nonaffine bilinear cells. Inputs and
    /// outputs use their respective mesh vertex order; nonfinite data or
    /// arithmetic are rejected.
    pub fn prolongate_q1(&self, parent_values: &[f64]) -> Result<Vec<f64>> {
        if parent_values.len() != self.original_vertices
            || parent_values.iter().any(|x| !x.is_finite())
        {
            return Err(error("invalid parent Q1 field"));
        }
        let mut result = parent_values.to_vec();
        for origin in &self.new_vertices {
            let value = match *origin {
                NewVertex::Edge([a, b]) => parent_values[a.0] * 0.5 + parent_values[b.0] * 0.5,
                NewVertex::Center(ids) => ids.iter().map(|v| parent_values[v.0] * 0.25).sum(),
            };
            if !value.is_finite() {
                return Err(error("nonfinite interpolated Q1 value"));
            }
            result.push(value);
        }
        Ok(result)
    }
}
