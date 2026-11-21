//! Uniform refinement of P1 tetrahedra.
use std::collections::BTreeMap;

use crate::base::Point3;
use crate::grid::{CellId, VertexId};
use crate::tetra::{AffineTet, BoundaryFaceId, TetMesh};
use crate::{error, Result};

fn edge_key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}
fn face_key(mut vertices: [VertexId; 3]) -> [VertexId; 3] {
    vertices.sort();
    vertices
}

/// One global conforming midpoint refinement of a tetrahedral mesh.
/// Each tetrahedron is split into four corner tetrahedra and a central
/// octahedron, which is split into four more tetrahedra. The octahedron's
/// shortest of its three diagonals is selected (ties follow local edge order).
/// This choice is deterministic but does not guarantee mesh quality. Shared
/// parent faces are split into the same four triangles regardless of that
/// choice. Conformity and nonoverlap are required of the input mesh.
///
/// Original vertex indices are retained, while cell and boundary face IDs
/// are mesh-specific. Only one refinement level is represented; no automatic
/// tag propagation, local adaptation or persistent hierarchy is provided.
#[derive(Debug, Clone)]
pub struct UniformRefinement {
    mesh: TetMesh,
    original_vertices: usize,
    midpoint_parents: Vec<[VertexId; 2]>,
    boundary_parents: Vec<BoundaryFaceId>,
}

impl UniformRefinement {
    /// Refine all cells, sharing edge midpoints across the whole mesh.
    /// Refined geometry is revalidated and positively oriented. Coordinates
    /// that cannot be represented, or child cells that fail the usual geometry
    /// validity threshold, cause an error. No quality improvement is guaranteed.
    pub fn new(parent: &TetMesh) -> Result<Self> {
        let count = parent
            .cells()
            .len()
            .checked_mul(8)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = parent.vertices().to_vec();
        let mut midpoint_parents = Vec::new();
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
                let point = Point3::new(
                    0.5 * p.x + 0.5 * q.x,
                    0.5 * p.y + 0.5 * q.y,
                    0.5 * p.z + 0.5 * q.z,
                )?;
                let id = VertexId(vertices.len());
                vertices.push(point);
                midpoint_parents.push([key.0, key.1]);
                midpoints.insert(key, id);
                Ok(id)
            };
            let ab = midpoint(a, b)?;
            let ac = midpoint(a, c)?;
            let ad = midpoint(a, d)?;
            let bc = midpoint(b, c)?;
            let bd = midpoint(b, d)?;
            let cd = midpoint(c, d)?;
            let diagonals = [
                (ab, cd, [ac, ad, bd, bc]),
                (ac, bd, [ab, ad, cd, bc]),
                (ad, bc, [ab, ac, cd, bd]),
            ];
            let length = |a: VertexId, b: VertexId| {
                let p = vertices[a.0];
                let q = vertices[b.0];
                (p.x - q.x).hypot(p.y - q.y).hypot(p.z - q.z)
            };
            let mut selected = diagonals[0];
            let mut shortest = length(selected.0, selected.1);
            for &choice in &diagonals[1..] {
                let candidate = length(choice.0, choice.1);
                if candidate < shortest {
                    selected = choice;
                    shortest = candidate;
                }
            }
            if !shortest.is_finite() || shortest <= 0.0 {
                return Err(error("invalid octahedron diagonal"));
            }
            let (e, f, [r0, r1, r2, r3]) = selected;
            let children = [
                [a, ab, ac, ad],
                [b, ab, bc, bd],
                [c, ac, bc, cd],
                [d, ad, bd, cd],
                [e, f, r0, r1],
                [e, f, r1, r2],
                [e, f, r2, r3],
                [e, f, r3, r0],
            ];
            for mut child in children {
                let points = child.map(|v| vertices[v.0]);
                if AffineTet::new(points).is_err() {
                    child.swap(0, 1);
                    AffineTet::new(child.map(|v| vertices[v.0]))?;
                }
                cells.push(child);
            }
        }
        let mesh = TetMesh::new(vertices, cells)?;
        let mut parent_of_face = BTreeMap::new();
        for (i, &[a, b, c]) in parent.boundary_faces().iter().enumerate() {
            let get = |a, b| {
                midpoints
                    .get(&edge_key(a, b))
                    .copied()
                    .ok_or_else(|| error("missing boundary midpoint"))
            };
            let ab = get(a, b)?;
            let bc = get(b, c)?;
            let ca = get(c, a)?;
            for face in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {
                parent_of_face.insert(face_key(face), BoundaryFaceId(i));
            }
        }
        let boundary_parents = mesh
            .boundary_faces()
            .iter()
            .map(|&face| {
                parent_of_face
                    .get(&face_key(face))
                    .copied()
                    .ok_or_else(|| error("refinement produced an unexpected boundary face"))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            mesh,
            original_vertices: parent.vertices().len(),
            midpoint_parents,
            boundary_parents,
        })
    }
    /// Refined mesh, whose cell and face IDs differ from the parent's.
    pub fn mesh(&self) -> &TetMesh {
        &self.mesh
    }
    /// Parent of a child cell, in the input mesh.
    pub fn parent_cell(&self, child: CellId) -> Result<CellId> {
        if child.0 >= self.mesh.cells().len() {
            return Err(error("cell index out of bounds"));
        }
        Ok(CellId(child.0 / 8))
    }
    /// Eight children in deterministic order: four corners (local vertices
    /// 0 to 3), then four around the chosen octahedron diagonal.
    pub fn children(&self, parent: CellId) -> Result<[CellId; 8]> {
        if parent.0 >= self.mesh.cells().len() / 8 {
            return Err(error("cell index out of bounds"));
        }
        let i = parent.0 * 8;
        Ok([
            CellId(i),
            CellId(i + 1),
            CellId(i + 2),
            CellId(i + 3),
            CellId(i + 4),
            CellId(i + 5),
            CellId(i + 6),
            CellId(i + 7),
        ])
    }
    /// Parent boundary face of a refined boundary face.
    pub fn parent_boundary_face(&self, child: BoundaryFaceId) -> Result<BoundaryFaceId> {
        self.boundary_parents
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary face index out of bounds"))
    }
    /// Prolong a continuous scalar P1 nodal field by edge-midpoint averages.
    /// Each child represents the same piecewise affine function as its parent
    /// in exact arithmetic. Inputs and outputs use their respective mesh's
    /// vertex order; nonfinite values or arithmetic are rejected.
    pub fn prolongate_p1(&self, parent_values: &[f64]) -> Result<Vec<f64>> {
        if parent_values.len() != self.original_vertices
            || parent_values.iter().any(|v| !v.is_finite())
        {
            return Err(error("invalid parent P1 field"));
        }
        let mut result = parent_values.to_vec();
        for &[a, b] in &self.midpoint_parents {
            let value = 0.5 * parent_values[a.0] + 0.5 * parent_values[b.0];
            if !value.is_finite() {
                return Err(error("nonfinite interpolated P1 value"));
            }
            result.push(value);
        }
        Ok(result)
    }
}
