//! Immutable, conforming two-dimensional triangle meshes.

use crate::base::{triangle_determinant, Point2};
pub use crate::local_refinement::LocalRefinement;
pub use crate::nonconforming_refinement::{HangingNode, NonconformingRefinement};
pub use crate::refinement_history::RefinementHistory;
use crate::{error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// Index into a mesh's vertices.
///
/// Typed IDs are index-based. A `VertexId` is valid only for the [`Mesh2`]
/// instance it came from: refinement, renumbering, or any other structural
/// change of that mesh invalidates it, and no automatic remapping is provided.
/// Persistent or path-based identifiers are out of scope.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VertexId(pub usize);
/// Index into a mesh's cells.
///
/// Typed IDs are index-based. A `CellId` is valid only for the [`Mesh2`]
/// instance it came from: refinement, renumbering, or any other structural
/// change of that mesh invalidates it, and no automatic remapping is provided.
/// In particular, after [`UniformRefinement`] the same numeric index addresses
/// a (different) child cell of the original mesh. Persistent or path-based
/// identifiers are out of scope.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CellId(pub usize);
/// Index into the boundary edges of a particular mesh.
///
/// Typed IDs are index-based. A `BoundaryEdgeId` is valid only for the
/// [`Mesh2`] instance it came from: refinement, renumbering, or any other
/// structural change of that mesh invalidates it, and no automatic remapping is
/// provided. In particular, after [`UniformRefinement`] the same numeric index
/// addresses a (different, half-length) boundary edge. Persistent or path-based
/// identifiers are out of scope.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BoundaryEdgeId(pub usize);

/// A boundary edge oriented as in its counterclockwise adjacent cell.
/// Its outward normal is the right-hand normal to that orientation.
#[derive(Debug, Copy, Clone)]
pub struct BoundaryEdgeGeometry {
    pub vertices: [VertexId; 2],
    pub points: [Point2; 2],
    pub outward_normal: [f64; 2],
    pub length: f64,
}

/// An immutable mesh whose triangles have counterclockwise orientation.
#[derive(Debug, Clone)]
pub struct Mesh2 {
    vertices: Vec<Point2>,
    cells: Vec<[VertexId; 3]>,
    boundary_edges: Vec<[VertexId; 2]>,
    oriented_boundary_edges: Vec<[VertexId; 2]>,
    boundary_vertices: Vec<VertexId>,
}

impl Mesh2 {
    /// Construct a mesh. Invalid indices, duplicate cells, degenerate or
    /// clockwise triangles, edges shared by more than two cells and shared
    /// edges with inconsistent orientations are rejected. Geometric
    /// intersections and disconnected vertex fans are not checked.
    pub fn new(vertices: Vec<Point2>, cells: Vec<[VertexId; 3]>) -> Result<Self> {
        if vertices.is_empty() || cells.is_empty() {
            return Err(error("mesh must contain vertices and cells"));
        }
        if vertices
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite())
        {
            return Err(error("mesh coordinates must be finite"));
        }
        // Keep the first cell's directed edge to recover an outward normal.
        let mut edges = BTreeMap::<(VertexId, VertexId), (usize, [VertexId; 2])>::new();
        let mut unique_cells = BTreeSet::new();
        let mut used = vec![false; vertices.len()];
        for cell in &cells {
            if cell.iter().any(|v| v.0 >= vertices.len()) {
                return Err(error("cell vertex index out of bounds"));
            }
            triangle_determinant([
                vertices[cell[0].0],
                vertices[cell[1].0],
                vertices[cell[2].0],
            ])?;
            let mut sorted = *cell;
            sorted.sort();
            if !unique_cells.insert(sorted) {
                return Err(error("duplicate cell"));
            }
            for v in cell {
                used[v.0] = true;
            }
            for &(a, b) in &[(cell[0], cell[1]), (cell[1], cell[2]), (cell[2], cell[0])] {
                let key = if a < b { (a, b) } else { (b, a) };
                let (count, first_orientation) = edges.entry(key).or_insert((0, [a, b]));
                *count += 1;
                if *count > 2 {
                    return Err(error("non-manifold edge"));
                }
                if *count == 2 && *first_orientation == [a, b] {
                    return Err(error("inconsistent orientation across shared edge"));
                }
            }
        }
        if used.iter().any(|v| !v) {
            return Err(error("mesh has unused vertices"));
        }
        let boundary_edges: Vec<_> = edges
            .iter()
            .filter(|(_, (count, _))| *count == 1)
            .map(|(&(a, b), _)| [a, b])
            .collect();
        let oriented_boundary_edges: Vec<_> = edges
            .values()
            .filter(|(count, _)| *count == 1)
            .map(|(_, vertices)| *vertices)
            .collect();
        let boundary_vertices: Vec<_> = boundary_edges
            .iter()
            .flat_map(|e| e.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            vertices,
            cells,
            boundary_edges,
            oriented_boundary_edges,
            boundary_vertices,
        })
    }

    /// Uniform unit square with `n` subdivisions on each axis, two triangles per square.
    pub fn unit_square(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(error("subdivisions must be positive"));
        }
        let side = n
            .checked_add(1)
            .ok_or_else(|| error("mesh size overflow"))?;
        let nv = side
            .checked_mul(side)
            .ok_or_else(|| error("mesh size overflow"))?;
        let nc = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(2))
            .ok_or_else(|| error("mesh size overflow"))?;
        let mut vertices = Vec::with_capacity(nv);
        for j in 0..side {
            for i in 0..side {
                vertices.push(Point2 {
                    x: i as f64 / n as f64,
                    y: j as f64 / n as f64,
                });
            }
        }
        let mut cells = Vec::with_capacity(nc);
        for j in 0..n {
            for i in 0..n {
                let a = j * side + i;
                let b = a + 1;
                let d = a + side;
                let c = d + 1;
                cells.push([VertexId(a), VertexId(b), VertexId(c)]);
                cells.push([VertexId(a), VertexId(c), VertexId(d)]);
            }
        }
        Self::new(vertices, cells)
    }

    pub fn vertices(&self) -> &[Point2] {
        &self.vertices
    }
    pub fn cells(&self) -> &[[VertexId; 3]] {
        &self.cells
    }
    /// Canonically sorted endpoint IDs. Their order is not an orientation.
    /// The position in this slice is the corresponding [`BoundaryEdgeId`].
    pub fn boundary_edges(&self) -> &[[VertexId; 2]] {
        &self.boundary_edges
    }
    /// Return the directed endpoints, length and unit outward normal.
    /// The edge ID is valid only for this mesh and must be in bounds.
    pub fn boundary_edge_geometry(&self, id: BoundaryEdgeId) -> Result<BoundaryEdgeGeometry> {
        let vertices = *self
            .oriented_boundary_edges
            .get(id.0)
            .ok_or_else(|| error("boundary edge index out of bounds"))?;
        let points = [self.vertices[vertices[0].0], self.vertices[vertices[1].0]];
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
    pub fn boundary_vertices(&self) -> &[VertexId] {
        &self.boundary_vertices
    }
    /// Return the three vertex coordinates of `id`.
    ///
    /// The id must be valid for this mesh; out-of-bounds ids are rejected.
    /// Like all typed IDs, `id` is invalidated by refinement or renumbering.
    pub fn cell_points(&self, id: CellId) -> Result<[Point2; 3]> {
        let cell = self
            .cells
            .get(id.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([
            self.vertices[cell[0].0],
            self.vertices[cell[1].0],
            self.vertices[cell[2].0],
        ])
    }
}

/// One uniform red refinement of an immutable triangle mesh.
///
/// Each parent produces four triangles, ordered by parent. Original vertices
/// retain their indices; each new vertex is the midpoint of one parent edge.
/// The two meshes are distinct: cell and boundary-edge IDs are not reusable.
/// No physical tags or other external metadata are transferred automatically.
#[derive(Debug, Clone)]
pub struct UniformRefinement {
    mesh: Mesh2,
    original_vertices: usize,
    midpoint_parents: Vec<[VertexId; 2]>,
    boundary_parents: Vec<BoundaryEdgeId>,
}

impl UniformRefinement {
    /// Construct a conforming midpoint refinement, returning an error if the
    /// resulting coordinates or triangles cannot be represented safely.
    pub fn new(parent: &Mesh2) -> Result<Self> {
        parent
            .cells
            .len()
            .checked_mul(4)
            .ok_or_else(|| error("refinement size overflow"))?;
        let mut vertices = parent.vertices.clone();
        let mut midpoint_parents = Vec::new();
        let mut midpoints = BTreeMap::new();
        let mut cells = Vec::new();
        for &[a, b, c] in &parent.cells {
            let mut midpoint = |a: VertexId, b: VertexId| -> Result<VertexId> {
                let key = if a < b { (a, b) } else { (b, a) };
                if let Some(&id) = midpoints.get(&key) {
                    return Ok(id);
                }
                let p = vertices[a.0];
                let q = vertices[b.0];
                // Weighted sum avoids overflowing when p and q are both large.
                let point = Point2::new(p.x * 0.5 + q.x * 0.5, p.y * 0.5 + q.y * 0.5)?;
                let id = VertexId(vertices.len());
                vertices.push(point);
                midpoint_parents.push([key.0, key.1]);
                midpoints.insert(key, id);
                Ok(id)
            };
            let ab = midpoint(a, b)?;
            let bc = midpoint(b, c)?;
            let ca = midpoint(c, a)?;
            cells.extend_from_slice(&[[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]]);
        }
        let mesh = Mesh2::new(vertices, cells)?;
        let mut parent_of_edge = BTreeMap::new();
        for (i, &[a, b]) in parent.boundary_edges.iter().enumerate() {
            let mid = midpoints[&(a, b)];
            for (x, y) in [(a, mid), (mid, b)] {
                let key = if x < y { (x, y) } else { (y, x) };
                parent_of_edge.insert(key, BoundaryEdgeId(i));
            }
        }
        let boundary_parents = mesh
            .boundary_edges
            .iter()
            .map(|edge| {
                parent_of_edge
                    .get(&(edge[0], edge[1]))
                    .copied()
                    .ok_or_else(|| error("refinement produced an unexpected boundary edge"))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            mesh,
            original_vertices: parent.vertices.len(),
            midpoint_parents,
            boundary_parents,
        })
    }
    /// The refined mesh. Its IDs refer to this mesh, not to its parent.
    pub fn mesh(&self) -> &Mesh2 {
        &self.mesh
    }
    /// Parent of a refined cell.
    pub fn parent_cell(&self, child: CellId) -> Result<CellId> {
        if child.0 >= self.mesh.cells.len() {
            return Err(error("cell index out of bounds"));
        }
        Ok(CellId(child.0 / 4))
    }
    /// Four children of a parent cell, in deterministic local order.
    pub fn children(&self, parent: CellId) -> Result<[CellId; 4]> {
        if parent.0 >= self.mesh.cells.len() / 4 {
            return Err(error("cell index out of bounds"));
        }
        let i = parent.0 * 4;
        Ok([CellId(i), CellId(i + 1), CellId(i + 2), CellId(i + 3)])
    }
    /// Parent boundary edge of a refined boundary edge.
    pub fn parent_boundary_edge(&self, child: BoundaryEdgeId) -> Result<BoundaryEdgeId> {
        self.boundary_parents
            .get(child.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Prolong a scalar continuous P1 nodal field by midpoint interpolation.
    /// The input must be in the original mesh's vertex order. The values of
    /// the same finite-element function are preserved exactly in real arithmetic.
    pub fn prolongate_p1(&self, parent_values: &[f64]) -> Result<Vec<f64>> {
        if parent_values.len() != self.original_vertices
            || parent_values.iter().any(|v| !v.is_finite())
        {
            return Err(error("invalid parent P1 field"));
        }
        let mut values = parent_values.to_vec();
        for &[a, b] in &self.midpoint_parents {
            let value = parent_values[a.0] * 0.5 + parent_values[b.0] * 0.5;
            if !value.is_finite() {
                return Err(error("nonfinite interpolated P1 value"));
            }
            values.push(value);
        }
        Ok(values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Point2;

    /// ADR-0003 §4: after uniform refinement the same numeric cell index
    /// addresses a different (child) cell, not the original.
    #[test]
    fn typed_cell_id_is_invalidated_by_uniform_refinement() {
        let coarse = Mesh2::unit_square(1).unwrap();
        let old = CellId(0);
        let old_points = coarse.cell_points(old).unwrap();

        let refinement = UniformRefinement::new(&coarse).unwrap();
        let fine = refinement.mesh();

        // The old numeric index still addresses a cell in the refined mesh...
        assert!(fine.cell_points(old).is_ok());
        // ...but that cell is a child of the original, hence a different cell.
        let new_points = fine.cell_points(old).unwrap();
        assert_ne!(new_points, old_points);
        assert_eq!(refinement.parent_cell(old).unwrap(), CellId(0));
    }

    /// ADR-0003 §4: boundary-edge ids are likewise not reusable across meshes.
    #[test]
    fn typed_boundary_edge_id_is_invalidated_by_uniform_refinement() {
        let coarse = Mesh2::unit_square(1).unwrap();
        let old = BoundaryEdgeId(0);
        let old_geom = coarse.boundary_edge_geometry(old).unwrap();

        let refinement = UniformRefinement::new(&coarse).unwrap();
        let fine = refinement.mesh();

        assert!(fine.boundary_edge_geometry(old).is_ok());
        let new_geom = fine.boundary_edge_geometry(old).unwrap();
        // The refined edge is half the length of its parent, so it is a
        // different edge even though the numeric index coincides.
        assert!((new_geom.length - old_geom.length).abs() > 1e-12);
    }

    /// Operations taking ids validate bounds and return `Result` errors.
    #[test]
    fn cell_points_rejects_out_of_bounds_id() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let n = mesh.cells().len();
        assert!(mesh.cell_points(CellId(n)).is_err());
        assert!(mesh.cell_points(CellId(usize::MAX)).is_err());
    }

    #[test]
    fn boundary_edge_geometry_rejects_out_of_bounds_id() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let n = mesh.boundary_edges().len();
        assert!(mesh.boundary_edge_geometry(BoundaryEdgeId(n)).is_err());
        assert!(mesh
            .boundary_edge_geometry(BoundaryEdgeId(usize::MAX))
            .is_err());
    }

    #[test]
    fn mesh_construction_rejects_out_of_bounds_vertex_id() {
        let p0 = Point2::new(0.0, 0.0).unwrap();
        let p1 = Point2::new(1.0, 0.0).unwrap();
        let p2 = Point2::new(0.0, 1.0).unwrap();
        let err = Mesh2::new(
            vec![p0, p1, p2],
            vec![[VertexId(0), VertexId(1), VertexId(99)]],
        );
        assert!(err.is_err());
    }

    #[test]
    fn uniform_refinement_helpers_reject_out_of_bounds_ids() {
        let coarse = Mesh2::unit_square(1).unwrap();
        let refinement = UniformRefinement::new(&coarse).unwrap();
        let fine_cells = refinement.mesh().cells().len();
        let coarse_cells = coarse.cells().len();

        assert!(refinement.parent_cell(CellId(fine_cells)).is_err());
        assert!(refinement.children(CellId(coarse_cells)).is_err());

        let fine_edges = refinement.mesh().boundary_edges().len();
        assert!(refinement
            .parent_boundary_edge(BoundaryEdgeId(fine_edges))
            .is_err());
    }
}
