//! P1 degree-of-freedom numbering.
use crate::grid::{CellId, Mesh2};
use crate::{error, Result};

/// Global degree-of-freedom index.
///
/// `DofId` is index-based and valid only for the [`DofHandler`] (and the
/// underlying [`Mesh2`]) it came from. Refinement, renumbering, or any other
/// structural change of that mesh invalidates it; no automatic remapping is
/// provided. Persistent or path-based identifiers are out of scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DofId(pub usize);

/// One DoF per mesh vertex. Borrowing prevents mesh mutation while in use.
///
/// The DoF identifiers produced by a handler are valid only for the mesh it
/// borrows; they are invalidated by refinement or renumbering of that mesh.
pub struct DofHandler<'a> {
    mesh: &'a Mesh2,
}
impl<'a> DofHandler<'a> {
    /// Borrow `mesh` for DoF numbering. The handler's ids are valid only for
    /// this mesh.
    pub fn new(mesh: &'a Mesh2) -> Self {
        Self { mesh }
    }
    /// Number of DoFs, one per mesh vertex.
    pub fn count(&self) -> usize {
        self.mesh.vertices().len()
    }
    /// DoFs of `cell`, in local vertex order.
    ///
    /// Returns an error if `cell` is not a valid cell id for this handler's
    /// mesh. The returned ids are valid only for this mesh and are invalidated
    /// by refinement or renumbering.
    pub fn cell_dofs(&self, cell: CellId) -> Result<[DofId; 3]> {
        let ids = self
            .mesh
            .cells()
            .get(cell.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([DofId(ids[0].0), DofId(ids[1].0), DofId(ids[2].0)])
    }
    /// Boundary DoFs, in boundary-vertex order.
    ///
    /// The returned ids are valid only for this mesh and are invalidated by
    /// refinement or renumbering.
    pub fn boundary_dofs(&self) -> Vec<DofId> {
        self.mesh
            .boundary_vertices()
            .iter()
            .map(|v| DofId(v.0))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{CellId, Mesh2};

    /// Shared vertices yield the same `DofId` in neighbouring cells.
    #[test]
    fn neighbouring_cells_share_dofs_at_common_vertices() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let dofs = DofHandler::new(&mesh);

        let first = dofs.cell_dofs(CellId(0)).unwrap();
        let second = dofs.cell_dofs(CellId(1)).unwrap();

        // The two triangles of the unit square share vertices 0 and 3.
        assert!(first.contains(&DofId(0)));
        assert!(second.contains(&DofId(0)));
        assert!(first.contains(&DofId(3)));
        assert!(second.contains(&DofId(3)));
    }

    /// Boundary DoFs are exactly the boundary vertices, in the same order.
    #[test]
    fn boundary_dofs_match_boundary_vertices() {
        let mesh = Mesh2::unit_square(2).unwrap();
        let dofs = DofHandler::new(&mesh);

        let boundary_dofs = dofs.boundary_dofs();
        let boundary_vertices = mesh.boundary_vertices();

        assert_eq!(boundary_dofs.len(), boundary_vertices.len());
        for (dof, vertex) in boundary_dofs.iter().zip(boundary_vertices) {
            assert_eq!(dof.0, vertex.0);
        }
    }

    /// Local vertex index maps to the same global `DofId` for every cell.
    #[test]
    fn cell_dofs_are_consistent_with_local_to_global_indices() {
        let mesh = Mesh2::unit_square(2).unwrap();
        let dofs = DofHandler::new(&mesh);

        assert_eq!(dofs.count(), mesh.vertices().len());

        for (i, cell) in mesh.cells().iter().enumerate() {
            let ids = dofs.cell_dofs(CellId(i)).unwrap();
            for (local, vertex) in cell.iter().enumerate() {
                assert_eq!(ids[local], DofId(vertex.0));
            }
        }
    }

    /// Operations taking ids validate bounds and return `Result` errors.
    #[test]
    fn cell_dofs_rejects_out_of_bounds_cell_id() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let dofs = DofHandler::new(&mesh);
        let n = mesh.cells().len();

        assert!(dofs.cell_dofs(CellId(n)).is_err());
        assert!(dofs.cell_dofs(CellId(usize::MAX)).is_err());
    }

    // No renumbering API exists in `DofHandler`: DoF ids are fixed to vertex
    // indices of the borrowed mesh. Invariance of an assembled solution under a
    // consistent renumbering therefore cannot be exercised here; the
    // local->global consistency test above documents the actual contract, and
    // refinement invalidates the ids as documented on `DofId`.
    #[test]
    fn dof_ids_are_fixed_to_vertex_indices_documenting_no_renumbering_api() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let dofs = DofHandler::new(&mesh);
        assert_eq!(dofs.cell_dofs(CellId(0)).unwrap()[0], DofId(0));
    }
}
