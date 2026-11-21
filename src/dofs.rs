//! P1 degree-of-freedom numbering.
use crate::grid::{CellId, Mesh2};
use crate::{error, Result};

/// Shared error message for connectivity accessors on a count-only handler.
const COUNT_ONLY_ERROR: &str = "no cell_dofs/boundary_dofs via this constructor";

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
///
/// # Count-only mode
///
/// [`DofHandler::from_vertex_count`] builds a handler for callers that know
/// only the vertex count of a nonconforming mesh. It supports [`count`](Self::count)
/// but carries no connectivity: [`cell_dofs`](Self::cell_dofs) and
/// [`try_boundary_dofs`](Self::try_boundary_dofs) return
/// `Err("no cell_dofs/boundary_dofs via this constructor")`, and
/// [`boundary_dofs`](Self::boundary_dofs) panics with the same message.
///
/// Nonconforming vertex-indexed DoFs come from the mesh itself: see
/// [`NonconformingMesh::vertices`](crate::grid::NonconformingMesh::vertices).
/// Hanging-node DoFs are numbered there as ordinary vertices and constrained
/// through `AffineConstraints`; no separate renumbering is needed.
pub struct DofHandler<'a> {
    /// `None` only for handlers made by [`from_vertex_count`](Self::from_vertex_count).
    mesh: Option<&'a Mesh2>,
    /// One DoF per vertex, stored so count-only handlers need no mesh.
    count: usize,
}
impl<'a> DofHandler<'a> {
    /// Borrow `mesh` for DoF numbering. The handler's ids are valid only for
    /// this mesh.
    pub fn new(mesh: &'a Mesh2) -> Self {
        Self {
            count: mesh.vertices().len(),
            mesh: Some(mesh),
        }
    }

    /// Build a count-only handler from a vertex count (RFC-0001 §6).
    ///
    /// Use this when only the number of vertices of a nonconforming mesh is
    /// available, for example at a history level whose mesh is not borrowed.
    /// The vertex-indexed DoFs themselves are the indices of
    /// [`NonconformingMesh::vertices`](crate::grid::NonconformingMesh::vertices);
    /// connectivity-based accessors are unavailable here.
    pub fn from_vertex_count(count: usize) -> Self {
        Self { mesh: None, count }
    }

    /// Number of DoFs, one per mesh vertex. Works in both modes.
    pub fn count(&self) -> usize {
        self.count
    }

    /// DoFs of `cell`, in local vertex order.
    ///
    /// Returns an error if `cell` is not a valid cell id for this handler's
    /// mesh. The returned ids are valid only for this mesh and are invalidated
    /// by refinement or renumbering.
    ///
    /// # Errors
    ///
    /// Returns `Err("no cell_dofs/boundary_dofs via this constructor")` when
    /// this handler came from [`from_vertex_count`](Self::from_vertex_count):
    /// a bare vertex count carries no cell connectivity.
    pub fn cell_dofs(&self, cell: CellId) -> Result<[DofId; 3]> {
        let mesh = self.mesh.ok_or_else(|| error(COUNT_ONLY_ERROR))?;
        let ids = mesh
            .cells()
            .get(cell.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([DofId(ids[0].0), DofId(ids[1].0), DofId(ids[2].0)])
    }

    /// Boundary DoFs, in boundary-vertex order.
    ///
    /// The returned ids are valid only for this mesh and are invalidated by
    /// refinement or renumbering.
    ///
    /// # Panics
    ///
    /// Panics with `"no cell_dofs/boundary_dofs via this constructor"` when
    /// this handler came from [`from_vertex_count`](Self::from_vertex_count).
    /// Use [`try_boundary_dofs`](Self::try_boundary_dofs) for a fallible result.
    pub fn boundary_dofs(&self) -> Vec<DofId> {
        self.try_boundary_dofs().expect(COUNT_ONLY_ERROR)
    }

    /// Fallible form of [`boundary_dofs`](Self::boundary_dofs).
    ///
    /// # Errors
    ///
    /// Returns `Err("no cell_dofs/boundary_dofs via this constructor")` for a
    /// handler built with [`from_vertex_count`](Self::from_vertex_count).
    pub fn try_boundary_dofs(&self) -> Result<Vec<DofId>> {
        let mesh = self.mesh.ok_or_else(|| error(COUNT_ONLY_ERROR))?;
        Ok(mesh
            .boundary_vertices()
            .iter()
            .map(|v| DofId(v.0))
            .collect())
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

    /// RFC-0001 §6: a count-only handler reports the stored count verbatim.
    #[test]
    fn from_vertex_count_reports_stored_count() {
        for n in [0, 1, 7, usize::MAX] {
            assert_eq!(DofHandler::from_vertex_count(n).count(), n);
        }
    }

    /// RFC-0001 §6: connectivity-based accessors are unavailable through the
    /// count-only constructor and report the documented error.
    #[test]
    fn from_vertex_count_rejects_connectivity_access() {
        let dofs = DofHandler::from_vertex_count(9);

        for cell in [CellId(0), CellId(usize::MAX)] {
            assert_eq!(
                dofs.cell_dofs(cell),
                Err(error(COUNT_ONLY_ERROR)),
                "cell {:?}",
                cell
            );
        }
        assert_eq!(dofs.try_boundary_dofs(), Err(error(COUNT_ONLY_ERROR)));
    }

    /// The infallible wrapper documents and enforces its count-only panic.
    #[test]
    #[should_panic(expected = "no cell_dofs/boundary_dofs via this constructor")]
    fn boundary_dofs_panics_for_count_only_handler() {
        DofHandler::from_vertex_count(3).boundary_dofs();
    }

    /// Mesh-backed handlers keep both boundary forms equivalent.
    #[test]
    fn mesh_backed_boundary_dofs_fallible_form_agrees() {
        let mesh = Mesh2::unit_square(1).unwrap();
        let dofs = DofHandler::new(&mesh);

        assert_eq!(dofs.try_boundary_dofs().unwrap(), dofs.boundary_dofs());
    }
}
