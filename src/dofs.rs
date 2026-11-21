//! P1 degree-of-freedom numbering.
use crate::grid::{CellId, Mesh2};
use crate::{error, Result};

/// Global degree-of-freedom index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct DofId(pub usize);

/// One DoF per mesh vertex. Borrowing prevents mesh mutation while in use.
pub struct DofHandler<'a> {
    mesh: &'a Mesh2,
}
impl<'a> DofHandler<'a> {
    pub fn new(mesh: &'a Mesh2) -> Self {
        Self { mesh }
    }
    pub fn count(&self) -> usize {
        self.mesh.vertices().len()
    }
    pub fn cell_dofs(&self, cell: CellId) -> Result<[DofId; 3]> {
        let ids = self
            .mesh
            .cells()
            .get(cell.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([DofId(ids[0].0), DofId(ids[1].0), DofId(ids[2].0)])
    }
    pub fn boundary_dofs(&self) -> Vec<DofId> {
        self.mesh
            .boundary_vertices()
            .iter()
            .map(|v| DofId(v.0))
            .collect()
    }
}
