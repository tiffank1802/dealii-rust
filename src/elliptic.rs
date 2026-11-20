//! Global P1 assembly for scalar symmetric diffusion-reaction problems.
use crate::base::Point2;
use crate::grid::{CellId, Mesh2};
use crate::integrators::{diffusion_reaction_p1, DiffusionTensor};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::{error, Error, Result};

/// Assemble the volume terms of `-div(A ∇u) + c u = f` on conforming
/// triangles. Coefficients receive both the cell ID and physical quadrature
/// point, allowing cellwise material data (e.g. Gmsh physical tags).
///
/// The local quadrature and sampled-data requirements are those of
/// [`diffusion_reaction_p1`]. No boundary condition is imposed. The result is
/// symmetric but is not automatically positive definite: pure diffusion with
/// only Neumann conditions has a constant nullspace per connected component.
/// Dirichlet constraints or a sufficiently coercive reaction are needed for
/// CG, together with suitable conditioning. For a Neumann edge, use
/// `numerics::assemble_neumann` with conormal flux `g = n · A ∇u`, add the
/// resulting vector to the volume load, then impose Dirichlet constraints.
/// This routine does not implement a nullspace-aware pure-Neumann solve.
pub fn assemble<D, C, F>(
    mesh: &Mesh2,
    diffusion: D,
    reaction: C,
    source: F,
) -> Result<(CsrMatrix, Vec<f64>)>
where
    D: Fn(CellId, Point2) -> DiffusionTensor,
    C: Fn(CellId, Point2) -> f64,
    F: Fn(CellId, Point2) -> f64,
{
    let n = mesh.vertices().len();
    let mut matrix = CooBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let id = CellId(i);
        let local = diffusion_reaction_p1(
            mesh.cell_points(id)?,
            |p| diffusion(id, p),
            |p| reaction(id, p),
            |p| source(id, p),
        )
        .map_err(|e| Error(format!("cell {}: {}", i, e)))?;
        for r in 0..3 {
            rhs[cell[r].0] += local.rhs[r];
            if !rhs[cell[r].0].is_finite() {
                return Err(error("nonfinite assembled right hand side"));
            }
            for c in 0..3 {
                matrix.add(cell[r].0, cell[c].0, local.matrix[r][c])?;
            }
        }
    }
    Ok((matrix.into_csr()?, rhs))
}
