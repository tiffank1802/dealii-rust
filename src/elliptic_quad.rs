//! Scalar anisotropic diffusion-reaction on conforming Q1 quadrilaterals.
use crate::base::Point2;
use crate::grid::CellId;
use crate::integrators::DiffusionTensor;
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quad::{self, BilinearMap, QuadMesh};
use crate::quadrature::square_degree_5;
use crate::{error, Error, Result};

/// Local volume matrix and load in the input Q1 vertex order.
#[derive(Debug, Clone)]
pub struct LocalSystem {
    pub matrix: [[f64; 4]; 4],
    pub rhs: [f64; 4],
}

/// Integrate `∫(∇φ_i · A∇φ_j + c φ_i φ_j)` and `∫f φ_i` over a convex
/// quadrilateral. Coefficients are sampled at 3×3 Gauss points: `A` must be
/// symmetric positive definite, `c` finite and nonnegative, and `f` finite
/// at these points. This does not guarantee positivity between samples,
/// uniform ellipticity, coercivity, or conditioning. On nonaffine cells, the
/// diffusion integrand is generally rational and quadrature is approximate.
/// Geometry is validated as in [`BilinearMap`]. No boundary terms are included.
pub fn diffusion_reaction_q1<D, C, F>(
    points: [Point2; 4],
    diffusion: D,
    reaction: C,
    source: F,
) -> Result<LocalSystem>
where
    D: Fn(Point2) -> DiffusionTensor,
    C: Fn(Point2) -> f64,
    F: Fn(Point2) -> f64,
{
    let map = BilinearMap::new(points)?;
    let mut matrix = [[0.0; 4]; 4];
    let mut rhs = [0.0; 4];
    for q in square_degree_5() {
        let mapped = map.evaluate(q.point)?;
        let tensor = diffusion(mapped.point);
        tensor.validate()?;
        let c = reaction(mapped.point);
        let f = source(mapped.point);
        if !c.is_finite() || c < 0.0 {
            return Err(error("reaction must be finite and nonnegative"));
        }
        if !f.is_finite() {
            return Err(error("source function returned a nonfinite value"));
        }
        let weight = mapped.determinant * q.weight;
        if !weight.is_finite() || weight <= 0.0 {
            return Err(error("invalid quadrature weight"));
        }
        let shape = quad::shape_values(q.point);
        for i in 0..4 {
            rhs[i] += f * shape[i] * weight;
            for j in i..4 {
                let [gx, gy] = mapped.gradients[j];
                let ag = [
                    tensor.xx * gx + tensor.xy * gy,
                    tensor.xy * gx + tensor.yy * gy,
                ];
                let term = mapped.gradients[i][0] * ag[0] + mapped.gradients[i][1] * ag[1];
                matrix[i][j] += (term + c * shape[i] * shape[j]) * weight;
            }
        }
    }
    for i in 0..4 {
        if !rhs[i].is_finite() {
            return Err(error("nonfinite local right hand side"));
        }
        for j in i..4 {
            if !matrix[i][j].is_finite() {
                return Err(error("nonfinite local matrix entry"));
            }
            matrix[j][i] = matrix[i][j];
        }
    }
    Ok(LocalSystem { matrix, rhs })
}

/// Assemble volume terms for `-div(A∇u) + cu = f`. Callbacks receive the
/// `CellId` (belonging to this mesh) and physical quadrature point. The matrix
/// is symmetric but not necessarily positive definite; no boundary conditions
/// are imposed. For mixed conditions, add `quad::assemble_neumann` with the
/// outward conormal flux `g=n·A∇u` before applying Dirichlet constraints.
/// Full Neumann with `c=0` has a constant nullspace per component; use
/// [`crate::neumann_quad::solve_diffusion`] with compatible data in this case.
/// Positive reaction can remove the nullspace, but this API does not guarantee it.
/// A conforming, nonoverlapping mesh is required. See [`diffusion_reaction_q1`]
/// for quadrature and coefficient validation limitations.
pub fn assemble<D, C, F>(
    mesh: &QuadMesh,
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
        let local = diffusion_reaction_q1(
            mesh.cell_points(id)?,
            |p| diffusion(id, p),
            |p| reaction(id, p),
            |p| source(id, p),
        )
        .map_err(|e| Error(format!("cell {}: {}", i, e)))?;
        for r in 0..4 {
            rhs[cell[r].0] += local.rhs[r];
            if !rhs[cell[r].0].is_finite() {
                return Err(error("nonfinite assembled right hand side"));
            }
            for c in 0..4 {
                matrix.add(cell[r].0, cell[c].0, local.matrix[r][c])?;
            }
        }
    }
    Ok((matrix.into_csr()?, rhs))
}
