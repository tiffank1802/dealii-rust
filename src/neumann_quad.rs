//! Pure Neumann Poisson on Q1 quadrilaterals, with integral mean zero per component.
use crate::base::Point2;
use crate::grid::{BoundaryEdgeId, CellId, VertexId};
use crate::lac::{CgReport, CsrMatrix};
pub use crate::neumann::NeumannOptions;
use crate::quad::{self, BilinearMap, QuadMesh};
use crate::quadrature::square_degree_5;
use crate::{elliptic_quad, integrators::DiffusionTensor};
use crate::{error, neumann_core, Result};

/// Diagnostics for one connected Q1 DoF component.
#[derive(Debug, Clone, Copy)]
pub struct ComponentReport {
    /// Temporarily pinned vertex, belonging to the input mesh.
    pub anchor: VertexId,
    /// Component area, equal to the sum of the Q1 integral weights.
    pub area: f64,
    /// Sum of the original assembled loads on this component.
    pub load_imbalance: f64,
    /// Allowed absolute imbalance based on the input tolerances.
    pub compatibility_tolerance: f64,
}

/// Normalized solution and diagnostics for the pinned and original systems.
#[derive(Debug)]
pub struct NeumannSolution {
    pub values: Vec<f64>,
    pub cg: CgReport,
    /// Euclidean norm of `b - A u` for the original, unpinned system.
    pub original_residual_norm: f64,
    pub components: Vec<ComponentReport>,
}

/// Solve `-Δu = source` with outward `∂u/∂n = flux` on the entire boundary.
///
/// Check the assembled compatibility `|sum(b_i)| <= abs_tol + rel_tol *
/// sum(|b_i|)` independently on each DoF-graph component. Incompatible loads
/// are rejected, not projected. Pin one vertex per component for CG, then
/// subtract the FEM integral mean, using `∫φ_i` computed with Gauss 3×3.
/// These Q1 weights need not be equal on a nonaffine quadrilateral. Components
/// are ordered by anchor; the original residual may exceed the pinned one
/// when a small imbalance is accepted. Compatibility includes quadrature error.
///
/// Requires a conforming, nonoverlapping mesh. Only scalar unit-coefficient
/// Laplacian and full Neumann boundary are supported; stiffness quadrature on
/// nonaffine cells is approximate. Tolerances are shared with [`NeumannOptions`].
pub fn solve<F, G>(
    mesh: &QuadMesh,
    source: F,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    F: Fn(Point2) -> f64,
    G: Fn(Point2, [f64; 2]) -> f64,
{
    neumann_core::validate(options)?;
    let (matrix, rhs) = quad::assemble_poisson(mesh, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

/// Solve `-div(A∇u)=source` with full-boundary conormal flux `n·A∇u`,
/// without reaction. Coefficient and source callbacks receive the cell ID
/// and physical quadrature point; the flux callback receives point and
/// outward normal and must supply the actual conormal flux, including material
/// data. It is not inferred from the cellwise diffusion tensor.
///
/// Uses the same discrete compatibility checks, pinning and Q1 integral-mean
/// normalization as [`solve`]. The caller must provide a conforming,
/// nonoverlapping mesh and compatible boundary/interface data. The checks in
/// [`elliptic_quad::assemble`] establish positivity only at quadrature points,
/// not global ellipticity or good conditioning. Nonaffine stiffness quadrature
/// and general data introduce errors. This is not a reaction solver.
pub fn solve_diffusion<D, F, G>(
    mesh: &QuadMesh,
    diffusion: D,
    source: F,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    D: Fn(CellId, Point2) -> DiffusionTensor,
    F: Fn(CellId, Point2) -> f64,
    G: Fn(Point2, [f64; 2]) -> f64,
{
    neumann_core::validate(options)?;
    let (matrix, rhs) = elliptic_quad::assemble(mesh, diffusion, |_, _| 0.0, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

fn finish<G>(
    mesh: &QuadMesh,
    matrix: CsrMatrix,
    mut rhs: Vec<f64>,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    G: Fn(Point2, [f64; 2]) -> f64,
{
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let load = quad::assemble_neumann(mesh, &edges, flux)?;
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
        if !b.is_finite() {
            return Err(error("nonfinite Neumann right hand side"));
        }
    }
    let mut weights = vec![0.0; mesh.vertices().len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i))?)?;
        for q in square_degree_5() {
            let mapped = map.evaluate(q.point)?;
            let factor = q.weight * mapped.determinant;
            if !factor.is_finite() || factor <= 0.0 {
                return Err(error("invalid finite-element integration weight"));
            }
            let shape = quad::shape_values(q.point);
            for a in 0..4 {
                weights[cell[a].0] += factor * shape[a];
            }
        }
    }
    let result = neumann_core::solve(matrix, rhs, weights, mesh.cells(), options)?;
    Ok(NeumannSolution {
        values: result.values,
        cg: result.cg,
        original_residual_norm: result.original_residual_norm,
        components: result
            .components
            .into_iter()
            .map(|c| ComponentReport {
                anchor: c.anchor,
                area: c.measure,
                load_imbalance: c.load_imbalance,
                compatibility_tolerance: c.compatibility_tolerance,
            })
            .collect(),
    })
}
