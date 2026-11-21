//! Pure Neumann Poisson on P1 tetrahedra, with integral mean zero per component.
use crate::base::Point3;
use crate::elliptic_tetra::{self, DiffusionTensor3};
use crate::grid::{CellId, VertexId};
use crate::lac::{CgReport, CsrMatrix};
pub use crate::neumann::NeumannOptions;
use crate::tetra::{self, AffineTet, BoundaryFaceId, TetMesh};
use crate::{error, neumann_core, Result};

/// Diagnostics for a connected tetrahedral component.
#[derive(Debug, Clone, Copy)]
pub struct ComponentReport {
    /// Temporarily pinned vertex ID in the original mesh.
    pub anchor: VertexId,
    /// Volume, equal to the sum of P1 integral weights in this component.
    pub volume: f64,
    /// Sum of the original assembled loads on the component.
    pub load_imbalance: f64,
    /// Allowed absolute imbalance based on the input tolerances.
    pub compatibility_tolerance: f64,
}
/// Normalized solution and residuals of both pinned and original systems.
#[derive(Debug)]
pub struct NeumannSolution {
    pub values: Vec<f64>,
    /// CG diagnostics for the temporarily pinned system.
    pub cg: CgReport,
    /// Euclidean norm of `b - A u` for the original, unpinned system.
    pub original_residual_norm: f64,
    pub components: Vec<ComponentReport>,
}

/// Solve `-Δu = source` with outward `∂u/∂n = flux` on the entire boundary.
/// On each connected component of the P1 DoF graph, checks the *assembled*
/// compatibility `|sum(b_i)| <= abs_tol + rel_tol * sum(|b_i|)`; quadrature
/// error is included. It rejects incompatibility rather than projecting data.
/// Temporarily pins one vertex per component, solves by CG, then subtracts the
/// FEM integral mean (not the nodal mean). Components are ordered by anchor.
/// The reported true residual of the original system can be larger than the
/// pinned residual when an imbalance is accepted within tolerance.
///
/// The mesh must be conforming/nonoverlapping; scalar Laplacian and full
/// Neumann boundary only. This is not a solver for general reaction/diffusion.
/// Defaults and tolerance semantics are shared with [`NeumannOptions`].
pub fn solve<F, G>(
    mesh: &TetMesh,
    source: F,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    F: Fn(Point3) -> f64,
    G: Fn(Point3, [f64; 3]) -> f64,
{
    neumann_core::validate(options)?;
    let (matrix, rhs) = tetra::assemble_poisson(mesh, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

/// Solve `-div(A∇u)=source` with full-boundary conormal flux `n·A∇u`,
/// without reaction. Coefficient and source callbacks receive the cell ID
/// and physical quadrature point; the flux callback receives point and outward
/// normal and must supply the actual conormal flux, including material data.
/// It is not inferred from the cellwise diffusion tensor.
///
/// Uses the same discrete compatibility checks, pinning and P1 integral-mean
/// normalization as [`solve`]. The caller must provide a conforming,
/// nonoverlapping mesh and compatible boundary/interface data. The checks in
/// [`elliptic_tetra::assemble`] establish positivity only at quadrature points,
/// not global ellipticity or good conditioning. Quadrature error affects
/// compatibility. This is not a reaction solver.
pub fn solve_diffusion<D, F, G>(
    mesh: &TetMesh,
    diffusion: D,
    source: F,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    D: Fn(CellId, Point3) -> DiffusionTensor3,
    F: Fn(CellId, Point3) -> f64,
    G: Fn(Point3, [f64; 3]) -> f64,
{
    neumann_core::validate(options)?;
    let (matrix, rhs) = elliptic_tetra::assemble(mesh, diffusion, |_, _| 0.0, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

fn finish<G>(
    mesh: &TetMesh,
    matrix: CsrMatrix,
    mut rhs: Vec<f64>,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    G: Fn(Point3, [f64; 3]) -> f64,
{
    let faces: Vec<_> = (0..mesh.boundary_faces().len())
        .map(BoundaryFaceId)
        .collect();
    let load = tetra::assemble_neumann(mesh, &faces, flux)?;
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
        if !b.is_finite() {
            return Err(error("nonfinite Neumann right hand side"));
        }
    }
    let mut weights = vec![0.0; mesh.vertices().len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        // det J = 6 * volume, ∫φ_i = volume/4.
        let weight = AffineTet::new(mesh.cell_points(CellId(i))?)?.determinant() / 24.0;
        if !weight.is_finite() || weight <= 0.0 {
            return Err(error("invalid finite-element integration weight"));
        }
        for v in cell {
            weights[v.0] += weight;
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
                volume: c.measure,
                load_imbalance: c.load_imbalance,
                compatibility_tolerance: c.compatibility_tolerance,
            })
            .collect(),
    })
}
