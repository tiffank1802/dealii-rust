//! Scalar P1 Poisson with pure Neumann boundary conditions.
//!
//! The constant nullspace is handled by temporarily pinning one vertex per
//! connected component of the discrete DoF graph, then subtracting the
//! finite-element integral mean on each component.

use crate::base::{triangle_determinant, Point2};
use crate::grid::{BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::lac::{CgOptions, CgReport, CsrMatrix};
use crate::numerics::{assemble_neumann, assemble_poisson};
use crate::{elliptic, integrators::DiffusionTensor};
use crate::{error, neumann_core, Result};

/// Options for checking discrete compatibility and solving the pinned system.
#[derive(Debug, Clone, Copy)]
pub struct NeumannOptions {
    /// Absolute tolerance for the sum of assembled loads on each component.
    pub compatibility_absolute_tolerance: f64,
    /// Relative tolerance multiplied by the sum of absolute nodal loads.
    pub compatibility_relative_tolerance: f64,
    /// Solver options apply to the temporarily pinned, nonsingular system.
    pub cg: CgOptions,
}

impl Default for NeumannOptions {
    fn default() -> Self {
        Self {
            compatibility_absolute_tolerance: 1e-12,
            compatibility_relative_tolerance: 1e-10,
            cg: CgOptions::default(),
        }
    }
}

/// Diagnostics for one connected component, ordered by increasing anchor ID.
#[derive(Debug, Clone, Copy)]
pub struct ComponentReport {
    /// Vertex temporarily fixed to zero before the mean is removed.
    pub anchor: VertexId,
    /// Area, also the sum of the P1 integral weights on this component.
    pub area: f64,
    /// Sum of the unmodified assembled loads on this component.
    pub load_imbalance: f64,
    /// Largest allowed absolute load imbalance for this component.
    pub compatibility_tolerance: f64,
}

/// Solution and diagnostics for the original and temporarily pinned systems.
#[derive(Debug)]
pub struct NeumannSolution {
    /// Nodal values, with zero finite-element integral mean per component.
    pub values: Vec<f64>,
    /// CG diagnostics for the pinned system, before mean normalization.
    pub cg: CgReport,
    /// Euclidean norm of `b - A u` for the original, unpinned system.
    pub original_residual_norm: f64,
    pub components: Vec<ComponentReport>,
}

/// Solve `-Δu = source` with `∂u/∂n = flux` on the whole boundary.
///
/// `flux` receives the physical point and outward unit normal. Zero flux may
/// be returned on any part of the boundary. A conforming, nonoverlapping
/// triangle mesh is required, just as for [`Mesh2::new`].
///
/// For each connected component of the P1 DoF graph, the assembled right-hand
/// side must satisfy `|sum(b_i)| <= abs_tol + rel_tol * sum(|b_i|)`. This is
/// a *discrete* test, including quadrature error; a continuously compatible
/// non-polynomial problem may need a different tolerance or quadrature.
/// Incompatible data are rejected, not silently projected. Data accepted within
/// the tolerance are left unchanged, so the residual of the original system
/// can exceed the pinned-system CG residual. Both quantities are reported.
///
/// The temporary pins make the P1 Laplace system positive definite, after
/// which the integral mean is removed separately on each component. No input
/// matrix or right-hand side is modified or exposed.
pub fn solve<F, G>(
    mesh: &Mesh2,
    source: F,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    F: Fn(Point2) -> f64,
    G: Fn(Point2, [f64; 2]) -> f64,
{
    neumann_core::validate(options)?;
    let (matrix, rhs) = assemble_poisson(mesh, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

/// Solve `-div(A∇u) = source` with conormal flux `g = n·A∇u` on the
/// entire boundary, no reaction term. `diffusion` and `source` receive the
/// cell ID and physical quadrature point; `flux` receives point and outward
/// normal and must provide the actual conormal flux (including material data).
/// In particular, boundary flux is not inferred from the cellwise tensor.
///
/// As in [`solve`], assembled compatibility is checked independently per
/// DoF component, then temporary pins and the FEM integral mean fix the
/// nullspace; no load projection is performed. Quadrature error is included
/// in the compatibility test. [`elliptic::assemble`] checks tensor positive
/// definiteness only at quadrature points, not global ellipticity or adequate
/// conditioning. The caller must provide a conforming, nonoverlapping mesh
/// and compatible boundary/interface data. This is not a reaction solver.
pub fn solve_diffusion<D, F, G>(
    mesh: &Mesh2,
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
    let (matrix, rhs) = elliptic::assemble(mesh, diffusion, |_, _| 0.0, source)?;
    finish(mesh, matrix, rhs, flux, options)
}

fn finish<G>(
    mesh: &Mesh2,
    matrix: CsrMatrix,
    mut rhs: Vec<f64>,
    flux: G,
    options: NeumannOptions,
) -> Result<NeumannSolution>
where
    G: Fn(Point2, [f64; 2]) -> f64,
{
    let boundary: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let flux_load = assemble_neumann(mesh, &boundary, flux)?;
    for (b, g) in rhs.iter_mut().zip(flux_load) {
        *b += g;
        if !b.is_finite() {
            return Err(error("nonfinite Neumann right hand side"));
        }
    }
    let mut weights = vec![0.0; mesh.vertices().len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let weight = triangle_determinant(mesh.cell_points(CellId(i))?)? / 6.0;
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
                area: c.measure,
                load_imbalance: c.load_imbalance,
                compatibility_tolerance: c.compatibility_tolerance,
            })
            .collect(),
    })
}
