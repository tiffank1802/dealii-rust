//! Scalar P1 Poisson with pure Neumann boundary conditions.
//!
//! The constant nullspace is handled by temporarily pinning one vertex per
//! connected component of the discrete DoF graph, then subtracting the
//! finite-element integral mean on each component.

use std::collections::BTreeMap;

use crate::base::{triangle_determinant, Point2};
use crate::grid::{BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::lac::{conjugate_gradient, CgOptions, CgReport};
use crate::numerics::{assemble_neumann, assemble_poisson};
use crate::{error, Error, Result};

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

fn root(parents: &mut [usize], mut v: usize) -> usize {
    while parents[v] != v {
        parents[v] = parents[parents[v]];
        v = parents[v];
    }
    v
}

fn component_vertices(mesh: &Mesh2) -> Vec<Vec<usize>> {
    let mut parents: Vec<usize> = (0..mesh.vertices().len()).collect();
    let mut sizes = vec![1usize; parents.len()];
    for cell in mesh.cells() {
        for i in 1..3 {
            let a = root(&mut parents, cell[0].0);
            let b = root(&mut parents, cell[i].0);
            if a != b {
                if sizes[a] >= sizes[b] {
                    parents[b] = a;
                    sizes[a] += sizes[b];
                } else {
                    parents[a] = b;
                    sizes[b] += sizes[a];
                }
            }
        }
    }
    let mut groups = BTreeMap::<usize, Vec<usize>>::new();
    for vertex in 0..parents.len() {
        groups
            .entry(root(&mut parents, vertex))
            .or_default()
            .push(vertex);
    }
    let mut components: Vec<Vec<usize>> = groups.into_values().collect();
    components.sort_by_key(|vertices| vertices[0]);
    components
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
    let abs_tol = options.compatibility_absolute_tolerance;
    let rel_tol = options.compatibility_relative_tolerance;
    if !abs_tol.is_finite() || abs_tol < 0.0 || !rel_tol.is_finite() || rel_tol < 0.0 {
        return Err(error("invalid Neumann compatibility tolerance"));
    }

    let (matrix, mut rhs) = assemble_poisson(mesh, source)?;
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

    let components = component_vertices(mesh);
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

    let mut reports = Vec::with_capacity(components.len());
    let mut anchors = Vec::with_capacity(components.len());
    for vertices in &components {
        let mut imbalance = 0.0;
        let mut load_scale = 0.0;
        let mut area = 0.0;
        for &i in vertices {
            imbalance += rhs[i];
            load_scale += rhs[i].abs();
            area += weights[i];
        }
        let tolerance = abs_tol + rel_tol * load_scale;
        if !imbalance.is_finite()
            || !load_scale.is_finite()
            || !area.is_finite()
            || area <= 0.0
            || !tolerance.is_finite()
        {
            return Err(error("nonfinite Neumann compatibility or area"));
        }
        let anchor = VertexId(vertices[0]);
        if imbalance.abs() > tolerance {
            return Err(Error(format!("incompatible Neumann load on component anchored at vertex {}: imbalance {} exceeds tolerance {}", anchor.0, imbalance, tolerance)));
        }
        reports.push(ComponentReport {
            anchor,
            area,
            load_imbalance: imbalance,
            compatibility_tolerance: tolerance,
        });
        anchors.push((anchor.0, 0.0));
    }

    let (pinned, pinned_rhs) = matrix.with_dirichlet(&rhs, &anchors)?;
    let (mut values, cg) = conjugate_gradient(&pinned, &pinned_rhs, options.cg)?;
    for (vertices, report) in components.iter().zip(&reports) {
        let integral: f64 = vertices.iter().map(|&i| weights[i] * values[i]).sum();
        let mean = integral / report.area;
        if !mean.is_finite() {
            return Err(error("nonfinite Neumann solution mean"));
        }
        for &i in vertices {
            values[i] -= mean;
        }
    }
    let ax = matrix.mul(&values)?;
    let original_residual_norm = rhs
        .iter()
        .zip(ax.iter())
        .fold(0.0_f64, |norm, (b, a)| norm.hypot(b - a));
    if !original_residual_norm.is_finite() {
        return Err(error("nonfinite original Neumann residual"));
    }
    Ok(NeumannSolution {
        values,
        cg,
        original_residual_norm,
        components: reports,
    })
}
