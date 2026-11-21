//! Residual indicators for homogeneous-Dirichlet scalar Poisson with P1 triangles.
use std::collections::BTreeMap;

use crate::base::Point2;
use crate::fe::AffineTriangle;
use crate::grid::{CellId, Mesh2, VertexId};
use crate::quadrature::triangle_degree_5;
use crate::{error, Result};

/// Cellwise indicators and their root-sum-square global value.
#[derive(Debug, Clone)]
pub struct ResidualEstimate {
    /// Square root of each nonnegative cell contribution, indexed by `CellId`.
    pub cell_indicators: Vec<f64>,
    /// Square root of the sum of squared cell indicators.
    pub global_indicator: f64,
}

/// Compute a residual indicator for `-Δu=f`, with homogeneous Dirichlet on
/// the entire boundary and constant unit diffusion. The caller must enforce
/// the boundary condition; it is not verified here.
///
/// For each cell K, `η_K² = h_K² ∫_K f² + 1/2 Σ_e h_e ∫_e [∂u_h/∂n]²`,
/// with the sum over *interior* edges, `h_K` the longest edge, and `h_e` the
/// edge length. For P1, Δu_h=0 inside each triangle and the jump is constant.
/// Degree-five triangle quadrature approximates the source integral. This is
/// an indicator, not a guaranteed bound or a deal.II-equivalent estimator.
/// Neumann boundary residuals, variable diffusion and nonconforming grids are
/// unsupported. As with `Mesh2`, conformity and non-intersection are assumed.
pub fn poisson_dirichlet<F>(mesh: &Mesh2, solution: &[f64], source: F) -> Result<ResidualEstimate>
where
    F: Fn(Point2) -> f64,
{
    if solution.len() != mesh.vertices().len() || solution.iter().any(|v| !v.is_finite()) {
        return Err(error("invalid solution vector"));
    }
    let mut squared = vec![0.0; mesh.cells().len()];
    let mut gradients = Vec::new();
    let mut edges: BTreeMap<(VertexId, VertexId), Vec<usize>> = BTreeMap::new();
    let quadrature = triangle_degree_5();
    for (i, cell) in mesh.cells().iter().enumerate() {
        let points = mesh.cell_points(CellId(i))?;
        let mapping = AffineTriangle::new(points)?;
        let basis_gradients = mapping.gradients();
        let mut grad = [0.0; 2];
        for j in 0..3 {
            for (d, g) in grad.iter_mut().enumerate() {
                *g += solution[cell[j].0] * basis_gradients[j][d];
            }
        }
        if grad.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite solution gradient"));
        }
        gradients.push(grad);
        let mut diameter: f64 = 0.0;
        for &(a, b) in &[(cell[0], cell[1]), (cell[1], cell[2]), (cell[2], cell[0])] {
            let p = mesh.vertices()[a.0];
            let q = mesh.vertices()[b.0];
            let length = (q.x - p.x).hypot(q.y - p.y);
            if !length.is_finite() || length <= 0.0 {
                return Err(error("invalid edge length"));
            }
            diameter = diameter.max(length);
            let key = if a < b { (a, b) } else { (b, a) };
            edges.entry(key).or_default().push(i);
        }
        let mut integral = 0.0;
        for q in &quadrature {
            let value = source(mapping.map(q.point));
            if !value.is_finite() {
                return Err(error("source function returned a nonfinite value"));
            }
            integral += value * value * mapping.determinant() * q.weight;
        }
        squared[i] = diameter * diameter * integral;
    }
    for ((a, b), adjacent) in edges {
        if adjacent.len() == 2 {
            let p = mesh.vertices()[a.0];
            let q = mesh.vertices()[b.0];
            let dx = q.x - p.x;
            let dy = q.y - p.y;
            let length = dx.hypot(dy);
            let [g0, g1] = [gradients[adjacent[0]], gradients[adjacent[1]]];
            let jump = (g0[0] - g1[0]) * (dy / length) - (g0[1] - g1[1]) * (dx / length);
            let contribution = 0.5 * length * length * jump * jump;
            squared[adjacent[0]] += contribution;
            squared[adjacent[1]] += contribution;
        }
    }
    if squared.iter().any(|x| !x.is_finite() || *x < 0.0) {
        return Err(error("nonfinite residual indicator"));
    }
    let sum: f64 = squared.iter().sum();
    if !sum.is_finite() {
        return Err(error("nonfinite global residual indicator"));
    }
    Ok(ResidualEstimate {
        cell_indicators: squared.into_iter().map(f64::sqrt).collect(),
        global_indicator: sum.sqrt(),
    })
}

/// Select a minimal prefix of cells whose squared indicators reach a fraction
/// `theta` of their sum (Dörfler/bulk marking). `theta` must be finite in
/// `[0,1]`; indicators must be finite and nonnegative, in cell order. Cells
/// are returned in decreasing indicator order, with ties broken by CellId.
/// Zero `theta` or all-zero indicators produce an empty selection. Scaling by
/// the largest indicator avoids overflow in squaring; very small contributions
/// can underflow and comparisons near the threshold depend on rounding. This
/// only selects cells: refinement, mesh conformity and error reduction are not
/// guaranteed by the marking rule.
pub fn mark_bulk(indicators: &[f64], theta: f64) -> Result<Vec<CellId>> {
    if !theta.is_finite()
        || !(0.0..=1.0).contains(&theta)
        || indicators.iter().any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(error("invalid indicators or bulk marking fraction"));
    }
    if theta == 0.0 || indicators.is_empty() {
        return Ok(Vec::new());
    }
    let mut indices: Vec<_> = (0..indicators.len()).collect();
    indices.sort_by(|&a, &b| indicators[b].total_cmp(&indicators[a]).then(a.cmp(&b)));
    let scale = indicators[indices[0]];
    if scale == 0.0 {
        return Ok(Vec::new());
    }
    let weights: Vec<f64> = indices
        .iter()
        .map(|&i| {
            let r = indicators[i] / scale;
            r * r
        })
        .collect();
    let total: f64 = weights.iter().sum();
    let target = theta * total;
    let mut sum = 0.0;
    let mut marked = Vec::new();
    for (&i, &weight) in indices.iter().zip(&weights) {
        if sum >= target {
            break;
        }
        sum += weight;
        marked.push(CellId(i));
    }
    Ok(marked)
}
