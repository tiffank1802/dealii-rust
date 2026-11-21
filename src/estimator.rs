//! Residual indicators for homogeneous-Dirichlet scalar Poisson with P1 triangles.
use std::collections::{BTreeMap, BTreeSet};

use crate::base::Point2;
use crate::fe::AffineTriangle;
use crate::grid::{CellId, Mesh2, NonconformingMesh, VertexId};
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

/// Canonical key for an undirected edge, matching `discover_hanging`.
fn edge_key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Collect the leaf sub-segments of a split master edge (RFC-0002 §3).
///
/// Mirrors `collect_subedges` (`nonconforming_refinement.rs`): the walk is
/// iterative and keeps a visited set, so malformed split maps cannot loop and
/// deep interfaces cannot exhaust the stack. A sub-segment is a leaf when its
/// key is absent from `edge_splits` or its midpoint coincides with an endpoint.
fn interface_leaves(
    edge_splits: &BTreeMap<(VertexId, VertexId), VertexId>,
    a: VertexId,
    b: VertexId,
) -> Vec<(VertexId, VertexId)> {
    let mut pending = vec![(a, b)];
    let mut seen = BTreeSet::new();
    let mut leaves = Vec::new();
    while let Some((x, y)) = pending.pop() {
        let key = edge_key(x, y);
        if !seen.insert(key) {
            continue;
        }
        match edge_splits.get(&key) {
            Some(&mid) if mid != key.0 && mid != key.1 => {
                pending.push((key.0, mid));
                pending.push((mid, key.1));
            }
            _ => leaves.push(key),
        }
    }
    leaves
}

/// Charge one interface leaf segment to its coarse and fine owners (RFC-0002 §3).
///
/// Sign convention: for a canonical leaf `(lo, hi)` with `lo < hi` and edge
/// vector `(dx, dy) = p_hi - p_lo`, the normal is `n = (dy, -dx)/l_leaf` (a
/// clockwise quarter turn, as in [`poisson_dirichlet`]). The jump is
/// `(∇u_coarse - ∇u_fine)·n`. Only `jump²` is charged, so reversing the normal
/// changes neither contribution. The coarse side uses the master length
/// `L_master`, the fine side the leaf length `l_leaf`; when `L_master == l_leaf`
/// both reduce to the conforming `0.5 · h_e² · jump²` term.
fn interface_leaf_contributions(
    leaf: (VertexId, VertexId),
    vertices: &[Point2],
    master_length: f64,
    grad_coarse: [f64; 2],
    grad_fine: [f64; 2],
) -> Result<(f64, f64)> {
    let (lo, hi) = leaf;
    let p = vertices[lo.0];
    let q = vertices[hi.0];
    let dx = q.x - p.x;
    let dy = q.y - p.y;
    let length = dx.hypot(dy);
    if !length.is_finite() || length <= 0.0 {
        return Err(error("invalid interface leaf length"));
    }
    let jump = (grad_coarse[0] - grad_fine[0]) * (dy / length)
        - (grad_coarse[1] - grad_fine[1]) * (dx / length);
    if !jump.is_finite() {
        return Err(error("nonfinite interface jump"));
    }
    let square = jump * jump;
    let coarse = 0.5 * master_length * length * square;
    let fine = 0.5 * length * length * square;
    if !coarse.is_finite() || !fine.is_finite() {
        return Err(error("nonfinite interface indicator"));
    }
    Ok((coarse, fine))
}

/// Residual indicator for the active cells of a nonconforming P1 mesh.
///
/// The cellwise formula, source quadrature and gradient evaluation match
/// [`poisson_dirichlet`] exactly. Edge incidence is built from the active cells
/// with canonical vertex keys, as `discover_hanging` does. An ordinary interior
/// edge (two incident active cells) keeps the frozen conforming jump term. A
/// split edge with one incident active cell is a coarse/fine interface root:
/// its sub-segments are walked through `edge_splits` to leaves, each leaf must
/// have exactly one incident fine cell, and the two-sided jump is charged with
/// the RFC-0002 §3 per-side `h_e` convention. Inconsistent incidence is an
/// error, never a skipped contribution.
///
/// The caller must pass the **reconstructed** solution vector: every
/// `mesh.hanging` rule is checked to a `1e-12 · (1 + max|u|)` tolerance over
/// the slave and both masters. The solution is never reconstructed here, so a
/// bad solve cannot be masked.
pub fn poisson_dirichlet_nonconforming<F>(
    mesh: &NonconformingMesh,
    solution: &[f64],
    source: F,
) -> Result<ResidualEstimate>
where
    F: Fn(Point2) -> f64,
{
    mesh.validate()?;
    if solution.len() != mesh.vertices.len() || solution.iter().any(|v| !v.is_finite()) {
        return Err(error("invalid solution vector"));
    }
    // The contract is the reconstructed vector: do not silently reconstruct,
    // which would hide a bad solve and interfere with Dirichlet rules. The
    // tolerance is `1e-12 · (1 + max|u|)` over the slave and both masters,
    // never an unscaled equality test.
    for h in &mesh.hanging {
        let slave = solution[h.vertex.0];
        let first = solution[h.masters[0].0];
        let second = solution[h.masters[1].0];
        let average = 0.5 * (first + second);
        if !average.is_finite() {
            return Err(error("nonfinite hanging node constraint value"));
        }
        let scale = 1.0 + slave.abs().max(first.abs()).max(second.abs());
        if (slave - average).abs() > 1e-12 * scale {
            return Err(error("solution does not satisfy hanging node constraints"));
        }
    }

    let mut squared = vec![0.0; mesh.cells.len()];
    let mut gradients = Vec::with_capacity(mesh.cells.len());
    let mut edges: BTreeMap<(VertexId, VertexId), Vec<usize>> = BTreeMap::new();
    let quadrature = triangle_degree_5();
    for (i, cell) in mesh.cells.iter().enumerate() {
        let points = [
            mesh.vertices[cell[0].0],
            mesh.vertices[cell[1].0],
            mesh.vertices[cell[2].0],
        ];
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
            let p = mesh.vertices[a.0];
            let q = mesh.vertices[b.0];
            let length = (q.x - p.x).hypot(q.y - p.y);
            if !length.is_finite() || length <= 0.0 {
                return Err(error("invalid edge length"));
            }
            diameter = diameter.max(length);
            edges.entry(edge_key(a, b)).or_default().push(i);
        }
        if !diameter.is_finite() || diameter <= 0.0 {
            return Err(error("invalid cell diameter"));
        }
        let mut integral = 0.0;
        for q in &quadrature {
            let value = source(mapping.map(q.point));
            if !value.is_finite() {
                return Err(error("source function returned a nonfinite value"));
            }
            integral += value * value * mapping.determinant() * q.weight;
            if !integral.is_finite() {
                return Err(error("nonfinite source integral"));
            }
        }
        let volume = diameter * diameter * integral;
        if !volume.is_finite() || volume < 0.0 {
            return Err(error("nonfinite volume indicator"));
        }
        squared[i] = volume;
    }
    // Borrow `edges` so leaf lookups can consult the same active-cell
    // incidence map used to find the roots.
    for (&(a, b), adjacent) in &edges {
        if adjacent.len() > 2 {
            return Err(error("edge has more than two incident active cells"));
        }
        if adjacent.len() == 2 {
            if mesh.edge_splits.contains_key(&(a, b)) {
                // A split master edge shared by two active cells has hanging
                // nodes on both sides with no fine cell to close the gap.
                return Err(error(&format!(
                    "split edge ({}, {}) has two incident active cells",
                    a.0, b.0
                )));
            }
            let p = mesh.vertices[a.0];
            let q = mesh.vertices[b.0];
            let dx = q.x - p.x;
            let dy = q.y - p.y;
            let length = dx.hypot(dy);
            if !length.is_finite() || length <= 0.0 {
                return Err(error("invalid edge length"));
            }
            let [g0, g1] = [gradients[adjacent[0]], gradients[adjacent[1]]];
            let jump = (g0[0] - g1[0]) * (dy / length) - (g0[1] - g1[1]) * (dx / length);
            if !jump.is_finite() {
                return Err(error("nonfinite edge jump"));
            }
            let contribution = 0.5 * length * length * jump * jump;
            if !contribution.is_finite() {
                return Err(error("nonfinite edge indicator"));
            }
            squared[adjacent[0]] += contribution;
            squared[adjacent[1]] += contribution;
            if !squared[adjacent[0]].is_finite() || !squared[adjacent[1]].is_finite() {
                return Err(error("nonfinite residual indicator"));
            }
            continue;
        }

        // One incident cell: either an exterior edge (no jump term, as in the
        // conforming estimator) or an interface root when the key is split.
        let &mid = match mesh.edge_splits.get(&(a, b)) {
            Some(mid) => mid,
            None => continue,
        };
        // `mesh.validate()` does not relate `edge_splits` to `hanging`, so a
        // hand-built mesh could charge interface jumps while the slave values
        // are unconstrained. Check the constraint behind this root only:
        // boundary split midpoints are Dirichlet DoFs and never appear as
        // roots, so they are intentionally outside this guard.
        if !mesh
            .hanging
            .iter()
            .any(|h| h.vertex == mid && edge_key(h.masters[0], h.masters[1]) == (a, b))
        {
            return Err(error(&format!(
                "interface edge ({}, {}) split at {} has no matching hanging node",
                a.0, b.0, mid.0
            )));
        }
        let coarse = adjacent[0];
        let p = mesh.vertices[a.0];
        let q = mesh.vertices[b.0];
        let dx = q.x - p.x;
        let dy = q.y - p.y;
        let master_length = dx.hypot(dy);
        if !master_length.is_finite() || master_length <= 0.0 {
            return Err(error("invalid interface master length"));
        }
        for leaf in interface_leaves(&mesh.edge_splits, a, b) {
            let (lo, hi) = leaf;
            let fine_cells = edges.get(&leaf).ok_or_else(|| {
                error(&format!(
                    "interface leaf ({}, {}) has no incident active cell",
                    lo.0, hi.0
                ))
            })?;
            if fine_cells.len() != 1 {
                return Err(error(&format!(
                    "interface leaf ({}, {}) must have exactly one incident active cell",
                    lo.0, hi.0
                )));
            }
            let fine = fine_cells[0];
            if fine == coarse {
                return Err(error("interface leaf is incident to the coarse owner"));
            }
            let (coarse_charge, fine_charge) = interface_leaf_contributions(
                leaf,
                &mesh.vertices,
                master_length,
                gradients[coarse],
                gradients[fine],
            )?;
            squared[coarse] += coarse_charge;
            squared[fine] += fine_charge;
            if !squared[coarse].is_finite() || !squared[fine].is_finite() {
                return Err(error("nonfinite residual indicator"));
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{CellId, Mesh2, NonconformingRefinement};

    fn close(x: f64, y: f64) {
        assert!((x - y).abs() < 1e-10, "{} != {}", x, y);
    }

    #[test]
    fn step1_matches_conforming_without_interfaces() {
        let parent = Mesh2::unit_square(1).unwrap();
        let marks: Vec<CellId> = (0..parent.cells().len()).map(CellId).collect();
        let transition = NonconformingRefinement::new(&parent, &marks).unwrap();
        let mesh = transition.mesh();
        let conforming = Mesh2::new(mesh.vertices.clone(), mesh.cells.clone()).unwrap();
        let solution: Vec<f64> = mesh
            .vertices
            .iter()
            .map(|p| p.x * p.x + p.y * p.y)
            .collect();
        let source = |_p: Point2| -4.0;

        let expected = poisson_dirichlet(&conforming, &solution, source).unwrap();
        let actual = poisson_dirichlet_nonconforming(mesh, &solution, source).unwrap();

        assert_eq!(actual.cell_indicators.len(), mesh.cells.len());
        for (got, want) in actual.cell_indicators.iter().zip(&expected.cell_indicators) {
            close(*got, *want);
        }
        close(actual.global_indicator, expected.global_indicator);
    }

    #[test]
    fn step1_rejects_invalid_input() {
        let parent = Mesh2::unit_square(1).unwrap();
        let transition = NonconformingRefinement::new(&parent, &[]).unwrap();
        let mesh = transition.mesh();

        assert!(poisson_dirichlet_nonconforming(mesh, &[0.0], |_| 0.0).is_err());
        assert!(poisson_dirichlet_nonconforming(mesh, &[f64::NAN; 4], |_| 0.0).is_err());
        assert!(poisson_dirichlet_nonconforming(mesh, &[0.0; 4], |_| f64::NAN).is_err());
        assert!(poisson_dirichlet_nonconforming(mesh, &[0.0; 4], |_| f64::MAX).is_err());

        let mut invalid = mesh.clone();
        invalid.parent_cells.clear();
        assert!(poisson_dirichlet_nonconforming(&invalid, &[0.0; 4], |_| 0.0).is_err());
    }
}
