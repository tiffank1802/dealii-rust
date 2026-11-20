//! Shared internal nullspace handling for nodal P1 and Q1 Laplacians.
use crate::grid::VertexId;
use crate::lac::{conjugate_gradient, CgReport, CsrMatrix};
use crate::neumann::NeumannOptions;
use crate::{error, Error, Result};
use std::collections::BTreeMap;

pub(crate) struct Report {
    pub anchor: VertexId,
    pub measure: f64,
    pub load_imbalance: f64,
    pub compatibility_tolerance: f64,
}
pub(crate) struct Solution {
    pub values: Vec<f64>,
    pub cg: CgReport,
    pub original_residual_norm: f64,
    pub components: Vec<Report>,
}

pub(crate) fn validate(options: NeumannOptions) -> Result<()> {
    let abs = options.compatibility_absolute_tolerance;
    let rel = options.compatibility_relative_tolerance;
    if !abs.is_finite() || abs < 0.0 || !rel.is_finite() || rel < 0.0 {
        return Err(error("invalid Neumann compatibility tolerance"));
    }
    Ok(())
}
fn root(parents: &mut [usize], mut vertex: usize) -> usize {
    while parents[vertex] != vertex {
        parents[vertex] = parents[parents[vertex]];
        vertex = parents[vertex];
    }
    vertex
}
fn components<const N: usize>(n: usize, cells: &[[VertexId; N]]) -> Vec<Vec<usize>> {
    let mut parents: Vec<_> = (0..n).collect();
    let mut sizes = vec![1usize; n];
    for cell in cells {
        for vertex in cell.iter().skip(1) {
            let a = root(&mut parents, cell[0].0);
            let b = root(&mut parents, vertex.0);
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
    for i in 0..n {
        groups.entry(root(&mut parents, i)).or_default().push(i);
    }
    let mut result: Vec<_> = groups.into_values().collect();
    result.sort_by_key(|vertices| vertices[0]);
    result
}

/// Requires a conforming nodal Laplacian, one weight per vertex and valid cell
/// indices. This is not a generic singular-matrix solver.
pub(crate) fn solve<const N: usize>(
    matrix: CsrMatrix,
    rhs: Vec<f64>,
    weights: Vec<f64>,
    cells: &[[VertexId; N]],
    options: NeumannOptions,
) -> Result<Solution> {
    validate(options)?;
    let n = rhs.len();
    if N == 0
        || matrix.rows() != n
        || matrix.cols() != n
        || weights.len() != n
        || rhs.iter().any(|x| !x.is_finite())
        || weights.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || cells.iter().flatten().any(|v| v.0 >= n)
    {
        return Err(error("invalid Neumann system or integration weights"));
    }
    let components = components(n, cells);
    let mut reports = Vec::with_capacity(components.len());
    let mut anchors = Vec::with_capacity(components.len());
    for vertices in &components {
        let mut imbalance = 0.0;
        let mut scale = 0.0;
        let mut measure = 0.0;
        for &i in vertices {
            imbalance += rhs[i];
            scale += rhs[i].abs();
            measure += weights[i];
        }
        let tolerance = options.compatibility_absolute_tolerance
            + options.compatibility_relative_tolerance * scale;
        if !imbalance.is_finite()
            || !scale.is_finite()
            || !measure.is_finite()
            || measure <= 0.0
            || !tolerance.is_finite()
        {
            return Err(error("nonfinite Neumann compatibility or measure"));
        }
        let anchor = VertexId(vertices[0]);
        if imbalance.abs() > tolerance {
            return Err(Error(format!("incompatible Neumann load on component anchored at vertex {}: imbalance {} exceeds tolerance {}",anchor.0,imbalance,tolerance)));
        }
        reports.push(Report {
            anchor,
            measure,
            load_imbalance: imbalance,
            compatibility_tolerance: tolerance,
        });
        anchors.push((anchor.0, 0.0));
    }
    let (pinned, pinned_rhs) = matrix.with_dirichlet(&rhs, &anchors)?;
    let (mut values, cg) = conjugate_gradient(&pinned, &pinned_rhs, options.cg)?;
    for (vertices, report) in components.iter().zip(&reports) {
        let integral: f64 = vertices.iter().map(|&i| weights[i] * values[i]).sum();
        let mean = integral / report.measure;
        if !mean.is_finite() {
            return Err(error("nonfinite Neumann solution mean"));
        }
        for &i in vertices {
            values[i] -= mean;
        }
    }
    let ax = matrix.mul(&values)?;
    let residual = rhs
        .iter()
        .zip(ax.iter())
        .fold(0.0_f64, |norm, (b, a)| norm.hypot(b - a));
    if !residual.is_finite() {
        return Err(error("nonfinite original Neumann residual"));
    }
    Ok(Solution {
        values,
        cg,
        original_residual_norm: residual,
        components: reports,
    })
}
