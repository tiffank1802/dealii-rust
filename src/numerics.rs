//! Assembly and error evaluation for scalar Poisson on triangles.
use crate::base::Point2;
use crate::dofs::DofHandler;
use crate::fe::{shape_values, AffineTriangle};
use crate::grid::{BoundaryEdgeId, CellId, Mesh2};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::{line_degree_5, triangle_degree_5};
use crate::{error, Result};
use std::collections::BTreeSet;

/// Assemble the unconstrained stiffness matrix and load vector for -Δu=f.
pub fn assemble_poisson<F: Fn(Point2) -> f64>(
    mesh: &Mesh2,
    source: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let dofs = DofHandler::new(mesh);
    let mut matrix = CooBuilder::new(dofs.count(), dofs.count());
    let mut rhs = vec![0.0; dofs.count()];
    let quadrature = triangle_degree_5();
    for cell in 0..mesh.cells().len() {
        let mapping = AffineTriangle::new(mesh.cell_points(CellId(cell))?)?;
        let ids = dofs.cell_dofs(CellId(cell))?;
        let gradients = mapping.gradients();
        let area = mapping.determinant() * 0.5;
        for i in 0..3 {
            for j in 0..3 {
                matrix.add(
                    ids[i].0,
                    ids[j].0,
                    area * (gradients[i][0] * gradients[j][0] + gradients[i][1] * gradients[j][1]),
                )?;
            }
        }
        for q in &quadrature {
            let value = source(mapping.map(q.point));
            if !value.is_finite() {
                return Err(error("source function returned a nonfinite value"));
            }
            let shapes = shape_values(q.point);
            for i in 0..3 {
                rhs[ids[i].0] += value * shapes[i] * mapping.determinant() * q.weight;
            }
        }
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        return Err(error("nonfinite assembled right hand side"));
    }
    Ok((matrix.into_csr()?, rhs))
}

/// Assemble the Neumann load on selected boundary edges.
///
/// `flux(point, normal)` is the prescribed outward flux: `∂u/∂n` for
/// `-Δu=f`, or `n · A ∇u` for `-div(A∇u) + cu = f`;
/// the weak form adds `∫g v ds` to the right-hand side. It returns a vector
/// to add to the volume load *before* Dirichlet elimination. Edge IDs must
/// refer to this mesh, be in bounds and be unique. The P1 shape functions on
/// an edge are integrated using a degree-five Gauss rule.
///
/// For pure diffusion (`c = 0`), a pure Neumann problem has a constant
/// nullspace on each connected component and requires load compatibility and
/// a normalization. This routine does not impose either; the existing CG
/// solver requires an SPD system. A positive reaction can remove this nullspace.
pub fn assemble_neumann<F>(mesh: &Mesh2, edges: &[BoundaryEdgeId], flux: F) -> Result<Vec<f64>>
where
    F: Fn(Point2, [f64; 2]) -> f64,
{
    let mut unique = BTreeSet::new();
    for &id in edges {
        if !unique.insert(id) {
            return Err(error("duplicate Neumann boundary edge"));
        }
        mesh.boundary_edge_geometry(id)?;
    }
    let mut load = vec![0.0; mesh.vertices().len()];
    let quadrature = line_degree_5();
    for &id in edges {
        let edge = mesh.boundary_edge_geometry(id)?;
        for q in &quadrature {
            let t = q.point;
            let point = Point2 {
                x: (1.0 - t) * edge.points[0].x + t * edge.points[1].x,
                y: (1.0 - t) * edge.points[0].y + t * edge.points[1].y,
            };
            let value = flux(point, edge.outward_normal);
            if !value.is_finite() {
                return Err(error("Neumann flux returned a nonfinite value"));
            }
            let weighted = value * edge.length * q.weight;
            load[edge.vertices[0].0] += weighted * (1.0 - t);
            load[edge.vertices[1].0] += weighted * t;
        }
    }
    if load.iter().any(|value| !value.is_finite()) {
        return Err(error("nonfinite Neumann load"));
    }
    Ok(load)
}

/// Integrated L2 error and H1 seminorm error, approximated by degree-5 quadrature.
pub fn errors<F, G>(
    mesh: &Mesh2,
    solution: &[f64],
    exact: F,
    exact_gradient: G,
) -> Result<(f64, f64)>
where
    F: Fn(Point2) -> f64,
    G: Fn(Point2) -> [f64; 2],
{
    if solution.len() != mesh.vertices().len() || solution.iter().any(|x| !x.is_finite()) {
        return Err(error("invalid solution vector"));
    }
    let dofs = DofHandler::new(mesh);
    let mut l2 = 0.0;
    let mut h1 = 0.0;
    let quadrature = triangle_degree_5();
    for cell in 0..mesh.cells().len() {
        let mapping = AffineTriangle::new(mesh.cell_points(CellId(cell))?)?;
        let ids = dofs.cell_dofs(CellId(cell))?;
        let gradients = mapping.gradients();
        let mut grad_h = [0.0; 2];
        for i in 0..3 {
            for (component, gradient) in grad_h.iter_mut().enumerate() {
                *gradient += solution[ids[i].0] * gradients[i][component];
            }
        }
        for q in &quadrature {
            let physical = mapping.map(q.point);
            let exact_value = exact(physical);
            let exact_grad = exact_gradient(physical);
            if !exact_value.is_finite() || exact_grad.iter().any(|x| !x.is_finite()) {
                return Err(error("exact solution returned a nonfinite value"));
            }
            let shapes = shape_values(q.point);
            let uh: f64 = (0..3).map(|i| shapes[i] * solution[ids[i].0]).sum();
            let diff = exact_value - uh;
            let dx = exact_grad[0] - grad_h[0];
            let dy = exact_grad[1] - grad_h[1];
            let jxw = mapping.determinant() * q.weight;
            l2 += diff * diff * jxw;
            h1 += (dx * dx + dy * dy) * jxw;
        }
    }
    if !l2.is_finite() || !h1.is_finite() {
        return Err(error("nonfinite integrated error"));
    }
    Ok((l2.sqrt(), h1.sqrt()))
}
