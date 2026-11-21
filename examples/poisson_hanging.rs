//! One-level nonconforming P1 Poisson with automatically generated midpoint constraints.
use dealii_rust::{
    base::Point2,
    fe::{self, AffineTriangle},
    grid::{CellId, Mesh2, NonconformingRefinement},
    lac::{conjugate_gradient, CgOptions},
    quadrature::triangle_degree_5,
};
fn exact(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y
}
fn errors(
    mesh: &NonconformingRefinement,
    u: &[f64],
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let (mut l2, mut h1) = (0., 0.);
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = AffineTriangle::new(mesh.cell_points(CellId(i))?)?;
        let grads = map.gradients();
        let mut grad = [0.; 2];
        for j in 0..3 {
            for (d, value) in grad.iter_mut().enumerate() {
                *value += u[cell[j].0] * grads[j][d];
            }
        }
        for q in triangle_degree_5() {
            let p = map.map(q.point);
            let shape = fe::shape_values(q.point);
            let value = (0..3).map(|j| shape[j] * u[cell[j].0]).sum::<f64>();
            let w = map.determinant() * q.weight;
            l2 += (value - exact(p)).powi(2) * w;
            h1 += ((grad[0] - 2. * p.x).powi(2) + (grad[1] - 2. * p.y).powi(2)) * w;
        }
    }
    Ok((f64::sqrt(l2), f64::sqrt(h1)))
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = std::env::args().nth(1).map_or(Ok(4), |s| s.parse())?;
    let parent = Mesh2::unit_square(n)?;
    // Deterministic nonuniform one-level pattern, not an adaptive estimator.
    let marked: Vec<_> = (0..parent.cells().len()).step_by(2).map(CellId).collect();
    let mesh = NonconformingRefinement::new(&parent, &marked)?;
    let (a, b) = mesh.assemble_poisson(|_| -4.)?;
    let mut constraints = mesh.hanging_constraints()?;
    for (id, value) in mesh.interpolate_boundary(exact)? {
        constraints.add(id, value, &[])?;
    }
    let constraints = constraints.close()?;
    let (a, b) = constraints.condense(&a, &b)?;
    let (values, iterations, residual) = if constraints.free_dofs().is_empty() {
        (constraints.reconstruct(&[])?, 0, 0.)
    } else {
        let (reduced, report) = conjugate_gradient(&a, &b, CgOptions::default())?;
        (
            constraints.reconstruct(&reduced)?,
            report.iterations,
            report.residual_norm,
        )
    };
    let (l2, h1) = errors(&mesh, &values)?;
    println!("n={} cells={} vertices={} hanging={} boundary_edges={} free={} iterations={} residual={:.3e} L2={:.6e} H1={:.6e}",n,mesh.cells().len(),mesh.vertices().len(),mesh.hanging_nodes().len(),mesh.boundary_edges().len(),constraints.free_dofs().len(),iterations,residual,l2,h1);
    Ok(())
}
