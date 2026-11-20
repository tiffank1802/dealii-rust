//! Uniform refinement and residual indicators for a manufactured Poisson problem.
use dealii_rust::{
    base::Point2,
    estimator,
    grid::{Mesh2, UniformRefinement},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};
use std::f64::consts::PI;

fn exact(p: Point2) -> f64 {
    (PI * p.x).sin() * (PI * p.y).sin()
}
fn source(p: Point2) -> f64 {
    2.0 * PI * PI * exact(p)
}
fn gradient(p: Point2) -> [f64; 2] {
    [
        PI * (PI * p.x).cos() * (PI * p.y).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).cos(),
    ]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = Mesh2::unit_square(4)?;
    for level in 0..4 {
        let (matrix, rhs) = numerics::assemble_poisson(&mesh, source)?;
        let constraints: Vec<_> = mesh
            .boundary_vertices()
            .iter()
            .map(|v| (v.0, 0.0))
            .collect();
        let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
        let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
        let (_, h1) = numerics::errors(&mesh, &solution, exact, gradient)?;
        let estimate = estimator::poisson_dirichlet(&mesh, &solution, source)?;
        println!("level={} vertices={} cells={} iterations={} residual={:.3e} H1={:.3e} indicator={:.3e}",
            level, mesh.vertices().len(), mesh.cells().len(), report.iterations, report.residual_norm, h1, estimate.global_indicator);
        if level < 3 {
            let refined = UniformRefinement::new(&mesh)?;
            // This transfer represents the same piecewise P1 function on the
            // fine mesh. The current CG always starts from zero, so the example
            // only demonstrates the transfer and solves anew at the next level.
            let transferred = refined.prolongate_p1(&solution)?;
            assert_eq!(transferred.len(), refined.mesh().vertices().len());
            mesh = refined.mesh().clone();
        }
    }
    Ok(())
}
