//! Uniform refinement and residual indicators for a manufactured Poisson problem.
use dealii_rust::{
    base::Point2,
    estimator,
    grid::{Mesh2, UniformRefinement},
    lac::{conjugate_gradient, conjugate_gradient_with_initial_guess, CgOptions},
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
    let mut initial: Option<Vec<f64>> = None;
    for level in 0..4 {
        let (matrix, rhs) = numerics::assemble_poisson(&mesh, source)?;
        let constraints: Vec<_> = mesh
            .boundary_vertices()
            .iter()
            .map(|v| (v.0, 0.0))
            .collect();
        let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
        let (solution, report) = match initial.take() {
            Some(mut guess) => {
                for &(i, value) in &constraints {
                    guess[i] = value;
                }
                conjugate_gradient_with_initial_guess(&matrix, &rhs, &guess, CgOptions::default())?
            }
            None => conjugate_gradient(&matrix, &rhs, CgOptions::default())?,
        };
        let (_, h1) = numerics::errors(&mesh, &solution, exact, gradient)?;
        let estimate = estimator::poisson_dirichlet(&mesh, &solution, source)?;
        println!("level={} vertices={} cells={} iterations={} residual={:.3e} H1={:.3e} indicator={:.3e}",
            level, mesh.vertices().len(), mesh.cells().len(), report.iterations, report.residual_norm, h1, estimate.global_indicator);
        if level < 3 {
            let refined = UniformRefinement::new(&mesh)?;
            // Prolong the coarse field as the initial guess for the next level.
            initial = Some(refined.prolongate_p1(&solution)?);
            mesh = refined.mesh().clone();
        }
    }
    Ok(())
}
