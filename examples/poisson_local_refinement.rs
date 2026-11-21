//! Indicator-driven, conforming local red-green refinement for P1 Poisson.
use dealii_rust::{
    base::Point2,
    estimator,
    grid::{Mesh2, RefinementHistory},
    lac::{conjugate_gradient, conjugate_gradient_with_initial_guess, CgOptions},
    numerics, vtk,
};
use std::f64::consts::PI;
use std::fs::File;
fn exact(p: Point2) -> f64 {
    (PI * p.x).sin() * (PI * p.y).sin()
}
fn source(p: Point2) -> f64 {
    2. * PI * PI * exact(p)
}
fn gradient(p: Point2) -> [f64; 2] {
    [
        PI * (PI * p.x).cos() * (PI * p.y).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).cos(),
    ]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut history = RefinementHistory::new(Mesh2::unit_square(4)?);
    let mut initial: Option<Vec<f64>> = None;
    for level in 0..5 {
        let mesh = history.mesh();
        let (a, b) = numerics::assemble_poisson(mesh, source)?;
        let boundary: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
        let (a, b) = a.with_dirichlet(&b, &boundary)?;
        let (u, report) = match initial.take() {
            Some(mut x) => {
                for &(i, g) in &boundary {
                    x[i] = g;
                }
                conjugate_gradient_with_initial_guess(&a, &b, &x, CgOptions::default())?
            }
            None => conjugate_gradient(&a, &b, CgOptions::default())?,
        };
        let (l2, h1) = numerics::errors(mesh, &u, exact, gradient)?;
        let estimate = estimator::poisson_dirichlet(mesh, &u, source)?;
        println!("level={} vertices={} cells={} iterations={} residual={:.3e} L2={:.3e} H1={:.3e} indicator={:.3e}",level,mesh.vertices().len(),mesh.cells().len(),report.iterations,report.residual_norm,l2,h1,estimate.global_indicator);
        if level == 4 {
            vtk::write_scalar_with_cell_scalar(
                File::create("poisson_local_refinement.vtk")?,
                mesh,
                &u,
                "indicator",
                &estimate.cell_indicators,
            )?;
            println!("Wrote poisson_local_refinement.vtk");
            // Demonstrate whole-level rollback with nodal injection. This is
            // not a new coarse solve and does not preserve child-only detail.
            let coarse = history.restrict_injection(level - 1, &u)?;
            history.coarsen_last()?;
            println!(
                "rolled_back_to_level={} restricted_values={}",
                history.level(),
                coarse.len()
            );
        }
        if level < 4 {
            // Mark cells accounting for half of the squared indicators.
            let marked = estimator::mark_bulk(&estimate.cell_indicators, 0.5)?;
            history.refine(&marked)?;
            initial = Some(history.transition(level)?.prolongate_p1(&u)?);
        }
    }
    Ok(())
}
