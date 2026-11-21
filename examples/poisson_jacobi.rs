//! Compare plain and Jacobi-preconditioned CG on a two-material P1 problem.
use dealii_rust::{
    elliptic,
    grid::Mesh2,
    integrators::DiffusionTensor,
    lac::{conjugate_gradient, conjugate_gradient_jacobi, CgOptions, JacobiPreconditioner},
    vtk,
};
use std::fs::File;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = Mesh2::unit_square(n)?;
    // Isotropic diffusion differs by 1000× between the two halves.
    let (a, b) = elliptic::assemble(
        &mesh,
        |_, p| {
            let k = if p.x < 0.5 { 1000.0 } else { 1.0 };
            DiffusionTensor {
                xx: k,
                xy: 0.0,
                yy: k,
            }
        },
        |_, _| 0.0,
        |_, _| 1.0,
    )?;
    let constraints: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, 0.0))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &constraints)?;
    // Construct the preconditioner from the already constrained matrix.
    let jacobi = JacobiPreconditioner::new(&a)?;
    let options = CgOptions {
        relative_tolerance: 1e-11,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (plain, plain_report) = conjugate_gradient(&a, &b, options)?;
    let (preconditioned, report) = conjugate_gradient_jacobi(&a, &b, &jacobi, options)?;
    let difference = plain
        .iter()
        .zip(&preconditioned)
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    println!("n={} vertices={} cells={} CG_iterations={} CG_residual={:.3e} Jacobi_iterations={} Jacobi_residual={:.3e} max_solution_difference={:.3e}", n,mesh.vertices().len(),mesh.cells().len(),plain_report.iterations,plain_report.residual_norm,report.iterations,report.residual_norm,difference);
    vtk::write_scalar(File::create("poisson_jacobi.vtk")?, &mesh, &preconditioned)?;
    println!("Wrote poisson_jacobi.vtk");
    Ok(())
}
