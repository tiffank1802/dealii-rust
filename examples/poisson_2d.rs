//! Solve -Δu=2π² sin(πx) sin(πy), u=0 on the unit square.
use dealii_rust::{
    dofs::DofHandler,
    grid::Mesh2,
    lac::{conjugate_gradient, CgOptions},
    numerics, vtk,
};
use std::f64::consts::PI;
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = Mesh2::unit_square(n)?;
    let (matrix, rhs) = numerics::assemble_poisson(&mesh, |p| {
        2.0 * PI * PI * (PI * p.x).sin() * (PI * p.y).sin()
    })?;
    let constraints: Vec<_> = DofHandler::new(&mesh)
        .boundary_dofs()
        .iter()
        .map(|d| (d.0, 0.0))
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
    let (l2, h1) = numerics::errors(
        &mesh,
        &solution,
        |p| (PI * p.x).sin() * (PI * p.y).sin(),
        |p| {
            [
                PI * (PI * p.x).cos() * (PI * p.y).sin(),
                PI * (PI * p.x).sin() * (PI * p.y).cos(),
            ]
        },
    )?;
    vtk::write_scalar(File::create("poisson_2d.vtk")?, &mesh, &solution)?;
    println!(
        "n={} vertices={} cells={} iterations={} residual={:.3e} L2={:.3e} H1_seminorm={:.3e}",
        n,
        mesh.vertices().len(),
        mesh.cells().len(),
        report.iterations,
        report.residual_norm,
        l2,
        h1
    );
    println!("Wrote poisson_2d.vtk");
    Ok(())
}
