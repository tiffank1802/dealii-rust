//! Poisson on Q1 quadrilaterals, homogeneous Dirichlet, sine solution.
use dealii_rust::{
    base::Point2,
    lac::{conjugate_gradient, CgOptions},
    quad::{self, QuadMesh},
};
use std::{f64::consts::PI, fs::File};

fn exact(p: Point2) -> f64 {
    (PI * p.x).sin() * (PI * p.y).sin()
}
fn gradient(p: Point2) -> [f64; 2] {
    [
        PI * (PI * p.x).cos() * (PI * p.y).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).cos(),
    ]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = QuadMesh::unit_square(n)?;
    let (matrix, rhs) = quad::assemble_poisson(&mesh, |p| 2.0 * PI * PI * exact(p))?;
    let constraints: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, 0.0))
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
    let (l2, h1) = quad::errors(&mesh, &solution, exact, gradient)?;
    quad::write_vtk(File::create("poisson_q1.vtk")?, &mesh, &solution)?;
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
    println!("Wrote poisson_q1.vtk");
    Ok(())
}
