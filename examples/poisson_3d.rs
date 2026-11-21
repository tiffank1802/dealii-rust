//! 3D Poisson P1 on tetrahedra, sine solution and homogeneous Dirichlet.
use dealii_rust::{
    base::Point3,
    lac::{conjugate_gradient, CgOptions},
    tetra::{self, TetMesh},
};
use std::{f64::consts::PI, fs::File};

fn exact(p: Point3) -> f64 {
    (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).sin()
}
fn gradient(p: Point3) -> [f64; 3] {
    [
        PI * (PI * p.x).cos() * (PI * p.y).sin() * (PI * p.z).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).cos() * (PI * p.z).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).cos(),
    ]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = TetMesh::unit_cube(n)?;
    let (matrix, rhs) = tetra::assemble_poisson(&mesh, |p| 3.0 * PI * PI * exact(p))?;
    let constraints: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, 0.0))
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
    let (l2, h1) = tetra::errors(&mesh, &solution, exact, gradient)?;
    tetra::write_vtk(File::create("poisson_3d.vtk")?, &mesh, &solution)?;
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
    println!("Wrote poisson_3d.vtk");
    Ok(())
}
