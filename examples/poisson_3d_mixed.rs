//! 3D Poisson P1: Dirichlet on x=0 and Neumann on the other faces.
//! Exact solution u=1+x+x²+2y+y²+3z+z², hence -Δu=-6.
use dealii_rust::{
    base::Point3,
    lac::{conjugate_gradient, CgOptions},
    tetra::{self, BoundaryFaceId, TetMesh},
};
use std::fs::File;

fn exact(p: Point3) -> f64 {
    1.0 + p.x + p.x * p.x + 2.0 * p.y + p.y * p.y + 3.0 * p.z + p.z * p.z
}
fn gradient(p: Point3) -> [f64; 3] {
    [1.0 + 2.0 * p.x, 2.0 + 2.0 * p.y, 3.0 + 2.0 * p.z]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = TetMesh::unit_cube(n)?;
    let mut dirichlet = Vec::new();
    let mut neumann = Vec::new();
    for i in 0..mesh.boundary_faces().len() {
        let id = BoundaryFaceId(i);
        let face = mesh.boundary_face_geometry(id)?;
        // This coordinate test is specific to the structured unit cube.
        if face.points.iter().all(|p| p.x == 0.0) {
            dirichlet.push(id);
        } else {
            neumann.push(id);
        }
    }
    let (matrix, mut rhs) = tetra::assemble_poisson(&mesh, |_| -6.0)?;
    let load = tetra::assemble_neumann(&mesh, &neumann, |p, n| {
        let g = gradient(p);
        g[0] * n[0] + g[1] * n[1] + g[2] * n[2]
    })?;
    for (b, value) in rhs.iter_mut().zip(load) {
        *b += value;
    }
    let constraints = tetra::interpolate_boundary(&mesh, &dirichlet, exact)?;
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (solution, report) = conjugate_gradient(&matrix, &rhs, options)?;
    let (l2, h1) = tetra::errors(&mesh, &solution, exact, gradient)?;
    tetra::write_vtk(File::create("poisson_3d_mixed.vtk")?, &mesh, &solution)?;
    println!("n={} vertices={} cells={} dirichlet_faces={} neumann_faces={} iterations={} residual={:.3e} L2={:.3e} H1_seminorm={:.3e}",n,mesh.vertices().len(),mesh.cells().len(),dirichlet.len(),neumann.len(),report.iterations,report.residual_norm,l2,h1);
    println!("Wrote poisson_3d_mixed.vtk");
    Ok(())
}
