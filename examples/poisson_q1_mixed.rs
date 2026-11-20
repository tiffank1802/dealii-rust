//! Q1 Poisson: Dirichlet at x=0, outward Neumann on the other sides.
//! Exact solution u=1+x+x²+2y+y², so -Δu=-4.
use dealii_rust::{
    base::Point2,
    grid::BoundaryEdgeId,
    lac::{conjugate_gradient, CgOptions},
    quad::{self, QuadMesh},
};
use std::fs::File;

fn exact(p: Point2) -> f64 {
    1.0 + p.x + p.x * p.x + 2.0 * p.y + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [1.0 + 2.0 * p.x, 2.0 + 2.0 * p.y]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = QuadMesh::unit_square(n)?;
    let mut dirichlet = Vec::new();
    let mut neumann = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        let edge = mesh.boundary_edge_geometry(id)?;
        // Selection by coordinates is specific to this structured square.
        if edge.points.iter().all(|p| p.x == 0.0) {
            dirichlet.push(id);
        } else {
            neumann.push(id);
        }
    }
    let (matrix, mut rhs) = quad::assemble_poisson(&mesh, |_| -4.0)?;
    let load = quad::assemble_neumann(&mesh, &neumann, |p, n| {
        let g = gradient(p);
        g[0] * n[0] + g[1] * n[1]
    })?;
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
    }
    let constraints = quad::interpolate_boundary(&mesh, &dirichlet, exact)?;
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (solution, report) = conjugate_gradient(&matrix, &rhs, options)?;
    let (l2, h1) = quad::errors(&mesh, &solution, exact, gradient)?;
    quad::write_vtk(File::create("poisson_q1_mixed.vtk")?, &mesh, &solution)?;
    println!("n={} vertices={} cells={} dirichlet_edges={} neumann_edges={} iterations={} residual={:.3e} L2={:.3e} H1_seminorm={:.3e}",n,mesh.vertices().len(),mesh.cells().len(),dirichlet.len(),neumann.len(),report.iterations,report.residual_norm,l2,h1);
    println!("Wrote poisson_q1_mixed.vtk");
    Ok(())
}
