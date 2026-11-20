//! Poisson with nonhomogeneous Dirichlet and Neumann conditions.
//!
//! u=1+x+x²+2y+y² on the unit square, so -Δu=-4. Dirichlet is
//! prescribed on x=0 and ∂u/∂n on the other three sides.
use dealii_rust::{
    base::Point2,
    dofs::DofHandler,
    grid::{BoundaryEdgeId, Mesh2},
    lac::{conjugate_gradient, CgOptions},
    numerics, vtk,
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
    let mesh = Mesh2::unit_square(n)?;
    let (matrix, mut rhs) = numerics::assemble_poisson(&mesh, |_| -4.0)?;
    let mut neumann_edges = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        let edge = mesh.boundary_edge_geometry(BoundaryEdgeId(i))?;
        // For this generated grid x=0 is exactly representable. Imported
        // meshes will need boundary markers instead of coordinate comparisons.
        if !(edge.points[0].x == 0.0 && edge.points[1].x == 0.0) {
            neumann_edges.push(BoundaryEdgeId(i));
        }
    }
    let load = numerics::assemble_neumann(&mesh, &neumann_edges, |p, normal| {
        let g = gradient(p);
        g[0] * normal[0] + g[1] * normal[1]
    })?;
    for (r, value) in rhs.iter_mut().zip(load) {
        *r += value;
    }
    let constraints: Vec<_> = DofHandler::new(&mesh)
        .boundary_dofs()
        .into_iter()
        .filter(|d| mesh.vertices()[d.0].x == 0.0)
        .map(|d| (d.0, exact(mesh.vertices()[d.0])))
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
    let (l2, h1) = numerics::errors(&mesh, &solution, exact, gradient)?;
    vtk::write_scalar(File::create("poisson_mixed.vtk")?, &mesh, &solution)?;
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
    println!("Wrote poisson_mixed.vtk");
    Ok(())
}
