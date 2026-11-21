//! Pure Neumann Poisson on the unit square: u=x²+y²-2/3.
use dealii_rust::{
    neumann_quad::{self, NeumannOptions},
    quad::{self, QuadMesh},
};
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = QuadMesh::unit_square(n)?;
    let result = neumann_quad::solve(
        &mesh,
        |_| -4.0,
        |p, normal| 2.0 * (p.x * normal[0] + p.y * normal[1]),
        NeumannOptions::default(),
    )?;
    let (l2, h1) = quad::errors(
        &mesh,
        &result.values,
        |p| p.x * p.x + p.y * p.y - 2.0 / 3.0,
        |p| [2.0 * p.x, 2.0 * p.y],
    )?;
    quad::write_vtk(
        File::create("poisson_q1_neumann.vtk")?,
        &mesh,
        &result.values,
    )?;
    println!("n={} vertices={} cells={} components={} iterations={} pinned_residual={:.3e} original_residual={:.3e} L2={:.3e} H1_seminorm={:.3e}", n, mesh.vertices().len(), mesh.cells().len(), result.components.len(), result.cg.iterations, result.cg.residual_norm, result.original_residual_norm, l2, h1);
    for (i, c) in result.components.iter().enumerate() {
        println!(
            "component={} anchor={} area={:.6} imbalance={:.3e} tolerance={:.3e}",
            i, c.anchor.0, c.area, c.load_imbalance, c.compatibility_tolerance
        );
    }
    println!("Wrote poisson_q1_neumann.vtk");
    Ok(())
}
