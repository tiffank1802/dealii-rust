//! Pure Neumann Poisson on the unit square, normalized to integral mean zero.
//!
//! u=x²+y²-2/3, -Δu=-4, ∂u/∂n=2x n_x+2y n_y.
use dealii_rust::{
    grid::Mesh2,
    neumann::{self, NeumannOptions},
    numerics, vtk,
};
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = Mesh2::unit_square(n)?;
    let result = neumann::solve(
        &mesh,
        |_| -4.0,
        |p, normal| 2.0 * p.x * normal[0] + 2.0 * p.y * normal[1],
        NeumannOptions::default(),
    )?;
    let (l2, h1) = numerics::errors(
        &mesh,
        &result.values,
        |p| p.x * p.x + p.y * p.y - 2.0 / 3.0,
        |p| [2.0 * p.x, 2.0 * p.y],
    )?;
    vtk::write_scalar(File::create("poisson_neumann.vtk")?, &mesh, &result.values)?;
    println!(
        "n={} vertices={} cells={} components={} iterations={} pinned_residual={:.3e} original_residual={:.3e} L2={:.3e} H1_seminorm={:.3e}",
        n, mesh.vertices().len(), mesh.cells().len(), result.components.len(), result.cg.iterations,
        result.cg.residual_norm, result.original_residual_norm, l2, h1
    );
    for (index, component) in result.components.iter().enumerate() {
        println!(
            "component={} anchor={} area={:.6} load_imbalance={:.3e} tolerance={:.3e}",
            index,
            component.anchor.0,
            component.area,
            component.load_imbalance,
            component.compatibility_tolerance
        );
    }
    println!("Wrote poisson_neumann.vtk");
    Ok(())
}
