//! Pure Neumann Poisson on the unit cube, with zero integral mean.
//! u = x²+y²+z²-1; -Δu=-6; g=2(x n_x+y n_y+z n_z).
use dealii_rust::{
    neumann::NeumannOptions,
    neumann3d,
    tetra::{self, TetMesh},
};
use std::fs::File;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = TetMesh::unit_cube(n)?;
    let result = neumann3d::solve(
        &mesh,
        |_| -6.0,
        |p, n| 2.0 * (p.x * n[0] + p.y * n[1] + p.z * n[2]),
        NeumannOptions::default(),
    )?;
    let (l2, h1) = tetra::errors(
        &mesh,
        &result.values,
        |p| p.x * p.x + p.y * p.y + p.z * p.z - 1.0,
        |p| [2.0 * p.x, 2.0 * p.y, 2.0 * p.z],
    )?;
    tetra::write_vtk(
        File::create("poisson_3d_neumann.vtk")?,
        &mesh,
        &result.values,
    )?;
    println!("n={} vertices={} cells={} components={} iterations={} pinned_residual={:.3e} original_residual={:.3e} L2={:.3e} H1_seminorm={:.3e}", n, mesh.vertices().len(), mesh.cells().len(), result.components.len(), result.cg.iterations, result.cg.residual_norm, result.original_residual_norm, l2, h1);
    for (i, c) in result.components.iter().enumerate() {
        println!(
            "component={} anchor={} volume={:.6} imbalance={:.3e} tolerance={:.3e}",
            i, c.anchor.0, c.volume, c.load_imbalance, c.compatibility_tolerance
        );
    }
    println!("Wrote poisson_3d_neumann.vtk");
    Ok(())
}
