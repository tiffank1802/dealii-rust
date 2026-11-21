//! Pure Neumann variable anisotropic diffusion on the unit cube.
//! u=x²+y²+z²-1, A has diagonal (2+x,3+y,4+z) and constant off-diagonals.
use dealii_rust::{
    base::Point3,
    elliptic_tetra::DiffusionTensor3,
    neumann3d::{self, NeumannOptions},
    tetra::{self, TetMesh},
};
use std::fs::File;
fn tensor(p: Point3) -> DiffusionTensor3 {
    DiffusionTensor3 {
        xx: 2.0 + p.x,
        xy: 0.5,
        xz: 0.25,
        yy: 3.0 + p.y,
        yz: 0.75,
        zz: 4.0 + p.z,
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = TetMesh::unit_cube(n)?;
    let r = neumann3d::solve_diffusion(
        &mesh,
        |_, p| tensor(p),
        |_, p| -18.0 - 4.0 * (p.x + p.y + p.z),
        |p, n| {
            let a = tensor(p);
            let [x, y, z] = [2.0 * p.x, 2.0 * p.y, 2.0 * p.z];
            n[0] * (a.xx * x + a.xy * y + a.xz * z)
                + n[1] * (a.xy * x + a.yy * y + a.yz * z)
                + n[2] * (a.xz * x + a.yz * y + a.zz * z)
        },
        NeumannOptions::default(),
    )?;
    let (l2, h1) = tetra::errors(
        &mesh,
        &r.values,
        |p| p.x * p.x + p.y * p.y + p.z * p.z - 1.0,
        |p| [2.0 * p.x, 2.0 * p.y, 2.0 * p.z],
    )?;
    tetra::write_vtk(File::create("elliptic_3d_neumann.vtk")?, &mesh, &r.values)?;
    println!("n={} vertices={} cells={} components={} iterations={} pinned_residual={:.3e} original_residual={:.3e} L2={:.3e} H1_seminorm={:.3e}", n,mesh.vertices().len(),mesh.cells().len(),r.components.len(),r.cg.iterations,r.cg.residual_norm,r.original_residual_norm,l2,h1);
    for (i, c) in r.components.iter().enumerate() {
        println!(
            "component={} anchor={} volume={:.6} imbalance={:.3e} tolerance={:.3e}",
            i, c.anchor.0, c.volume, c.load_imbalance, c.compatibility_tolerance
        );
    }
    println!("Wrote elliptic_3d_neumann.vtk");
    Ok(())
}
