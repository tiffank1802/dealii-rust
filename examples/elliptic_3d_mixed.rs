//! Mixed 3D diffusion-reaction with variable anisotropic tensor.
use dealii_rust::{
    base::Point3,
    elliptic_tetra::{self, DiffusionTensor3},
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
    let (matrix, mut rhs) = elliptic_tetra::assemble(
        &mesh,
        |_, p| tensor(p),
        |_, p| 1.0 + p.x + p.y + p.z,
        |_, p| -(24.0 + 4.0 * (p.x + p.y + p.z)) + (1.0 + p.x + p.y + p.z) * exact(p),
    )?;
    let mut dirichlet = Vec::new();
    let mut neumann = Vec::new();
    for i in 0..mesh.boundary_faces().len() {
        let id = BoundaryFaceId(i);
        if mesh
            .boundary_face_geometry(id)?
            .points
            .iter()
            .all(|p| p.x == 0.0)
        {
            dirichlet.push(id);
        } else {
            neumann.push(id);
        }
    }
    let load = tetra::assemble_neumann(&mesh, &neumann, |p, n| {
        let a = tensor(p);
        let [x, y, z] = gradient(p);
        n[0] * (a.xx * x + a.xy * y + a.xz * z)
            + n[1] * (a.xy * x + a.yy * y + a.yz * z)
            + n[2] * (a.xz * x + a.yz * y + a.zz * z)
    })?;
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
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
    tetra::write_vtk(File::create("elliptic_3d_mixed.vtk")?, &mesh, &solution)?;
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
    println!("Wrote elliptic_3d_mixed.vtk");
    Ok(())
}
