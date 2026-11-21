//! Variable anisotropic diffusion-reaction with mixed boundary conditions.
//! u = 1+x+x²+2y+y²; A = [[2+x, 1/2],[1/2,3+y]], c=1+x+y.
use dealii_rust::{
    base::Point2,
    elliptic, field,
    grid::{BoundaryEdgeId, Mesh2},
    integrators::DiffusionTensor,
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
fn diffusion(p: Point2) -> DiffusionTensor {
    DiffusionTensor {
        xx: 2.0 + p.x,
        xy: 0.5,
        yy: 3.0 + p.y,
    }
}
fn source(p: Point2) -> f64 {
    -(13.0 + 4.0 * p.x + 4.0 * p.y) + (1.0 + p.x + p.y) * exact(p)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 16,
    };
    let mesh = Mesh2::unit_square(n)?;
    let (a, mut rhs) = elliptic::assemble(
        &mesh,
        |_, p| diffusion(p),
        |_, p| 1.0 + p.x + p.y,
        |_, p| source(p),
    )?;
    let mut dirichlet = Vec::new();
    let mut neumann = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        let edge = mesh.boundary_edge_geometry(id)?;
        if edge.points.iter().all(|p| p.x == 0.0) {
            dirichlet.push(id);
        } else {
            neumann.push(id);
        }
    }
    let load = numerics::assemble_neumann(&mesh, &neumann, |p, n| {
        let a = diffusion(p);
        let g = gradient(p);
        n[0] * (a.xx * g[0] + a.xy * g[1]) + n[1] * (a.xy * g[0] + a.yy * g[1])
    })?;
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
    }
    let constraints = field::interpolate_boundary(&mesh, &dirichlet, exact)?;
    let (a, rhs) = a.with_dirichlet(&rhs, &constraints)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (solution, report) = conjugate_gradient(&a, &rhs, options)?;
    let (l2, h1) = numerics::errors(&mesh, &solution, exact, gradient)?;
    vtk::write_scalar(File::create("elliptic_mixed.vtk")?, &mesh, &solution)?;
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
    println!("Wrote elliptic_mixed.vtk");
    Ok(())
}
