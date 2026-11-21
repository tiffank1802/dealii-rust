//! Uniform Q1 refinement with explicit boundary-label inheritance.
use dealii_rust::{
    base::Point2,
    grid::BoundaryEdgeId,
    lac::{conjugate_gradient, conjugate_gradient_with_initial_guess, CgOptions},
    quad::{self, QuadMesh, UniformRefinement},
};
use std::fs::File;
fn exact(p: Point2) -> f64 {
    1.0 + p.x + p.x * p.x + 2.0 * p.y + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [1.0 + 2.0 * p.x, 2.0 + 2.0 * p.y]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = QuadMesh::unit_square(4)?;
    // Initially classify x=0 as Dirichlet, then inherit the classification
    // from boundary parents rather than repeating the coordinate test.
    let mut dirichlet = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        dirichlet.push(
            mesh.boundary_edge_geometry(BoundaryEdgeId(i))?
                .points
                .iter()
                .all(|p| p.x == 0.0),
        );
    }
    let mut initial: Option<Vec<f64>> = None;
    for level in 0..4 {
        let (a, mut b) = quad::assemble_poisson(&mesh, |_| -4.0)?;
        let mut d = Vec::new();
        let mut n = Vec::new();
        for (i, &is_d) in dirichlet.iter().enumerate() {
            if is_d {
                d.push(BoundaryEdgeId(i));
            } else {
                n.push(BoundaryEdgeId(i));
            }
        }
        let flux = quad::assemble_neumann(&mesh, &n, |p, n| {
            let g = gradient(p);
            g[0] * n[0] + g[1] * n[1]
        })?;
        for (b, f) in b.iter_mut().zip(flux) {
            *b += f;
        }
        let boundary = quad::interpolate_boundary(&mesh, &d, exact)?;
        let (a, b) = a.with_dirichlet(&b, &boundary)?;
        let options = CgOptions {
            relative_tolerance: 1e-12,
            absolute_tolerance: 1e-13,
            ..CgOptions::default()
        };
        let (solution, report) = match initial.take() {
            Some(mut guess) => {
                // Newly created boundary nodes must satisfy the fine-grid data.
                for &(i, value) in &boundary {
                    guess[i] = value;
                }
                conjugate_gradient_with_initial_guess(&a, &b, &guess, options)?
            }
            None => conjugate_gradient(&a, &b, options)?,
        };
        let (l2, h1) = quad::errors(&mesh, &solution, exact, gradient)?;
        println!("level={} vertices={} cells={} dirichlet_edges={} neumann_edges={} iterations={} residual={:.3e} L2={:.3e} H1={:.3e}",level,mesh.vertices().len(),mesh.cells().len(),d.len(),n.len(),report.iterations,report.residual_norm,l2,h1);
        if level == 3 {
            quad::write_vtk(File::create("poisson_q1_refinement.vtk")?, &mesh, &solution)?;
            println!("Wrote poisson_q1_refinement.vtk");
        } else {
            let r = UniformRefinement::new(&mesh)?;
            initial = Some(r.prolongate_q1(&solution)?);
            let mut child_dirichlet = Vec::new();
            for i in 0..r.mesh().boundary_edges().len() {
                child_dirichlet.push(dirichlet[r.parent_boundary_edge(BoundaryEdgeId(i))?.0]);
            }
            dirichlet = child_dirichlet;
            mesh = r.mesh().clone();
        }
    }
    Ok(())
}
