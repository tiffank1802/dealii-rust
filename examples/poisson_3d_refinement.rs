//! Global tetrahedral refinement for a mixed Poisson problem.
use dealii_rust::{
    base::Point3,
    lac::{conjugate_gradient, conjugate_gradient_with_initial_guess, CgOptions},
    tetra::{self, BoundaryFaceId, TetMesh, UniformRefinement},
};
use std::fs::File;
fn exact(p: Point3) -> f64 {
    1.0 + p.x + p.x * p.x + 2.0 * p.y + p.y * p.y + 3.0 * p.z + p.z * p.z
}
fn gradient(p: Point3) -> [f64; 3] {
    [1.0 + 2.0 * p.x, 2.0 + 2.0 * p.y, 3.0 + 2.0 * p.z]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut mesh = TetMesh::unit_cube(2)?;
    let mut dirichlet = Vec::new();
    for i in 0..mesh.boundary_faces().len() {
        dirichlet.push(
            mesh.boundary_face_geometry(BoundaryFaceId(i))?
                .points
                .iter()
                .all(|p| p.x == 0.0),
        );
    }
    let mut initial: Option<Vec<f64>> = None;
    for level in 0..4 {
        let (a, mut b) = tetra::assemble_poisson(&mesh, |_| -6.0)?;
        let mut d = Vec::new();
        let mut n = Vec::new();
        for (i, &is_d) in dirichlet.iter().enumerate() {
            if is_d {
                d.push(BoundaryFaceId(i));
            } else {
                n.push(BoundaryFaceId(i));
            }
        }
        let load = tetra::assemble_neumann(&mesh, &n, |p, n| {
            let g = gradient(p);
            g[0] * n[0] + g[1] * n[1] + g[2] * n[2]
        })?;
        for (b, g) in b.iter_mut().zip(load) {
            *b += g;
        }
        let constraints = tetra::interpolate_boundary(&mesh, &d, exact)?;
        let (a, b) = a.with_dirichlet(&b, &constraints)?;
        let options = CgOptions {
            relative_tolerance: 1e-12,
            absolute_tolerance: 1e-13,
            ..CgOptions::default()
        };
        let (u, report) = match initial.take() {
            Some(mut guess) => {
                // Newly created boundary nodes must satisfy the fine-grid data.
                for &(i, value) in &constraints {
                    guess[i] = value;
                }
                conjugate_gradient_with_initial_guess(&a, &b, &guess, options)?
            }
            None => conjugate_gradient(&a, &b, options)?,
        };
        let (l2, h1) = tetra::errors(&mesh, &u, exact, gradient)?;
        println!("level={} vertices={} cells={} dirichlet_faces={} neumann_faces={} iterations={} residual={:.3e} L2={:.3e} H1={:.3e}",level,mesh.vertices().len(),mesh.cells().len(),d.len(),n.len(),report.iterations,report.residual_norm,l2,h1);
        if level == 3 {
            tetra::write_vtk(File::create("poisson_3d_refinement.vtk")?, &mesh, &u)?;
            println!("Wrote poisson_3d_refinement.vtk");
        } else {
            let r = UniformRefinement::new(&mesh)?;
            initial = Some(r.prolongate_p1(&u)?);
            let mut labels = Vec::new();
            for i in 0..r.mesh().boundary_faces().len() {
                labels.push(dirichlet[r.parent_boundary_face(BoundaryFaceId(i))?.0]);
            }
            dirichlet = labels;
            mesh = r.mesh().clone();
        }
    }
    Ok(())
}
