//! Read a tagged Gmsh 4.1 or 2.2 ASCII mesh and solve a mixed Poisson problem.
//!
//! This example expects physical curve 40 on the Dirichlet boundary and
//! curves 10, 20, 30 on the Neumann boundary, as in the provided fixture.
use std::fs::File;

use dealii_rust::{
    base::Point2,
    gmsh,
    grid::BoundaryEdgeId,
    lac::{conjugate_gradient, CgOptions},
    numerics, vtk,
};

fn exact(p: Point2) -> f64 {
    1.0 + p.x + 2.0 * p.y
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let first = args.next();
    let use_22 = first.as_deref() == Some("--msh22");
    let path = if use_22 {
        args.next()
            .unwrap_or_else(|| "tests/fixtures/square22.msh".to_string())
    } else {
        first.unwrap_or_else(|| "tests/fixtures/square41.msh".to_string())
    };
    if args.next().is_some() {
        return Err("usage: poisson_gmsh [--msh22] [path]".into());
    }
    let file = File::open(path)?;
    let imported = if use_22 {
        gmsh::read_msh22(file)?
    } else {
        gmsh::read(file)?
    };
    let mesh = imported.mesh();
    let left_edges = imported.boundary_edges_with_tag(40);
    if left_edges.is_empty() {
        return Err("missing Dirichlet physical curve 40".into());
    }
    let mut neumann_edges = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        let tags = imported.boundary_physical_tags(id)?;
        let dirichlet = tags.contains(&40);
        let neumann = [10, 20, 30].iter().any(|tag| tags.contains(tag));
        if dirichlet == neumann {
            return Err("every boundary edge must have exactly one boundary type".into());
        }
        if neumann {
            neumann_edges.push(id);
        }
    }
    let (matrix, mut rhs) = numerics::assemble_poisson(mesh, |_| 0.0)?;
    let load = numerics::assemble_neumann(mesh, &neumann_edges, |_, normal| {
        normal[0] + 2.0 * normal[1]
    })?;
    for (value, contribution) in rhs.iter_mut().zip(load) {
        *value += contribution;
    }
    let constraints: Vec<_> = imported
        .boundary_vertices_with_tag(40)
        .into_iter()
        .map(|v| (v.0, exact(mesh.vertices()[v.0])))
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints)?;
    let (solution, report) = conjugate_gradient(&matrix, &rhs, CgOptions::default())?;
    let (l2, h1) = numerics::errors(mesh, &solution, exact, |_| [1.0, 2.0])?;
    vtk::write_scalar(File::create("poisson_gmsh.vtk")?, mesh, &solution)?;
    println!("vertices={} triangles={} dirichlet_edges={} neumann_edges={} iterations={} residual={:.3e} L2={:.3e} H1={:.3e}",
        mesh.vertices().len(), mesh.cells().len(), left_edges.len(), neumann_edges.len(), report.iterations, report.residual_norm, l2, h1);
    println!("Wrote poisson_gmsh.vtk");
    Ok(())
}
