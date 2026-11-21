//! Interpolate, L2-project and evaluate a quadratic field on P1 tetrahedra.
use dealii_rust::{
    base::Point3,
    grid::CellId,
    lac::CgOptions,
    tetra::{self, TetMesh, UniformRefinement},
    tetra_field,
};
use std::fs::File;
fn f(p: Point3) -> f64 {
    p.x * p.x + p.y * p.y + p.z * p.z
}
fn gradient(p: Point3) -> [f64; 3] {
    [2.0 * p.x, 2.0 * p.y, 2.0 * p.z]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 4,
    };
    let mesh = TetMesh::unit_cube(n)?;
    let interpolated = tetra_field::interpolate(&mesh, f)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    };
    let (projected, report) = tetra_field::project_l2(&mesh, f, options)?;
    let l2_interp = tetra::errors(&mesh, &interpolated, f, gradient)?.0;
    let l2_projection = tetra::errors(&mesh, &projected, f, gradient)?.0;
    let point = Point3::new(0.31, 0.42, 0.21)?;
    let location = tetra_field::locate(&mesh, point, 1e-12)?.ok_or("point outside mesh")?;
    let value =
        tetra_field::evaluate(&mesh, &projected, point, 1e-12)?.ok_or("point outside mesh")?;
    println!(
        "point=({:.2},{:.2},{:.2}) cell={} barycentric={:?} projected_value={:.6e} exact={:.6e}",
        point.x,
        point.y,
        point.z,
        location.cell.0,
        location.barycentric,
        value,
        f(point)
    );
    let outside = Point3::new(-0.0001 / n as f64, 0.42, 0.21)?;
    let location = tetra_field::locate(&mesh, outside, 0.001)?.ok_or("point outside tolerance")?;
    println!(
        "outside_point_extrapolated={} barycentric={:?}",
        location.barycentric.iter().any(|&v| v < 0.0),
        location.barycentric
    );
    let refined = UniformRefinement::new(&mesh)?;
    let prolonged = refined.prolongate_p1(&projected)?;
    println!("n={} vertices={} cells={} iterations={} residual={:.3e} L2_interpolation={:.6e} L2_projection={:.6e}",n,mesh.vertices().len(),mesh.cells().len(),report.iterations,report.residual_norm,l2_interp,l2_projection);
    println!(
        "prolonged_vertices={} first_parent_children={:?}",
        prolonged.len(),
        refined.children(CellId(0))?
    );
    tetra::write_vtk(
        File::create("tetra_field_projection.vtk")?,
        &mesh,
        &projected,
    )?;
    println!("Wrote tetra_field_projection.vtk");
    Ok(())
}
