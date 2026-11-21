//! Interpolate, project and evaluate a quadratic field on Q1 quadrilaterals.
use dealii_rust::{
    base::Point2,
    grid::CellId,
    lac::CgOptions,
    quad::{self, QuadMesh, UniformRefinement},
    quad_field,
};
use std::fs::File;
fn f(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [2.0 * p.x, 2.0 * p.y]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n = match std::env::args().nth(1) {
        Some(s) => s.parse()?,
        None => 8,
    };
    let mesh = QuadMesh::unit_square(n)?;
    let interpolated = quad_field::interpolate(&mesh, f)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    };
    let (projected, report) = quad_field::project_l2(&mesh, f, options)?;
    let l2_interpolated = quad::errors(&mesh, &interpolated, f, gradient)?.0;
    let l2_projected = quad::errors(&mesh, &projected, f, gradient)?.0;
    let point = Point2::new(0.31, 0.42)?;
    let loc = quad_field::locate(&mesh, point, 1e-12)?.ok_or("point outside mesh")?;
    let value =
        quad_field::evaluate(&mesh, &projected, point, 1e-12)?.ok_or("point outside mesh")?;
    println!(
        "point=({:.2},{:.2}) cell={} reference=({:.5},{:.5}) projected_value={:.6e} exact={:.6e}",
        point.x,
        point.y,
        loc.cell.0,
        loc.reference.x,
        loc.reference.y,
        value,
        f(point)
    );
    let outside = Point2::new(-0.0001 / n as f64, 0.42)?;
    let boundary = quad_field::locate(&mesh, outside, 0.001)?.ok_or("point outside tolerance")?;
    println!(
        "outside_point_projected={} distance={:.3e} boundary=({:.6},{:.6})",
        boundary.projected, boundary.distance, boundary.point.x, boundary.point.y
    );
    let refined = UniformRefinement::new(&mesh)?;
    let prolonged = refined.prolongate_q1(&projected)?;
    // Refinement preserves the projected field; this is not a new projection.
    let child = refined.children(CellId(0))?;
    println!("n={} vertices={} cells={} iterations={} residual={:.3e} L2_interpolation={:.6e} L2_projection={:.6e}",n,mesh.vertices().len(),mesh.cells().len(),report.iterations,report.residual_norm,l2_interpolated,l2_projected);
    println!(
        "prolonged_vertices={} first_parent_children={:?}",
        prolonged.len(),
        child
    );
    quad::write_vtk(
        File::create("quad_field_projection.vtk")?,
        &mesh,
        &projected,
    )?;
    println!("Wrote quad_field_projection.vtk");
    Ok(())
}
