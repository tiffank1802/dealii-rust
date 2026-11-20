//! Interpolate, project and evaluate a quadratic field on a P1 mesh.
use dealii_rust::{base::Point2, field, grid::Mesh2, lac::CgOptions, numerics};

fn target(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [2.0 * p.x, 2.0 * p.y]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = Mesh2::unit_square(8)?;
    let interpolated = field::interpolate(&mesh, target)?;
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    };
    let (projected, report) = field::project_l2(&mesh, target, options)?;
    let (interpolation_error, _) = numerics::errors(&mesh, &interpolated, target, gradient)?;
    let (projection_error, _) = numerics::errors(&mesh, &projected, target, gradient)?;
    let point = Point2::new(0.31, 0.42)?;
    let value = field::evaluate(&mesh, &projected, point, 1e-12)?.ok_or("point outside mesh")?;
    println!(
        "vertices={} cells={} iterations={} residual={:.3e}",
        mesh.vertices().len(),
        mesh.cells().len(),
        report.iterations,
        report.residual_norm
    );
    println!(
        "L2 interpolation={:.6e} projection={:.6e}",
        interpolation_error, projection_error
    );
    println!(
        "point=({:.2}, {:.2}) projected={:.6e} exact={:.6e}",
        point.x,
        point.y,
        value,
        target(point)
    );
    Ok(())
}
