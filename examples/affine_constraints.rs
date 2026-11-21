//! Manually assembled P1 hanging-node patch: algebraic constraints only.
//! No nonconforming mesh topology or automatic hanging-node discovery is used.
use dealii_rust::{
    base::Point2,
    constraints::AffineConstraints,
    integrators::{diffusion_reaction_p1, DiffusionTensor},
    lac::{conjugate_gradient, CgOptions, CooBuilder},
};
fn exact(p: Point2) -> f64 {
    1. + p.x + 2. * p.y
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let p = |x, y| Point2::new(x, y).unwrap();
    let points = [
        p(0., 0.),
        p(0.5, 0.),
        p(0.5, 1.),
        p(0., 1.),
        p(1., 0.),
        p(1., 1.),
        p(0.5, 0.5),
        p(0.75, 0.5),
    ];
    // Coarse edge (1,2), fine edges (1,6) and (6,2); 7 is interior.
    let cells = [
        [0, 1, 2],
        [0, 2, 3],
        [1, 4, 7],
        [4, 5, 7],
        [5, 2, 7],
        [2, 6, 7],
        [6, 1, 7],
    ];
    let mut coo = CooBuilder::new(points.len(), points.len());
    for cell in cells {
        let local = diffusion_reaction_p1(
            cell.map(|i| points[i]),
            |_| DiffusionTensor::identity(),
            |_| 0.,
            |_| 0.,
        )?;
        for i in 0..3 {
            for j in 0..3 {
                coo.add(cell[i], cell[j], local.matrix[i][j])?;
            }
        }
    }
    let a = coo.into_csr()?;
    let mut c = AffineConstraints::new(points.len());
    // Boundary conditions and the hanging-node rule are explicit inputs.
    for (i, &point) in points.iter().take(6).enumerate() {
        c.add(i, exact(point), &[])?;
    }
    c.add(6, 0., &[(1, 0.5), (2, 0.5)])?;
    let c = c.close()?;
    let (reduced, b) = c.condense(&a, &[0.; 8])?;
    let (y, report) = conjugate_gradient(&reduced, &b, CgOptions::default())?;
    let x = c.reconstruct(&y)?;
    let max_error = x
        .iter()
        .zip(points)
        .map(|(&v, p)| (v - exact(p)).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "full_dofs={} free_dofs={} iterations={} residual={:.3e} max_error={:.3e}",
        points.len(),
        c.free_dofs().len(),
        report.iterations,
        report.residual_norm,
        max_error
    );
    println!(
        "hanging_value={:.6} endpoint_mean={:.6}",
        x[6],
        (x[1] + x[2]) * 0.5
    );
    Ok(())
}
