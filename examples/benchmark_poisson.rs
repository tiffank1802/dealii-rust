//! Small serial timing harness, not a statistically rigorous benchmark.
use dealii_rust::{
    base::Point2,
    grid::Mesh2,
    lac::{conjugate_gradient, CgOptions},
    numerics,
};
use std::{f64::consts::PI, time::Instant};
fn source(p: Point2) -> f64 {
    2. * PI * PI * (PI * p.x).sin() * (PI * p.y).sin()
}
fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = std::env::args().nth(1).map_or(Ok(32), |s| s.parse())?;
    let repetitions: usize = std::env::args().nth(2).map_or(Ok(3), |s| s.parse())?;
    if n == 0 || repetitions == 0 {
        return Err("n and repetitions must be positive".into());
    }
    let mut mesh_times = Vec::new();
    let mut assembly_times = Vec::new();
    let mut constraint_times = Vec::new();
    let mut solve_times = Vec::new();
    let mut last = None;
    for _ in 0..repetitions {
        let start = Instant::now();
        let mesh = Mesh2::unit_square(n)?;
        mesh_times.push(start.elapsed().as_secs_f64() * 1000.);
        let start = Instant::now();
        let (a, b) = numerics::assemble_poisson(&mesh, source)?;
        assembly_times.push(start.elapsed().as_secs_f64() * 1000.);
        let start = Instant::now();
        let constraints: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
        let (a, b) = a.with_dirichlet(&b, &constraints)?;
        constraint_times.push(start.elapsed().as_secs_f64() * 1000.);
        let start = Instant::now();
        let (u, report) = conjugate_gradient(&a, &b, CgOptions::default())?;
        solve_times.push(start.elapsed().as_secs_f64() * 1000.);
        let residual = a
            .mul(&u)?
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y) * (x - y))
            .sum::<f64>()
            .sqrt();
        if !residual.is_finite() {
            return Err("nonfinite true residual".into());
        }
        last = Some((
            mesh.vertices().len(),
            mesh.cells().len(),
            report.iterations,
            residual,
            a.stored_entries(),
            a.logical_storage_bytes(),
        ));
    }
    let (vertices, cells, iterations, residual, stored_entries, csr_bytes) =
        last.ok_or("no runs")?;
    println!(
        "n={} repeats={} vertices={} cells={} iterations={} true_residual={:.3e} csr_entries={} csr_logical_bytes={}",
        n, repetitions, vertices, cells, iterations, residual, stored_entries, csr_bytes
    );
    println!(
        "median_ms mesh={:.3} assembly={:.3} dirichlet={:.3} solve={:.3}",
        median(&mut mesh_times),
        median(&mut assembly_times),
        median(&mut constraint_times),
        median(&mut solve_times)
    );
    Ok(())
}
