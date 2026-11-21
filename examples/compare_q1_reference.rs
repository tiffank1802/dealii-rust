//! Compare a deal.II Q1 reference file against the Rust Q1 implementation.
use dealii_rust::{
    base::Point2,
    lac::CgOptions,
    quad::{self, QuadMesh},
    quad_field,
};
use std::{env, fs};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = env::args()
        .nth(1)
        .ok_or("usage: compare_q1_reference <reference-file>")?;
    let text = fs::read_to_string(path)?;
    let numbers: Vec<f64> = text
        .split_whitespace()
        .map(str::parse::<f64>)
        .collect::<Result<_, _>>()?;
    if numbers.len() != 25 || numbers.iter().any(|x| !x.is_finite()) {
        return Err("expected 25 finite numbers".into());
    }
    let mesh = QuadMesh::unit_square(1)?;
    // Mesh vertex order is row-major: (0,0), (1,0), (0,1), (1,1).
    // Reference order is (0,0), (1,0), (1,1), (0,1).
    let order = [0, 1, 3, 2];
    let f = |p: Point2| p.x * p.x + p.y * p.y;
    let (stiffness, _) = quad::assemble_poisson(&mesh, |_| 0.)?;
    let (_, load) = quad_field::projection_system(&mesh, f)?;
    let opts = CgOptions {
        relative_tolerance: 1e-13,
        absolute_tolerance: 1e-15,
        ..CgOptions::default()
    };
    let (solution, _) = quad_field::project_l2(&mesh, f, opts)?;
    let l2 = quad::errors(&mesh, &solution, f, |p| [2. * p.x, 2. * p.y])?.0;
    let mut matrix_diff: f64 = 0.;
    let mut load_diff: f64 = 0.;
    let mut solution_diff: f64 = 0.;
    for i in 0..4 {
        for j in 0..4 {
            matrix_diff =
                matrix_diff.max((stiffness.get(order[i], order[j])? - numbers[5 * i + j]).abs());
        }
        load_diff = load_diff.max((load[order[i]] - numbers[5 * i + 4]).abs());
        solution_diff = solution_diff.max((solution[order[i]] - numbers[20 + i]).abs());
    }
    let error_diff = (l2 - numbers[24]).abs();
    println!(
        "max_abs: stiffness={:.3e} load={:.3e} projection={:.3e} L2_error={:.3e}",
        matrix_diff, load_diff, solution_diff, error_diff
    );
    if [matrix_diff, load_diff, solution_diff, error_diff]
        .iter()
        .any(|&v| v > 1e-10)
    {
        return Err("reference comparison exceeds 1e-10".into());
    }
    Ok(())
}
