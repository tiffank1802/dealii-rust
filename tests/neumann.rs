use dealii_rust::{
    base::{triangle_determinant, Point2},
    grid::{CellId, Mesh2, VertexId},
    lac::CgOptions,
    neumann::{self, NeumannOptions},
    numerics,
};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{} != {}",
        actual,
        expected
    );
}
fn integral(mesh: &Mesh2, values: &[f64]) -> f64 {
    mesh.cells()
        .iter()
        .enumerate()
        .map(|(i, cell)| {
            let det = triangle_determinant(mesh.cell_points(CellId(i)).unwrap()).unwrap();
            det / 6.0 * (values[cell[0].0] + values[cell[1].0] + values[cell[2].0])
        })
        .sum()
}
fn two_components() -> Mesh2 {
    Mesh2::new(
        vec![
            p(0., 0.),
            p(1., 0.),
            p(0., 1.),
            p(3., 0.),
            p(4., 0.),
            p(3., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(3), VertexId(4), VertexId(5)],
        ],
    )
    .unwrap()
}

#[test]
fn zero_data_returns_zero_solution_and_diagnostics() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let result = neumann::solve(&mesh, |_| 0.0, |_, _| 0.0, NeumannOptions::default()).unwrap();
    assert!(result.values.iter().all(|x| *x == 0.0));
    assert_eq!(result.cg.iterations, 0);
    assert_eq!(result.original_residual_norm, 0.0);
    assert_eq!(result.components.len(), 1);
    near(result.components[0].area, 1.0, 1e-14);
    assert_eq!(result.components[0].anchor, VertexId(0));
}

#[test]
fn affine_solution_is_recovered_up_to_integral_mean() {
    let mesh = Mesh2::unit_square(8).unwrap();
    let result = neumann::solve(
        &mesh,
        |_| 0.0,
        |_, n| n[0] + 2.0 * n[1],
        NeumannOptions::default(),
    )
    .unwrap();
    for (value, point) in result.values.iter().zip(mesh.vertices()) {
        near(*value, point.x + 2.0 * point.y - 1.5, 2e-9);
    }
    near(integral(&mesh, &result.values), 0.0, 1e-12);
    assert!(result.original_residual_norm < 1e-8);
    assert!(
        result.components[0].load_imbalance.abs() < result.components[0].compatibility_tolerance
    );
}

fn quadratic(n: usize) -> (f64, f64) {
    let mesh = Mesh2::unit_square(n).unwrap();
    let result = neumann::solve(
        &mesh,
        |_| -4.0,
        |p, normal| 2.0 * p.x * normal[0] + 2.0 * p.y * normal[1],
        NeumannOptions::default(),
    )
    .unwrap();
    near(integral(&mesh, &result.values), 0.0, 1e-11);
    assert!(result.original_residual_norm < 1e-7);
    numerics::errors(
        &mesh,
        &result.values,
        |p| p.x * p.x + p.y * p.y - 2.0 / 3.0,
        |p| [2.0 * p.x, 2.0 * p.y],
    )
    .unwrap()
}

#[test]
fn quadratic_solution_converges_in_l2_and_h1() {
    let a = quadratic(8);
    let b = quadratic(16);
    let c = quadratic(32);
    for (old, new) in [(a, b), (b, c)] {
        let l2 = (old.0 / new.0).log2();
        let h1 = (old.1 / new.1).log2();
        assert!((1.8..2.2).contains(&l2), "L2 rate {}", l2);
        assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
    }
}

#[test]
fn mean_is_integral_weighted_not_a_nodal_average() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let result = neumann::solve(
        &mesh,
        |_| -4.0,
        |p, n| 2.0 * p.x * n[0] + 2.0 * p.y * n[1],
        NeumannOptions::default(),
    )
    .unwrap();
    near(integral(&mesh, &result.values), 0.0, 1e-12);
    let nodal_mean = result.values.iter().sum::<f64>() / result.values.len() as f64;
    assert!(
        nodal_mean.abs() > 1e-3,
        "fixture must distinguish integral and nodal means"
    );
}

#[test]
fn independent_components_are_anchored_and_normalized_separately() {
    let mesh = two_components();
    let result = neumann::solve(
        &mesh,
        |_| 0.0,
        |_, n| n[0] + 2.0 * n[1],
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(result.components.len(), 2);
    assert_eq!(result.components[0].anchor, VertexId(0));
    assert_eq!(result.components[1].anchor, VertexId(3));
    for report in &result.components {
        near(report.area, 0.5, 1e-14);
    }
    for (i, value) in result.values.iter().enumerate() {
        near(*value, [-1.0, 0.0, 1.0][i % 3], 1e-12);
    }
    near(result.original_residual_norm, 0.0, 1e-12);
}

#[test]
fn global_balance_is_not_enough_for_disconnected_components() {
    let mesh = two_components();
    // Each triangle has area 1/2: their loads cancel globally, not locally.
    let err = neumann::solve(
        &mesh,
        |p| if p.x < 2.0 { 1.0 } else { -1.0 },
        |_, _| 0.0,
        NeumannOptions::default(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("incompatible Neumann load"));
}

#[test]
fn near_incompatible_data_is_not_silently_projected() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let options = NeumannOptions {
        compatibility_absolute_tolerance: 1e-5,
        compatibility_relative_tolerance: 0.0,
        ..NeumannOptions::default()
    };
    let result = neumann::solve(&mesh, |_| 1e-6, |_, _| 0.0, options).unwrap();
    near(result.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(result.original_residual_norm > 9e-7);
    assert!(result.cg.residual_norm < 1e-12);
    near(integral(&mesh, &result.values), 0.0, 1e-18);
}

#[test]
fn invalid_inputs_and_nonconvergence_are_reported() {
    let mesh = Mesh2::unit_square(4).unwrap();
    assert!(neumann::solve(&mesh, |_| 1.0, |_, _| 0.0, NeumannOptions::default()).is_err());
    assert!(neumann::solve(&mesh, |_| f64::NAN, |_, _| 0.0, NeumannOptions::default()).is_err());
    assert!(neumann::solve(
        &mesh,
        |_| 0.0,
        |_, _| f64::INFINITY,
        NeumannOptions::default()
    )
    .is_err());
    let invalid = NeumannOptions {
        compatibility_relative_tolerance: -1.0,
        ..NeumannOptions::default()
    };
    assert!(neumann::solve(&mesh, |_| 0.0, |_, _| 0.0, invalid).is_err());
    let invalid = NeumannOptions {
        compatibility_absolute_tolerance: f64::NAN,
        ..NeumannOptions::default()
    };
    assert!(neumann::solve(&mesh, |_| 0.0, |_, _| 0.0, invalid).is_err());
    let no_iterations = NeumannOptions {
        cg: CgOptions {
            max_iterations: 0,
            ..CgOptions::default()
        },
        ..NeumannOptions::default()
    };
    assert!(neumann::solve(&mesh, |_| 0.0, |_, n| n[0] + n[1], no_iterations).is_err());
}

#[test]
fn compatibility_relative_tolerance_uses_absolute_load_scale() {
    let mesh = Mesh2::unit_square(4).unwrap();
    // A small constant flux offset creates a known nonzero global imbalance,
    // while the linear flux supplies cancelling positive and negative loads.
    let loose = NeumannOptions {
        compatibility_absolute_tolerance: 0.0,
        compatibility_relative_tolerance: 1e-5,
        ..NeumannOptions::default()
    };
    let result = neumann::solve(&mesh, |_| 0.0, |_, n| n[0] + 1e-6, loose).unwrap();
    near(result.components[0].load_imbalance, 4e-6, 1e-13);
    assert!(result.components[0].compatibility_tolerance > 4e-6);
    assert!(result.original_residual_norm > 3e-6);
    let strict = NeumannOptions {
        compatibility_relative_tolerance: 1e-8,
        ..loose
    };
    assert!(neumann::solve(&mesh, |_| 0.0, |_, n| n[0] + 1e-6, strict).is_err());
}
