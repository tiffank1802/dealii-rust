use dealii_rust::{
    base::Point3,
    grid::{CellId, VertexId},
    lac::CgOptions,
    neumann::NeumannOptions,
    neumann3d,
    tetra::{self, AffineTet, TetMesh},
};
fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn integral(mesh: &TetMesh, u: &[f64]) -> f64 {
    mesh.cells()
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let w = AffineTet::new(mesh.cell_points(CellId(i)).unwrap())
                .unwrap()
                .determinant()
                / 24.0;
            w * c.iter().map(|v| u[v.0]).sum::<f64>()
        })
        .sum()
}
fn disconnected() -> TetMesh {
    TetMesh::new(
        vec![
            p(0., 0., 0.),
            p(1., 0., 0.),
            p(0., 1., 0.),
            p(0., 0., 1.),
            p(3., 0., 0.),
            p(4., 0., 0.),
            p(3., 1., 0.),
            p(3., 0., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap()
}

#[test]
fn zero_data_returns_zero_solution_and_diagnostics() {
    let mesh = TetMesh::unit_cube(3).unwrap();
    let result = neumann3d::solve(&mesh, |_| 0., |_, _| 0., NeumannOptions::default()).unwrap();
    assert!(result.values.iter().all(|x| *x == 0.));
    assert_eq!(result.cg.iterations, 0);
    assert_eq!(result.original_residual_norm, 0.);
    assert_eq!(result.components.len(), 1);
    near(result.components[0].volume, 1., 1e-13);
    assert_eq!(result.components[0].anchor, VertexId(0));
}

#[test]
fn affine_solution_reproduced_with_zero_integral_mean() {
    let mesh = TetMesh::unit_cube(4).unwrap();
    let result = neumann3d::solve(
        &mesh,
        |_| 0.,
        |_, n| n[0] + 2. * n[1] + 3. * n[2],
        NeumannOptions::default(),
    )
    .unwrap();
    for (value, pt) in result.values.iter().zip(mesh.vertices()) {
        near(*value, pt.x + 2. * pt.y + 3. * pt.z - 3., 3e-9);
    }
    near(integral(&mesh, &result.values), 0., 1e-12);
    assert!(result.original_residual_norm < 1e-8);
    assert!(
        result.components[0].load_imbalance.abs() < result.components[0].compatibility_tolerance
    );
}

fn quadratic(n: usize) -> (f64, f64) {
    let mesh = TetMesh::unit_cube(n).unwrap();
    let result = neumann3d::solve(
        &mesh,
        |_| -6.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1] + p.z * n[2]),
        NeumannOptions::default(),
    )
    .unwrap();
    near(integral(&mesh, &result.values), 0., 1e-11);
    assert!(result.original_residual_norm < 1e-7);
    tetra::errors(
        &mesh,
        &result.values,
        |p| p.x * p.x + p.y * p.y + p.z * p.z - 1.,
        |p| [2. * p.x, 2. * p.y, 2. * p.z],
    )
    .unwrap()
}
#[test]
fn quadratic_solution_converges_in_l2_and_h1() {
    let results: Vec<_> = [4, 8, 16].iter().map(|&n| quadratic(n)).collect();
    for pair in results.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
        assert!((0.8..1.2).contains(&h1), "H1 rate {}", h1);
    }
}

#[test]
fn integral_mean_differs_from_arithmetic_nodal_average() {
    let mesh = TetMesh::unit_cube(4).unwrap();
    let result = neumann3d::solve(
        &mesh,
        |_| -6.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1] + p.z * n[2]),
        NeumannOptions::default(),
    )
    .unwrap();
    near(integral(&mesh, &result.values), 0., 1e-12);
    let nodal_mean = result.values.iter().sum::<f64>() / result.values.len() as f64;
    assert!(nodal_mean.abs() > 1e-3);
}

#[test]
fn disconnected_components_are_normalized_independently() {
    let mesh = disconnected();
    let result = neumann3d::solve(
        &mesh,
        |_| 0.,
        |_, n| n[0] + 2. * n[1] + 3. * n[2],
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(result.components.len(), 2);
    assert_eq!(result.components[0].anchor, VertexId(0));
    assert_eq!(result.components[1].anchor, VertexId(4));
    for report in result.components {
        near(report.volume, 1. / 6., 1e-14);
    }
    for (i, v) in result.values.iter().enumerate() {
        near(*v, [-1.5, -0.5, 0.5, 1.5][i % 4], 1e-12);
    }
    near(result.original_residual_norm, 0., 1e-12);
}

#[test]
fn globally_balanced_but_locally_incompatible_load_is_rejected() {
    let mesh = disconnected();
    let err = neumann3d::solve(
        &mesh,
        |p| if p.x < 2. { 1. } else { -1. },
        |_, _| 0.,
        NeumannOptions::default(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("incompatible Neumann load"));
}

#[test]
fn accepted_imbalance_is_not_projected_and_original_residual_is_reported() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let options = NeumannOptions {
        compatibility_absolute_tolerance: 1e-5,
        compatibility_relative_tolerance: 0.,
        ..NeumannOptions::default()
    };
    let result = neumann3d::solve(&mesh, |_| 1e-6, |_, _| 0., options).unwrap();
    near(result.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(result.original_residual_norm > 9e-7);
    assert!(result.cg.residual_norm < 1e-12);
    near(integral(&mesh, &result.values), 0., 1e-18);
}

#[test]
fn relative_tolerance_scales_with_absolute_load() {
    let mesh = TetMesh::unit_cube(3).unwrap();
    let opts = NeumannOptions {
        compatibility_absolute_tolerance: 0.,
        compatibility_relative_tolerance: 1e-5,
        ..NeumannOptions::default()
    };
    let result = neumann3d::solve(&mesh, |_| 0., |_, n| n[0] + 1e-6, opts).unwrap();
    near(result.components[0].load_imbalance, 6e-6, 1e-13);
    assert!(result.components[0].compatibility_tolerance > 6e-6);
    assert!(result.original_residual_norm > 5e-6);
    let strict = NeumannOptions {
        compatibility_relative_tolerance: 1e-8,
        ..opts
    };
    assert!(neumann3d::solve(&mesh, |_| 0., |_, n| n[0] + 1e-6, strict).is_err());
}

#[test]
fn invalid_data_tolerances_and_nonconvergence_are_reported() {
    let mesh = TetMesh::unit_cube(3).unwrap();
    assert!(neumann3d::solve(&mesh, |_| 1., |_, _| 0., NeumannOptions::default()).is_err());
    assert!(neumann3d::solve(&mesh, |_| f64::NAN, |_, _| 0., NeumannOptions::default()).is_err());
    assert!(neumann3d::solve(
        &mesh,
        |_| 0.,
        |_, _| f64::INFINITY,
        NeumannOptions::default()
    )
    .is_err());
    let invalid = NeumannOptions {
        compatibility_absolute_tolerance: f64::NAN,
        ..NeumannOptions::default()
    };
    assert!(neumann3d::solve(&mesh, |_| 0., |_, _| 0., invalid).is_err());
    let invalid = NeumannOptions {
        compatibility_relative_tolerance: -1.,
        ..NeumannOptions::default()
    };
    assert!(neumann3d::solve(&mesh, |_| 0., |_, _| 0., invalid).is_err());
    let opts = NeumannOptions {
        cg: CgOptions {
            max_iterations: 0,
            ..CgOptions::default()
        },
        ..NeumannOptions::default()
    };
    assert!(neumann3d::solve(&mesh, |_| 0., |_, n| n[0] + n[1], opts).is_err());
}
