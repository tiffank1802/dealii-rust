use dealii_rust::{
    base::Point2,
    grid::{CellId, VertexId},
    lac::CgOptions,
    neumann_quad::{self, NeumannOptions},
    quad::{self, BilinearMap, QuadMesh},
    quadrature::square_degree_5,
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn integral(mesh: &QuadMesh, u: &[f64]) -> f64 {
    let mut sum = 0.;
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i)).unwrap()).unwrap();
        for q in square_degree_5() {
            let values = quad::shape_values(q.point);
            let uh: f64 = (0..4).map(|j| values[j] * u[cell[j].0]).sum();
            sum += q.weight * map.evaluate(q.point).unwrap().determinant * uh;
        }
    }
    sum
}
fn disconnected() -> QuadMesh {
    QuadMesh::new(
        vec![
            p(0., 0.),
            p(2., 0.),
            p(1., 1.),
            p(0., 1.),
            p(3., 0.),
            p(4., 0.),
            p(4., 1.),
            p(3., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap()
}
#[test]
fn zero_problem_and_diagnostics() {
    let mesh = QuadMesh::unit_square(4).unwrap();
    let result = neumann_quad::solve(&mesh, |_| 0., |_, _| 0., NeumannOptions::default()).unwrap();
    assert!(result.values.iter().all(|v| *v == 0.));
    assert_eq!(result.cg.iterations, 0);
    assert_eq!(result.original_residual_norm, 0.);
    assert_eq!(result.components.len(), 1);
    near(result.components[0].area, 1., 1e-13);
    assert_eq!(result.components[0].anchor, VertexId(0));
}
#[test]
fn nonaffine_trapezoid_uses_true_integral_mean() {
    let mesh = disconnected();
    let result = neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, n| n[0] + 2. * n[1],
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(result.components.len(), 2);
    assert_eq!(result.components[0].anchor, VertexId(0));
    assert_eq!(result.components[1].anchor, VertexId(4));
    near(result.components[0].area, 1.5, 1e-13);
    near(result.components[1].area, 1., 1e-13);
    // Trapezoid centroid (7/9,4/9), integral mean x+2y=5/3;
    // analytic Q1 shape-integral weights 5/12,5/12,1/3,1/3.
    for (i, pt) in mesh.vertices().iter().enumerate() {
        let mean = if i < 4 { 5. / 3. } else { 4.5 };
        near(result.values[i], pt.x + 2. * pt.y - mean, 2e-12);
    }
    let u = &result.values;
    near(
        (5. * u[0] + 5. * u[1] + 4. * u[2] + 4. * u[3]) / 12.,
        0.,
        1e-12,
    );
    near(u[4..].iter().sum::<f64>() / 4., 0., 1e-12);
    assert!(u[..4].iter().sum::<f64>().abs() > 0.1); // not a nodal mean
    assert!(result.original_residual_norm < 1e-12);
}
#[test]
fn global_balance_cannot_hide_incompatible_components() {
    let mesh = disconnected();
    let err = neumann_quad::solve(
        &mesh,
        |p| if p.x < 2.5 { 2. } else { -3. },
        |_, _| 0.,
        NeumannOptions::default(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("incompatible Neumann load"));
}
#[test]
fn accepted_imbalance_is_not_projected() {
    let mesh = QuadMesh::unit_square(3).unwrap();
    let options = NeumannOptions {
        compatibility_absolute_tolerance: 1e-5,
        compatibility_relative_tolerance: 0.,
        ..NeumannOptions::default()
    };
    let r = neumann_quad::solve(&mesh, |_| 1e-6, |_, _| 0., options).unwrap();
    near(r.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(r.original_residual_norm > 9e-7);
    assert!(r.cg.residual_norm < 1e-12);
    near(integral(&mesh, &r.values), 0., 1e-17);
}
#[test]
fn relative_tolerance_scales_with_load() {
    let mesh = QuadMesh::unit_square(4).unwrap();
    let options = NeumannOptions {
        compatibility_absolute_tolerance: 0.,
        compatibility_relative_tolerance: 1e-5,
        ..NeumannOptions::default()
    };
    let r = neumann_quad::solve(&mesh, |_| 0., |_, n| n[0] + 1e-6, options).unwrap();
    near(r.components[0].load_imbalance, 4e-6, 1e-13);
    assert!(r.components[0].compatibility_tolerance > 4e-6);
    assert!(r.original_residual_norm > 3e-6);
    assert!(neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, n| n[0] + 1e-6,
        NeumannOptions {
            compatibility_relative_tolerance: 1e-8,
            ..options
        }
    )
    .is_err());
}
fn quadratic(n: usize, distorted: bool) -> (f64, f64) {
    let base = QuadMesh::unit_square(n).unwrap();
    let mesh = if distorted {
        let points = base
            .vertices()
            .iter()
            .map(|p| Point2 {
                x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
                y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
            })
            .collect();
        QuadMesh::new(points, base.cells().to_vec()).unwrap()
    } else {
        base
    };
    let r = neumann_quad::solve(
        &mesh,
        |_| -4.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1]),
        NeumannOptions::default(),
    )
    .unwrap();
    near(integral(&mesh, &r.values), 0., 1e-11);
    assert!(r.original_residual_norm < 1e-7);
    // Domain remains the unit square under the deformation.
    quad::errors(
        &mesh,
        &r.values,
        |p| p.x * p.x + p.y * p.y - 2. / 3.,
        |p| [2. * p.x, 2. * p.y],
    )
    .unwrap()
}
#[test]
fn quadratic_convergence_on_regular_and_deformed_meshes() {
    for distorted in [false, true] {
        let v: Vec<_> = [8, 16, 32]
            .iter()
            .map(|&n| quadratic(n, distorted))
            .collect();
        for pair in v.windows(2) {
            let l2 = (pair[0].0 / pair[1].0).log2();
            let h1 = (pair[0].1 / pair[1].1).log2();
            assert!((1.8..2.2).contains(&l2), "L2 rate {}", l2);
            assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
        }
    }
}
#[test]
fn rejects_bad_data_options_and_nonconvergence() {
    let mesh = QuadMesh::unit_square(4).unwrap();
    assert!(neumann_quad::solve(&mesh, |_| 1., |_, _| 0., NeumannOptions::default()).is_err());
    assert!(
        neumann_quad::solve(&mesh, |_| f64::NAN, |_, _| 0., NeumannOptions::default()).is_err()
    );
    assert!(neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, _| f64::INFINITY,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, _| 0.,
        NeumannOptions {
            compatibility_absolute_tolerance: f64::NAN,
            ..NeumannOptions::default()
        }
    )
    .is_err());
    assert!(neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, _| 0.,
        NeumannOptions {
            compatibility_relative_tolerance: -1.,
            ..NeumannOptions::default()
        }
    )
    .is_err());
    assert!(neumann_quad::solve(
        &mesh,
        |_| 0.,
        |_, n| n[0] + n[1],
        NeumannOptions {
            cg: CgOptions {
                max_iterations: 0,
                ..CgOptions::default()
            },
            ..NeumannOptions::default()
        }
    )
    .is_err());
}
