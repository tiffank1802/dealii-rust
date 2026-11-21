use dealii_rust::{
    base::Point2,
    grid::{CellId, VertexId},
    quad::{self, BilinearMap, QuadMesh, UniformRefinement},
    quad_field,
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn one(points: [Point2; 4]) -> QuadMesh {
    QuadMesh::new(
        points.to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}
fn trapezoid() -> QuadMesh {
    one([p(0., 0.), p(2., 0.), p(1., 1.), p(0., 1.)])
}
#[test]
fn inverse_recovers_bilinear_reference_points_and_q1_field() {
    let mesh = trapezoid();
    let values = [0., 1., 1., 0.];
    for i in 1..20 {
        for j in 1..20 {
            let xi = i as f64 / 20.;
            let eta = j as f64 / 20.;
            let point = p(xi * (2. - eta), eta);
            let found = quad_field::locate(&mesh, point, 0.).unwrap().unwrap();
            assert_eq!(found.cell, CellId(0));
            assert!(!found.projected);
            assert_eq!(found.distance, 0.);
            near(found.reference.x, xi, 1e-13);
            near(found.reference.y, eta, 1e-13);
            near(
                quad_field::evaluate(&mesh, &values, point, 0.)
                    .unwrap()
                    .unwrap(),
                xi,
                1e-13,
            );
        }
    }
}
#[test]
fn skew_nonaffine_cell_round_trips_reference_grid_and_nodal_values() {
    let points = [p(1., 2.), p(4., 2.5), p(3.5, 5.), p(0.5, 4.)];
    let mesh = one(points);
    let map = BilinearMap::new(points).unwrap();
    let values = [2., -1., 7., 3.];
    for i in 0..=10 {
        for j in 0..=10 {
            let r = p(i as f64 / 10., j as f64 / 10.);
            let point = map.evaluate(r).unwrap().point;
            let loc = quad_field::locate(&mesh, point, 1e-12).unwrap().unwrap();
            near(loc.reference.x, r.x, 2e-12);
            near(loc.reference.y, r.y, 2e-12);
            let shape = quad::shape_values(r);
            let expected: f64 = (0..4).map(|k| shape[k] * values[k]).sum();
            near(
                quad_field::evaluate(&mesh, &values, point, 1e-12)
                    .unwrap()
                    .unwrap(),
                expected,
                2e-11,
            );
        }
    }
}
#[test]
fn shared_edges_use_first_exact_cell_and_actual_cell_beats_tolerance() {
    let mesh = QuadMesh::unit_square(2).unwrap();
    let values: Vec<_> = mesh.vertices().iter().map(|p| p.x + 2. * p.y).collect();
    let on_edge = p(0.5, 0.25);
    let location = quad_field::locate(&mesh, on_edge, 0.).unwrap().unwrap();
    assert_eq!(location.cell, CellId(0));
    assert!(!location.projected);
    near(
        quad_field::evaluate(&mesh, &values, on_edge, 0.)
            .unwrap()
            .unwrap(),
        1.,
        1e-13,
    );
    let inside_right = p(0.51, 0.25);
    let location = quad_field::locate(&mesh, inside_right, 0.1)
        .unwrap()
        .unwrap();
    assert_eq!(location.cell, CellId(1));
    assert!(!location.projected);
    near(
        quad_field::evaluate(&mesh, &values, inside_right, 0.1)
            .unwrap()
            .unwrap(),
        1.01,
        1e-13,
    );
}
#[test]
fn tolerance_projects_onto_nearest_skew_boundary_not_extrapolating() {
    let mesh = trapezoid();
    let delta = 0.01;
    let component = delta / 2.0_f64.sqrt();
    let point = p(1.5 + component, 0.5 + component);
    assert!(quad_field::locate(&mesh, point, 0.).unwrap().is_none());
    assert!(quad_field::locate(&mesh, point, 0.0049).unwrap().is_none()); // longest edge=2
    let loc = quad_field::locate(&mesh, point, 0.0051).unwrap().unwrap();
    assert!(loc.projected);
    near(loc.distance, delta, 1e-13);
    near(loc.point.x, 1.5, 1e-13);
    near(loc.point.y, 0.5, 1e-13);
    near(loc.reference.x, 1., 1e-13);
    near(loc.reference.y, 0.5, 1e-13);
    near(
        quad_field::evaluate(&mesh, &[0., 1., 1., 0.], point, 0.0051)
            .unwrap()
            .unwrap(),
        1.,
        1e-13,
    );
}
#[test]
fn tolerance_chooses_first_eligible_cell_but_prefers_any_exact_cell() {
    let mesh = QuadMesh::new(
        vec![
            p(0., 0.),
            p(1., 0.),
            p(1., 1.),
            p(0., 1.),
            p(1.3, 0.),
            p(2.3, 0.),
            p(2.3, 1.),
            p(1.3, 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap();
    let loc = quad_field::locate(&mesh, p(1.2, 0.5), 0.3)
        .unwrap()
        .unwrap();
    assert_eq!(loc.cell, CellId(0));
    assert!(loc.projected);
    near(loc.distance, 0.2, 1e-13);
    let loc = quad_field::locate(&mesh, p(1.31, 0.5), 0.4)
        .unwrap()
        .unwrap();
    assert_eq!(loc.cell, CellId(1));
    assert!(!loc.projected);
    assert!(quad_field::locate(&mesh, p(10., 10.), 0.01)
        .unwrap()
        .is_none());
}
#[test]
fn projected_corner_and_axis_aligned_edges_have_consistent_reference_coordinates() {
    let mesh = QuadMesh::unit_square(1).unwrap();
    let values: Vec<_> = mesh.vertices().iter().map(|p| p.x + 2. * p.y).collect();
    for (query, expected, value) in [
        (p(-0.01, 0.2), p(0., 0.2), 0.4),
        (p(1.01, 0.2), p(1., 0.2), 1.4),
        (p(0.4, -0.01), p(0.4, 0.), 0.4),
        (p(0.4, 1.01), p(0.4, 1.), 2.4),
        (p(-0.01, -0.01), p(0., 0.), 0.),
    ] {
        let loc = quad_field::locate(&mesh, query, 0.02).unwrap().unwrap();
        assert!(loc.projected);
        near(loc.point.x, expected.x, 1e-13);
        near(loc.point.y, expected.y, 1e-13);
        near(
            quad_field::evaluate(&mesh, &values, query, 0.02)
                .unwrap()
                .unwrap(),
            value,
            1e-13,
        );
    }
}
#[test]
fn tolerance_scales_with_cell_size() {
    for scale in [0.01, 1., 100.] {
        let mesh = one([p(0., 0.), p(scale, 0.), p(scale, scale), p(0., scale)]);
        let point = p(-0.01 * scale, 0.5 * scale);
        assert!(quad_field::locate(&mesh, point, 0.009).unwrap().is_none());
        let loc = quad_field::locate(&mesh, point, 0.011).unwrap().unwrap();
        near(loc.distance / scale, 0.01, 1e-13);
        assert!(loc.projected);
    }
}
#[test]
fn evaluation_of_prolonged_field_matches_parent_inside_cells() {
    let base = QuadMesh::unit_square(3).unwrap();
    let points = base
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
        })
        .collect();
    let mesh = QuadMesh::new(points, base.cells().to_vec()).unwrap();
    let values: Vec<_> = (0..mesh.vertices().len())
        .map(|i| ((i * i + 3 * i) % 19) as f64)
        .collect();
    let r = UniformRefinement::new(&mesh).unwrap();
    let fine = r.prolongate_q1(&values).unwrap();
    for i in 0..mesh.cells().len() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i)).unwrap()).unwrap();
        for q in [p(0.1, 0.2), p(0.5, 0.5), p(0.8, 0.7), p(0., 0.3)] {
            let point = map.evaluate(q).unwrap().point;
            let old = quad_field::evaluate(&mesh, &values, point, 1e-11)
                .unwrap()
                .unwrap();
            let new = quad_field::evaluate(r.mesh(), &fine, point, 1e-11)
                .unwrap()
                .unwrap();
            near(old, new, 1e-9);
        }
    }
}
#[test]
fn translated_cell_and_deformed_affine_function() {
    let origin = 1e9;
    let mesh = one([
        p(origin, origin),
        p(origin + 2., origin),
        p(origin + 1.5, origin + 1.),
        p(origin, origin + 1.),
    ]);
    let map = BilinearMap::new(mesh.cell_points(CellId(0)).unwrap()).unwrap();
    let r = p(0.3, 0.7);
    let point = map.evaluate(r).unwrap().point;
    let loc = quad_field::locate(&mesh, point, 1e-8).unwrap().unwrap();
    near(loc.reference.x, r.x, 2e-7);
    near(loc.reference.y, r.y, 2e-7);
}
#[test]
fn rejects_invalid_inputs_and_extreme_arithmetic() {
    let mesh = trapezoid();
    let point = p(1., 0.5);
    assert!(quad_field::locate(&mesh, Point2 { x: f64::NAN, y: 0. }, 0.).is_err());
    assert!(quad_field::locate(&mesh, point, -1.).is_err());
    assert!(quad_field::locate(&mesh, point, f64::INFINITY).is_err());
    assert!(quad_field::locate(&mesh, p(10., 10.), f64::MAX).is_err());
    assert!(quad_field::evaluate(&mesh, &[1.], point, 0.).is_err());
    assert!(quad_field::evaluate(&mesh, &[f64::NAN; 4], point, 0.).is_err());
    assert!(quad_field::evaluate(
        &mesh,
        &[1.; 4],
        Point2 {
            x: 0.,
            y: f64::INFINITY
        },
        0.
    )
    .is_err());
    assert!(quad_field::evaluate(&mesh, &[1.; 4], p(10., 10.), 0.)
        .unwrap()
        .is_none());
}
