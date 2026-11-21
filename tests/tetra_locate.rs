use dealii_rust::{
    base::Point3,
    grid::{CellId, VertexId},
    tetra::{AffineTet, TetMesh, UniformRefinement},
    tetra_field,
};
fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}, tol {}", a, b, tol);
}
fn single(points: [Point3; 4]) -> TetMesh {
    TetMesh::new(
        points.to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}
fn reference() -> TetMesh {
    single([p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)])
}
fn joined() -> TetMesh {
    TetMesh::new(
        vec![
            p(0., 0., 0.),
            p(1., 0., 0.),
            p(0., 1., 0.),
            p(0., 0., 1.),
            p(1., 1., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(1), VertexId(2), VertexId(3), VertexId(4)],
        ],
    )
    .unwrap()
}
#[test]
fn reference_vertices_edges_faces_and_interior_are_evaluated() {
    let mesh = reference();
    let values = [2., 4., 8., 16.];
    for (point, weights, expected) in [
        (p(0., 0., 0.), [1., 0., 0., 0.], 2.),
        (p(1., 0., 0.), [0., 1., 0., 0.], 4.),
        (p(0., 1., 0.), [0., 0., 1., 0.], 8.),
        (p(0., 0., 1.), [0., 0., 0., 1.], 16.),
        (p(0.5, 0.5, 0.), [0., 0.5, 0.5, 0.], 6.),
        (p(0.25, 0.25, 0.5), [0., 0.25, 0.25, 0.5], 11.),
        (p(0.25, 0.25, 0.25), [0.25; 4], 7.5),
    ] {
        let loc = tetra_field::locate(&mesh, point, 0.).unwrap().unwrap();
        assert_eq!(loc.cell, CellId(0));
        for (a, b) in loc.barycentric.iter().zip(weights) {
            near(*a, b, 1e-14);
        }
        near(
            tetra_field::evaluate(&mesh, &values, point, 0.)
                .unwrap()
                .unwrap(),
            expected,
            1e-14,
        );
    }
}
#[test]
fn skew_tetrahedron_recovers_reference_coordinates_and_field() {
    let points = [p(1., 2., 3.), p(3., 2., 3.), p(2., 5., 3.), p(1., 3., 7.)];
    let mesh = single(points);
    let map = AffineTet::new(points).unwrap();
    let values = [-3., 1., 7., 4.];
    for i in 0..=8 {
        for j in 0..=8 - i {
            for k in 0..=8 - i - j {
                let (x, y, z) = (i as f64 / 8., j as f64 / 8., k as f64 / 8.);
                let point = map.map(p(x, y, z)).unwrap();
                let l = tetra_field::locate(&mesh, point, 1e-12).unwrap().unwrap();
                for (a, b) in l.barycentric.iter().zip([1. - x - y - z, x, y, z]) {
                    near(*a, b, 1e-12);
                }
                let expected =
                    (1. - x - y - z) * values[0] + x * values[1] + y * values[2] + z * values[3];
                near(
                    tetra_field::evaluate(&mesh, &values, point, 1e-12)
                        .unwrap()
                        .unwrap(),
                    expected,
                    1e-11,
                );
            }
        }
    }
}
#[test]
fn shared_face_and_exact_containment_have_priority() {
    let mesh = joined();
    let values = [2., 4., 8., 16., 20.];
    let face = p(0.25, 0.25, 0.5);
    let loc = tetra_field::locate(&mesh, face, 0.).unwrap().unwrap();
    assert_eq!(loc.cell, CellId(0));
    near(
        tetra_field::evaluate(&mesh, &values, face, 0.)
            .unwrap()
            .unwrap(),
        11.,
        1e-13,
    );
    let point = p(0.34, 0.34, 0.34);
    let loc = tetra_field::locate(&mesh, point, 0.1).unwrap().unwrap();
    assert_eq!(loc.cell, CellId(1));
    assert!(loc.barycentric.iter().all(|v| *v >= 0.));
}
#[test]
fn barycentric_tolerance_allows_extrapolation_not_projection() {
    let mesh = reference();
    let point = p(-0.001, 0.25, 0.25);
    assert!(tetra_field::locate(&mesh, point, 0.0009).unwrap().is_none());
    let loc = tetra_field::locate(&mesh, point, 0.0011).unwrap().unwrap();
    near(loc.barycentric[1], -0.001, 1e-15);
    near(
        tetra_field::evaluate(&mesh, &[0., 1., 0., 0.], point, 0.0011)
            .unwrap()
            .unwrap(),
        -0.001,
        1e-15,
    );
    assert!(
        tetra_field::evaluate(&mesh, &[0., 1., 0., 0.], p(10., 10., 10.), 0.1)
            .unwrap()
            .is_none()
    );
}
#[test]
fn disconnected_mesh_selects_first_tolerance_match_but_prefers_exact() {
    let mesh = TetMesh::new(
        vec![
            p(0., 0., 0.),
            p(1., 0., 0.),
            p(0., 1., 0.),
            p(0., 0., 1.),
            p(1.2, 0., 0.),
            p(2.2, 0., 0.),
            p(1.2, 1., 0.),
            p(1.2, 0., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap();
    let loc = tetra_field::locate(&mesh, p(1.1, 0.01, 0.01), 0.2)
        .unwrap()
        .unwrap();
    assert_eq!(loc.cell, CellId(0));
    assert!(loc.barycentric.iter().any(|x| *x < 0.));
    let loc = tetra_field::locate(&mesh, p(1.21, 0.01, 0.01), 0.3)
        .unwrap()
        .unwrap();
    assert_eq!(loc.cell, CellId(1));
    assert!(loc.barycentric.iter().all(|x| *x >= 0.));
    assert!(tetra_field::locate(&mesh, p(5., 5., 5.), 0.01)
        .unwrap()
        .is_none());
}
#[test]
fn barycentric_tolerance_is_invariant_under_scaling() {
    for scale in [0.01, 1., 100.] {
        let mesh = single([
            p(0., 0., 0.),
            p(scale, 0., 0.),
            p(0., scale, 0.),
            p(0., 0., scale),
        ]);
        let pt = p(-0.01 * scale, 0.25 * scale, 0.25 * scale);
        assert!(tetra_field::locate(&mesh, pt, 0.009).unwrap().is_none());
        let loc = tetra_field::locate(&mesh, pt, 0.011).unwrap().unwrap();
        near(loc.barycentric[1], -0.01, 1e-14);
    }
}
#[test]
fn large_translation_preserves_interior_coordinates_to_representable_accuracy() {
    let t = 1e9;
    let points = [
        p(t, t, t),
        p(t + 2., t, t),
        p(t + 1., t + 3., t),
        p(t, t + 1., t + 4.),
    ];
    let mesh = single(points);
    let point = AffineTet::new(points)
        .unwrap()
        .map(p(0.3, 0.2, 0.1))
        .unwrap();
    let loc = tetra_field::locate(&mesh, point, 1e-10).unwrap().unwrap();
    for (a, b) in loc.barycentric.iter().zip([0.4, 0.3, 0.2, 0.1]) {
        near(*a, b, 2e-7);
    }
}
#[test]
fn prolongated_field_matches_piecewise_parent_field_including_faces() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let values: Vec<_> = (0..mesh.vertices().len())
        .map(|i| ((i * i + 3 * i) % 19) as f64)
        .collect();
    let r = UniformRefinement::new(&mesh).unwrap();
    let fine = r.prolongate_p1(&values).unwrap();
    for i in 0..mesh.cells().len() {
        let map = AffineTet::new(mesh.cell_points(CellId(i)).unwrap()).unwrap();
        for q in [
            p(0.1, 0.2, 0.3),
            p(0.5, 0.25, 0.25),
            p(0., 0., 0.),
            p(0.2, 0.3, 0.),
        ] {
            let point = map.map(q).unwrap();
            let a = tetra_field::evaluate(&mesh, &values, point, 1e-12)
                .unwrap()
                .unwrap();
            let b = tetra_field::evaluate(r.mesh(), &fine, point, 1e-12)
                .unwrap()
                .unwrap();
            near(a, b, 1e-10);
        }
    }
}
#[test]
fn detects_bad_inputs_and_unrepresentable_arithmetic() {
    let mesh = reference();
    let point = p(0.1, 0.2, 0.3);
    for bad in [
        Point3 {
            x: f64::NAN,
            y: 0.,
            z: 0.,
        },
        Point3 {
            x: 0.,
            y: f64::INFINITY,
            z: 0.,
        },
    ] {
        assert!(tetra_field::locate(&mesh, bad, 0.).is_err());
    }
    for t in [-1., f64::NAN, f64::INFINITY] {
        assert!(tetra_field::locate(&mesh, point, t).is_err());
    }
    assert!(tetra_field::evaluate(&mesh, &[1.], point, 0.).is_err());
    assert!(tetra_field::evaluate(&mesh, &[f64::NAN; 4], point, 0.).is_err());
    assert!(tetra_field::evaluate(&mesh, &[0., f64::MAX, 0., 0.], p(2., 0., 0.), 1.).is_err());
    assert!(tetra_field::locate(&mesh, p(f64::MAX, f64::MAX, f64::MAX), 0.).is_err());
}
