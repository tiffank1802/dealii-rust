#![allow(clippy::needless_range_loop)] // Readable analytical matrix fixtures.
use dealii_rust::{
    base::Point2,
    grid::{CellId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    quad::{self, BilinearMap, QuadMesh},
    quadrature::square_degree_5,
};
use std::f64::consts::PI;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn square() -> [Point2; 4] {
    [p(0., 0.), p(1., 0.), p(1., 1.), p(0., 1.)]
}

#[test]
fn tensor_quadrature_integrates_polynomial_degrees_five_per_axis() {
    for a in 0..=5 {
        for b in 0..=5 {
            let actual: f64 = square_degree_5()
                .iter()
                .map(|q| q.weight * q.point.x.powi(a) * q.point.y.powi(b))
                .sum();
            near(actual, 1. / ((a + 1) * (b + 1)) as f64, 1e-14);
        }
    }
}

#[test]
fn shape_functions_and_reference_gradients() {
    for (i, &pt) in square().iter().enumerate() {
        let vals = quad::shape_values(pt);
        for j in 0..4 {
            near(vals[j], if i == j { 1. } else { 0. }, 1e-15);
        }
    }
    for pt in [p(0., 0.), p(1., 1.), p(0.25, 0.72)] {
        near(quad::shape_values(pt).iter().sum(), 1., 1e-15);
        let grad = quad::reference_gradients(pt);
        for d in 0..2 {
            near(grad.iter().map(|g| g[d]).sum(), 0., 1e-15);
        }
    }
}

#[test]
fn bilinear_map_area_and_physical_gradients_on_trapezoid() {
    let points = [p(0., 0.), p(2., 0.), p(3., 1.), p(0., 1.)];
    let map = BilinearMap::new(points).unwrap();
    for (i, &q) in square().iter().enumerate() {
        assert_eq!(map.evaluate(q).unwrap().point, points[i]);
    }
    let q = p(0.4, 0.6);
    let mapped = map.evaluate(q).unwrap();
    near(mapped.point.x, 1.04, 1e-15);
    near(mapped.point.y, 0.6, 1e-15);
    near(mapped.determinant, 2.6, 1e-15);
    for d in 0..2 {
        near(mapped.gradients.iter().map(|g| g[d]).sum(), 0., 1e-14);
        near(
            (0..4).map(|i| points[i].x * mapped.gradients[i][d]).sum(),
            if d == 0 { 1. } else { 0. },
            1e-14,
        );
        near(
            (0..4).map(|i| points[i].y * mapped.gradients[i][d]).sum(),
            if d == 1 { 1. } else { 0. },
            1e-14,
        );
    }
    let area: f64 = square_degree_5()
        .iter()
        .map(|q| map.evaluate(q.point).unwrap().determinant * q.weight)
        .sum();
    near(area, 2.5, 1e-14);
    assert!(map.evaluate(p(-0.01, 0.5)).is_err());
    assert!(map
        .evaluate(Point2 {
            x: f64::NAN,
            y: 0.5
        })
        .is_err());
}

#[test]
fn q1_mesh_counts_boundaries_and_rejects_bad_cells() {
    for n in [1, 2, 5] {
        let mesh = QuadMesh::unit_square(n).unwrap();
        assert_eq!(mesh.vertices().len(), (n + 1) * (n + 1));
        assert_eq!(mesh.cells().len(), n * n);
        assert_eq!(mesh.boundary_edges().len(), 4 * n);
        assert_eq!(mesh.boundary_vertices().len(), 4 * n);
        assert!(mesh.cell_points(CellId(n * n)).is_err());
        let area: f64 = (0..mesh.cells().len())
            .map(|i| {
                let map = BilinearMap::new(mesh.cell_points(CellId(i)).unwrap()).unwrap();
                square_degree_5()
                    .iter()
                    .map(|q| map.evaluate(q.point).unwrap().determinant * q.weight)
                    .sum::<f64>()
            })
            .sum();
        near(area, 1., 1e-13);
    }
    assert!(QuadMesh::unit_square(0).is_err());
    assert!(QuadMesh::unit_square(usize::MAX).is_err());
    for points in [
        [p(0., 0.), p(0., 1.), p(1., 1.), p(1., 0.)], // clockwise
        [p(0., 0.), p(1., 1.), p(1., 0.), p(0., 1.)], // bow tie
        [p(0., 0.), p(1., 0.), p(0.2, 0.2), p(0., 1.)], // concave
        [p(0., 0.), p(1., 0.), p(1., 1e-16), p(0., 1e-16)], // nearly degenerate
    ] {
        assert!(BilinearMap::new(points).is_err());
    }
    let v = square().to_vec();
    assert!(QuadMesh::new(
        v.clone(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(4)]]
    )
    .is_err());
    assert!(QuadMesh::new(
        v.clone(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(2)]]
    )
    .is_err());
    assert!(QuadMesh::new(
        v.clone(),
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(1), VertexId(2), VertexId(3), VertexId(0)]
        ]
    )
    .is_err());
    let mut v = v;
    v.push(p(2., 2.));
    assert!(QuadMesh::new(
        v,
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]]
    )
    .is_err());
}

#[test]
fn non_manifold_and_inconsistent_shared_edges_are_rejected() {
    let v = vec![
        p(0., 0.),
        p(1., 0.),
        p(1., 1.),
        p(0., 1.),
        p(0.1, 0.1),
        p(0.9, 0.1),
        p(0.8, 1.),
        p(0.2, 1.),
    ];
    // Two overlapping CCW cells sharing edge 0->1 in the same direction.
    assert!(QuadMesh::new(
        v,
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(0), VertexId(1), VertexId(5), VertexId(4)]
        ]
    )
    .is_err());
    let v = vec![
        p(0., 0.),
        p(1., 0.),
        p(1., 1.),
        p(0., 1.),
        p(1., -1.),
        p(0., -1.),
        p(0.9, 0.5),
        p(0.1, 0.5),
    ];
    assert!(QuadMesh::new(
        v,
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(1), VertexId(0), VertexId(5), VertexId(4)],
            [VertexId(0), VertexId(1), VertexId(6), VertexId(7)]
        ]
    )
    .is_err());
}

#[test]
fn reference_square_local_matrix_and_rhs_match_analytical_values() {
    let mesh = QuadMesh::new(
        square().to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let (a, b) = quad::assemble_poisson(&mesh, |_| 2.).unwrap();
    let expected = [
        [4., -1., -2., -1.],
        [-1., 4., -1., -2.],
        [-2., -1., 4., -1.],
        [-1., -2., -1., 4.],
    ];
    for i in 0..4 {
        near(b[i], 0.5, 1e-14);
        for j in 0..4 {
            near(a.get(i, j).unwrap(), expected[i][j] / 6., 1e-14);
        }
    }
    for value in a.mul(&[1., 1., 1., 1.]).unwrap() {
        near(value, 0., 1e-14);
    }
}

fn solve(
    mesh: &QuadMesh,
    source: impl Fn(Point2) -> f64,
    boundary: impl Fn(Point2) -> f64,
) -> Vec<f64> {
    let (a, b) = quad::assemble_poisson(mesh, source).unwrap();
    let constraints: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, boundary(mesh.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    conjugate_gradient(&a, &b, options).unwrap().0
}
#[test]
fn distorted_mesh_reproduces_affine_patch_with_nonzero_dirichlet() {
    let base = QuadMesh::unit_square(4).unwrap();
    let vertices: Vec<_> = base
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
        })
        .collect();
    let mesh = QuadMesh::new(vertices, base.cells().to_vec()).unwrap();
    let exact = |p: Point2| 1. + 2. * p.x - 3. * p.y;
    let values = solve(&mesh, |_| 0., exact);
    for (&u, &pt) in values.iter().zip(mesh.vertices()) {
        near(u, exact(pt), 1e-10);
    }
    let (l2, h1) = quad::errors(&mesh, &values, exact, |_| [2., -3.]).unwrap();
    assert!(l2 < 1e-10 && h1 < 1e-9);
}

#[test]
fn manufactured_solution_has_expected_convergence() {
    let exact = |p: Point2| (PI * p.x).sin() * (PI * p.y).sin();
    let grad = |p: Point2| {
        [
            PI * (PI * p.x).cos() * (PI * p.y).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).cos(),
        ]
    };
    let mut results = Vec::new();
    for n in [8, 16, 32] {
        let mesh = QuadMesh::unit_square(n).unwrap();
        let values = solve(&mesh, |p| 2. * PI * PI * exact(p), |_| 0.);
        results.push(quad::errors(&mesh, &values, exact, grad).unwrap());
    }
    for pair in results.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.8..2.2).contains(&l2), "L2 rate {}", l2);
        assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
    }
}

#[test]
fn invalid_numerics_and_vtk_are_reported() {
    let mesh = QuadMesh::unit_square(1).unwrap();
    assert!(quad::assemble_poisson(&mesh, |_| f64::NAN).is_err());
    assert!(quad::errors(&mesh, &[0.], |_| 0., |_| [0., 0.]).is_err());
    assert!(quad::errors(&mesh, &[f64::NAN; 4], |_| 0., |_| [0., 0.]).is_err());
    assert!(quad::errors(&mesh, &[0.; 4], |_| f64::INFINITY, |_| [0., 0.]).is_err());
    assert!(quad::write_vtk(Vec::new(), &mesh, &[0.]).is_err());
    let mut output = Vec::new();
    quad::write_vtk(&mut output, &mesh, &[0., 1., 2., 3.]).unwrap();
    let text = String::from_utf8(output).unwrap();
    assert!(text.contains("CELLS 1 5\n4 0 1 3 2\n"));
    assert!(text.contains("CELL_TYPES 1\n9\n"));
    assert!(text.contains("POINT_DATA 4"));
}

#[test]
fn manufactured_solution_converges_on_smooth_distorted_quadrilaterals() {
    let exact = |p: Point2| (PI * p.x).sin() * (PI * p.y).sin();
    let grad = |p: Point2| {
        [
            PI * (PI * p.x).cos() * (PI * p.y).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).cos(),
        ]
    };
    let mut results = Vec::new();
    for n in [8, 16, 32] {
        let base = QuadMesh::unit_square(n).unwrap();
        let vertices: Vec<_> = base
            .vertices()
            .iter()
            .map(|p| Point2 {
                x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
                y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
            })
            .collect();
        let mesh = QuadMesh::new(vertices, base.cells().to_vec()).unwrap();
        let values = solve(&mesh, |p| 2. * PI * PI * exact(p), |_| 0.);
        results.push(quad::errors(&mesh, &values, exact, grad).unwrap());
    }
    for pair in results.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.8..2.2).contains(&l2), "L2 rate {}", l2);
        assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
    }
}

#[test]
fn physical_gradients_satisfy_chain_rule_on_skew_nonaffine_cell() {
    let map = BilinearMap::new([p(1., 2.), p(4., 3.), p(4.5, 6.), p(0., 5.)]).unwrap();
    for q in [p(0.2, 0.3), p(0.75, 0.7)] {
        let result = map.evaluate(q).unwrap();
        let step = 1e-6;
        let plus_x = map.evaluate(p(q.x + step, q.y)).unwrap().point;
        let minus_x = map.evaluate(p(q.x - step, q.y)).unwrap().point;
        let plus_y = map.evaluate(p(q.x, q.y + step)).unwrap().point;
        let minus_y = map.evaluate(p(q.x, q.y - step)).unwrap().point;
        let dxi = [
            (plus_x.x - minus_x.x) / (2. * step),
            (plus_x.y - minus_x.y) / (2. * step),
        ];
        let deta = [
            (plus_y.x - minus_y.x) / (2. * step),
            (plus_y.y - minus_y.y) / (2. * step),
        ];
        let grad = quad::reference_gradients(q);
        for i in 0..4 {
            near(
                result.gradients[i][0] * dxi[0] + result.gradients[i][1] * dxi[1],
                grad[i][0],
                1e-9,
            );
            near(
                result.gradients[i][0] * deta[0] + result.gradients[i][1] * deta[1],
                grad[i][1],
                1e-9,
            );
        }
    }
}
