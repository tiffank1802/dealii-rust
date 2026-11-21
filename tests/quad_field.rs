#![allow(clippy::needless_range_loop)] // Explicit analytical matrices.
use dealii_rust::{
    base::Point2,
    grid::{BoundaryEdgeId, VertexId},
    lac::{conjugate_gradient, conjugate_gradient_jacobi, CgOptions, JacobiPreconditioner},
    quad::{self, QuadMesh, UniformRefinement},
    quad_field,
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn opts() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-13,
        absolute_tolerance: 1e-14,
        max_iterations: 1000,
    }
}
fn square() -> QuadMesh {
    QuadMesh::new(
        vec![p(0., 0.), p(1., 0.), p(1., 1.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}
fn trapezoid() -> QuadMesh {
    QuadMesh::new(
        vec![p(0., 0.), p(2., 0.), p(1., 1.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}
fn distorted(n: usize) -> QuadMesh {
    let mesh = QuadMesh::unit_square(n).unwrap();
    let points = mesh
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
        })
        .collect();
    QuadMesh::new(points, mesh.cells().to_vec()).unwrap()
}
#[test]
fn reference_square_has_analytical_consistent_mass_matrix() {
    let m = quad_field::mass_matrix(&square()).unwrap();
    let expected = [
        [4., 2., 1., 2.],
        [2., 4., 2., 1.],
        [1., 2., 4., 2.],
        [2., 1., 2., 4.],
    ];
    for i in 0..4 {
        for j in 0..4 {
            near(m.get(i, j).unwrap(), expected[i][j] / 36., 1e-14);
            assert_eq!(m.get(i, j).unwrap(), m.get(j, i).unwrap());
        }
    }
}
#[test]
fn trapezoid_mass_is_analytically_integrated_and_positive() {
    let mesh = trapezoid();
    let m = quad_field::mass_matrix(&mesh).unwrap();
    // Map x=ξ(2−η), y=η, det J=2−η. Factor x and y integrals.
    let x = [[1. / 3., 1. / 6.], [1. / 6., 1. / 3.]];
    let y = [[7. / 12., 1. / 4.], [1. / 4., 5. / 12.]];
    let nodes = [(0, 0), (1, 0), (1, 1), (0, 1)];
    let mut total = 0.;
    for i in 0..4 {
        for j in 0..4 {
            let (xi, yi) = nodes[i];
            let (xj, yj) = nodes[j];
            let v = m.get(i, j).unwrap();
            near(v, x[xi][xj] * y[yi][yj], 1e-14);
            total += v;
        }
    }
    near(total, 1.5, 1e-13);
    let v = [1., -2., 3., -4.];
    let mv = m.mul(&v).unwrap();
    assert!(v.iter().zip(mv).map(|(a, b)| a * b).sum::<f64>() > 0.);
}
#[test]
fn cubic_source_integrates_exactly_on_nonaffine_trapezoid() {
    let (_, b) = quad_field::projection_system(&trapezoid(), |p| p.x.powi(3)).unwrap();
    // ∫ξ³ l_i(ξ) (2−η)^4 l_j(η) dξdη, independent analytic result.
    for (v, e) in b.iter().zip([43. / 200., 43. / 50., 19. / 50., 19. / 200.]) {
        near(*v, e, 2e-14);
    }
}
#[test]
fn quartic_source_is_exact_on_affine_square() {
    let (_, b) = quad_field::projection_system(&square(), |p| p.x.powi(4)).unwrap();
    for (v, e) in b.iter().zip([1. / 60., 1. / 12., 1. / 12., 1. / 60.]) {
        near(*v, e, 2e-14);
    }
}
#[test]
fn cyclic_vertex_order_does_not_change_global_mass_or_load() {
    let first = trapezoid();
    let rotated = QuadMesh::new(
        first.vertices().to_vec(),
        vec![[VertexId(1), VertexId(2), VertexId(3), VertexId(0)]],
    )
    .unwrap();
    let (a, x) = quad_field::projection_system(&first, |p| p.x.powi(3) + p.y).unwrap();
    let (b, y) = quad_field::projection_system(&rotated, |p| p.x.powi(3) + p.y).unwrap();
    for i in 0..4 {
        near(x[i], y[i], 1e-13);
        for j in 0..4 {
            near(a.get(i, j).unwrap(), b.get(i, j).unwrap(), 1e-13);
        }
    }
}
#[test]
fn interpolates_and_projects_affine_physical_fields() {
    for mesh in [square(), trapezoid(), distorted(5)] {
        for f in [|_: Point2| 1., |p: Point2| 3. + 2. * p.x - 5. * p.y] {
            let interpolated = quad_field::interpolate(&mesh, f).unwrap();
            let (u, r) = quad_field::project_l2(&mesh, f, opts()).unwrap();
            assert!(r.residual_norm < 1e-11);
            for (a, b) in u.iter().zip(interpolated) {
                near(*a, b, 1e-9);
            }
        }
    }
}
#[test]
fn nonaffine_q1_reference_coordinate_is_reproduced_and_transferred() {
    let mesh = trapezoid();
    let f = |p: Point2| p.x / (2. - p.y);
    let (u, _) = quad_field::project_l2(&mesh, f, opts()).unwrap();
    for (u, expected) in u.iter().zip([0., 1., 1., 0.]) {
        near(*u, expected, 1e-12);
    }
    let refined = UniformRefinement::new(&mesh).unwrap();
    let transferred = refined.prolongate_q1(&u).unwrap();
    let (fine, _) = quad_field::project_l2(refined.mesh(), f, opts()).unwrap();
    for (a, b) in fine.iter().zip(transferred) {
        near(*a, b, 1e-11);
    }
}
#[test]
fn quadratic_projection_square_has_analytic_coefficients_and_orthogonality() {
    let mesh = square();
    let f = |p: Point2| p.x * p.x;
    let (matrix, b) = quad_field::projection_system(&mesh, f).unwrap();
    let (u, r) = conjugate_gradient(&matrix, &b, opts()).unwrap();
    assert!(r.residual_norm < 1e-13);
    for (u, e) in u.iter().zip([-1. / 6., 5. / 6., 5. / 6., -1. / 6.]) {
        near(*u, e, 1e-13);
    }
    for (a, b) in matrix.mul(&u).unwrap().iter().zip(b) {
        near(*a, b, 1e-13);
    }
    let interp = quad_field::interpolate(&mesh, f).unwrap();
    let l2 = quad::errors(&mesh, &u, f, |p| [2. * p.x, 0.]).unwrap().0;
    let i2 = quad::errors(&mesh, &interp, f, |p| [2. * p.x, 0.])
        .unwrap()
        .0;
    near(l2, 1. / 180.0_f64.sqrt(), 1e-13);
    assert!(l2 < i2 * 0.5);
}
#[test]
fn constrained_projection_and_jacobi_agree() {
    let mesh = distorted(4);
    let f = |p: Point2| 1. + p.x + 2. * p.y;
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let constraints = quad::interpolate_boundary(&mesh, &edges, f).unwrap();
    let (m, b) = quad_field::projection_system(&mesh, f).unwrap();
    let (m, b) = m.with_dirichlet(&b, &constraints).unwrap();
    let j = JacobiPreconditioner::new(&m).unwrap();
    let (u, _) = conjugate_gradient_jacobi(&m, &b, &j, opts()).unwrap();
    for (u, p) in u.iter().zip(mesh.vertices()) {
        near(*u, f(*p), 1e-10);
    }
}
#[test]
fn disconnected_mesh_projects_each_component() {
    let mesh = QuadMesh::new(
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
    .unwrap();
    let f = |p: Point2| 1. + 2. * p.x - p.y;
    let (u, _) = quad_field::project_l2(&mesh, f, opts()).unwrap();
    for (u, p) in u.iter().zip(mesh.vertices()) {
        near(*u, f(*p), 1e-11);
    }
    let m = quad_field::mass_matrix(&mesh).unwrap();
    for i in 0..4 {
        for j in 4..8 {
            near(m.get(i, j).unwrap(), 0., 1e-14);
        }
    }
}
#[test]
fn quadratic_projection_converges_on_regular_and_deformed_meshes() {
    let f = |p: Point2| p.x * p.x + p.y * p.y;
    let grad = |p: Point2| [2. * p.x, 2. * p.y];
    for deform in [false, true] {
        let mut old = f64::INFINITY;
        for n in [4, 8, 16, 32] {
            let mesh = if deform {
                distorted(n)
            } else {
                QuadMesh::unit_square(n).unwrap()
            };
            let interpolation = quad_field::interpolate(&mesh, f).unwrap();
            let (projection, _) = quad_field::project_l2(&mesh, f, opts()).unwrap();
            let l2 = quad::errors(&mesh, &projection, f, grad).unwrap().0;
            let l2_i = quad::errors(&mesh, &interpolation, f, grad).unwrap().0;
            assert!(l2 < l2_i, "projection {} interpolation {}", l2, l2_i);
            assert!(l2 < old * 0.4, "current {} previous {}", l2, old);
            old = l2;
        }
    }
}
#[test]
fn invalid_data_unrepresentable_values_and_solver_failure_are_rejected() {
    let mesh = square();
    assert!(quad_field::interpolate(&mesh, |_| f64::NAN).is_err());
    assert!(quad_field::projection_system(&mesh, |_| f64::INFINITY).is_err());
    let large = QuadMesh::new(
        vec![p(0., 0.), p(10., 0.), p(10., 10.), p(0., 10.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(quad_field::projection_system(&large, |_| f64::MAX).is_err());
    assert!(quad_field::project_l2(
        &mesh,
        |_| 1.,
        CgOptions {
            max_iterations: 0,
            ..opts()
        }
    )
    .is_err());
    assert!(quad_field::project_l2(
        &mesh,
        |_| 1.,
        CgOptions {
            relative_tolerance: f64::NAN,
            ..opts()
        }
    )
    .is_err());
    let e = 4e-162;
    let tiny = QuadMesh::new(
        vec![p(0., 0.), p(e, 0.), p(e, e), p(0., e)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(quad_field::mass_matrix(&tiny).is_err());
}
