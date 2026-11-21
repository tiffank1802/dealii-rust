#![allow(clippy::needless_range_loop)] // Analytical local arrays use indices.
use dealii_rust::{
    base::Point3,
    grid::VertexId,
    lac::{conjugate_gradient, conjugate_gradient_jacobi, CgOptions, JacobiPreconditioner},
    tetra::{self, BoundaryFaceId, TetMesh, UniformRefinement},
    tetra_field,
};
fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t, "{} != {} (tol {})", a, b, t);
}
fn reference() -> TetMesh {
    TetMesh::new(
        vec![p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}
fn opts() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    }
}
fn affine(p: Point3) -> f64 {
    1. + 2. * p.x - 3. * p.y + 4. * p.z
}
#[test]
fn reference_mass_is_analytical_symmetric_and_positive() {
    let m = tetra_field::mass_matrix(&reference()).unwrap();
    let mut sum = 0.;
    for i in 0..4 {
        for j in 0..4 {
            let x = m.get(i, j).unwrap();
            near(x, if i == j { 1. / 60. } else { 1. / 120. }, 1e-15);
            assert_eq!(x, m.get(j, i).unwrap());
            sum += x;
        }
    }
    near(sum, 1. / 6., 1e-15);
    let v = [1., -2., 3., -4.];
    assert!(
        v.iter()
            .zip(m.mul(&v).unwrap())
            .map(|(a, b)| a * b)
            .sum::<f64>()
            > 0.
    );
}
#[test]
fn global_mass_sums_to_cube_volume_and_uses_shared_vertices() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let m = tetra_field::mass_matrix(&mesh).unwrap();
    let ones = vec![1.; mesh.vertices().len()];
    near(m.mul(&ones).unwrap().iter().sum(), 1., 1e-13);
    for i in 0..ones.len() {
        for j in 0..ones.len() {
            assert_eq!(m.get(i, j).unwrap(), m.get(j, i).unwrap());
        }
    }
}
#[test]
fn degree_four_source_load_is_exact() {
    let (_, b) = tetra_field::projection_system(&reference(), |p| p.x.powi(4)).unwrap();
    for (x, y) in b
        .iter()
        .zip([1. / 1680., 1. / 336., 1. / 1680., 1. / 1680.])
    {
        near(*x, y, 1e-15);
    }
}
#[test]
fn quadratic_projection_has_analytic_solution_and_orthogonality() {
    let mesh = reference();
    let f = |p: Point3| p.x * p.x;
    let (m, b) = tetra_field::projection_system(&mesh, f).unwrap();
    let (u, report) = conjugate_gradient(&m, &b, opts()).unwrap();
    assert!(report.residual_norm < 1e-13);
    for (a, b) in u.iter().zip([-1. / 15., 3. / 5., -1. / 15., -1. / 15.]) {
        near(*a, b, 1e-13);
    }
    for (a, b) in m.mul(&u).unwrap().iter().zip(b) {
        near(*a, b, 1e-14);
    }
    let interp = tetra_field::interpolate(&mesh, f).unwrap();
    let e = tetra::errors(&mesh, &u, f, |p| [2. * p.x, 0., 0.])
        .unwrap()
        .0;
    let ei = tetra::errors(&mesh, &interp, f, |p| [2. * p.x, 0., 0.])
        .unwrap()
        .0;
    assert!(e < ei);
}
#[test]
fn affine_functions_are_reproduced_on_cube_and_skew_mesh() {
    let skew = TetMesh::new(
        vec![p(1., 2., 3.), p(3., 2., 3.), p(2., 5., 3.), p(1., 3., 7.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    for mesh in [reference(), skew, TetMesh::unit_cube(3).unwrap()] {
        for f in [|_: Point3| 1., affine] {
            let nodal = tetra_field::interpolate(&mesh, f).unwrap();
            let (u, report) = tetra_field::project_l2(&mesh, f, opts()).unwrap();
            assert!(report.residual_norm < 1e-11);
            for (a, b) in u.iter().zip(nodal) {
                near(*a, b, 1e-9);
            }
        }
    }
}
#[test]
fn vertex_permutations_do_not_change_global_system() {
    let m = reference();
    let rotated = TetMesh::new(
        m.vertices().to_vec(),
        vec![[VertexId(1), VertexId(2), VertexId(0), VertexId(3)]],
    )
    .unwrap();
    let (a, x) = tetra_field::projection_system(&m, |p| p.x.powi(3) + p.y).unwrap();
    let (b, y) = tetra_field::projection_system(&rotated, |p| p.x.powi(3) + p.y).unwrap();
    for i in 0..4 {
        near(x[i], y[i], 1e-14);
        for j in 0..4 {
            near(a.get(i, j).unwrap(), b.get(i, j).unwrap(), 1e-15);
        }
    }
}
#[test]
fn constrained_projection_and_jacobi_reproduce_affine() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let faces: Vec<_> = (0..mesh.boundary_faces().len())
        .map(BoundaryFaceId)
        .collect();
    let constraints = tetra::interpolate_boundary(&mesh, &faces, affine).unwrap();
    let (m, b) = tetra_field::projection_system(&mesh, affine).unwrap();
    let (m, b) = m.with_dirichlet(&b, &constraints).unwrap();
    let j = JacobiPreconditioner::new(&m).unwrap();
    let (u, _) = conjugate_gradient_jacobi(&m, &b, &j, opts()).unwrap();
    for (a, p) in u.iter().zip(mesh.vertices()) {
        near(*a, affine(*p), 1e-10);
    }
}
#[test]
fn disconnected_components_and_prolongation() {
    let mesh = TetMesh::new(
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
    .unwrap();
    let (u, _) = tetra_field::project_l2(&mesh, affine, opts()).unwrap();
    let r = UniformRefinement::new(&mesh).unwrap();
    let fine = r.prolongate_p1(&u).unwrap();
    for (a, p) in fine.iter().zip(r.mesh().vertices()) {
        near(*a, affine(*p), 1e-10);
    }
}
#[test]
fn quadratic_projection_converges_under_refinement() {
    let f = |p: Point3| p.x * p.x + p.y * p.y + p.z * p.z;
    let grad = |p: Point3| [2. * p.x, 2. * p.y, 2. * p.z];
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let mesh = TetMesh::unit_cube(n).unwrap();
        let (u, _) = tetra_field::project_l2(&mesh, f, opts()).unwrap();
        let interp = tetra_field::interpolate(&mesh, f).unwrap();
        let e = tetra::errors(&mesh, &u, f, grad).unwrap().0;
        let i = tetra::errors(&mesh, &interp, f, grad).unwrap().0;
        assert!(e < i);
        errors.push(e);
    }
    assert!(errors[0] / errors[1] > 3.5);
    assert!(errors[1] / errors[2] > 3.5);
}
#[test]
fn invalid_values_overflow_and_solver_failure_are_rejected() {
    let mesh = reference();
    assert!(tetra_field::interpolate(&mesh, |_| f64::NAN).is_err());
    assert!(tetra_field::projection_system(&mesh, |_| f64::INFINITY).is_err());
    assert!(tetra_field::projection_system(&mesh, |_| f64::MAX).is_ok()); // bounded domain, finite integrals
    let large = TetMesh::new(
        vec![
            p(0., 0., 0.),
            p(100., 0., 0.),
            p(0., 100., 0.),
            p(0., 0., 100.),
        ],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(tetra_field::projection_system(&large, |_| f64::MAX).is_err());
    let options = CgOptions {
        max_iterations: 0,
        ..opts()
    };
    assert!(tetra_field::project_l2(&mesh, affine, options).is_err());
}
