use dealii_rust::{
    base::{Point2, Point3},
    grid::{Mesh2, UniformRefinement as TriRefinement},
    lac::{
        conjugate_gradient, conjugate_gradient_jacobi_with_initial_guess as warm_pcg,
        conjugate_gradient_with_initial_guess as warm_cg, CgOptions, CooBuilder, CsrMatrix,
        JacobiPreconditioner,
    },
    numerics,
    quad::{self, QuadMesh, UniformRefinement as QuadRefinement},
    tetra::{self, TetMesh, UniformRefinement as TetRefinement},
};
fn diagonal(values: &[f64]) -> CsrMatrix {
    let mut coo = CooBuilder::new(values.len(), values.len());
    for (i, &value) in values.iter().enumerate() {
        coo.add(i, i, value).unwrap();
    }
    coo.into_csr().unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn strict() -> CgOptions {
    CgOptions {
        absolute_tolerance: 1e-13,
        relative_tolerance: 1e-12,
        ..CgOptions::default()
    }
}
#[test]
fn exact_initial_guess_returns_without_iterations_or_mutating_input() {
    let a = diagonal(&[2., 3., 4.]);
    let rhs = [2., 6., 12.];
    let initial = [1., 2., 3.];
    let (u, report) = warm_cg(
        &a,
        &rhs,
        &initial,
        CgOptions {
            max_iterations: 0,
            ..strict()
        },
    )
    .unwrap();
    assert_eq!(u, initial);
    assert_eq!(initial, [1., 2., 3.]);
    assert_eq!(report.iterations, 0);
    assert_eq!(report.residual_norm, 0.);
}
#[test]
fn nonzero_guess_zero_rhs_must_solve_not_return_early() {
    let a = diagonal(&[2., 3., 4.]);
    let initial = [1., 2., 3.];
    assert!(warm_cg(
        &a,
        &[0.; 3],
        &initial,
        CgOptions {
            max_iterations: 0,
            ..strict()
        }
    )
    .is_err());
    let (u, r) = warm_cg(&a, &[0.; 3], &initial, strict()).unwrap();
    assert!(r.iterations > 0);
    assert!(r.residual_norm <= strict().absolute_tolerance);
    for v in u {
        near(v, 0., 1e-12);
    }
}
#[test]
fn relative_threshold_is_based_on_rhs_not_initial_residual() {
    let a = diagonal(&[1.]);
    let opts = CgOptions {
        max_iterations: 0,
        relative_tolerance: 0.01,
        absolute_tolerance: 0.,
    };
    // Residual 0.005 is acceptable relative to b=1 despite nonzero x0.
    let (_, r) = warm_cg(&a, &[1.], &[0.995], opts).unwrap();
    assert_eq!(r.iterations, 0);
    near(r.residual_norm, 0.005, 1e-14);
    // With zero RHS, relative tolerance alone is zero even for a small residual.
    assert!(warm_cg(&a, &[0.], &[0.005], opts).is_err());
    // A bad guess cannot relax the threshold by inflating its initial residual.
    assert!(warm_cg(&a, &[1.], &[10.], opts).is_err());
}
#[test]
fn zero_guess_matches_original_api() {
    let a = diagonal(&[1., 2., 3., 4.]);
    let b = [1., 1., 1., 1.];
    let (u, r) = conjugate_gradient(&a, &b, strict()).unwrap();
    let (v, s) = warm_cg(&a, &b, &[0.; 4], strict()).unwrap();
    assert_eq!(u, v);
    assert_eq!(r.iterations, s.iterations);
    assert_eq!(r.residual_norm, s.residual_norm);
}
#[test]
fn warm_start_can_reduce_iterations_but_not_a_general_guarantee() {
    let a = diagonal(&[1., 2., 3., 4.]);
    let b = [1., 2., 3., 4.];
    let (_, cold) = conjugate_gradient(&a, &b, strict()).unwrap();
    let (u, warm) = warm_cg(&a, &b, &[1., 1., 1., 0.], strict()).unwrap();
    assert_eq!(warm.iterations, 1);
    assert!(cold.iterations > warm.iterations);
    for v in u {
        near(v, 1., 1e-13);
    }
}
#[test]
fn rejects_invalid_inputs_nonconvergence_and_numerical_breakdown() {
    let a = diagonal(&[1., 2.]);
    assert!(warm_cg(&a, &[1., 1.], &[0.], strict()).is_err());
    assert!(warm_cg(&a, &[1.], &[0.], strict()).is_err());
    assert!(warm_cg(&a, &[1., 1.], &[0., f64::NAN], strict()).is_err());
    assert!(warm_cg(&a, &[1., 1.], &[0., f64::INFINITY], strict()).is_err());
    assert!(warm_cg(&a, &[f64::NAN, 1.], &[0., 0.], strict()).is_err());
    assert!(warm_cg(
        &a,
        &[1., 1.],
        &[0., 0.],
        CgOptions {
            relative_tolerance: -1.,
            ..strict()
        }
    )
    .is_err());
    assert!(warm_cg(
        &a,
        &[1., 1.],
        &[0., 0.],
        CgOptions {
            absolute_tolerance: f64::INFINITY,
            ..strict()
        }
    )
    .is_err());
    assert!(warm_cg(
        &a,
        &[1., 1.],
        &[0., 0.],
        CgOptions {
            max_iterations: 0,
            ..strict()
        }
    )
    .is_err());
    assert!(warm_cg(&diagonal(&[-1.]), &[1.], &[0.], strict()).is_err());
    assert!(warm_cg(&diagonal(&[f64::MAX]), &[0.], &[2.], strict()).is_err());
    let exact = CgOptions {
        absolute_tolerance: 0.,
        relative_tolerance: 0.,
        ..strict()
    };
    assert!(warm_cg(&diagonal(&[1.]), &[1e-200], &[0.], exact).is_err());
}
#[test]
fn residual_is_true_residual_including_zero_iteration_case() {
    let a = diagonal(&[1., 2.]);
    let b = [1., 2.];
    let options = CgOptions {
        absolute_tolerance: 0.5,
        relative_tolerance: 0.,
        ..strict()
    };
    let (u, r) = warm_cg(&a, &b, &[1., 0.8], options).unwrap();
    assert_eq!(r.iterations, 0);
    let ax = a.mul(&u).unwrap();
    let actual = (b[0] - ax[0]).hypot(b[1] - ax[1]);
    near(r.residual_norm, actual, 1e-15);
    let (u, r) = warm_cg(&a, &b, &[1., 0.8], strict()).unwrap();
    let ax = a.mul(&u).unwrap();
    let actual = (b[0] - ax[0]).hypot(b[1] - ax[1]);
    near(r.residual_norm, actual, 1e-15);
}
fn exact2(p: Point2) -> f64 {
    1. + p.x * p.x + p.y * p.y
}
fn exact3(p: Point3) -> f64 {
    1. + p.x * p.x + p.y * p.y + p.z * p.z
}
fn compare(a: &CsrMatrix, b: &[f64], constraints: &[(usize, f64)], mut initial: Vec<f64>) {
    // New boundary values differ from an interpolated nonlinear coarse trace.
    assert!(constraints
        .iter()
        .any(|&(i, v)| (initial[i] - v).abs() > 1e-7));
    for &(i, value) in constraints {
        initial[i] = value;
    }
    let (warm, wr) = warm_cg(a, b, &initial, strict()).unwrap();
    let jacobi = JacobiPreconditioner::new(a).unwrap();
    let (preconditioned, pr) = warm_pcg(a, b, &initial, &jacobi, strict()).unwrap();
    assert!(pr.residual_norm < 1e-8);
    for (x, y) in preconditioned.iter().zip(&warm) {
        near(*x, *y, 1e-8);
    }
    let (cold, cr) = conjugate_gradient(a, b, strict()).unwrap();
    assert!(wr.residual_norm < 1e-8 && cr.residual_norm < 1e-8);
    for (x, y) in warm.iter().zip(&cold) {
        near(*x, *y, 1e-8);
    }
    for &(i, value) in constraints {
        near(warm[i], value, 1e-9);
    }
}
#[test]
fn triangle_refinement_warm_start_matches_cold_with_nonzero_dirichlet() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let (a, b) = numerics::assemble_poisson(&mesh, |_| -4.).unwrap();
    let d: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact2(mesh.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, strict()).unwrap();
    let r = TriRefinement::new(&mesh).unwrap();
    let fine = r.mesh();
    let x = r.prolongate_p1(&u).unwrap();
    let (a, b) = numerics::assemble_poisson(fine, |_| -4.).unwrap();
    let d: Vec<_> = fine
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact2(fine.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    compare(&a, &b, &d, x);
}
#[test]
fn quad_refinement_warm_start_matches_cold_with_nonzero_dirichlet() {
    let mesh = QuadMesh::unit_square(3).unwrap();
    let (a, b) = quad::assemble_poisson(&mesh, |_| -4.).unwrap();
    let d: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact2(mesh.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, strict()).unwrap();
    let r = QuadRefinement::new(&mesh).unwrap();
    let fine = r.mesh();
    let x = r.prolongate_q1(&u).unwrap();
    let (a, b) = quad::assemble_poisson(fine, |_| -4.).unwrap();
    let d: Vec<_> = fine
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact2(fine.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    compare(&a, &b, &d, x);
}
#[test]
fn tetra_refinement_warm_start_matches_cold_with_nonzero_dirichlet() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let (a, b) = tetra::assemble_poisson(&mesh, |_| -6.).unwrap();
    let d: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact3(mesh.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, strict()).unwrap();
    let r = TetRefinement::new(&mesh).unwrap();
    let fine = r.mesh();
    let x = r.prolongate_p1(&u).unwrap();
    let (a, b) = tetra::assemble_poisson(fine, |_| -6.).unwrap();
    let d: Vec<_> = fine
        .boundary_vertices()
        .iter()
        .map(|v| (v.0, exact3(fine.vertices()[v.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    compare(&a, &b, &d, x);
}
