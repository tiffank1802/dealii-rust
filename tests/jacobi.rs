use dealii_rust::{
    elliptic,
    grid::Mesh2,
    integrators::DiffusionTensor,
    lac::{
        conjugate_gradient, conjugate_gradient_jacobi as pcg,
        conjugate_gradient_jacobi_with_initial_guess as warm_pcg, CgOptions, CooBuilder, CsrMatrix,
        JacobiPreconditioner,
    },
};
fn matrix(n: usize, entries: &[(usize, usize, f64)]) -> CsrMatrix {
    let mut a = CooBuilder::new(n, n);
    for &(i, j, x) in entries {
        a.add(i, j, x).unwrap();
    }
    a.into_csr().unwrap()
}
fn diagonal(d: &[f64]) -> CsrMatrix {
    let entries: Vec<_> = d.iter().enumerate().map(|(i, &x)| (i, i, x)).collect();
    matrix(d.len(), &entries)
}
fn opts() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-11,
        absolute_tolerance: 1e-13,
        max_iterations: 10000,
    }
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t, "{} != {}", a, b);
}
fn true_residual(a: &CsrMatrix, b: &[f64], x: &[f64]) -> f64 {
    let ax = a.mul(x).unwrap();
    b.iter().zip(ax).fold(0.0_f64, |r, (b, ax)| r.hypot(b - ax))
}
#[test]
fn builds_and_applies_jacobi_from_csr_diagonal() {
    let a = matrix(
        3,
        &[
            (0, 0, 2.),
            (0, 1, -1.),
            (1, 0, -1.),
            (1, 1, 4.),
            (2, 2, 0.5),
        ],
    );
    let j = JacobiPreconditioner::new(&a).unwrap();
    assert_eq!(j.apply(&[4., -2., 3.]).unwrap(), vec![2., -0.5, 6.]);
    assert!(j.apply(&[]).is_err());
    assert!(j.apply(&[1., f64::NAN, 1.]).is_err());
    assert!(j.apply(&[1., f64::INFINITY, 1.]).is_err());
    let tiny = f64::from_bits(1);
    let j = JacobiPreconditioner::new(&diagonal(&[tiny])).unwrap();
    assert_eq!(j.apply(&[tiny]).unwrap(), vec![1.]); // no explicit overflowing reciprocal
    assert!(j.apply(&[1.]).is_err());
}
#[test]
fn rejects_missing_zero_negative_diagonal_or_nonsquare() {
    assert!(JacobiPreconditioner::new(&CooBuilder::new(2, 3).into_csr().unwrap()).is_err());
    assert!(JacobiPreconditioner::new(&matrix(2, &[(0, 0, 1.)])).is_err());
    assert!(JacobiPreconditioner::new(&diagonal(&[1., 0.])).is_err());
    assert!(JacobiPreconditioner::new(&diagonal(&[1., -2.])).is_err());
}
#[test]
fn diagonal_system_converges_in_one_jacobi_step() {
    let a = diagonal(&[1., 2., 5., 10., 100.]);
    let b = [1., 2., 5., 10., 100.];
    let j = JacobiPreconditioner::new(&a).unwrap();
    let (u, r) = pcg(&a, &b, &j, opts()).unwrap();
    let (_, cold) = conjugate_gradient(&a, &b, opts()).unwrap();
    assert_eq!(r.iterations, 1);
    assert!(cold.iterations > r.iterations);
    for v in u {
        near(v, 1., 1e-13);
    }
}
#[test]
fn coupled_spd_solution_and_initial_guess_agree_with_cg() {
    let a = matrix(
        3,
        &[
            (0, 0, 4.),
            (0, 1, -1.),
            (1, 0, -1.),
            (1, 1, 3.),
            (1, 2, -0.5),
            (2, 1, -0.5),
            (2, 2, 2.),
        ],
    );
    let solution = [1., -2., 3.];
    let b = a.mul(&solution).unwrap();
    let j = JacobiPreconditioner::new(&a).unwrap();
    let guess = [9., -7., 2.];
    let (u, r) = warm_pcg(&a, &b, &guess, &j, opts()).unwrap();
    let (v, _) = conjugate_gradient(&a, &b, opts()).unwrap();
    assert_eq!(guess, [9., -7., 2.]);
    for i in 0..3 {
        near(u[i], solution[i], 1e-11);
        near(u[i], v[i], 1e-11);
    }
    near(r.residual_norm, true_residual(&a, &b, &u), 1e-15);
}
#[test]
fn identity_jacobi_matches_unpreconditioned_cg() {
    let a = matrix(
        3,
        &[
            (0, 0, 1.),
            (0, 1, -0.2),
            (1, 0, -0.2),
            (1, 1, 1.),
            (1, 2, -0.1),
            (2, 1, -0.1),
            (2, 2, 1.),
        ],
    );
    let b = [1., 2., 3.];
    let j = JacobiPreconditioner::new(&a).unwrap();
    let (u, r) = pcg(&a, &b, &j, opts()).unwrap();
    let (v, s) = conjugate_gradient(&a, &b, opts()).unwrap();
    assert_eq!(u, v);
    assert_eq!(r.iterations, s.iterations);
    assert_eq!(r.residual_norm, s.residual_norm);
}
#[test]
fn residual_and_stopping_test_use_unpreconditioned_norm() {
    let a = diagonal(&[1e8, 1.]);
    let b = [1., 0.];
    let j = JacobiPreconditioner::new(&a).unwrap();
    let o = CgOptions {
        absolute_tolerance: 0.1,
        relative_tolerance: 0.,
        ..opts()
    };
    // Preconditioned residual sqrt(rᵀM⁻¹r)=1e-4, but the actual residual is 1.
    assert!(pcg(
        &a,
        &b,
        &j,
        CgOptions {
            max_iterations: 0,
            ..o
        }
    )
    .is_err());
    let (u, r) = pcg(&a, &b, &j, o).unwrap();
    assert_eq!(r.iterations, 1);
    near(r.residual_norm, true_residual(&a, &b, &u), 1e-15);
    near(u[0], 1e-8, 1e-20);
}
#[test]
fn warm_start_uses_rhs_relative_tolerance_and_checks_zero_rhs() {
    let a = diagonal(&[2., 3.]);
    let j = JacobiPreconditioner::new(&a).unwrap();
    let o = CgOptions {
        max_iterations: 0,
        relative_tolerance: 0.01,
        absolute_tolerance: 0.,
    };
    let (x, r) = warm_pcg(&a, &[2., 3.], &[1., 1.], &j, o).unwrap();
    assert_eq!(x, vec![1., 1.]);
    assert_eq!(r.iterations, 0);
    assert!(warm_pcg(&a, &[0., 0.], &[0.01, 0.], &j, o).is_err());
    let (u, r) = warm_pcg(&a, &[0., 0.], &[0.01, 0.], &j, opts()).unwrap();
    assert!(r.iterations > 0);
    assert!(r.residual_norm < opts().absolute_tolerance);
    near(u[0], 0., 1e-14);
}
#[test]
fn dimensions_inputs_breakdown_and_nonconvergence() {
    let a = diagonal(&[1., 2.]);
    let j = JacobiPreconditioner::new(&a).unwrap();
    let wrong = JacobiPreconditioner::new(&diagonal(&[1.])).unwrap();
    assert!(pcg(&a, &[0., 0.], &wrong, opts()).is_err()); // dimension checked even on zero RHS
    assert!(warm_pcg(&a, &[1., 2.], &[0.], &j, opts()).is_err());
    assert!(warm_pcg(&a, &[1., 2.], &[0., f64::NAN], &j, opts()).is_err());
    assert!(pcg(&a, &[1., f64::INFINITY], &j, opts()).is_err());
    assert!(pcg(
        &a,
        &[1., 2.],
        &j,
        CgOptions {
            max_iterations: 0,
            ..opts()
        }
    )
    .is_err());
    assert!(pcg(
        &a,
        &[1., 2.],
        &j,
        CgOptions {
            relative_tolerance: -1.,
            ..opts()
        }
    )
    .is_err());
    let indefinite = matrix(2, &[(0, 0, 1.), (0, 1, 2.), (1, 0, 2.), (1, 1, 1.)]);
    let j = JacobiPreconditioner::new(&indefinite).unwrap();
    assert!(pcg(&indefinite, &[1., -1.], &j, opts())
        .unwrap_err()
        .to_string()
        .contains("breakdown"));
    let a = diagonal(&[1e-300]);
    let j = JacobiPreconditioner::new(&a).unwrap();
    assert!(pcg(&a, &[1e10], &j, opts()).is_err());
}
#[test]
fn diagonal_is_snapshot_and_should_be_rebuilt_after_constraints() {
    let a = matrix(2, &[(0, 0, 4.), (0, 1, -1.), (1, 0, -1.), (1, 1, 2.)]);
    let before = JacobiPreconditioner::new(&a).unwrap();
    let (constrained, _) = a.with_dirichlet(&[0., 0.], &[(0, 2.)]).unwrap();
    let after = JacobiPreconditioner::new(&constrained).unwrap();
    assert_eq!(before.apply(&[4., 2.]).unwrap(), vec![1., 1.]);
    assert_eq!(after.apply(&[4., 2.]).unwrap(), vec![4., 1.]);
}
#[test]
fn fem_discontinuous_material_pcg_agrees_with_plain_cg() {
    let mesh = Mesh2::unit_square(16).unwrap();
    let (a, b) = elliptic::assemble(
        &mesh,
        |_, p| {
            let k = if p.x < 0.5 { 1000. } else { 1. };
            DiffusionTensor {
                xx: k,
                xy: 0.,
                yy: k,
            }
        },
        |_, _| 0.,
        |_, _| 1.,
    )
    .unwrap();
    let d: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    let j = JacobiPreconditioner::new(&a).unwrap();
    let (u, r) = pcg(&a, &b, &j, opts()).unwrap();
    let (v, s) = conjugate_gradient(&a, &b, opts()).unwrap();
    for (a, b) in u.iter().zip(v) {
        near(*a, b, 1e-8);
    }
    near(r.residual_norm, true_residual(&a, &b, &u), 1e-15);
    assert!(
        r.residual_norm
            <= (opts().relative_tolerance * b.iter().fold(0.0_f64, |v, x| v.hypot(*x)))
                .max(opts().absolute_tolerance)
    );
    assert!(s.residual_norm < 1e-9);
}
