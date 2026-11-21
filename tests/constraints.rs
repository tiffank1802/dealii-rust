#![allow(clippy::needless_range_loop)] // Dense reference formulas are deliberately indexed.
use dealii_rust::{
    base::Point2,
    constraints::AffineConstraints,
    grid::Mesh2,
    lac::{conjugate_gradient, CgOptions, CooBuilder, CsrMatrix},
    numerics,
};
fn matrix(values: &[&[f64]]) -> CsrMatrix {
    let rows = values.len();
    let cols = values.first().map_or(0, |r| r.len());
    let mut b = CooBuilder::new(rows, cols);
    for (i, row) in values.iter().enumerate() {
        for (j, &x) in row.iter().enumerate() {
            b.add(i, j, x).unwrap();
        }
    }
    b.into_csr().unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{} != {}", a, b);
}
fn options() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    }
}
#[test]
fn identity_constraints_preserve_system_and_vector() {
    let a = matrix(&[&[2., -1.], &[-1., 2.]]);
    let c = AffineConstraints::new(2).close().unwrap();
    assert_eq!(c.size(), 2);
    assert_eq!(c.free_dofs(), &[0, 1]);
    let (b, r) = c.condense(&a, &[1., 3.]).unwrap();
    for i in 0..2 {
        near(r[i], [1., 3.][i]);
        for j in 0..2 {
            near(a.get(i, j).unwrap(), b.get(i, j).unwrap());
        }
    }
    assert_eq!(c.reconstruct(&[7., 8.]).unwrap(), vec![7., 8.]);
}
#[test]
fn analytical_midpoint_constraint_reduces_system_and_reconstructs_solution() {
    let a = matrix(&[&[2., -1., 0.], &[-1., 2., -1.], &[0., -1., 2.]]);
    let rhs = [1., 0., 3.];
    let mut c = AffineConstraints::new(3);
    c.add(1, 1., &[(0, 0.5), (2, 0.5)]).unwrap();
    let c = c.close().unwrap();
    assert_eq!(c.free_dofs(), &[0, 2]);
    let (ar, br) = c.condense(&a, &rhs).unwrap();
    for i in 0..2 {
        near(br[i], [1., 3.][i]);
        for j in 0..2 {
            near(ar.get(i, j).unwrap(), [[1.5, -0.5], [-0.5, 1.5]][i][j]);
        }
    }
    let (y, _) = conjugate_gradient(&ar, &br, options()).unwrap();
    let x = c.reconstruct(&y).unwrap();
    for (a, b) in x.iter().zip([1.5, 3., 2.5]) {
        near(*a, b);
    }
    let ax = a.mul(&x).unwrap();
    let residual: Vec<_> = rhs.iter().zip(ax).map(|(b, v)| b - v).collect();
    near(residual[0] + 0.5 * residual[1], 0.);
    near(residual[2] + 0.5 * residual[1], 0.);
}
#[test]
fn nested_rules_expand_constants_coefficients_and_preserve_index_order() {
    let mut c = AffineConstraints::new(5);
    c.add(2, 3., &[(1, 0.5)]).unwrap();
    c.add(1, 2., &[(0, 2.)]).unwrap();
    c.add(4, -1., &[(2, 1.), (1, -1.), (3, 3.)]).unwrap();
    let c = c.close().unwrap();
    assert_eq!(c.free_dofs(), &[0, 3]);
    assert_eq!(
        c.reconstruct(&[4., 5.]).unwrap(),
        vec![4., 10., 8., 5., 12.]
    );
}
#[test]
fn dense_reference_agrees_for_nested_nonsymmetric_system() {
    let dense = [
        [4., 1., -2., 0.5],
        [-3., 5., 1., 2.],
        [1., 2., 3., -4.],
        [0.25, -1., 2., 7.],
    ];
    let a = matrix(&[&dense[0], &dense[1], &dense[2], &dense[3]]);
    let rhs = [2., -1., 3., 4.];
    let mut b = AffineConstraints::new(4);
    b.add(0, 1., &[(1, 2.), (3, -1.)]).unwrap();
    b.add(2, -2., &[(0, 0.5)]).unwrap();
    let c = b.close().unwrap();
    let transform = [[2., -1.], [1., 0.], [1., -0.5], [0., 1.]];
    let g = [1., 0., -1.5, 0.];
    let (r, q) = c.condense(&a, &rhs).unwrap();
    for k in 0..2 {
        let mut expected_rhs = 0.;
        for i in 0..4 {
            let mut shift = rhs[i];
            for j in 0..4 {
                shift -= dense[i][j] * g[j];
            }
            expected_rhs += transform[i][k] * shift;
        }
        near(q[k], expected_rhs);
        for l in 0..2 {
            let mut expected = 0.;
            for i in 0..4 {
                for j in 0..4 {
                    expected += transform[i][k] * dense[i][j] * transform[j][l];
                }
            }
            near(r.get(k, l).unwrap(), expected);
        }
    }
    let v = c.reconstruct(&[3., -2.]).unwrap();
    for i in 0..4 {
        near(v[i], g[i] + transform[i][0] * 3. - transform[i][1] * 2.);
    }
}
#[test]
fn fully_constrained_and_zero_dimensional_systems_are_supported() {
    let mut b = AffineConstraints::new(3);
    b.add(0, 3., &[]).unwrap();
    b.add(1, 1., &[(0, 2.)]).unwrap();
    b.add(2, -1., &[(1, 0.5)]).unwrap();
    let c = b.close().unwrap();
    assert!(c.free_dofs().is_empty());
    assert_eq!(c.reconstruct(&[]).unwrap(), vec![3., 7., 2.5]);
    let a = matrix(&[&[1., 0., 0.], &[0., 1., 0.], &[0., 0., 1.]]);
    let (r, b) = c.condense(&a, &[1., 2., 3.]).unwrap();
    assert_eq!(r.rows(), 0);
    assert_eq!(r.cols(), 0);
    assert!(b.is_empty());
    let c = AffineConstraints::new(0).close().unwrap();
    let a = CooBuilder::new(0, 0).into_csr().unwrap();
    assert!(c.reconstruct(&[]).unwrap().is_empty());
    assert!(c.condense(&a, &[]).unwrap().1.is_empty());
}
#[test]
fn rejects_cycles_duplicates_invalid_bounds_and_nonfinite_data() {
    let mut b = AffineConstraints::new(3);
    assert!(b.add(3, 0., &[]).is_err());
    assert!(b.add(0, f64::NAN, &[]).is_err());
    assert!(b.add(0, 0., &[(3, 1.)]).is_err());
    assert!(b.add(0, 0., &[(1, f64::INFINITY)]).is_err());
    assert!(b.add(0, 0., &[(1, 0.), (1, 1.)]).is_err());
    b.add(0, 0., &[(1, 1.)]).unwrap();
    assert!(b.add(0, 0., &[]).is_err());
    b.add(1, 0., &[(2, 1.)]).unwrap();
    b.add(2, 0., &[(0, 1.)]).unwrap();
    assert!(b.close().is_err());
    let mut b = AffineConstraints::new(1);
    b.add(0, 0., &[(0, 1.)]).unwrap();
    assert!(b.close().is_err());
    let mut b = AffineConstraints::new(1);
    b.add(0, 2., &[(0, 0.)]).unwrap();
    assert_eq!(b.close().unwrap().reconstruct(&[]).unwrap(), vec![2.]);
    let c = AffineConstraints::new(2).close().unwrap();
    assert!(c.reconstruct(&[1.]).is_err());
    assert!(c.reconstruct(&[f64::NAN, 0.]).is_err());
    let a = matrix(&[&[1., 0.], &[0., 1.]]);
    assert!(c.condense(&a, &[1.]).is_err());
    assert!(c.condense(&a, &[1., f64::NAN]).is_err());
    assert!(c.condense(&matrix(&[&[1., 2.]]), &[1., 2.]).is_err());
}
#[test]
fn rejects_overflow_during_expansion_reconstruction_and_reduction() {
    let mut b = AffineConstraints::new(3);
    b.add(1, 0., &[(0, f64::MAX)]).unwrap();
    b.add(2, 0., &[(1, 2.)]).unwrap();
    assert!(b.close().is_err());
    let mut b = AffineConstraints::new(2);
    b.add(0, f64::MAX, &[]).unwrap();
    b.add(1, 0., &[(0, 2.)]).unwrap();
    assert!(b.close().is_err());
    let mut b = AffineConstraints::new(2);
    b.add(1, 0., &[(0, f64::MAX)]).unwrap();
    let c = b.close().unwrap();
    assert!(c.reconstruct(&[2.]).is_err());
    assert!(c
        .condense(&matrix(&[&[1., 0.], &[0., 2.]]), &[0., 0.])
        .is_err());
    let mut b = AffineConstraints::new(1);
    b.add(0, f64::MAX, &[]).unwrap();
    let c = b.close().unwrap();
    assert!(c.condense(&matrix(&[&[2.]]), &[0.]).is_err());
}
#[test]
fn deep_dependencies_are_resolved_without_recursion() {
    let n = 20_000;
    let mut b = AffineConstraints::new(n);
    for i in 1..n {
        b.add(i, 1., &[(i - 1, 1.)]).unwrap();
    }
    let c = b.close().unwrap();
    assert_eq!(c.free_dofs(), &[0]);
    let x = c.reconstruct(&[2.]).unwrap();
    near(x[n - 1], n as f64 + 1.);
}
#[test]
fn constrained_fem_affine_patch_with_dirichlet() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let f = |p: Point2| 1. + p.x + 2. * p.y;
    let (a, b) = numerics::assemble_poisson(&mesh, |_| 0.).unwrap();
    let mut constraints = AffineConstraints::new(mesh.vertices().len());
    for v in mesh.boundary_vertices() {
        constraints.add(v.0, f(mesh.vertices()[v.0]), &[]).unwrap();
    }
    // Interior node (1/3,1/3) is the midpoint of (0,1/3) and (2/3,1/3).
    constraints.add(5, 0., &[(4, 0.5), (6, 0.5)]).unwrap();
    let c = constraints.close().unwrap();
    let (reduced, rhs) = c.condense(&a, &b).unwrap();
    assert_eq!(c.free_dofs().len(), 3);
    for i in 0..reduced.rows() {
        for j in 0..reduced.cols() {
            near(reduced.get(i, j).unwrap(), reduced.get(j, i).unwrap());
        }
    }
    let (y, _) = conjugate_gradient(&reduced, &rhs, options()).unwrap();
    let u = c.reconstruct(&y).unwrap();
    for (x, p) in u.iter().zip(mesh.vertices()) {
        near(*x, f(*p));
    }
    let (l2, h1) = numerics::errors(&mesh, &u, f, |_| [1., 2.]).unwrap();
    assert!(l2 < 1e-11 && h1 < 1e-11);
}

#[test]
fn fixed_constraints_agree_with_existing_dirichlet_elimination() {
    let a = matrix(&[&[2., -1., 0.], &[-1., 2., -1.], &[0., -1., 2.]]);
    let b = [0., 1., 2.];
    let (full, rhs) = a.with_dirichlet(&b, &[(0, 3.), (2, 4.)]).unwrap();
    let (expected, _) = conjugate_gradient(&full, &rhs, options()).unwrap();
    let mut builder = AffineConstraints::new(3);
    builder.add(0, 3., &[]).unwrap();
    builder.add(2, 4., &[]).unwrap();
    let c = builder.close().unwrap();
    let (r, rhs) = c.condense(&a, &b).unwrap();
    let (y, _) = conjugate_gradient(&r, &rhs, options()).unwrap();
    let x = c.reconstruct(&y).unwrap();
    for (a, b) in x.iter().zip(expected) {
        near(*a, b);
    }
}
#[test]
fn manually_assembled_nonconforming_patch_has_continuous_trace() {
    // A coarse edge (1,2) meets two fine edges (1,6),(6,2). Assemble
    // local cells directly: Mesh2 cannot identify nonconforming topology.
    use dealii_rust::integrators::{diffusion_reaction_p1, DiffusionTensor};
    let p = |x, y| Point2::new(x, y).unwrap();
    let vertices = [
        p(0., 0.),
        p(0.5, 0.),
        p(0.5, 1.),
        p(0., 1.),
        p(1., 0.),
        p(1., 1.),
        p(0.5, 0.5),
        p(0.75, 0.5),
    ];
    let cells = [
        [0, 1, 2],
        [0, 2, 3],
        [1, 4, 7],
        [4, 5, 7],
        [5, 2, 7],
        [2, 6, 7],
        [6, 1, 7],
    ];
    let mut coo = CooBuilder::new(vertices.len(), vertices.len());
    for cell in cells {
        let local = diffusion_reaction_p1(
            cell.map(|i| vertices[i]),
            |_| DiffusionTensor::identity(),
            |_| 0.,
            |_| 0.,
        )
        .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                coo.add(cell[i], cell[j], local.matrix[i][j]).unwrap();
            }
        }
    }
    let a = coo.into_csr().unwrap();
    let exact = |p: Point2| 1. + p.x + 2. * p.y;
    let mut builder = AffineConstraints::new(vertices.len());
    for i in 0..6 {
        builder.add(i, exact(vertices[i]), &[]).unwrap();
    }
    builder.add(6, 0., &[(1, 0.5), (2, 0.5)]).unwrap();
    let c = builder.close().unwrap();
    assert_eq!(c.free_dofs(), &[7]);
    let (r, b) = c.condense(&a, &[0.; 8]).unwrap();
    let (y, _) = conjugate_gradient(&r, &b, options()).unwrap();
    let x = c.reconstruct(&y).unwrap();
    for (value, p) in x.iter().zip(vertices) {
        near(*value, exact(p));
    }
    for k in 0..=10 {
        let t = k as f64 / 10.;
        let coarse = (1. - t) * x[1] + t * x[2];
        let fine = if t <= 0.5 {
            (1. - 2. * t) * x[1] + 2. * t * x[6]
        } else {
            (2. - 2. * t) * x[6] + (2. * t - 1.) * x[2]
        };
        near(coarse, fine);
    }
}
