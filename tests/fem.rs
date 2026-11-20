#![allow(clippy::needless_range_loop)] // Indexed matrices keep the analytical fixtures readable.

use dealii_rust::{
    base::{triangle_determinant, Point2},
    dofs::DofHandler,
    fe::{shape_values, AffineTriangle},
    grid::{CellId, Mesh2, VertexId},
    lac::{conjugate_gradient, CgOptions, CooBuilder},
    numerics,
    quadrature::triangle_degree_5,
    vtk,
};
use std::f64::consts::PI;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |a, b| a * b as f64)
}

#[test]
fn mesh_geometry_boundary_and_dofs() {
    for n in &[1, 2, 5] {
        let mesh = Mesh2::unit_square(*n).unwrap();
        assert_eq!(mesh.vertices().len(), (n + 1) * (n + 1));
        assert_eq!(mesh.cells().len(), 2 * n * n);
        assert_eq!(mesh.boundary_edges().len(), 4 * n);
        assert_eq!(mesh.boundary_vertices().len(), 4 * n);
        let area: f64 = (0..mesh.cells().len())
            .map(|i| triangle_determinant(mesh.cell_points(CellId(i)).unwrap()).unwrap() / 2.0)
            .sum();
        near(area, 1.0, 1e-12);
    }
    let mesh = Mesh2::unit_square(1).unwrap();
    let dofs = DofHandler::new(&mesh);
    assert_eq!(dofs.count(), 4);
    let a = dofs.cell_dofs(CellId(0)).unwrap();
    let b = dofs.cell_dofs(CellId(1)).unwrap();
    assert_eq!(a.iter().filter(|x| b.contains(x)).count(), 2);
    assert_eq!(dofs.boundary_dofs().len(), 4);
    assert!(dofs.cell_dofs(CellId(2)).is_err());
}

#[test]
fn invalid_meshes_and_degenerate_triangles_are_rejected() {
    assert!(Point2::new(f64::NAN, 0.0).is_err());
    assert!(Mesh2::unit_square(0).is_err());
    assert!(Mesh2::unit_square(usize::MAX).is_err());
    let vertices = vec![p(0., 0.), p(1., 0.), p(0., 1.)];
    assert!(Mesh2::new(
        vertices.clone(),
        vec![[VertexId(0), VertexId(2), VertexId(1)]]
    )
    .is_err());
    assert!(Mesh2::new(
        vertices.clone(),
        vec![[VertexId(0), VertexId(1), VertexId(3)]]
    )
    .is_err());
    assert!(Mesh2::new(
        vertices.clone(),
        vec![[VertexId(0), VertexId(1), VertexId(1)]]
    )
    .is_err());
    assert!(Mesh2::new(
        vertices,
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(0), VertexId(1), VertexId(2)]
        ]
    )
    .is_err());
    let v = vec![p(0., 0.), p(1., 0.), p(0., 1.), p(0.3, 1.), p(0.6, 1.)];
    assert!(Mesh2::new(
        v,
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(0), VertexId(1), VertexId(3)],
            [VertexId(0), VertexId(1), VertexId(4)]
        ]
    )
    .is_err());
    assert!(triangle_determinant([p(0., 0.), p(1., 0.), p(0., 1e-16)]).is_err());
}

#[test]
fn quadrature_integrates_monomials_through_degree_five() {
    let q = triangle_degree_5();
    for a in 0..=5 {
        for b in 0..=5 - a {
            let computed: f64 = q
                .iter()
                .map(|v| v.weight * v.point.x.powi(a as i32) * v.point.y.powi(b as i32))
                .sum();
            near(
                computed,
                factorial(a) * factorial(b) / factorial(a + b + 2),
                2e-14,
            );
        }
    }
}

#[test]
fn p1_shapes_and_physical_gradients() {
    for (i, point) in [p(0., 0.), p(1., 0.), p(0., 1.)].iter().enumerate() {
        let values = shape_values(*point);
        for j in 0..3 {
            near(values[j], if i == j { 1.0 } else { 0.0 }, 1e-14);
        }
    }
    let map = AffineTriangle::new([p(1., 2.), p(3., 2.), p(1., 5.)]).unwrap();
    near(map.determinant(), 6., 1e-14);
    assert_eq!(map.map(p(0.5, 0.25)), p(2., 2.75));
    let g = map.gradients();
    for i in 0..3 {
        near(g[i][0], [-0.5, 0.5, 0.0][i], 1e-14);
        near(g[i][1], [-1.0 / 3.0, 0.0, 1.0 / 3.0][i], 1e-14);
    }
}

#[test]
fn local_stiffness_and_load_are_correct() {
    let mesh = Mesh2::new(
        vec![p(0., 0.), p(1., 0.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    let (a, b) = numerics::assemble_poisson(&mesh, |_| 1.0).unwrap();
    let expected = [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]];
    for i in 0..3 {
        for j in 0..3 {
            near(a.get(i, j).unwrap(), expected[i][j], 1e-13);
        }
        near(b[i], 1.0 / 6.0, 1e-13);
    }
    assert!(numerics::assemble_poisson(&mesh, |_| f64::NAN).is_err());
}

#[test]
fn coo_duplicates_empty_rows_and_invalid_values() {
    let mut b = CooBuilder::new(3, 3);
    b.add(0, 0, 2.0).unwrap();
    b.add(0, 0, 3.0).unwrap();
    b.add(2, 1, -2.0).unwrap();
    assert!(b.add(3, 0, 1.0).is_err());
    assert!(b.add(0, 0, f64::NAN).is_err());
    let a = b.into_csr().unwrap();
    assert_eq!(a.mul(&[1., 2., 3.]).unwrap(), vec![5., 0., -4.]);
    assert!(a.mul(&[1.]).is_err());
    assert!(a.get(0, 3).is_err());
    assert!(CooBuilder::new(usize::MAX, 1).into_csr().is_err());
    let mut b = CooBuilder::new(1, 1);
    b.add(0, 0, f64::MAX).unwrap();
    b.add(0, 0, f64::MAX).unwrap();
    assert!(b.into_csr().is_err());
}

#[test]
fn nonzero_dirichlet_is_symmetric_and_cg_solves() {
    let mut coo = CooBuilder::new(3, 3);
    for &(i, j, v) in &[
        (0, 0, 2.),
        (0, 1, -1.),
        (1, 0, -1.),
        (1, 1, 2.),
        (1, 2, -1.),
        (2, 1, -1.),
        (2, 2, 2.),
    ] {
        coo.add(i, j, v).unwrap();
    }
    let matrix = coo.into_csr().unwrap();
    let (a, b) = matrix
        .with_dirichlet(&[0., 0., 0.], &[(0, 2.), (2, 3.)])
        .unwrap();
    assert_eq!(b, vec![2., 5., 3.]);
    for i in 0..3 {
        for j in 0..3 {
            near(a.get(i, j).unwrap(), a.get(j, i).unwrap(), 1e-14);
        }
    }
    let (x, report) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    for (v, expected) in x.iter().zip([2., 2.5, 3.].iter()) {
        near(*v, *expected, 1e-10);
    }
    assert!(report.residual_norm < 1e-10);
    assert!(matrix
        .with_dirichlet(&[0., 0., 0.], &[(0, 1.), (0, 1.)])
        .is_err());
    assert!(matrix.with_dirichlet(&[0., 0., 0.], &[(3, 1.)]).is_err());
    assert_eq!(matrix.get(0, 1).unwrap(), -1.0);
}

#[test]
fn cg_errors_and_zero_rhs() {
    let mut b = CooBuilder::new(1, 1);
    b.add(0, 0, 1.0).unwrap();
    let a = b.into_csr().unwrap();
    let (x, report) = conjugate_gradient(&a, &[0.0], CgOptions::default()).unwrap();
    assert_eq!(x, vec![0.0]);
    assert_eq!(report.iterations, 0);
    let mut opts = CgOptions {
        max_iterations: 0,
        ..CgOptions::default()
    };
    assert!(conjugate_gradient(&a, &[1.0], opts).is_err());
    assert!(conjugate_gradient(&a, &[], CgOptions::default()).is_err());
    opts.relative_tolerance = -1.0;
    assert!(conjugate_gradient(&a, &[1.0], opts).is_err());
    let mut neg = CooBuilder::new(1, 1);
    neg.add(0, 0, -1.0).unwrap();
    assert!(conjugate_gradient(&neg.into_csr().unwrap(), &[1.0], CgOptions::default()).is_err());
}

fn solve(n: usize) -> (f64, f64) {
    let mesh = Mesh2::unit_square(n).unwrap();
    let (a, b) = numerics::assemble_poisson(&mesh, |p| {
        2.0 * PI * PI * (PI * p.x).sin() * (PI * p.y).sin()
    })
    .unwrap();
    let boundary: Vec<_> = DofHandler::new(&mesh)
        .boundary_dofs()
        .iter()
        .map(|d| (d.0, 0.0))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &boundary).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    numerics::errors(
        &mesh,
        &u,
        |p| (PI * p.x).sin() * (PI * p.y).sin(),
        |p| {
            [
                PI * (PI * p.x).cos() * (PI * p.y).sin(),
                PI * (PI * p.x).sin() * (PI * p.y).cos(),
            ]
        },
    )
    .unwrap()
}

#[test]
fn poisson_manufactured_solution_converges() {
    let e8 = solve(8);
    let e16 = solve(16);
    let e32 = solve(32);
    for (old, new) in [(e8, e16), (e16, e32)].iter() {
        let l2_rate = (old.0 / new.0).log2();
        let h1_rate = (old.1 / new.1).log2();
        assert!(l2_rate > 1.8 && l2_rate < 2.2, "L2 rate {}", l2_rate);
        assert!(h1_rate > 0.9 && h1_rate < 1.1, "H1 rate {}", h1_rate);
    }
}

#[test]
fn vtk_contains_mesh_and_point_data() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let mut data = Vec::new();
    vtk::write_scalar(&mut data, &mesh, &[0., 1., 2., 3.]).unwrap();
    let s = String::from_utf8(data).unwrap();
    assert!(s.contains("POINTS 4 double"));
    assert!(s.contains("CELLS 2 8"));
    assert!(s.contains("CELL_TYPES 2\n5\n5"));
    assert!(s.contains("POINT_DATA 4"));
    assert!(vtk::write_scalar(Vec::new(), &mesh, &[0.]).is_err());
}

#[test]
fn global_stiffness_on_two_triangles() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let (a, b) = numerics::assemble_poisson(&mesh, |_| 0.0).unwrap();
    let expected = [
        [1.0, -0.5, -0.5, 0.0],
        [-0.5, 1.0, 0.0, -0.5],
        [-0.5, 0.0, 1.0, -0.5],
        [0.0, -0.5, -0.5, 1.0],
    ];
    for i in 0..4 {
        near(b[i], 0.0, 1e-14);
        for j in 0..4 {
            near(a.get(i, j).unwrap(), expected[i][j], 1e-14);
        }
    }
}

#[test]
fn affine_field_has_zero_integrated_error() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let u: Vec<_> = mesh
        .vertices()
        .iter()
        .map(|p| 1.0 + 2.0 * p.x - 3.0 * p.y)
        .collect();
    let (l2, h1) =
        numerics::errors(&mesh, &u, |p| 1.0 + 2.0 * p.x - 3.0 * p.y, |_| [2.0, -3.0]).unwrap();
    assert!(l2 < 1e-12 && h1 < 1e-12);
    assert!(numerics::errors(&mesh, &u[..1], |_| 0.0, |_| [0.0, 0.0]).is_err());
    assert!(numerics::errors(&mesh, &u, |_| f64::NAN, |_| [0.0, 0.0]).is_err());
}
