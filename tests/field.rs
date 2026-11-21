use dealii_rust::{
    base::Point2,
    field,
    grid::{BoundaryEdgeId, CellId, Mesh2, UniformRefinement, VertexId},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11, "{} != {}", a, b);
}
fn triangle() -> Mesh2 {
    Mesh2::new(
        vec![p(0., 0.), p(1., 0.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap()
}
fn strict_cg() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-13,
        absolute_tolerance: 1e-14,
        max_iterations: 1000,
    }
}

#[test]
fn consistent_mass_matrix_matches_analytic_and_global_assembly() {
    let mesh = triangle();
    let mass = field::mass_matrix(&mesh).unwrap();
    for i in 0..3 {
        for j in 0..3 {
            close(
                mass.get(i, j).unwrap(),
                if i == j { 1.0 / 12.0 } else { 1.0 / 24.0 },
            );
        }
    }
    let square = Mesh2::unit_square(1).unwrap();
    let m = field::mass_matrix(&square).unwrap();
    let mut sum = 0.0;
    for i in 0..4 {
        for j in 0..4 {
            close(m.get(i, j).unwrap(), m.get(j, i).unwrap());
            sum += m.get(i, j).unwrap();
        }
    }
    close(sum, 1.0);
    close(m.get(0, 3).unwrap(), 1.0 / 12.0); // shared diagonal in both triangles
    close(m.get(1, 2).unwrap(), 0.0);
}

#[test]
fn interpolation_and_projection_reproduce_affine_fields() {
    let mesh = Mesh2::unit_square(4).unwrap();
    for f in [|_: Point2| 1.0, |p: Point2| 3.0 + 2.0 * p.x - 5.0 * p.y] {
        let interpolated = field::interpolate(&mesh, f).unwrap();
        let (projected, report) = field::project_l2(&mesh, f, strict_cg()).unwrap();
        assert!(report.residual_norm < 1e-12);
        for (x, y) in interpolated.iter().zip(projected) {
            close(*x, y);
        }
    }
}

#[test]
fn quadratic_projection_has_analytic_coefficients_and_orthogonality() {
    let mesh = triangle();
    let f = |p: Point2| p.x * p.x;
    let interpolated = field::interpolate(&mesh, f).unwrap();
    let (matrix, rhs) = field::projection_system(&mesh, f).unwrap();
    for (actual, expected) in rhs.iter().zip([1.0 / 60.0, 1.0 / 20.0, 1.0 / 60.0]) {
        close(*actual, expected);
    }
    let (u, _) = conjugate_gradient(&matrix, &rhs, strict_cg()).unwrap();
    for (actual, expected) in u.iter().zip([-0.1, 0.7, -0.1]) {
        close(*actual, expected);
    }
    let product = matrix.mul(&u).unwrap();
    for (x, y) in product.iter().zip(rhs) {
        close(*x, y);
    }
    let exact = |p: Point2| p.x * p.x;
    let grad = |p: Point2| [2.0 * p.x, 0.0];
    let (projected_error, _) = numerics::errors(&mesh, &u, exact, grad).unwrap();
    let (interpolated_error, _) = numerics::errors(&mesh, &interpolated, exact, grad).unwrap();
    assert!(projected_error < interpolated_error * 0.6);
}

#[test]
fn projection_system_can_impose_dirichlet_values() {
    let mesh = Mesh2::unit_square(2).unwrap();
    let f = |p: Point2| 1.0 + p.x + 2.0 * p.y;
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let boundary = field::interpolate_boundary(&mesh, &edges, f).unwrap();
    let (m, b) = field::projection_system(&mesh, f).unwrap();
    let (m, b) = m.with_dirichlet(&b, &boundary).unwrap();
    let (u, _) = conjugate_gradient(&m, &b, strict_cg()).unwrap();
    for (i, p) in mesh.vertices().iter().enumerate() {
        close(u[i], f(*p));
    }
}

#[test]
fn boundary_interpolation_deduplicates_vertices_and_checks_inputs() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let b = field::interpolate_boundary(&mesh, &edges, |p| p.x + 2.0 * p.y).unwrap();
    assert_eq!(b, vec![(0, 0.0), (1, 1.0), (2, 2.0), (3, 3.0)]);
    assert!(field::interpolate_boundary(&mesh, &[], |_| f64::NAN)
        .unwrap()
        .is_empty());
    assert!(field::interpolate_boundary(&mesh, &[edges[0], edges[0]], |_| 1.0).is_err());
    assert!(field::interpolate_boundary(&mesh, &[BoundaryEdgeId(4)], |_| 1.0).is_err());
    assert!(field::interpolate_boundary(&mesh, &edges, |_| f64::INFINITY).is_err());
}

#[test]
fn evaluate_skew_triangle_vertices_edges_inside_outside() {
    let mesh = Mesh2::new(
        vec![p(2., 1.), p(5., 2.), p(3., 4.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    let values = [2.0, 8.0, -1.0];
    let location = field::locate(&mesh, p(3.0, 2.0), 0.0).unwrap().unwrap();
    assert_eq!(location.cell, CellId(0));
    for (actual, expected) in location.barycentric.iter().zip([0.5, 0.25, 0.25]) {
        close(*actual, expected);
    }
    close(
        field::evaluate(&mesh, &values, p(3., 2.), 0.0)
            .unwrap()
            .unwrap(),
        2.75,
    );
    for (i, &pt) in mesh.vertices().iter().enumerate() {
        close(
            field::evaluate(&mesh, &values, pt, 1e-14).unwrap().unwrap(),
            values[i],
        );
    }
    close(
        field::evaluate(&mesh, &values, p(3.5, 1.5), 1e-14)
            .unwrap()
            .unwrap(),
        5.0,
    );
    assert!(field::locate(&mesh, p(1., 1.), 0.0).unwrap().is_none());
    assert!(field::evaluate(&mesh, &values, p(10., 10.), 0.0)
        .unwrap()
        .is_none());
}

#[test]
fn tolerance_is_explicit_and_exact_cell_is_preferred() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let just_outside = p(-1e-6, 0.5);
    assert!(field::locate(&mesh, just_outside, 0.0).unwrap().is_none());
    let located = field::locate(&mesh, just_outside, 1e-5).unwrap().unwrap();
    assert!(located.barycentric.iter().any(|&v| v < 0.0));
    let values = field::interpolate(&mesh, |p| p.x).unwrap();
    close(
        field::evaluate(&mesh, &values, just_outside, 1e-5)
            .unwrap()
            .unwrap(),
        -1e-6,
    );
    // (0.49, 0.51) is just across the diagonal: cell 0 is within the
    // tolerance but cell 1 contains it. No extrapolation should be selected.
    assert_eq!(
        field::locate(&mesh, p(0.49, 0.51), 0.1)
            .unwrap()
            .unwrap()
            .cell,
        CellId(1)
    );
    assert_eq!(
        field::locate(&mesh, p(0.5, 0.5), 0.0)
            .unwrap()
            .unwrap()
            .cell,
        CellId(0)
    );
}

#[test]
fn evaluation_matches_prolonged_piecewise_p1_including_shared_edges() {
    let coarse = Mesh2::unit_square(2).unwrap();
    let values: Vec<_> = (0..coarse.vertices().len())
        .map(|i| (i * i) as f64)
        .collect();
    let refinement = UniformRefinement::new(&coarse).unwrap();
    let fine_values = refinement.prolongate_p1(&values).unwrap();
    for x in 0..=16 {
        for y in 0..=16 {
            let point = p(x as f64 / 16.0, y as f64 / 16.0);
            let a = field::evaluate(&coarse, &values, point, 1e-12)
                .unwrap()
                .unwrap();
            let b = field::evaluate(refinement.mesh(), &fine_values, point, 1e-12)
                .unwrap()
                .unwrap();
            close(a, b);
        }
    }
}

#[test]
fn rejects_nonfinite_and_invalid_inputs_and_solver_failure() {
    let mesh = triangle();
    assert!(field::interpolate(&mesh, |_| f64::NAN).is_err());
    assert!(field::projection_system(&mesh, |_| f64::INFINITY).is_err());
    let large = Mesh2::new(
        vec![p(0., 0.), p(10., 0.), p(0., 10.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    assert!(field::projection_system(&large, |_| f64::MAX).is_err());
    assert!(field::evaluate(&mesh, &[1.0], p(0.1, 0.1), 0.0).is_err());
    assert!(field::evaluate(&mesh, &[f64::NAN; 3], p(0.1, 0.1), 0.0).is_err());
    assert!(field::locate(
        &mesh,
        Point2 {
            x: f64::NAN,
            y: 0.0
        },
        0.0
    )
    .is_err());
    assert!(field::locate(&mesh, p(0.1, 0.1), -0.1).is_err());
    assert!(field::locate(&mesh, p(0.1, 0.1), f64::INFINITY).is_err());
    let opts = CgOptions {
        max_iterations: 0,
        ..strict_cg()
    };
    assert!(field::project_l2(&mesh, |_| 1.0, opts).is_err());
    let opts = CgOptions {
        relative_tolerance: f64::NAN,
        ..strict_cg()
    };
    assert!(field::project_l2(&mesh, |_| 1.0, opts).is_err());
}

#[test]
fn l2_projection_converges_and_is_more_accurate_than_interpolation() {
    let f = |p: Point2| p.x * p.x + p.y * p.y;
    let grad = |p: Point2| [2.0 * p.x, 2.0 * p.y];
    let mut prev = f64::INFINITY;
    for n in [2, 4, 8, 16] {
        let mesh = Mesh2::unit_square(n).unwrap();
        let interpolated = field::interpolate(&mesh, f).unwrap();
        let (projected, _) = field::project_l2(&mesh, f, strict_cg()).unwrap();
        let (pi, _) = numerics::errors(&mesh, &interpolated, f, grad).unwrap();
        let (pr, _) = numerics::errors(&mesh, &projected, f, grad).unwrap();
        assert!(pr < pi, "projection={}, interpolation={}", pr, pi);
        assert!(pr < prev * 0.4, "projection={}, previous={}", pr, prev);
        prev = pr;
    }
}

#[test]
fn disconnected_mesh_projection_and_lookup() {
    let mesh = Mesh2::new(
        vec![
            p(0., 0.),
            p(1., 0.),
            p(0., 1.),
            p(3., 0.),
            p(4., 0.),
            p(3., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(3), VertexId(4), VertexId(5)],
        ],
    )
    .unwrap();
    let f = |p: Point2| 1.0 + 2.0 * p.x - p.y;
    let (u, _) = field::project_l2(&mesh, f, strict_cg()).unwrap();
    for (value, pt) in u.iter().zip(mesh.vertices()) {
        close(*value, f(*pt));
    }
    assert!(field::evaluate(&mesh, &u, p(2., 0.2), 1e-12)
        .unwrap()
        .is_none());
    close(
        field::evaluate(&mesh, &u, p(3.2, 0.3), 1e-12)
            .unwrap()
            .unwrap(),
        f(p(3.2, 0.3)),
    );
}

#[test]
fn gmsh_boundary_interpolation_drives_mixed_poisson() {
    let imported = dealii_rust::gmsh::parse(include_str!("fixtures/square41.msh")).unwrap();
    let mesh = imported.mesh();
    let f = |p: Point2| 1.0 + p.x + 2.0 * p.y;
    let dirichlet =
        field::interpolate_boundary(mesh, &imported.boundary_edges_with_tag(40), f).unwrap();
    let mut neumann = Vec::new();
    for tag in [10, 20, 30] {
        neumann.extend(imported.boundary_edges_with_tag(tag));
    }
    let (matrix, mut rhs) = numerics::assemble_poisson(mesh, |_| 0.0).unwrap();
    let load = numerics::assemble_neumann(mesh, &neumann, |_, n| n[0] + 2.0 * n[1]).unwrap();
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
    }
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &dirichlet).unwrap();
    let (solution, _) = conjugate_gradient(&matrix, &rhs, strict_cg()).unwrap();
    let point = p(0.73, 0.24);
    close(
        field::evaluate(mesh, &solution, point, 1e-12)
            .unwrap()
            .unwrap(),
        f(point),
    );
}

#[test]
fn tiny_unrepresentable_mass_entries_are_rejected() {
    let e = 4e-162;
    let mesh = Mesh2::new(
        vec![p(0., 0.), p(e, 0.), p(0., e)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    assert!(field::mass_matrix(&mesh).is_err());
}
