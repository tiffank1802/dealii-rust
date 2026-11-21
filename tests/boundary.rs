use dealii_rust::{
    base::Point2,
    dofs::DofHandler,
    grid::{BoundaryEdgeId, Mesh2, VertexId},
    lac::{conjugate_gradient, CgOptions},
    numerics,
    quadrature::line_degree_5,
};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{} != {}", a, b);
}
fn all_edges(mesh: &Mesh2) -> Vec<BoundaryEdgeId> {
    (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect()
}

#[test]
fn line_gauss_integrates_polynomials_to_degree_five() {
    for degree in 0..=5 {
        let result: f64 = line_degree_5()
            .iter()
            .map(|q| q.weight * q.point.powi(degree))
            .sum();
        near(result, 1.0 / (degree as f64 + 1.0), 1e-14);
    }
}

#[test]
fn boundary_normals_follow_the_adjacent_cell_not_sorted_vertex_ids() {
    let mesh = Mesh2::new(
        vec![p(0.0, 0.0), p(2.0, 0.0), p(0.0, 1.0)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    for id in all_edges(&mesh) {
        let edge = mesh.boundary_edge_geometry(id).unwrap();
        let midpoint = p(
            (edge.points[0].x + edge.points[1].x) / 2.0,
            (edge.points[0].y + edge.points[1].y) / 2.0,
        );
        let center = p(2.0 / 3.0, 1.0 / 3.0);
        assert!(
            edge.outward_normal[0] * (midpoint.x - center.x)
                + edge.outward_normal[1] * (midpoint.y - center.y)
                > 0.0
        );
        near(
            edge.outward_normal[0].hypot(edge.outward_normal[1]),
            1.0,
            1e-14,
        );
        let canonical = mesh.boundary_edges()[id.0];
        assert!(canonical.contains(&edge.vertices[0]) && canonical.contains(&edge.vertices[1]));
        if canonical == [VertexId(0), VertexId(2)] {
            assert_eq!(edge.vertices, [VertexId(2), VertexId(0)]);
            near(edge.outward_normal[0], -1.0, 1e-14);
            near(edge.outward_normal[1], 0.0, 1e-14);
        }
    }
    assert!(mesh.boundary_edge_geometry(BoundaryEdgeId(3)).is_err());
}

#[test]
fn shared_edges_must_be_oriented_in_opposite_directions() {
    // Both CCW triangles lie on the same side of the shared edge.
    let mesh = Mesh2::new(
        vec![p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0), p(0.5, 1.0)],
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(0), VertexId(1), VertexId(3)],
        ],
    );
    assert!(mesh.is_err());

    // The first two cells are oppositely oriented on the shared edge; a
    // third incident cell is still non-manifold.
    let mesh = Mesh2::new(
        vec![
            p(0.0, 0.0),
            p(1.0, 0.0),
            p(0.0, 1.0),
            p(0.0, -1.0),
            p(0.5, 1.0),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(1), VertexId(0), VertexId(3)],
            [VertexId(0), VertexId(1), VertexId(4)],
        ],
    );
    assert!(mesh.is_err());
}

#[test]
fn neumann_integrates_linear_flux_on_skew_edge() {
    let mesh = Mesh2::new(
        vec![p(0.0, 0.0), p(2.0, 0.0), p(0.0, 1.0)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    let id = BoundaryEdgeId(
        mesh.boundary_edges()
            .iter()
            .position(|e| *e == [VertexId(1), VertexId(2)])
            .unwrap(),
    );
    let load = numerics::assemble_neumann(&mesh, &[id], |point, normal| {
        near(normal[0], 1.0 / 5.0_f64.sqrt(), 1e-14);
        near(normal[1], 2.0 / 5.0_f64.sqrt(), 1e-14);
        point.x
    })
    .unwrap();
    near(load[0], 0.0, 1e-14);
    near(load[1], 2.0 * 5.0_f64.sqrt() / 3.0, 1e-14);
    near(load[2], 5.0_f64.sqrt() / 3.0, 1e-14);
    assert!(numerics::assemble_neumann(&mesh, &[id], |_, _| f64::MAX).is_err());
}

#[test]
fn neumann_load_validates_edges_and_flux() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let ids = all_edges(&mesh);
    assert_eq!(
        numerics::assemble_neumann(&mesh, &[], |_, _| 1.0).unwrap(),
        vec![0.0; 4]
    );
    assert!(numerics::assemble_neumann(&mesh, &[ids[0], ids[0]], |_, _| 1.0).is_err());
    assert!(numerics::assemble_neumann(&mesh, &[BoundaryEdgeId(4)], |_, _| 1.0).is_err());
    assert!(numerics::assemble_neumann(&mesh, &ids, |_, _| f64::NAN).is_err());
    let load = numerics::assemble_neumann(&mesh, &ids, |_, _| 1.0).unwrap();
    near(load.iter().sum(), 4.0, 1e-14);
}

fn mixed_solution(n: usize, quadratic: bool) -> (f64, f64, Vec<f64>, Mesh2) {
    let mesh = Mesh2::unit_square(n).unwrap();
    let exact = |p: Point2| {
        if quadratic {
            1.0 + p.x + p.x * p.x + 2.0 * p.y + p.y * p.y
        } else {
            1.0 + p.x + 2.0 * p.y
        }
    };
    let gradient = |p: Point2| {
        if quadratic {
            [1.0 + 2.0 * p.x, 2.0 + 2.0 * p.y]
        } else {
            [1.0, 2.0]
        }
    };
    let (a, mut b) =
        numerics::assemble_poisson(&mesh, |_| if quadratic { -4.0 } else { 0.0 }).unwrap();
    let ids: Vec<_> = all_edges(&mesh)
        .into_iter()
        .filter(|id| {
            let edge = mesh.boundary_edge_geometry(*id).unwrap();
            !(edge.points[0].x == 0.0 && edge.points[1].x == 0.0)
        })
        .collect();
    let load = numerics::assemble_neumann(&mesh, &ids, |p, normal| {
        let g = gradient(p);
        g[0] * normal[0] + g[1] * normal[1]
    })
    .unwrap();
    for (b, value) in b.iter_mut().zip(load) {
        *b += value;
    }
    let constraints: Vec<_> = DofHandler::new(&mesh)
        .boundary_dofs()
        .into_iter()
        .filter(|d| mesh.vertices()[d.0].x == 0.0)
        .map(|d| (d.0, exact(mesh.vertices()[d.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    let (l2, h1) = numerics::errors(&mesh, &u, exact, gradient).unwrap();
    (l2, h1, u, mesh)
}

#[test]
fn mixed_boundaries_reproduce_affine_solution() {
    for n in [2, 4, 8] {
        let (l2, h1, u, mesh) = mixed_solution(n, false);
        assert!(l2 < 1e-9 && h1 < 1e-9, "L2={}, H1={}", l2, h1);
        for (value, point) in u.iter().zip(mesh.vertices()) {
            near(*value, 1.0 + point.x + 2.0 * point.y, 1e-9);
        }
    }
}

#[test]
fn mixed_quadratic_solution_converges_in_l2_and_h1() {
    let a = mixed_solution(8, true);
    let b = mixed_solution(16, true);
    let c = mixed_solution(32, true);
    for (old, new) in [(&a, &b), (&b, &c)] {
        let l2_rate = (old.0 / new.0).log2();
        let h1_rate = (old.1 / new.1).log2();
        assert!((1.8..2.2).contains(&l2_rate), "L2 rate {}", l2_rate);
        assert!((0.9..1.1).contains(&h1_rate), "H1 rate {}", h1_rate);
    }
}

#[test]
fn pure_neumann_matrix_has_constant_kernel_and_load_must_balance() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let (a, _) = numerics::assemble_poisson(&mesh, |_| 0.0).unwrap();
    for value in a.mul(&vec![1.0; mesh.vertices().len()]).unwrap() {
        near(value, 0.0, 1e-14);
    }
    let compatible =
        numerics::assemble_neumann(&mesh, &all_edges(&mesh), |_, n| n[0] + 2.0 * n[1]).unwrap();
    near(compatible.iter().sum(), 0.0, 1e-14);
    let incompatible = numerics::assemble_neumann(&mesh, &all_edges(&mesh), |_, _| 1.0).unwrap();
    near(incompatible.iter().sum(), 4.0, 1e-14);
    // A normalization or nullspace-aware solver is required; no CG solve here.
}
