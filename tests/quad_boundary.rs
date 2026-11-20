use dealii_rust::{
    base::Point2,
    grid::{BoundaryEdgeId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    quad::{self, QuadMesh},
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn edges(mesh: &QuadMesh) -> Vec<BoundaryEdgeId> {
    (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect()
}
fn trapezoid() -> QuadMesh {
    QuadMesh::new(
        vec![p(0., 0.), p(2., 0.), p(3., 1.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap()
}

#[test]
fn normals_use_cell_orientation_not_canonical_vertex_order() {
    let mesh = trapezoid();
    for id in edges(&mesh) {
        let geometry = mesh.boundary_edge_geometry(id).unwrap();
        let midpoint = p(
            (geometry.points[0].x + geometry.points[1].x) / 2.,
            (geometry.points[0].y + geometry.points[1].y) / 2.,
        );
        let center = p(1.25, 0.5);
        assert!(
            geometry.outward_normal[0] * (midpoint.x - center.x)
                + geometry.outward_normal[1] * (midpoint.y - center.y)
                > 0.
        );
        near(
            geometry.outward_normal[0].hypot(geometry.outward_normal[1]),
            1.,
            1e-14,
        );
        let canonical = mesh.boundary_edges()[id.0];
        assert!(
            canonical.contains(&geometry.vertices[0]) && canonical.contains(&geometry.vertices[1])
        );
        if canonical == [VertexId(0), VertexId(3)] {
            assert_eq!(geometry.vertices, [VertexId(3), VertexId(0)]);
            near(geometry.outward_normal[0], -1., 1e-14);
        }
    }
    assert!(mesh.boundary_edge_geometry(BoundaryEdgeId(4)).is_err());
}

#[test]
fn linear_flux_integral_on_skew_edge_is_analytic() {
    let mesh = trapezoid();
    let id = edges(&mesh)
        .into_iter()
        .find(|id| mesh.boundary_edges()[id.0] == [VertexId(1), VertexId(2)])
        .unwrap();
    let length = 2.0_f64.sqrt();
    let b = quad::assemble_neumann(&mesh, &[id], |p, n| {
        near(n[0], 1. / length, 1e-14);
        near(n[1], -1. / length, 1e-14);
        p.x
    })
    .unwrap();
    near(b[0], 0., 1e-14);
    near(b[1], 7. * length / 6., 1e-14);
    near(b[2], 8. * length / 6., 1e-14);
    near(b[3], 0., 1e-14);
}

#[test]
fn boundary_load_and_data_validate_inputs() {
    let mesh = QuadMesh::unit_square(1).unwrap();
    let ids = edges(&mesh);
    assert_eq!(
        quad::assemble_neumann(&mesh, &[], |_, _| f64::NAN).unwrap(),
        vec![0.; 4]
    );
    assert!(quad::assemble_neumann(&mesh, &[ids[0], ids[0]], |_, _| 1.).is_err());
    assert!(quad::assemble_neumann(&mesh, &[BoundaryEdgeId(ids.len())], |_, _| 1.).is_err());
    assert!(quad::assemble_neumann(&mesh, &ids, |_, _| f64::NAN).is_err());
    assert!(quad::assemble_neumann(&mesh, &ids, |_, _| f64::INFINITY).is_err());
    assert!(quad::interpolate_boundary(&mesh, &[ids[0], ids[0]], |_| 1.).is_err());
    assert!(quad::interpolate_boundary(&mesh, &[BoundaryEdgeId(ids.len())], |_| 1.).is_err());
    assert!(quad::interpolate_boundary(&mesh, &ids, |_| f64::NAN).is_err());
    assert!(quad::interpolate_boundary(&mesh, &[], |_| f64::NAN)
        .unwrap()
        .is_empty());
    let constraints = quad::interpolate_boundary(&mesh, &ids, |p| 1. + p.x + 2. * p.y).unwrap();
    assert_eq!(constraints, vec![(0, 1.), (1, 2.), (2, 3.), (3, 4.)]);
    let big = QuadMesh::new(
        vec![p(0., 0.), p(10., 0.), p(10., 10.), p(0., 10.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(quad::assemble_neumann(&big, &edges(&big), |_, _| f64::MAX).is_err());
}

#[test]
fn neumann_load_has_correct_perimeter_balance_and_constant_kernel() {
    let mesh = QuadMesh::unit_square(4).unwrap();
    let ids = edges(&mesh);
    let b = quad::assemble_neumann(&mesh, &ids, |_, _| 1.).unwrap();
    near(b.iter().sum(), 4., 1e-13);
    let (matrix, source) = quad::assemble_poisson(&mesh, |_| -4.).unwrap();
    let g = quad::assemble_neumann(&mesh, &ids, |p, n| 2. * p.x * n[0] + 2. * p.y * n[1]).unwrap();
    near(
        source.iter().sum::<f64>() + g.iter().sum::<f64>(),
        0.,
        1e-12,
    );
    for value in matrix.mul(&vec![1.; mesh.vertices().len()]).unwrap() {
        near(value, 0., 1e-13);
    }
    // Without normalization, the pure-Neumann matrix remains singular.
}

fn exact(p: Point2, quadratic: bool) -> f64 {
    1. + p.x + 2. * p.y + if quadratic { p.x * p.x + p.y * p.y } else { 0. }
}
fn gradient(p: Point2, quadratic: bool) -> [f64; 2] {
    if quadratic {
        [1. + 2. * p.x, 2. + 2. * p.y]
    } else {
        [1., 2.]
    }
}
fn mixed(n: usize, quadratic: bool, distorted: bool) -> (f64, f64) {
    let mut mesh = QuadMesh::unit_square(n).unwrap();
    if distorted {
        let vertices = mesh
            .vertices()
            .iter()
            .map(|p| Point2 {
                x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
                y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
            })
            .collect();
        mesh = QuadMesh::new(vertices, mesh.cells().to_vec()).unwrap();
    }
    let mut dirichlet = Vec::new();
    let mut neumann = Vec::new();
    for id in edges(&mesh) {
        let edge = mesh.boundary_edge_geometry(id).unwrap();
        if edge.points.iter().all(|p| p.x == 0.) {
            dirichlet.push(id)
        } else {
            neumann.push(id)
        }
    }
    let (a, mut b) = quad::assemble_poisson(&mesh, |_| if quadratic { -4. } else { 0. }).unwrap();
    let load = quad::assemble_neumann(&mesh, &neumann, |p, n| {
        let g = gradient(p, quadratic);
        g[0] * n[0] + g[1] * n[1]
    })
    .unwrap();
    for (value, g) in b.iter_mut().zip(load) {
        *value += g;
    }
    let constraints =
        quad::interpolate_boundary(&mesh, &dirichlet, |p| exact(p, quadratic)).unwrap();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (u, _) = conjugate_gradient(&a, &b, options).unwrap();
    quad::errors(
        &mesh,
        &u,
        |p| exact(p, quadratic),
        |p| gradient(p, quadratic),
    )
    .unwrap()
}
#[test]
fn mixed_affine_solution_is_reproduced_even_on_distorted_mesh() {
    for n in [2, 4, 8] {
        for distorted in [false, true] {
            let (l2, h1) = mixed(n, false, distorted);
            assert!(l2 < 1e-9 && h1 < 1e-9, "L2={} H1={}", l2, h1);
        }
    }
}
#[test]
fn mixed_quadratic_solution_converges_on_square_and_distorted_meshes() {
    for distorted in [false, true] {
        let values: Vec<_> = [8, 16, 32]
            .iter()
            .map(|&n| mixed(n, true, distorted))
            .collect();
        for pair in values.windows(2) {
            let l2 = (pair[0].0 / pair[1].0).log2();
            let h1 = (pair[0].1 / pair[1].1).log2();
            assert!((1.8..2.2).contains(&l2), "L2 rate {}", l2);
            assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
        }
    }
}

#[test]
fn mixed_affine_solution_on_skew_nonaffine_domain() {
    let base = QuadMesh::unit_square(5).unwrap();
    let points = base
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.3 * p.y + 0.1 * p.x * (1. - p.x) * p.y * (1. - p.y),
            y: p.y + 0.2 * p.x + 0.1 * p.x * (1. - p.x) * p.y * (1. - p.y),
        })
        .collect();
    let mesh = QuadMesh::new(points, base.cells().to_vec()).unwrap();
    let mut d = Vec::new();
    let mut n = Vec::new();
    for id in edges(&mesh) {
        // IDs correspond to the original mesh because its connectivity is unchanged.
        if mesh.boundary_edges()[id.0]
            .iter()
            .all(|v| base.vertices()[v.0].x == 0.)
        {
            d.push(id);
        } else {
            n.push(id);
        }
    }
    let exact = |p: Point2| 1. + 2. * p.x - 3. * p.y;
    let (a, mut rhs) = quad::assemble_poisson(&mesh, |_| 0.).unwrap();
    let g = quad::assemble_neumann(&mesh, &n, |_, normal| 2. * normal[0] - 3. * normal[1]).unwrap();
    for (b, value) in rhs.iter_mut().zip(g) {
        *b += value;
    }
    let constraints = quad::interpolate_boundary(&mesh, &d, exact).unwrap();
    let (a, b) = a.with_dirichlet(&rhs, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (u, _) = conjugate_gradient(&a, &b, options).unwrap();
    let (l2, h1) = quad::errors(&mesh, &u, exact, |_| [2., -3.]).unwrap();
    assert!(l2 < 1e-9 && h1 < 1e-9, "L2={} H1={}", l2, h1);
}
