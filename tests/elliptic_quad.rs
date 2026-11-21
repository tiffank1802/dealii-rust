#![allow(clippy::needless_range_loop)] // Indexed reference matrices.
use dealii_rust::{
    base::Point2,
    elliptic_quad::{self, diffusion_reaction_q1},
    grid::{BoundaryEdgeId, CellId, VertexId},
    integrators::DiffusionTensor,
    lac::{conjugate_gradient, CgOptions},
    quad::{self, QuadMesh},
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn tensor(xx: f64, xy: f64, yy: f64) -> DiffusionTensor {
    DiffusionTensor { xx, xy, yy }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn square() -> [Point2; 4] {
    [p(0., 0.), p(1., 0.), p(1., 1.), p(0., 1.)]
}

#[test]
fn analytical_anisotropic_diffusion_reaction_and_load_on_square() {
    let s = diffusion_reaction_q1(square(), |_| tensor(2., 1., 3.), |_| 4., |_| 6.).unwrap();
    let expected = [
        [47. / 18., 1. / 18., -11. / 9., -4. / 9.],
        [1. / 18., 29. / 18., -4. / 9., -2. / 9.],
        [-11. / 9., -4. / 9., 47. / 18., 1. / 18.],
        [-4. / 9., -2. / 9., 1. / 18., 29. / 18.],
    ];
    for i in 0..4 {
        near(s.rhs[i], 1.5, 1e-13);
        for j in 0..4 {
            near(s.matrix[i][j], expected[i][j], 1e-13);
            assert_eq!(s.matrix[i][j], s.matrix[j][i]);
        }
    }
}
#[test]
fn poisson_identity_matches_existing_on_nonaffine_mesh() {
    let base = QuadMesh::unit_square(3).unwrap();
    let points = base
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
        })
        .collect();
    let mesh = QuadMesh::new(points, base.cells().to_vec()).unwrap();
    let source = |p: Point2| 2. + p.x * p.x - p.y;
    let (a, b) = elliptic_quad::assemble(
        &mesh,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, p| source(p),
    )
    .unwrap();
    let (expected, load) = quad::assemble_poisson(&mesh, source).unwrap();
    for i in 0..mesh.vertices().len() {
        near(b[i], load[i], 1e-13);
        for j in 0..mesh.vertices().len() {
            near(a.get(i, j).unwrap(), expected.get(i, j).unwrap(), 1e-13);
        }
    }
}
#[test]
fn nonaffine_cyclic_permutation_preserves_local_system() {
    let points = [p(0., 0.), p(2., 0.), p(2.5, 1.), p(-0.2, 1.2)];
    let d = |p: Point2| tensor(2. + p.x, 0.2, 3. + p.y);
    let c = |p: Point2| 1. + p.y;
    let f = |p: Point2| p.x * p.y;
    let a = diffusion_reaction_q1(points, d, c, f).unwrap();
    let b = diffusion_reaction_q1([points[1], points[2], points[3], points[0]], d, c, f).unwrap();
    for i in 0..4 {
        near(b.rhs[i], a.rhs[(i + 1) % 4], 1e-12);
        for j in 0..4 {
            near(b.matrix[i][j], a.matrix[(i + 1) % 4][(j + 1) % 4], 1e-12);
        }
    }
}
#[test]
fn cell_ids_drive_coefficients_and_assembled_matrix() {
    let mesh = QuadMesh::unit_square(2).unwrap();
    let (a, b) = elliptic_quad::assemble(
        &mesh,
        |id, _| tensor(1. + id.0 as f64, 0., 1. + id.0 as f64),
        |id, _| id.0 as f64,
        |id, _| 1. + id.0 as f64,
    )
    .unwrap();
    let mut expected = vec![vec![0.; mesh.vertices().len()]; mesh.vertices().len()];
    let mut load = vec![0.; mesh.vertices().len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let v = i as f64;
        let local = diffusion_reaction_q1(
            mesh.cell_points(CellId(i)).unwrap(),
            |_| tensor(1. + v, 0., 1. + v),
            |_| v,
            |_| 1. + v,
        )
        .unwrap();
        for j in 0..4 {
            load[cell[j].0] += local.rhs[j];
            for k in 0..4 {
                expected[cell[j].0][cell[k].0] += local.matrix[j][k];
            }
        }
    }
    for i in 0..b.len() {
        near(b[i], load[i], 1e-13);
        for j in 0..b.len() {
            near(a.get(i, j).unwrap(), expected[i][j], 1e-13);
        }
    }
}
fn exact(p: Point2) -> f64 {
    1. + p.x + p.x * p.x + 2. * p.y + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [1. + 2. * p.x, 2. + 2. * p.y]
}
fn diffusion(p: Point2) -> DiffusionTensor {
    tensor(2. + p.x, 0.5, 3. + p.y)
}
fn flux(p: Point2, n: [f64; 2]) -> f64 {
    let a = diffusion(p);
    let g = gradient(p);
    n[0] * (a.xx * g[0] + a.xy * g[1]) + n[1] * (a.xy * g[0] + a.yy * g[1])
}
fn mixed(n: usize, distorted: bool) -> (f64, f64) {
    let base = QuadMesh::unit_square(n).unwrap();
    let mesh = if distorted {
        let pts = base
            .vertices()
            .iter()
            .map(|p| Point2 {
                x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
                y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
            })
            .collect();
        QuadMesh::new(pts, base.cells().to_vec()).unwrap()
    } else {
        base
    };
    let (a, mut b) = elliptic_quad::assemble(
        &mesh,
        |_, p| diffusion(p),
        |_, p| 1. + p.x + p.y,
        |_, p| -(13. + 4. * p.x + 4. * p.y) + (1. + p.x + p.y) * exact(p),
    )
    .unwrap();
    let mut d = vec![];
    let mut neumann = vec![];
    for i in 0..mesh.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        if mesh
            .boundary_edge_geometry(id)
            .unwrap()
            .points
            .iter()
            .all(|p| p.x == 0.)
        {
            d.push(id)
        } else {
            neumann.push(id)
        }
    }
    let g = quad::assemble_neumann(&mesh, &neumann, flux).unwrap();
    for (b, g) in b.iter_mut().zip(g) {
        *b += g;
    }
    let constraints = quad::interpolate_boundary(&mesh, &d, exact).unwrap();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let (u, _) = conjugate_gradient(
        &a,
        &b,
        CgOptions {
            relative_tolerance: 1e-12,
            absolute_tolerance: 1e-13,
            ..CgOptions::default()
        },
    )
    .unwrap();
    quad::errors(&mesh, &u, exact, gradient).unwrap()
}
#[test]
fn anisotropic_variable_mixed_problem_converges_on_both_meshes() {
    for distorted in [false, true] {
        let errors: Vec<_> = [8, 16, 32].iter().map(|&n| mixed(n, distorted)).collect();
        for pair in errors.windows(2) {
            let l2 = (pair[0].0 / pair[1].0).log2();
            let h1 = (pair[0].1 / pair[1].1).log2();
            assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
            assert!((0.9..1.1).contains(&h1), "H1 rate {}", h1);
        }
    }
}
#[test]
fn positive_reaction_all_neumann_recovers_affine_solution_on_trapezoid() {
    let mesh = QuadMesh::new(
        vec![p(0., 0.), p(2., 0.), p(2.5, 1.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let u = |p: Point2| 1. + p.x + 2. * p.y;
    let (a, mut b) = elliptic_quad::assemble(
        &mesh,
        |_, _| tensor(2., 0.5, 3.),
        |_, _| 2.,
        |_, p| 2. * u(p),
    )
    .unwrap();
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let g = quad::assemble_neumann(&mesh, &edges, |_, n| 3. * n[0] + 6.5 * n[1]).unwrap();
    for (b, g) in b.iter_mut().zip(g) {
        *b += g;
    }
    let (result, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    for (v, p) in result.iter().zip(mesh.vertices()) {
        near(*v, u(*p), 1e-10);
    }
}
#[test]
fn invalid_sampled_data_geometry_overflow_and_cell_error() {
    for a in [
        tensor(0., 0., 1.),
        tensor(1., 1., 1.),
        tensor(1., 2., 1.),
        tensor(f64::NAN, 0., 1.),
        tensor(1., 0., f64::INFINITY),
    ] {
        assert!(diffusion_reaction_q1(square(), |_| a, |_| 0., |_| 0.).is_err());
    }
    for c in [-1., f64::NAN, f64::INFINITY] {
        assert!(
            diffusion_reaction_q1(square(), |_| DiffusionTensor::identity(), |_| c, |_| 0.)
                .is_err()
        );
    }
    assert!(diffusion_reaction_q1(
        square(),
        |_| DiffusionTensor::identity(),
        |_| 0.,
        |_| f64::NAN
    )
    .is_err());
    assert!(diffusion_reaction_q1(
        [square()[0], square()[3], square()[2], square()[1]],
        |_| DiffusionTensor::identity(),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    assert!(
        diffusion_reaction_q1(square(), |_| tensor(f64::MAX, 0., f64::MAX), |_| 0., |_| 0.)
            .is_err()
    );
    let mesh = QuadMesh::unit_square(2).unwrap();
    let err = elliptic_quad::assemble(
        &mesh,
        |id, _| {
            if id.0 == 2 {
                tensor(-1., 0., 1.)
            } else {
                DiffusionTensor::identity()
            }
        },
        |_, _| 0.,
        |_, _| 0.,
    )
    .unwrap_err();
    assert!(err.to_string().contains("cell 2"));
}
#[test]
fn sampled_validation_does_not_guarantee_positivity_between_points() {
    assert!(diffusion_reaction_q1(
        square(),
        |p| tensor(if p.x == 0. && p.y == 0. { -1. } else { 1. }, 0., 1.),
        |_| 0.,
        |_| 0.
    )
    .is_ok());
    assert!(diffusion_reaction_q1(
        square(),
        |p| tensor(if p.x == 0.5 && p.y == 0.5 { 0. } else { 1. }, 0., 1.),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    assert!(
        diffusion_reaction_q1(square(), |_| tensor(1e-200, 0.5, 1e200), |_| 0., |_| 0.).is_ok()
    );
}
