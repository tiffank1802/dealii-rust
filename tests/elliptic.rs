#![allow(clippy::needless_range_loop)] // Indexed analytical local matrices.

use dealii_rust::{
    base::Point2,
    elliptic, field, gmsh,
    grid::{BoundaryEdgeId, CellId, Mesh2},
    integrators::{diffusion_reaction_p1, DiffusionTensor},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn tensor(xx: f64, xy: f64, yy: f64) -> DiffusionTensor {
    DiffusionTensor { xx, xy, yy }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn reference() -> [Point2; 3] {
    [p(0., 0.), p(1., 0.), p(0., 1.)]
}

#[test]
fn local_constant_tensor_reaction_source_matches_analytical_matrix() {
    let local = diffusion_reaction_p1(reference(), |_| tensor(2., 1., 3.), |_| 4., |_| 6.).unwrap();
    let stiffness = [[3.5, -1.5, -2.], [-1.5, 1., 0.5], [-2., 0.5, 1.5]];
    for i in 0..3 {
        close(local.rhs[i], 1., 1e-13);
        for j in 0..3 {
            let mass = if i == j { 1. / 3. } else { 1. / 6. };
            close(local.matrix[i][j], stiffness[i][j] + mass, 1e-13);
            assert_eq!(local.matrix[i][j], local.matrix[j][i]);
        }
    }
}

#[test]
fn variable_polynomial_data_integrates_at_expected_order() {
    let local = diffusion_reaction_p1(
        reference(),
        |p| tensor(1. + p.x, 0., 1. + p.y),
        |_| 0.,
        |p| p.x + p.y,
    )
    .unwrap();
    // ∫(1+x) = ∫(1+y) = 2/3 on the reference triangle.
    let expected = [
        [4. / 3., -2. / 3., -2. / 3.],
        [-2. / 3., 2. / 3., 0.],
        [-2. / 3., 0., 2. / 3.],
    ];
    for i in 0..3 {
        close(local.rhs[i], [1. / 12., 1. / 8., 1. / 8.][i], 1e-13);
        for j in 0..3 {
            close(local.matrix[i][j], expected[i][j], 1e-13);
        }
    }
}

#[test]
fn elliptic_identity_recovers_existing_poisson_assembly() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let source = |p: Point2| 2. + p.x * p.x - p.y;
    let (a, b) = elliptic::assemble(
        &mesh,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, p| source(p),
    )
    .unwrap();
    let (reference, reference_b) = numerics::assemble_poisson(&mesh, source).unwrap();
    for i in 0..mesh.vertices().len() {
        close(b[i], reference_b[i], 1e-13);
        for j in 0..mesh.vertices().len() {
            close(a.get(i, j).unwrap(), reference.get(i, j).unwrap(), 1e-13);
        }
    }
}

#[test]
fn physical_groups_drive_cellwise_materials() {
    let mesh = gmsh::parse(include_str!("fixtures/square41.msh")).unwrap();
    let m = mesh.mesh();
    let (a, _) = elliptic::assemble(
        m,
        |id, _| {
            assert_eq!(mesh.cell_physical_tags(id).unwrap(), &[100, 101]);
            // Distinct material per cell exercises the ID passed to coefficients.
            if id.0 % 2 == 0 {
                tensor(2., 0., 2.)
            } else {
                tensor(3., 0., 3.)
            }
        },
        |_, _| 0.,
        |_, _| 0.,
    )
    .unwrap();
    let mut expected = vec![vec![0.; m.vertices().len()]; m.vertices().len()];
    for (i, cell) in m.cells().iter().enumerate() {
        let d = if i % 2 == 0 { 2. } else { 3. };
        let local = diffusion_reaction_p1(
            m.cell_points(CellId(i)).unwrap(),
            |_| tensor(d, 0., d),
            |_| 0.,
            |_| 0.,
        )
        .unwrap();
        for j in 0..3 {
            for k in 0..3 {
                expected[cell[j].0][cell[k].0] += local.matrix[j][k];
            }
        }
    }
    for (i, row) in expected.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            close(a.get(i, j).unwrap(), *value, 1e-13);
        }
    }
    for x in a.mul(&vec![1.; m.vertices().len()]).unwrap() {
        close(x, 0., 1e-13);
    }
}

fn solution(p: Point2) -> f64 {
    1. + p.x + p.x * p.x + 2. * p.y + p.y * p.y
}
fn gradient(p: Point2) -> [f64; 2] {
    [1. + 2. * p.x, 2. + 2. * p.y]
}
fn diffusion(p: Point2) -> DiffusionTensor {
    tensor(2. + p.x, 0.5, 3. + p.y)
}
fn reaction(p: Point2) -> f64 {
    1. + p.x + p.y
}
fn flux(p: Point2, n: [f64; 2]) -> f64 {
    let g = gradient(p);
    let a = diffusion(p);
    n[0] * (a.xx * g[0] + a.xy * g[1]) + n[1] * (a.xy * g[0] + a.yy * g[1])
}
fn source(p: Point2) -> f64 {
    -(13. + 4. * p.x + 4. * p.y) + reaction(p) * solution(p)
}
fn mixed(n: usize) -> (f64, f64) {
    let mesh = Mesh2::unit_square(n).unwrap();
    let (a, mut b) = elliptic::assemble(
        &mesh,
        |_, p| diffusion(p),
        |_, p| reaction(p),
        |_, p| source(p),
    )
    .unwrap();
    let mut left = Vec::new();
    let mut other = Vec::new();
    for i in 0..mesh.boundary_edges().len() {
        let edge = mesh.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
        if edge.points.iter().all(|p| p.x == 0.) {
            left.push(BoundaryEdgeId(i));
        } else {
            other.push(BoundaryEdgeId(i));
        }
    }
    let load = numerics::assemble_neumann(&mesh, &other, flux).unwrap();
    for (value, g) in b.iter_mut().zip(load) {
        *value += g;
    }
    let constraints = field::interpolate_boundary(&mesh, &left, solution).unwrap();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (u, _) = conjugate_gradient(&a, &b, options).unwrap();
    numerics::errors(&mesh, &u, solution, gradient).unwrap()
}

#[test]
fn anisotropic_variable_diffusion_reaction_mixed_problem_converges() {
    let results: Vec<_> = [8, 16, 32].iter().map(|&n| mixed(n)).collect();
    for pair in results.windows(2) {
        let l2_rate = (pair[0].0 / pair[1].0).log2();
        let h1_rate = (pair[0].1 / pair[1].1).log2();
        assert!((1.8..2.2).contains(&l2_rate), "L2 rate {}", l2_rate);
        assert!((0.9..1.1).contains(&h1_rate), "H1 rate {}", h1_rate);
    }
}

#[test]
fn positive_reaction_all_neumann_does_not_require_poisson_normalization() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let u = |p: Point2| 1. + p.x + 2. * p.y;
    let a = tensor(2., 0.5, 3.);
    let (matrix, mut rhs) =
        elliptic::assemble(&mesh, |_, _| a, |_, _| 2., |_, p| 2. * u(p)).unwrap();
    let edges: Vec<_> = (0..mesh.boundary_edges().len())
        .map(BoundaryEdgeId)
        .collect();
    let load = numerics::assemble_neumann(&mesh, &edges, |_, n| n[0] * 3. + n[1] * 6.5).unwrap();
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
    }
    let (solution, _) = conjugate_gradient(&matrix, &rhs, CgOptions::default()).unwrap();
    for (v, point) in solution.iter().zip(mesh.vertices()) {
        close(*v, u(*point), 1e-8);
    }
}

#[test]
fn invalid_coefficients_and_arithmetic_are_rejected() {
    let bad = [
        tensor(0., 0., 1.),
        tensor(-1., 0., 1.),
        tensor(1., 1., 1.),
        tensor(1., 2., 1.),
        tensor(1., f64::NAN, 1.),
        tensor(f64::INFINITY, 0., 1.),
    ];
    for a in bad {
        assert!(diffusion_reaction_p1(reference(), |_| a, |_| 0., |_| 0.).is_err());
    }
    for c in [-1., f64::NAN, f64::INFINITY] {
        assert!(
            diffusion_reaction_p1(reference(), |_| DiffusionTensor::identity(), |_| c, |_| 0.)
                .is_err()
        );
    }
    assert!(diffusion_reaction_p1(
        reference(),
        |_| DiffusionTensor::identity(),
        |_| 0.,
        |_| f64::NAN
    )
    .is_err());
    assert!(diffusion_reaction_p1(
        [reference()[0], reference()[2], reference()[1]],
        |_| DiffusionTensor::identity(),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    assert!(diffusion_reaction_p1(
        reference(),
        |_| tensor(f64::MAX, 0., f64::MAX),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    let mesh = Mesh2::unit_square(1).unwrap();
    let err = elliptic::assemble(
        &mesh,
        |id, _| {
            if id.0 == 1 {
                tensor(-1., 0., 1.)
            } else {
                DiffusionTensor::identity()
            }
        },
        |_, _| 0.,
        |_, _| 0.,
    )
    .unwrap_err();
    assert!(err.to_string().contains("cell 1"));
}

#[test]
fn linear_reaction_contributes_a_consistent_weighted_mass_matrix() {
    let base = diffusion_reaction_p1(reference(), |_| DiffusionTensor::identity(), |_| 0., |_| 0.)
        .unwrap();
    let with_reaction = diffusion_reaction_p1(
        reference(),
        |_| DiffusionTensor::identity(),
        |p| p.x + p.y,
        |_| 0.,
    )
    .unwrap();
    let expected = [
        [1. / 30., 1. / 40., 1. / 40.],
        [1. / 40., 1. / 15., 1. / 30.],
        [1. / 40., 1. / 30., 1. / 15.],
    ];
    for i in 0..3 {
        for j in 0..3 {
            close(
                with_reaction.matrix[i][j] - base.matrix[i][j],
                expected[i][j],
                1e-13,
            );
        }
    }
}

#[test]
fn positivity_is_checked_at_quadrature_samples_and_not_elsewhere() {
    // The first quadrature point is the barycenter, where xx=0; other
    // points can be positive. A negative value only at a vertex is not sampled.
    assert!(diffusion_reaction_p1(
        reference(),
        |p| tensor((p.x - 1. / 3.).abs() + (p.y - 1. / 3.).abs(), 0., 1.),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    assert!(diffusion_reaction_p1(
        reference(),
        |p| tensor(if p.x == 0. && p.y == 0. { -1. } else { 1. }, 0., 1.),
        |_| 0.,
        |_| 0.
    )
    .is_ok());
}

#[test]
fn local_assembly_respects_cyclic_vertex_permutation() {
    let points = [p(1., 2.), p(3., 2.), p(2., 5.)];
    let diffusion = |p: Point2| tensor(2. + p.x, 0.25, 3. + p.y);
    let reaction = |p: Point2| 1. + p.x;
    let source = |p: Point2| p.x * p.y;
    let a = diffusion_reaction_p1(points, diffusion, reaction, source).unwrap();
    let b = diffusion_reaction_p1(
        [points[1], points[2], points[0]],
        diffusion,
        reaction,
        source,
    )
    .unwrap();
    for i in 0..3 {
        close(b.rhs[i], a.rhs[(i + 1) % 3], 1e-12);
        for j in 0..3 {
            close(b.matrix[i][j], a.matrix[(i + 1) % 3][(j + 1) % 3], 1e-12);
        }
    }
}

#[test]
fn strongly_scaled_positive_definite_tensor_is_not_rejected_by_test_overflow() {
    let local =
        diffusion_reaction_p1(reference(), |_| tensor(1e-200, 0.5, 1e200), |_| 0., |_| 0.).unwrap();
    assert!(local.matrix.iter().flatten().all(|x| x.is_finite()));
    assert!(
        diffusion_reaction_p1(reference(), |_| tensor(1e-200, 2., 1e200), |_| 0., |_| 0.).is_err()
    );
}
