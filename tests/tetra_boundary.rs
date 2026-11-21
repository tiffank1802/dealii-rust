use dealii_rust::{
    base::Point3,
    grid::VertexId,
    lac::{conjugate_gradient, CgOptions},
    tetra::{self, BoundaryFaceId, TetMesh},
};

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn faces(mesh: &TetMesh) -> Vec<BoundaryFaceId> {
    (0..mesh.boundary_faces().len())
        .map(BoundaryFaceId)
        .collect()
}

#[test]
fn linear_flux_on_skew_face_has_analytic_nodal_loads() {
    let mesh = TetMesh::new(
        vec![p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let id = faces(&mesh)
        .into_iter()
        .find(|&id| !mesh.boundary_faces()[id.0].contains(&VertexId(0)))
        .unwrap();
    let area = 3.0_f64.sqrt() / 2.0;
    let load = tetra::assemble_neumann(&mesh, &[id], |p, n| {
        for &normal in &n {
            near(normal, 1.0 / 3.0_f64.sqrt(), 1e-14);
        }
        p.x
    })
    .unwrap();
    near(load[0], 0., 1e-14);
    near(load[1], area / 6., 1e-14);
    near(load[2], area / 12., 1e-14);
    near(load[3], area / 12., 1e-14);
    near(load.iter().sum(), area / 3., 1e-14);
}

#[test]
fn face_flux_and_boundary_data_validate_inputs() {
    let mesh = TetMesh::unit_cube(1).unwrap();
    let all = faces(&mesh);
    assert_eq!(
        tetra::assemble_neumann(&mesh, &[], |_, _| f64::NAN).unwrap(),
        vec![0.; 8]
    );
    assert!(tetra::assemble_neumann(&mesh, &[all[0], all[0]], |_, _| 1.).is_err());
    assert!(tetra::assemble_neumann(&mesh, &[BoundaryFaceId(all.len())], |_, _| 1.).is_err());
    assert!(tetra::assemble_neumann(&mesh, &all, |_, _| f64::NAN).is_err());
    assert!(tetra::assemble_neumann(&mesh, &all, |_, _| f64::INFINITY).is_err());
    assert_eq!(
        tetra::interpolate_boundary(&mesh, &[], |_| f64::NAN).unwrap(),
        Vec::new()
    );
    assert!(tetra::interpolate_boundary(&mesh, &[all[0], all[0]], |_| 1.).is_err());
    assert!(tetra::interpolate_boundary(&mesh, &[BoundaryFaceId(all.len())], |_| 1.).is_err());
    assert!(tetra::interpolate_boundary(&mesh, &all, |_| f64::NAN).is_err());
    let constraints =
        tetra::interpolate_boundary(&mesh, &all, |p| 1. + p.x + 2. * p.y + 3. * p.z).unwrap();
    assert_eq!(constraints.len(), 8);
    for &(id, value) in &constraints {
        let pt = mesh.vertices()[id];
        near(value, 1. + pt.x + 2. * pt.y + 3. * pt.z, 1e-14);
    }
    assert!(constraints.windows(2).all(|w| w[0].0 < w[1].0));
    let large = TetMesh::new(
        vec![
            p(0., 0., 0.),
            p(10., 0., 0.),
            p(0., 10., 0.),
            p(0., 0., 10.),
        ],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(tetra::assemble_neumann(&large, &faces(&large), |_, _| f64::MAX).is_err());
}

#[test]
fn all_face_integrals_have_expected_balance_and_laplace_kernel() {
    let mesh = TetMesh::unit_cube(3).unwrap();
    let all = faces(&mesh);
    let boundary = tetra::assemble_neumann(&mesh, &all, |_, _| 1.).unwrap();
    near(boundary.iter().sum(), 6., 1e-12);
    let (a, source) = tetra::assemble_poisson(&mesh, |_| -6.).unwrap();
    let flux = tetra::assemble_neumann(&mesh, &all, |p, n| {
        2. * (p.x * n[0] + p.y * n[1] + p.z * n[2])
    })
    .unwrap();
    near(
        source.iter().sum::<f64>() + flux.iter().sum::<f64>(),
        0.,
        2e-12,
    );
    for x in a.mul(&vec![1.; mesh.vertices().len()]).unwrap() {
        near(x, 0., 1e-13);
    }
    // Pure Neumann remains singular; no direct CG solve is attempted.
}

fn exact(p: Point3, quadratic: bool) -> f64 {
    1. + p.x
        + 2. * p.y
        + 3. * p.z
        + if quadratic {
            p.x * p.x + p.y * p.y + p.z * p.z
        } else {
            0.
        }
}
fn grad(p: Point3, quadratic: bool) -> [f64; 3] {
    if quadratic {
        [1. + 2. * p.x, 2. + 2. * p.y, 3. + 2. * p.z]
    } else {
        [1., 2., 3.]
    }
}
fn mixed(n: usize, quadratic: bool) -> (f64, f64) {
    let mesh = TetMesh::unit_cube(n).unwrap();
    let mut left = Vec::new();
    let mut other = Vec::new();
    for id in faces(&mesh) {
        let face = mesh.boundary_face_geometry(id).unwrap();
        if face.points.iter().all(|p| p.x == 0.) {
            left.push(id)
        } else {
            other.push(id)
        }
    }
    let (a, mut rhs) =
        tetra::assemble_poisson(&mesh, |_| if quadratic { -6. } else { 0. }).unwrap();
    let neumann = tetra::assemble_neumann(&mesh, &other, |p, n| {
        let g = grad(p, quadratic);
        g[0] * n[0] + g[1] * n[1] + g[2] * n[2]
    })
    .unwrap();
    for (value, load) in rhs.iter_mut().zip(neumann) {
        *value += load;
    }
    let constraints = tetra::interpolate_boundary(&mesh, &left, |p| exact(p, quadratic)).unwrap();
    let (a, rhs) = a.with_dirichlet(&rhs, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    let (u, _) = conjugate_gradient(&a, &rhs, options).unwrap();
    tetra::errors(&mesh, &u, |p| exact(p, quadratic), |p| grad(p, quadratic)).unwrap()
}
#[test]
fn mixed_problem_recovers_affine_solution() {
    for n in [1, 2, 4] {
        let (l2, h1) = mixed(n, false);
        assert!(l2 < 1e-9 && h1 < 1e-9, "L2={} H1={}", l2, h1);
    }
}
#[test]
fn mixed_quadratic_problem_converges() {
    let results: Vec<_> = [4, 8, 16].iter().map(|&n| mixed(n, true)).collect();
    for pair in results.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
        assert!((0.8..1.2).contains(&h1), "H1 rate {}", h1);
    }
}
