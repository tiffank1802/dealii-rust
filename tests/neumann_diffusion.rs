use dealii_rust::{
    base::{Point2, Point3},
    elliptic_tetra::DiffusionTensor3,
    grid::{CellId, Mesh2, VertexId},
    integrators::DiffusionTensor,
    lac::CgOptions,
    neumann::{self, NeumannOptions},
    neumann3d, neumann_quad, numerics,
    quad::{self, QuadMesh},
    tetra::{self, TetMesh},
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn tensor2(p: Point2) -> DiffusionTensor {
    DiffusionTensor {
        xx: 2. + p.x,
        xy: 0.5,
        yy: 3. + p.y,
    }
}
fn source2(p: Point2) -> f64 {
    -10. - 4. * (p.x + p.y)
}
fn flux2(p: Point2, n: [f64; 2]) -> f64 {
    let a = tensor2(p);
    let x = 2. * p.x;
    let y = 2. * p.y;
    n[0] * (a.xx * x + a.xy * y) + n[1] * (a.xy * x + a.yy * y)
}
fn exact2(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y - 2. / 3.
}
fn grad2(p: Point2) -> [f64; 2] {
    [2. * p.x, 2. * p.y]
}
fn tensor3(p: Point3) -> DiffusionTensor3 {
    DiffusionTensor3 {
        xx: 2. + p.x,
        xy: 0.5,
        xz: 0.25,
        yy: 3. + p.y,
        yz: 0.75,
        zz: 4. + p.z,
    }
}
fn source3(p: Point3) -> f64 {
    -18. - 4. * (p.x + p.y + p.z)
}
fn flux3(p: Point3, n: [f64; 3]) -> f64 {
    let a = tensor3(p);
    let x = 2. * p.x;
    let y = 2. * p.y;
    let z = 2. * p.z;
    n[0] * (a.xx * x + a.xy * y + a.xz * z)
        + n[1] * (a.xy * x + a.yy * y + a.yz * z)
        + n[2] * (a.xz * x + a.yz * y + a.zz * z)
}
fn exact3(p: Point3) -> f64 {
    p.x * p.x + p.y * p.y + p.z * p.z - 1.
}
fn grad3(p: Point3) -> [f64; 3] {
    [2. * p.x, 2. * p.y, 2. * p.z]
}
fn check_rates(values: &[(f64, f64)]) {
    for pair in values.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
        assert!((0.8..1.2).contains(&h1), "H1 rate {}", h1);
    }
}
#[test]
fn triangle_variable_anisotropic_problem_converges() {
    let errors: Vec<_> = [8, 16, 32]
        .iter()
        .map(|&n| {
            let mesh = Mesh2::unit_square(n).unwrap();
            let r = neumann::solve_diffusion(
                &mesh,
                |_, p| tensor2(p),
                |_, p| source2(p),
                flux2,
                NeumannOptions::default(),
            )
            .unwrap();
            assert!(r.original_residual_norm < 1e-7);
            numerics::errors(&mesh, &r.values, exact2, grad2).unwrap()
        })
        .collect();
    check_rates(&errors);
}
#[test]
fn quadrilateral_variable_anisotropic_problem_converges_on_deformed_mesh() {
    let errors: Vec<_> = [8, 16, 32]
        .iter()
        .map(|&n| {
            let base = QuadMesh::unit_square(n).unwrap();
            let points = base
                .vertices()
                .iter()
                .map(|p| Point2 {
                    x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
                    y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
                })
                .collect();
            let mesh = QuadMesh::new(points, base.cells().to_vec()).unwrap();
            let r = neumann_quad::solve_diffusion(
                &mesh,
                |_, p| tensor2(p),
                |_, p| source2(p),
                flux2,
                NeumannOptions::default(),
            )
            .unwrap();
            assert!(r.original_residual_norm < 1e-7);
            quad::errors(&mesh, &r.values, exact2, grad2).unwrap()
        })
        .collect();
    check_rates(&errors);
}
#[test]
fn tetrahedron_variable_anisotropic_problem_converges() {
    let errors: Vec<_> = [4, 8, 16]
        .iter()
        .map(|&n| {
            let mesh = TetMesh::unit_cube(n).unwrap();
            let r = neumann3d::solve_diffusion(
                &mesh,
                |_, p| tensor3(p),
                |_, p| source3(p),
                flux3,
                NeumannOptions::default(),
            )
            .unwrap();
            assert!(r.original_residual_norm < 1e-7);
            tetra::errors(&mesh, &r.values, exact3, grad3).unwrap()
        })
        .collect();
    check_rates(&errors);
}
fn two_triangles() -> Mesh2 {
    Mesh2::new(
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
    .unwrap()
}
fn two_quads() -> QuadMesh {
    QuadMesh::new(
        vec![
            p(0., 0.),
            p(2., 0.),
            p(1., 1.),
            p(0., 1.),
            p(3., 0.),
            p(4., 0.),
            p(4., 1.),
            p(3., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap()
}
fn two_tets() -> TetMesh {
    TetMesh::new(
        vec![
            p3(0., 0., 0.),
            p3(1., 0., 0.),
            p3(0., 1., 0.),
            p3(0., 0., 1.),
            p3(3., 0., 0.),
            p3(4., 0., 0.),
            p3(3., 1., 0.),
            p3(3., 0., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap()
}
fn cell_tensor2(id: CellId) -> DiffusionTensor {
    let scale = 1. + id.0 as f64;
    DiffusionTensor {
        xx: 2. * scale,
        xy: 0.5 * scale,
        yy: 3. * scale,
    }
}
fn cell_flux2(pt: Point2, n: [f64; 2]) -> f64 {
    let s = if pt.x < 2.5 { 1. } else { 2. };
    s * (3. * n[0] + 6.5 * n[1])
}
#[test]
fn cellwise_materials_and_disconnected_components_triangle_and_quad() {
    let tri = two_triangles();
    let r = neumann::solve_diffusion(
        &tri,
        |id, _| cell_tensor2(id),
        |_, _| 0.,
        cell_flux2,
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(r.components.len(), 2);
    for (i, u) in r.values.iter().enumerate() {
        near(*u, [-1., 0., 1.][i % 3], 1e-12);
    }
    assert!(r.original_residual_norm < 1e-12);
    let quad = two_quads();
    let r = neumann_quad::solve_diffusion(
        &quad,
        |id, _| cell_tensor2(id),
        |_, _| 0.,
        cell_flux2,
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(r.components.len(), 2);
    for (i, pt) in quad.vertices().iter().enumerate() {
        near(
            r.values[i],
            pt.x + 2. * pt.y - if i < 4 { 5. / 3. } else { 4.5 },
            2e-12,
        );
    }
    near(
        (5. * r.values[0] + 5. * r.values[1] + 4. * r.values[2] + 4. * r.values[3]) / 12.,
        0.,
        1e-12,
    );
    assert!(r.original_residual_norm < 1e-12);
}
#[test]
fn discontinuous_coefficients_with_continuous_conormal_interface() {
    let mesh = QuadMesh::new(
        vec![
            p(0., 0.),
            p(1., 0.),
            p(2., 0.),
            p(0., 1.),
            p(1., 1.),
            p(2., 1.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(4), VertexId(3)],
            [VertexId(1), VertexId(2), VertexId(5), VertexId(4)],
        ],
    )
    .unwrap();
    // u=y has the same normal flux at the interface although A_xx jumps.
    let r = neumann_quad::solve_diffusion(
        &mesh,
        |id, _| DiffusionTensor {
            xx: 2. + id.0 as f64,
            xy: 0.,
            yy: 3.,
        },
        |_, _| 0.,
        |_, n| 3. * n[1],
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(r.components.len(), 1);
    for (u, point) in r.values.iter().zip(mesh.vertices()) {
        near(*u, point.y - 0.5, 1e-12);
    }
    assert!(r.original_residual_norm < 1e-12);
}
#[test]
fn cellwise_materials_and_disconnected_components_tetra() {
    let mesh = two_tets();
    let r = neumann3d::solve_diffusion(
        &mesh,
        |id, _| {
            let s = 1. + id.0 as f64;
            DiffusionTensor3 {
                xx: 2. * s,
                xy: 0.5 * s,
                xz: 0.25 * s,
                yy: 3. * s,
                yz: 0.75 * s,
                zz: 4. * s,
            }
        },
        |_, _| 0.,
        |p, n| {
            let s = if p.x < 2. { 1. } else { 2. };
            s * (3.75 * n[0] + 8.75 * n[1] + 13.75 * n[2])
        },
        NeumannOptions::default(),
    )
    .unwrap();
    assert_eq!(r.components.len(), 2);
    for (i, u) in r.values.iter().enumerate() {
        near(*u, [-1.5, -0.5, 0.5, 1.5][i % 4], 1e-12);
    }
    assert!(r.original_residual_norm < 1e-12);
}
#[test]
fn rejects_locally_incompatible_load_even_when_global_balance_is_zero() {
    let tri = two_triangles();
    assert!(neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |id, _| if id.0 == 0 { 1. } else { -1. },
        |_, _| 0.,
        NeumannOptions::default()
    )
    .unwrap_err()
    .to_string()
    .contains("incompatible Neumann"));
    let quad = two_quads();
    assert!(neumann_quad::solve_diffusion(
        &quad,
        |_, _| DiffusionTensor::identity(),
        |id, _| if id.0 == 0 { 2. } else { -3. },
        |_, _| 0.,
        NeumannOptions::default()
    )
    .unwrap_err()
    .to_string()
    .contains("incompatible Neumann"));
    let tet = two_tets();
    assert!(neumann3d::solve_diffusion(
        &tet,
        |_, _| DiffusionTensor3::identity(),
        |id, _| if id.0 == 0 { 1. } else { -1. },
        |_, _| 0.,
        NeumannOptions::default()
    )
    .unwrap_err()
    .to_string()
    .contains("incompatible Neumann"));
}
#[test]
fn accepted_imbalance_is_not_projected() {
    let options = NeumannOptions {
        compatibility_absolute_tolerance: 1e-5,
        compatibility_relative_tolerance: 0.,
        ..NeumannOptions::default()
    };
    let tri = Mesh2::unit_square(2).unwrap();
    let a = neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |_, _| 1e-6,
        |_, _| 0.,
        options,
    )
    .unwrap();
    near(a.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(a.original_residual_norm > 9e-7 && a.cg.residual_norm < 1e-12);
    let quad = QuadMesh::unit_square(2).unwrap();
    let b = neumann_quad::solve_diffusion(
        &quad,
        |_, _| DiffusionTensor::identity(),
        |_, _| 1e-6,
        |_, _| 0.,
        options,
    )
    .unwrap();
    near(b.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(b.original_residual_norm > 9e-7 && b.cg.residual_norm < 1e-12);
    let tet = TetMesh::unit_cube(2).unwrap();
    let c = neumann3d::solve_diffusion(
        &tet,
        |_, _| DiffusionTensor3::identity(),
        |_, _| 1e-6,
        |_, _| 0.,
        options,
    )
    .unwrap();
    near(c.components[0].load_imbalance, 1e-6, 1e-16);
    assert!(c.original_residual_norm > 9e-7 && c.cg.residual_norm < 1e-12);
}
#[test]
fn invalid_tensor_flux_source_tolerance_and_failed_cg_are_rejected() {
    let tri = Mesh2::unit_square(2).unwrap();
    let quad = QuadMesh::unit_square(2).unwrap();
    let tet = TetMesh::unit_cube(2).unwrap();
    let invalid = DiffusionTensor {
        xx: 1.,
        xy: 2.,
        yy: 1.,
    };
    let invalid3 = DiffusionTensor3 {
        xx: 1.,
        yy: 1.,
        zz: 1.,
        xy: 0.9,
        xz: 0.9,
        yz: -0.9,
    };
    assert!(neumann::solve_diffusion(
        &tri,
        |_, _| invalid,
        |_, _| 0.,
        |_, _| 0.,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann_quad::solve_diffusion(
        &quad,
        |_, _| invalid,
        |_, _| 0.,
        |_, _| 0.,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann3d::solve_diffusion(
        &tet,
        |_, _| invalid3,
        |_, _| 0.,
        |_, _| 0.,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |_, _| f64::NAN,
        |_, _| 0.,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann_quad::solve_diffusion(
        &quad,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, _| f64::NAN,
        NeumannOptions::default()
    )
    .is_err());
    assert!(neumann3d::solve_diffusion(
        &tet,
        |_, _| DiffusionTensor3::identity(),
        |_, _| 0.,
        |_, _| f64::INFINITY,
        NeumannOptions::default()
    )
    .is_err());
    let bad = NeumannOptions {
        compatibility_relative_tolerance: -1.,
        ..NeumannOptions::default()
    };
    assert!(neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, _| 0.,
        bad
    )
    .is_err());
    assert!(neumann_quad::solve_diffusion(
        &quad,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, _| 0.,
        bad
    )
    .is_err());
    assert!(neumann3d::solve_diffusion(
        &tet,
        |_, _| DiffusionTensor3::identity(),
        |_, _| 0.,
        |_, _| 0.,
        bad
    )
    .is_err());
    let no = NeumannOptions {
        cg: CgOptions {
            max_iterations: 0,
            ..CgOptions::default()
        },
        ..NeumannOptions::default()
    };
    assert!(neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, n| n[0],
        no
    )
    .is_err());
    assert!(neumann_quad::solve_diffusion(
        &quad,
        |_, _| DiffusionTensor::identity(),
        |_, _| 0.,
        |_, n| n[0],
        no
    )
    .is_err());
    assert!(neumann3d::solve_diffusion(
        &tet,
        |_, _| DiffusionTensor3::identity(),
        |_, _| 0.,
        |_, n| n[0],
        no
    )
    .is_err());
}
#[test]
fn identity_diffusion_reproduces_original_solvers() {
    let tri = Mesh2::unit_square(3).unwrap();
    let a = neumann::solve(
        &tri,
        |_| -4.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1]),
        NeumannOptions::default(),
    )
    .unwrap();
    let b = neumann::solve_diffusion(
        &tri,
        |_, _| DiffusionTensor::identity(),
        |_, _| -4.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1]),
        NeumannOptions::default(),
    )
    .unwrap();
    for (x, y) in a.values.iter().zip(b.values) {
        near(*x, y, 1e-9);
    }
    let q = QuadMesh::unit_square(3).unwrap();
    let a = neumann_quad::solve(
        &q,
        |_| -4.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1]),
        NeumannOptions::default(),
    )
    .unwrap();
    let b = neumann_quad::solve_diffusion(
        &q,
        |_, _| DiffusionTensor::identity(),
        |_, _| -4.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1]),
        NeumannOptions::default(),
    )
    .unwrap();
    for (x, y) in a.values.iter().zip(b.values) {
        near(*x, y, 1e-9);
    }
    let t = TetMesh::unit_cube(2).unwrap();
    let a = neumann3d::solve(
        &t,
        |_| -6.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1] + p.z * n[2]),
        NeumannOptions::default(),
    )
    .unwrap();
    let b = neumann3d::solve_diffusion(
        &t,
        |_, _| DiffusionTensor3::identity(),
        |_, _| -6.,
        |p, n| 2. * (p.x * n[0] + p.y * n[1] + p.z * n[2]),
        NeumannOptions::default(),
    )
    .unwrap();
    for (x, y) in a.values.iter().zip(b.values) {
        near(*x, y, 1e-9);
    }
}
