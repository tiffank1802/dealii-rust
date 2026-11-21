#![allow(clippy::needless_range_loop)] // Indexed analytic matrices.
use dealii_rust::{
    base::Point3,
    elliptic_tetra::{self, diffusion_reaction_p1, DiffusionTensor3},
    grid::{CellId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    quadrature::tetra_degree_5,
    tetra::{self, BoundaryFaceId, TetMesh},
};
fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn reference() -> [Point3; 4] {
    [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]
}
fn tensor() -> DiffusionTensor3 {
    DiffusionTensor3 {
        xx: 2.,
        xy: 0.5,
        xz: 0.25,
        yy: 3.,
        yz: 0.75,
        zz: 4.,
    }
}
#[test]
fn analytic_anisotropic_stiffness_reaction_and_load() {
    let s = diffusion_reaction_p1(reference(), |_| tensor(), |_| 4., |_| 6.).unwrap();
    // Gradients on the reference tetrahedron and A are integrated over volume 1/6.
    let grad = [[-1., -1., -1.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let a = [[2., 0.5, 0.25], [0.5, 3., 0.75], [0.25, 0.75, 4.]];
    for i in 0..4 {
        near(s.rhs[i], 0.25, 1e-14);
        for j in 0..4 {
            let mut value = 0.;
            for k in 0..3 {
                for l in 0..3 {
                    value += grad[i][k] * a[k][l] * grad[j][l] / 6.;
                }
            }
            value += if i == j { 1. / 15. } else { 1. / 30. };
            near(s.matrix[i][j], value, 2e-14);
            assert_eq!(s.matrix[i][j], s.matrix[j][i]);
        }
    }
}
#[test]
fn fourth_degree_source_is_integrated_exactly() {
    let s = diffusion_reaction_p1(
        reference(),
        |_| DiffusionTensor3::identity(),
        |_| 0.,
        |p| p.x.powi(4),
    )
    .unwrap();
    for i in 0..4 {
        near(s.rhs[i], if i == 1 { 1. / 336. } else { 1. / 1680. }, 1e-15);
    }
}
#[test]
fn polynomial_diffusion_and_reaction_reach_documented_degree() {
    let local = diffusion_reaction_p1(
        reference(),
        |p| DiffusionTensor3 {
            xx: 1. + p.x.powi(5),
            ..DiffusionTensor3::identity()
        },
        |p| p.x.powi(3),
        |p| p.x.powi(4),
    )
    .unwrap();
    // ∫x^5 = 1/336 on the reference tetrahedron. For vertex 1, φ=x and ∇φ=e_x.
    near(local.matrix[1][1], 1. / 6. + 2. / 336., 1e-14);
    near(local.rhs[1], 1. / 336., 1e-15);
}
#[test]
fn identity_recovers_existing_poisson_assembly() {
    let mesh = TetMesh::unit_cube(2).unwrap();
    let source = |p: Point3| 1. + p.x * p.x + p.y - p.z;
    let (a, b) = elliptic_tetra::assemble(
        &mesh,
        |_, _| DiffusionTensor3::identity(),
        |_, _| 0.,
        |_, p| source(p),
    )
    .unwrap();
    let (r, rb) = tetra::assemble_poisson(&mesh, source).unwrap();
    for i in 0..b.len() {
        near(b[i], rb[i], 2e-14);
        for j in 0..b.len() {
            near(a.get(i, j).unwrap(), r.get(i, j).unwrap(), 2e-14);
        }
    }
}
#[test]
fn cellwise_coefficients_are_assembled_and_errors_identify_cell() {
    let mesh = TetMesh::unit_cube(1).unwrap();
    let (a, b) = elliptic_tetra::assemble(
        &mesh,
        |id, _| DiffusionTensor3 {
            xx: 1. + id.0 as f64,
            yy: 1. + id.0 as f64,
            zz: 1. + id.0 as f64,
            ..DiffusionTensor3::identity()
        },
        |id, _| id.0 as f64,
        |id, _| 1. + id.0 as f64,
    )
    .unwrap();
    let mut expected = vec![vec![0.; b.len()]; b.len()];
    let mut load = vec![0.; b.len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let v = i as f64;
        let local = diffusion_reaction_p1(
            mesh.cell_points(CellId(i)).unwrap(),
            |_| DiffusionTensor3 {
                xx: 1. + v,
                yy: 1. + v,
                zz: 1. + v,
                ..DiffusionTensor3::identity()
            },
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
        near(b[i], load[i], 1e-14);
        for j in 0..b.len() {
            near(a.get(i, j).unwrap(), expected[i][j], 1e-14);
        }
    }
    let err = elliptic_tetra::assemble(
        &mesh,
        |id, _| {
            if id.0 == 3 {
                DiffusionTensor3 {
                    xx: -1.,
                    ..tensor()
                }
            } else {
                tensor()
            }
        },
        |_, _| 0.,
        |_, _| 0.,
    )
    .unwrap_err();
    assert!(err.to_string().contains("cell 3"));
}
#[test]
fn vertex_permutation_preserves_local_system() {
    let points = [
        p(0., 0., 0.),
        p(2., 0., 0.),
        p(0.3, 1., 0.),
        p(0.2, 0.5, 1.),
    ];
    let d = |p: Point3| DiffusionTensor3 {
        xx: 2. + p.x,
        yy: 3. + p.y,
        zz: 4. + p.z,
        ..tensor()
    };
    let c = |p: Point3| 1. + p.x;
    let f = |p: Point3| p.x + p.y * p.z;
    let a = diffusion_reaction_p1(points, d, c, f).unwrap();
    // Even permutation preserves positive orientation.
    let b = diffusion_reaction_p1([points[1], points[2], points[0], points[3]], d, c, f).unwrap();
    let perm = [1, 2, 0, 3];
    for i in 0..4 {
        near(b.rhs[i], a.rhs[perm[i]], 1e-13);
        for j in 0..4 {
            near(b.matrix[i][j], a.matrix[perm[i]][perm[j]], 1e-13);
        }
    }
}
#[test]
fn detects_indefiniteness_and_invalid_data_including_full_3d_condition() {
    let mut cases = vec![
        DiffusionTensor3 { xx: 0., ..tensor() },
        DiffusionTensor3 {
            xy: f64::NAN,
            ..tensor()
        },
        DiffusionTensor3 {
            zz: f64::INFINITY,
            ..tensor()
        },
    ];
    // Every 2x2 principal minor is positive, but the 3x3 determinant is negative.
    cases.push(DiffusionTensor3 {
        xx: 1.,
        yy: 1.,
        zz: 1.,
        xy: 0.9,
        xz: 0.9,
        yz: -0.9,
    });
    cases.push(DiffusionTensor3 {
        xx: 1.,
        yy: 1.,
        zz: 1.,
        xy: 0.,
        xz: 0.,
        yz: 1.,
    });
    for a in cases {
        assert!(diffusion_reaction_p1(reference(), |_| a, |_| 0., |_| 0.).is_err());
    }
    for c in [-1., f64::NAN, f64::INFINITY] {
        assert!(diffusion_reaction_p1(reference(), |_| tensor(), |_| c, |_| 0.).is_err());
    }
    assert!(diffusion_reaction_p1(reference(), |_| tensor(), |_| 0., |_| f64::NAN).is_err());
    assert!(diffusion_reaction_p1(
        [
            reference()[0],
            reference()[2],
            reference()[1],
            reference()[3]
        ],
        |_| tensor(),
        |_| 0.,
        |_| 0.
    )
    .is_err());
    assert!(diffusion_reaction_p1(
        reference(),
        |_| DiffusionTensor3 {
            xx: f64::MAX,
            yy: f64::MAX,
            zz: f64::MAX,
            ..DiffusionTensor3::identity()
        },
        |_| 0.,
        |_| 0.
    )
    .is_err());
}
#[test]
fn sampled_positivity_and_strongly_scaled_tensor() {
    // Nonpositive at an unsampled vertex is not detected.
    assert!(diffusion_reaction_p1(
        reference(),
        |p| DiffusionTensor3 {
            xx: if p.x == 0. && p.y == 0. && p.z == 0. {
                -1.
            } else {
                1.
            },
            ..DiffusionTensor3::identity()
        },
        |_| 0.,
        |_| 0.
    )
    .is_ok());
    let sample = tetra_degree_5()[0].point;
    assert!(diffusion_reaction_p1(
        reference(),
        |p| DiffusionTensor3 {
            xx: if p == sample { 0. } else { 1. },
            ..DiffusionTensor3::identity()
        },
        |_| 0.,
        |_| 0.
    )
    .is_err());
    let extreme = DiffusionTensor3 {
        xx: 1e-200,
        xy: 0.5,
        xz: 0.,
        yy: 1e200,
        yz: 0.,
        zz: 1.,
    };
    assert!(diffusion_reaction_p1(reference(), |_| extreme, |_| 0., |_| 0.).is_ok());
    assert!(diffusion_reaction_p1(
        reference(),
        |_| DiffusionTensor3 { xy: 2., ..extreme },
        |_| 0.,
        |_| 0.
    )
    .is_err());
}
fn exact(p: Point3) -> f64 {
    1. + p.x + p.x * p.x + 2. * p.y + p.y * p.y + 3. * p.z + p.z * p.z
}
fn grad(p: Point3) -> [f64; 3] {
    [1. + 2. * p.x, 2. + 2. * p.y, 3. + 2. * p.z]
}
fn variable(p: Point3) -> DiffusionTensor3 {
    DiffusionTensor3 {
        xx: 2. + p.x,
        yy: 3. + p.y,
        zz: 4. + p.z,
        ..tensor()
    }
}
fn conormal(p: Point3, n: [f64; 3]) -> f64 {
    let a = variable(p);
    let [x, y, z] = grad(p);
    n[0] * (a.xx * x + a.xy * y + a.xz * z)
        + n[1] * (a.xy * x + a.yy * y + a.yz * z)
        + n[2] * (a.xz * x + a.yz * y + a.zz * z)
}
fn mixed(n: usize) -> (f64, f64) {
    let mesh = TetMesh::unit_cube(n).unwrap();
    let (a, mut b) = elliptic_tetra::assemble(
        &mesh,
        |_, p| variable(p),
        |_, p| 1. + p.x + p.y + p.z,
        |_, p| -(24. + 4. * (p.x + p.y + p.z)) + (1. + p.x + p.y + p.z) * exact(p),
    )
    .unwrap();
    let mut d = vec![];
    let mut neumann = vec![];
    for i in 0..mesh.boundary_faces().len() {
        let id = BoundaryFaceId(i);
        if mesh
            .boundary_face_geometry(id)
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
    let g = tetra::assemble_neumann(&mesh, &neumann, conormal).unwrap();
    for (b, g) in b.iter_mut().zip(g) {
        *b += g;
    }
    let values = tetra::interpolate_boundary(&mesh, &d, exact).unwrap();
    let (a, b) = a.with_dirichlet(&b, &values).unwrap();
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
    tetra::errors(&mesh, &u, exact, grad).unwrap()
}
#[test]
fn anisotropic_mixed_problem_converges_in_3d() {
    let errors: Vec<_> = [4, 8, 16].iter().map(|&n| mixed(n)).collect();
    for pair in errors.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
        assert!((0.8..1.2).contains(&h1), "H1 rate {}", h1);
    }
}
#[test]
fn positive_reaction_all_neumann_recovers_affine_solution() {
    let mesh = TetMesh::new(
        reference().to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let u = |p: Point3| 1. + p.x + 2. * p.y + 3. * p.z;
    let (a, mut b) =
        elliptic_tetra::assemble(&mesh, |_, _| tensor(), |_, _| 2., |_, p| 2. * u(p)).unwrap();
    let faces: Vec<_> = (0..mesh.boundary_faces().len())
        .map(BoundaryFaceId)
        .collect();
    // A*[1,2,3] = [3.75,8.75,13.75].
    let g = tetra::assemble_neumann(&mesh, &faces, |_, n| {
        3.75 * n[0] + 8.75 * n[1] + 13.75 * n[2]
    })
    .unwrap();
    for (b, g) in b.iter_mut().zip(g) {
        *b += g;
    }
    let (values, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    for (value, p) in values.iter().zip(mesh.vertices()) {
        near(*value, u(*p), 1e-10);
    }
}
