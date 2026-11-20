#![allow(clippy::needless_range_loop)] // Analytical local arrays use indices.
use dealii_rust::{
    base::Point3,
    grid::{CellId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    quadrature::tetra_degree_5,
    tetra::{self, AffineTet, BoundaryFaceId, TetMesh},
};
use std::f64::consts::PI;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn reference() -> [Point3; 4] {
    [p(0., 0., 0.), p(1., 0., 0.), p(0., 1., 0.), p(0., 0., 1.)]
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn fact(n: usize) -> f64 {
    (1..=n).fold(1., |a, x| a * x as f64)
}

#[test]
fn quadrature_integrates_all_monomials_of_total_degree_five() {
    let q = tetra_degree_5();
    assert_eq!(q.len(), 64);
    for a in 0..=5 {
        for b in 0..=5 - a {
            for c in 0..=5 - a - b {
                let value: f64 = q
                    .iter()
                    .map(|q| {
                        q.weight
                            * q.point.x.powi(a as i32)
                            * q.point.y.powi(b as i32)
                            * q.point.z.powi(c as i32)
                    })
                    .sum();
                near(
                    value,
                    fact(a) * fact(b) * fact(c) / fact(a + b + c + 3),
                    1e-14,
                );
            }
        }
    }
}

#[test]
fn affine_map_reproduces_coordinates_and_transforms_gradients() {
    let points = [p(1., 2., 3.), p(3., 2., 3.), p(2., 5., 3.), p(1., 3., 7.)];
    let map = AffineTet::new(points).unwrap();
    near(map.determinant(), 24., 1e-14);
    for (i, pt) in reference().iter().enumerate() {
        assert_eq!(map.map(*pt).unwrap(), points[i]);
    }
    for (i, pt) in reference().iter().enumerate() {
        let shape = tetra::shape_values(*pt);
        for j in 0..4 {
            near(shape[j], if i == j { 1. } else { 0. }, 1e-14);
        }
    }
    let g = map.gradients();
    for d in 0..3 {
        near(g.iter().map(|v| v[d]).sum(), 0., 1e-14);
        for coord in 0..3 {
            let value: f64 = (0..4)
                .map(|i| [points[i].x, points[i].y, points[i].z][coord] * g[i][d])
                .sum();
            near(value, if coord == d { 1. } else { 0. }, 1e-14);
        }
    }
    assert!(map
        .map(Point3 {
            x: f64::NAN,
            y: 0.,
            z: 0.
        })
        .is_err());
}

#[test]
fn cube_mesh_volumes_faces_and_outward_normals() {
    for n in [1, 2, 3] {
        let mesh = TetMesh::unit_cube(n).unwrap();
        assert_eq!(mesh.vertices().len(), (n + 1).pow(3));
        assert_eq!(mesh.cells().len(), 6 * n.pow(3));
        assert_eq!(mesh.boundary_faces().len(), 12 * n * n);
        assert_eq!(
            mesh.boundary_vertices().len(),
            (n + 1).pow(3) - n.saturating_sub(1).pow(3)
        );
        let volume: f64 = (0..mesh.cells().len())
            .map(|i| {
                AffineTet::new(mesh.cell_points(CellId(i)).unwrap())
                    .unwrap()
                    .determinant()
                    / 6.
            })
            .sum();
        near(volume, 1., 1e-12);
        let mut area = 0.;
        for (i, face) in mesh.boundary_faces().iter().enumerate() {
            let geometry = mesh.boundary_face_geometry(BoundaryFaceId(i)).unwrap();
            area += geometry.area;
            near(
                geometry.outward_normal.iter().map(|x| x * x).sum(),
                1.,
                1e-14,
            );
            let touching: Vec<_> = mesh
                .cells()
                .iter()
                .filter(|cell| face.iter().all(|v| cell.contains(v)))
                .collect();
            assert_eq!(touching.len(), 1);
            let opposite = touching[0].iter().find(|v| !face.contains(v)).unwrap();
            let v = mesh.vertices()[opposite.0];
            let a = geometry.points[0];
            let normal = geometry.outward_normal;
            assert!(
                normal[0] * (v.x - a.x) + normal[1] * (v.y - a.y) + normal[2] * (v.z - a.z) < 0.
            );
        }
        near(area, 6., 1e-12);
        assert!(mesh
            .boundary_face_geometry(BoundaryFaceId(mesh.boundary_faces().len()))
            .is_err());
    }
    assert!(TetMesh::unit_cube(0).is_err());
    assert!(TetMesh::unit_cube(usize::MAX).is_err());
}

#[test]
fn invalid_meshes_and_orientations_are_rejected() {
    let r = reference();
    assert!(AffineTet::new([r[0], r[2], r[1], r[3]]).is_err());
    assert!(AffineTet::new([r[0], r[1], r[2], p(0., 0., 1e-16)]).is_err());
    assert!(AffineTet::new([
        r[0],
        r[1],
        r[2],
        Point3 {
            x: 0.,
            y: 0.,
            z: f64::NAN
        }
    ])
    .is_err());
    assert!(TetMesh::new(
        r.to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(4)]]
    )
    .is_err());
    assert!(TetMesh::new(
        r.to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(2)]]
    )
    .is_err());
    assert!(TetMesh::new(
        r.to_vec(),
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(1), VertexId(2), VertexId(0), VertexId(3)]
        ]
    )
    .is_err());
    let mut v = r.to_vec();
    v.push(p(0., 0., 2.));
    assert!(TetMesh::new(
        v,
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(0), VertexId(1), VertexId(2), VertexId(4)]
        ]
    )
    .is_err());
}

#[test]
fn reference_tetra_local_stiffness_and_load_are_analytic() {
    let mesh = TetMesh::new(
        reference().to_vec(),
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let (a, b) = tetra::assemble_poisson(&mesh, |_| 6.).unwrap();
    let expected = [
        [3., -1., -1., -1.],
        [-1., 1., 0., 0.],
        [-1., 0., 1., 0.],
        [-1., 0., 0., 1.],
    ];
    for i in 0..4 {
        near(b[i], 0.25, 1e-14);
        for j in 0..4 {
            near(a.get(i, j).unwrap(), expected[i][j] / 6., 1e-14);
        }
    }
    for value in a.mul(&[1.; 4]).unwrap() {
        near(value, 0., 1e-14);
    }
}

fn solve(
    mesh: &TetMesh,
    source: impl Fn(Point3) -> f64,
    boundary: impl Fn(Point3) -> f64,
) -> Vec<f64> {
    let (a, b) = tetra::assemble_poisson(mesh, source).unwrap();
    let constraints: Vec<_> = mesh
        .boundary_vertices()
        .iter()
        .map(|id| (id.0, boundary(mesh.vertices()[id.0])))
        .collect();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let options = CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-13,
        ..CgOptions::default()
    };
    conjugate_gradient(&a, &b, options).unwrap().0
}

#[test]
fn nonzero_dirichlet_affine_patch_on_distorted_mesh() {
    let base = TetMesh::unit_cube(3).unwrap();
    let points: Vec<_> = base
        .vertices()
        .iter()
        .map(|p| Point3 {
            x: p.x + 0.08 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.06 * p.y * (1. - p.y) * p.z,
            z: p.z + 0.04 * p.z * (1. - p.z) * p.x,
        })
        .collect();
    let mesh = TetMesh::new(points, base.cells().to_vec()).unwrap();
    let exact = |p: Point3| 1. + 2. * p.x - 3. * p.y + 4. * p.z;
    let u = solve(&mesh, |_| 0., exact);
    for (&value, &pt) in u.iter().zip(mesh.vertices()) {
        near(value, exact(pt), 1e-10);
    }
    let (l2, h1) = tetra::errors(&mesh, &u, exact, |_| [2., -3., 4.]).unwrap();
    assert!(l2 < 1e-10 && h1 < 1e-9);
}

fn exact(p: Point3) -> f64 {
    (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).sin()
}
fn gradient(p: Point3) -> [f64; 3] {
    [
        PI * (PI * p.x).cos() * (PI * p.y).sin() * (PI * p.z).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).cos() * (PI * p.z).sin(),
        PI * (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).cos(),
    ]
}
#[test]
fn manufactured_solution_converges_under_uniform_refinement() {
    let mut results = Vec::new();
    for n in [4, 8, 16] {
        let mesh = TetMesh::unit_cube(n).unwrap();
        let u = solve(&mesh, |p| 3. * PI * PI * exact(p), |_| 0.);
        results.push(tetra::errors(&mesh, &u, exact, gradient).unwrap());
    }
    for pair in results.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.7..2.3).contains(&l2), "L2 rate {}", l2);
        assert!((0.8..1.2).contains(&h1), "H1 rate {}", h1);
    }
}

#[test]
fn rejects_bad_numerical_values_and_exports_tetrahedra() {
    let mesh = TetMesh::unit_cube(1).unwrap();
    assert!(tetra::assemble_poisson(&mesh, |_| f64::NAN).is_err());
    assert!(tetra::errors(&mesh, &[0.], |_| 0., |_| [0.; 3]).is_err());
    assert!(tetra::errors(&mesh, &[f64::NAN; 8], |_| 0., |_| [0.; 3]).is_err());
    assert!(tetra::errors(&mesh, &[0.; 8], |_| f64::INFINITY, |_| [0.; 3]).is_err());
    assert!(tetra::write_vtk(Vec::new(), &mesh, &[0.]).is_err());
    let mut out = Vec::new();
    tetra::write_vtk(&mut out, &mesh, &[0.; 8]).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text.contains("CELLS 6 30"));
    assert!(text.contains("CELL_TYPES 6\n10\n10\n"));
    assert!(text.contains("POINT_DATA 8"));
}

#[test]
fn shared_face_orientation_and_nonmanifold_faces_are_checked() {
    let points = vec![
        p(0., 0., 0.),
        p(1., 0., 0.),
        p(0., 1., 0.),
        p(0., 0., 1.),
        p(0., 0., -1.),
        p(0.2, 0.2, 1.),
    ];
    let a = [VertexId(0), VertexId(1), VertexId(2), VertexId(3)];
    let b = [VertexId(0), VertexId(2), VertexId(1), VertexId(4)];
    let c = [VertexId(0), VertexId(1), VertexId(2), VertexId(5)];
    let err = TetMesh::new(points.clone(), vec![a, c]).unwrap_err();
    assert!(
        err.to_string().contains("inconsistent orientation"),
        "{}",
        err
    );
    let err = TetMesh::new(points, vec![a, b, c]).unwrap_err();
    assert!(err.to_string().contains("non-manifold face"), "{}", err);
}
