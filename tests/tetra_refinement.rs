use dealii_rust::{
    base::Point3,
    grid::{CellId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    tetra::{self, AffineTet, BoundaryFaceId, TetMesh, UniformRefinement},
};
fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn volume(mesh: &TetMesh, id: CellId) -> f64 {
    AffineTet::new(mesh.cell_points(id).unwrap())
        .unwrap()
        .determinant()
        / 6.
}
fn integral(mesh: &TetMesh, id: CellId, v: &[f64]) -> f64 {
    volume(mesh, id) * mesh.cells()[id.0].iter().map(|id| v[id.0]).sum::<f64>() / 4.
}
fn skew(n: usize) -> TetMesh {
    let mesh = TetMesh::unit_cube(n).unwrap();
    let vertices = mesh
        .vertices()
        .iter()
        .map(|p| Point3 {
            x: 1. + 2. * p.x + 0.3 * p.y + 0.2 * p.z,
            y: -2. + 0.2 * p.x + 1.5 * p.y + 0.1 * p.z,
            z: 0.5 + 0.1 * p.x + 0.2 * p.y + 1.2 * p.z,
        })
        .collect();
    TetMesh::new(vertices, mesh.cells().to_vec()).unwrap()
}
#[test]
fn topology_volume_orientation_and_face_lineage() {
    for n in [1, 2] {
        let parent = skew(n);
        let r = UniformRefinement::new(&parent).unwrap();
        let mesh = r.mesh();
        assert_eq!(mesh.cells().len(), 8 * parent.cells().len());
        assert_eq!(mesh.vertices().len(), (2 * n + 1).pow(3));
        assert_eq!(
            mesh.boundary_faces().len(),
            4 * parent.boundary_faces().len()
        );
        assert_eq!(
            &mesh.vertices()[..parent.vertices().len()],
            parent.vertices()
        );
        for i in 0..parent.cells().len() {
            let children = r.children(CellId(i)).unwrap();
            for &c in &children {
                assert_eq!(r.parent_cell(c).unwrap(), CellId(i));
                near(volume(mesh, c), volume(&parent, CellId(i)) / 8., 1e-12);
            }
            near(
                children.iter().map(|&c| volume(mesh, c)).sum(),
                volume(&parent, CellId(i)),
                1e-12,
            );
        }
        let mut counts = vec![0; parent.boundary_faces().len()];
        for i in 0..mesh.boundary_faces().len() {
            let id = BoundaryFaceId(i);
            let parent_id = r.parent_boundary_face(id).unwrap();
            counts[parent_id.0] += 1;
            let f = mesh.boundary_face_geometry(id).unwrap();
            let original = parent.boundary_face_geometry(parent_id).unwrap();
            near(f.area * 4., original.area, 1e-12);
            for k in 0..3 {
                near(f.outward_normal[k], original.outward_normal[k], 1e-12);
            }
        }
        assert!(counts.iter().all(|&n| n == 4));
        assert!(r.children(CellId(parent.cells().len())).is_err());
        assert!(r.parent_cell(CellId(mesh.cells().len())).is_err());
        assert!(r
            .parent_boundary_face(BoundaryFaceId(mesh.boundary_faces().len()))
            .is_err());
    }
}
#[test]
fn skew_disconnected_tetrahedra_preserve_affine_data() {
    let m = TetMesh::new(
        vec![
            p(1., 2., 3.),
            p(3., 2., 3.),
            p(1.5, 4., 3.),
            p(1.2, 2.5, 5.),
            p(7., 2., 3.),
            p(9., 2., 3.),
            p(7.5, 4., 3.),
            p(7.2, 2.5, 5.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2), VertexId(3)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap();
    let r = UniformRefinement::new(&m).unwrap();
    assert_eq!(r.mesh().vertices().len(), 20);
    assert_eq!(r.mesh().boundary_faces().len(), 32);
    let f = |p: Point3| 2. + 3. * p.x - 4. * p.y + 0.5 * p.z;
    let values: Vec<_> = m.vertices().iter().copied().map(f).collect();
    for (u, &pt) in r
        .prolongate_p1(&values)
        .unwrap()
        .iter()
        .zip(r.mesh().vertices())
    {
        near(*u, f(pt), 1e-12);
    }
    for i in 0..2 {
        near(
            r.children(CellId(i))
                .unwrap()
                .iter()
                .map(|&c| volume(r.mesh(), c))
                .sum(),
            volume(&m, CellId(i)),
            1e-12,
        );
    }
}
#[test]
fn every_child_preserves_piecewise_affine_values_gradients_and_integral() {
    let parent = skew(2);
    let r = UniformRefinement::new(&parent).unwrap();
    let fine = r.mesh();
    let values: Vec<_> = (0..parent.vertices().len())
        .map(|i| ((i * i + 3 * i) % 19) as f64 - 7.)
        .collect();
    let transferred = r.prolongate_p1(&values).unwrap();
    assert_eq!(&transferred[..values.len()], values.as_slice());
    for (i, cell) in parent.cells().iter().enumerate() {
        let map = AffineTet::new(parent.cell_points(CellId(i)).unwrap()).unwrap();
        let grads = map.gradients();
        let mut old_grad = [0.; 3];
        for j in 0..4 {
            for (k, value) in old_grad.iter_mut().enumerate() {
                *value += values[cell[j].0] * grads[j][k];
            }
        }
        let origin = parent.vertices()[cell[0].0];
        for child in r.children(CellId(i)).unwrap() {
            let child_map = AffineTet::new(fine.cell_points(child).unwrap()).unwrap();
            let shape_grads = child_map.gradients();
            let child_ids = fine.cells()[child.0];
            let mut new_grad = [0.; 3];
            for j in 0..4 {
                for (k, value) in new_grad.iter_mut().enumerate() {
                    *value += transferred[child_ids[j].0] * shape_grads[j][k];
                }
            }
            for k in 0..3 {
                near(old_grad[k], new_grad[k], 1e-11);
            }
            for q in [p(0.1, 0.2, 0.3), p(0.25, 0.25, 0.25), p(0.4, 0.1, 0.2)] {
                let point = child_map.map(q).unwrap();
                let shapes = tetra::shape_values(q);
                let eval: f64 = (0..4)
                    .map(|j| shapes[j] * transferred[child_ids[j].0])
                    .sum();
                let expected = values[cell[0].0]
                    + old_grad[0] * (point.x - origin.x)
                    + old_grad[1] * (point.y - origin.y)
                    + old_grad[2] * (point.z - origin.z);
                near(eval, expected, 1e-11);
            }
        }
        let total: f64 = r
            .children(CellId(i))
            .unwrap()
            .iter()
            .map(|&c| integral(fine, c, &transferred))
            .sum();
        near(total, integral(&parent, CellId(i), &values), 1e-12);
    }
}
#[test]
fn shared_edge_midpoints_are_unique_and_averaged() {
    let parent = TetMesh::unit_cube(2).unwrap();
    let r = UniformRefinement::new(&parent).unwrap();
    let values: Vec<_> = (0..parent.vertices().len())
        .map(|i| (i * i) as f64)
        .collect();
    let u = r.prolongate_p1(&values).unwrap();
    for cell in parent.cells() {
        for a in 0..4 {
            for b in a + 1..4 {
                let va = cell[a];
                let vb = cell[b];
                let a = parent.vertices()[va.0];
                let b = parent.vertices()[vb.0];
                let pt = p(
                    a.x * 0.5 + b.x * 0.5,
                    a.y * 0.5 + b.y * 0.5,
                    a.z * 0.5 + b.z * 0.5,
                );
                let ids: Vec<_> = r
                    .mesh()
                    .vertices()
                    .iter()
                    .enumerate()
                    .filter(|(_, q)| **q == pt)
                    .collect();
                assert_eq!(ids.len(), 1);
                near(u[ids[0].0], values[va.0] * 0.5 + values[vb.0] * 0.5, 1e-12);
            }
        }
    }
}
#[test]
fn vertex_permutations_all_split_into_eight_equal_volumes() {
    let v = vec![
        p(0., 0., 0.),
        p(2., 0., 0.),
        p(0.5, 2., 0.),
        p(0.2, 0.4, 1.5),
    ];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let ids = [a, b, c, d];
                    if (0..4).any(|i| (i + 1..4).any(|j| ids[i] == ids[j])) {
                        continue;
                    }
                    let coords = ids.map(|i| v[i]);
                    if AffineTet::new(coords).is_err() {
                        continue;
                    }
                    let mesh = TetMesh::new(v.clone(), vec![ids.map(VertexId)]).unwrap();
                    let r = UniformRefinement::new(&mesh).unwrap();
                    assert_eq!(r.mesh().cells().len(), 8);
                    assert_eq!(r.mesh().boundary_faces().len(), 16);
                    for child in r.children(CellId(0)).unwrap() {
                        near(
                            volume(r.mesh(), child),
                            volume(&mesh, CellId(0)) / 8.,
                            1e-14,
                        );
                    }
                }
            }
        }
    }
}
#[test]
fn central_children_share_a_shortest_octahedron_diagonal() {
    let points = [
        p(0., 0., 0.),
        p(2., 0.3, 0.),
        p(0.5, 1.3, 0.2),
        p(0.3, 0.2, 2.),
    ];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let ids = [a, b, c, d];
                    if (0..4).any(|i| (i + 1..4).any(|j| ids[i] == ids[j]))
                        || AffineTet::new(ids.map(|i| points[i])).is_err()
                    {
                        continue;
                    }
                    let parent = TetMesh::new(points.to_vec(), vec![ids.map(VertexId)]).unwrap();
                    let r = UniformRefinement::new(&parent).unwrap();
                    let children = r.children(CellId(0)).unwrap();
                    let first = r.mesh().cells()[children[4].0];
                    let common: Vec<_> = first
                        .iter()
                        .copied()
                        .filter(|id| {
                            children[5..]
                                .iter()
                                .all(|c| r.mesh().cells()[c.0].contains(id))
                        })
                        .collect();
                    assert_eq!(common.len(), 2);
                    let len = |a: Point3, b: Point3| (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z);
                    let mid = |a: Point3, b: Point3| {
                        p((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
                    };
                    let [a, b, c, d] = ids.map(|i| points[i]);
                    let candidates = [
                        len(mid(a, b), mid(c, d)),
                        len(mid(a, c), mid(b, d)),
                        len(mid(a, d), mid(b, c)),
                    ];
                    let shortest = candidates.iter().copied().fold(f64::INFINITY, f64::min);
                    near(
                        len(
                            r.mesh().vertices()[common[0].0],
                            r.mesh().vertices()[common[1].0],
                        ),
                        shortest,
                        1e-14,
                    );
                    let again = UniformRefinement::new(&parent).unwrap();
                    assert_eq!(r.mesh().cells(), again.mesh().cells());
                }
            }
        }
    }
}
#[test]
fn invalid_fields_and_unrepresentable_children_are_rejected() {
    let mesh = TetMesh::unit_cube(1).unwrap();
    let r = UniformRefinement::new(&mesh).unwrap();
    assert!(r.prolongate_p1(&[]).is_err());
    assert!(r
        .prolongate_p1(&vec![f64::NAN; mesh.vertices().len()])
        .is_err());
    assert!(r
        .prolongate_p1(&vec![f64::INFINITY; mesh.vertices().len()])
        .is_err());
    assert!(r
        .prolongate_p1(&vec![f64::MAX; mesh.vertices().len()])
        .unwrap()
        .iter()
        .all(|v| *v == f64::MAX));
    let large = TetMesh::new(
        vec![
            p(1e16, 0., 0.),
            p(1e16 + 2., 0., 0.),
            p(1e16, 2., 0.),
            p(1e16, 0., 2.),
        ],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(UniformRefinement::new(&large).is_err());
}
#[test]
fn repeated_refinement_composes_lineage_and_conserves_integral() {
    let coarse = skew(1);
    let r1 = UniformRefinement::new(&coarse).unwrap();
    let r2 = UniformRefinement::new(r1.mesh()).unwrap();
    let values: Vec<_> = (0..coarse.vertices().len())
        .map(|i| (i * i) as f64)
        .collect();
    let fine = r2
        .prolongate_p1(&r1.prolongate_p1(&values).unwrap())
        .unwrap();
    for i in 0..r2.mesh().cells().len() {
        let mid = r2.parent_cell(CellId(i)).unwrap();
        assert!(r1.parent_cell(mid).unwrap().0 < coarse.cells().len());
    }
    let old: f64 = (0..coarse.cells().len())
        .map(|i| integral(&coarse, CellId(i), &values))
        .sum();
    let new: f64 = (0..r2.mesh().cells().len())
        .map(|i| integral(r2.mesh(), CellId(i), &fine))
        .sum();
    near(old, new, 1e-11);
    for i in 0..r2.mesh().boundary_faces().len() {
        let mid = r2.parent_boundary_face(BoundaryFaceId(i)).unwrap();
        assert!(r1.parent_boundary_face(mid).unwrap().0 < coarse.boundary_faces().len());
    }
}
#[test]
fn boundary_labels_follow_parent_faces() {
    let mut mesh = TetMesh::unit_cube(2).unwrap();
    let mut labels: Vec<_> = (0..mesh.boundary_faces().len())
        .map(|i| {
            mesh.boundary_face_geometry(BoundaryFaceId(i))
                .unwrap()
                .points
                .iter()
                .all(|p| p.x == 0.)
        })
        .collect();
    for _ in 0..2 {
        let r = UniformRefinement::new(&mesh).unwrap();
        let next: Vec<_> = (0..r.mesh().boundary_faces().len())
            .map(|i| labels[r.parent_boundary_face(BoundaryFaceId(i)).unwrap().0])
            .collect();
        assert_eq!(
            next.iter().filter(|&&b| b).count(),
            4 * labels.iter().filter(|&&b| b).count()
        );
        for (i, &b) in next.iter().enumerate() {
            assert_eq!(
                b,
                r.mesh()
                    .boundary_face_geometry(BoundaryFaceId(i))
                    .unwrap()
                    .points
                    .iter()
                    .all(|p| p.x == 0.)
            );
        }
        mesh = r.mesh().clone();
        labels = next;
    }
}
fn solve(mesh: &TetMesh) -> (f64, f64) {
    use std::f64::consts::PI;
    let exact = |p: Point3| (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).sin();
    let (a, b) = tetra::assemble_poisson(mesh, |p| 3. * PI * PI * exact(p)).unwrap();
    let d: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
    let (a, b) = a.with_dirichlet(&b, &d).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    tetra::errors(mesh, &u, exact, |p| {
        [
            PI * (PI * p.x).cos() * (PI * p.y).sin() * (PI * p.z).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).cos() * (PI * p.z).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).sin() * (PI * p.z).cos(),
        ]
    })
    .unwrap()
}
#[test]
fn repeated_refinement_solves_poisson_and_converges() {
    let mut mesh = TetMesh::unit_cube(4).unwrap();
    let mut errors = vec![];
    for level in 0..3 {
        errors.push(solve(&mesh));
        if level < 2 {
            mesh = UniformRefinement::new(&mesh).unwrap().mesh().clone();
        }
    }
    for pair in errors.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.4..2.5).contains(&l2), "L2 {}", l2);
        assert!((0.7..1.3).contains(&h1), "H1 {}", h1);
    }
}
