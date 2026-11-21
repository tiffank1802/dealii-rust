use dealii_rust::{
    base::Point2,
    grid::{BoundaryEdgeId, CellId, VertexId},
    lac::{conjugate_gradient, CgOptions},
    quad::{self, BilinearMap, QuadMesh, UniformRefinement},
    quadrature::square_degree_5,
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}", a, b);
}
fn area(mesh: &QuadMesh, id: CellId) -> f64 {
    let map = BilinearMap::new(mesh.cell_points(id).unwrap()).unwrap();
    square_degree_5()
        .iter()
        .map(|q| q.weight * map.evaluate(q.point).unwrap().determinant)
        .sum()
}
fn cell_integral(mesh: &QuadMesh, id: CellId, values: &[f64]) -> f64 {
    let map = BilinearMap::new(mesh.cell_points(id).unwrap()).unwrap();
    let cell = mesh.cells()[id.0];
    square_degree_5()
        .iter()
        .map(|q| {
            let shape = quad::shape_values(q.point);
            let u: f64 = (0..4).map(|i| shape[i] * values[cell[i].0]).sum();
            q.weight * map.evaluate(q.point).unwrap().determinant * u
        })
        .sum()
}
fn distorted(n: usize) -> QuadMesh {
    let base = QuadMesh::unit_square(n).unwrap();
    let pts = base
        .vertices()
        .iter()
        .map(|p| Point2 {
            x: p.x + 0.2 * p.x * (1. - p.x) * p.y,
            y: p.y + 0.1 * p.y * (1. - p.y) * p.x,
        })
        .collect();
    QuadMesh::new(pts, base.cells().to_vec()).unwrap()
}
#[test]
fn topology_counts_areas_and_boundary_lineage() {
    for n in [1, 2, 4] {
        let parent = distorted(n);
        let r = UniformRefinement::new(&parent).unwrap();
        let mesh = r.mesh();
        assert_eq!(mesh.cells().len(), 4 * parent.cells().len());
        assert_eq!(mesh.vertices().len(), (2 * n + 1) * (2 * n + 1));
        assert_eq!(
            mesh.boundary_edges().len(),
            2 * parent.boundary_edges().len()
        );
        assert_eq!(
            &mesh.vertices()[..parent.vertices().len()],
            parent.vertices()
        );
        for i in 0..parent.cells().len() {
            let children = r.children(CellId(i)).unwrap();
            for &child in &children {
                assert_eq!(r.parent_cell(child).unwrap(), CellId(i));
            }
            near(
                children.iter().map(|&child| area(mesh, child)).sum(),
                area(&parent, CellId(i)),
                1e-13,
            );
        }
        let mut counts = vec![0; parent.boundary_edges().len()];
        for i in 0..mesh.boundary_edges().len() {
            let id = BoundaryEdgeId(i);
            let orig = r.parent_boundary_edge(id).unwrap();
            counts[orig.0] += 1;
            let old = parent.boundary_edge_geometry(orig).unwrap();
            let new = mesh.boundary_edge_geometry(id).unwrap();
            near(new.length * 2., old.length, 1e-13);
            near(new.outward_normal[0], old.outward_normal[0], 1e-13);
            near(new.outward_normal[1], old.outward_normal[1], 1e-13);
        }
        assert!(counts.iter().all(|&n| n == 2));
        assert!(r.children(CellId(parent.cells().len())).is_err());
        assert!(r.parent_cell(CellId(mesh.cells().len())).is_err());
        assert!(r
            .parent_boundary_edge(BoundaryEdgeId(mesh.boundary_edges().len()))
            .is_err());
    }
}
#[test]
fn field_and_mapping_are_preserved_inside_each_nonaffine_child() {
    let parent = distorted(3);
    let r = UniformRefinement::new(&parent).unwrap();
    let fine = r.mesh();
    // A general continuous piecewise Q1 field, not a global affine function.
    let values: Vec<_> = (0..parent.vertices().len())
        .map(|i| ((i * i + 3 * i) % 19) as f64 - 7.)
        .collect();
    let transferred = r.prolongate_q1(&values).unwrap();
    assert_eq!(&transferred[..values.len()], values.as_slice());
    let offsets = [(0., 0.), (0.5, 0.), (0.5, 0.5), (0., 0.5)];
    for i in 0..parent.cells().len() {
        let old = BilinearMap::new(parent.cell_points(CellId(i)).unwrap()).unwrap();
        let old_cell = parent.cells()[i];
        for (k, child) in r.children(CellId(i)).unwrap().iter().enumerate() {
            let new = BilinearMap::new(fine.cell_points(*child).unwrap()).unwrap();
            let new_cell = fine.cells()[child.0];
            for local in [p(0.17, 0.22), p(0.5, 0.5), p(0.8, 0.7)] {
                let parent_point = p(offsets[k].0 + 0.5 * local.x, offsets[k].1 + 0.5 * local.y);
                let x = old.evaluate(parent_point).unwrap();
                let y = new.evaluate(local).unwrap();
                near(x.point.x, y.point.x, 1e-13);
                near(x.point.y, y.point.y, 1e-13);
                near(x.determinant * 0.25, y.determinant, 1e-13);
                let os = quad::shape_values(parent_point);
                let ns = quad::shape_values(local);
                let mut u = 0.;
                let mut v = 0.;
                let mut ug = [0.; 2];
                let mut vg = [0.; 2];
                for j in 0..4 {
                    u += values[old_cell[j].0] * os[j];
                    v += transferred[new_cell[j].0] * ns[j];
                    for d in 0..2 {
                        ug[d] += values[old_cell[j].0] * x.gradients[j][d];
                        vg[d] += transferred[new_cell[j].0] * y.gradients[j][d];
                    }
                }
                near(u, v, 1e-12);
                near(ug[0], vg[0], 1e-11);
                near(ug[1], vg[1], 1e-11);
            }
        }
        let total: f64 = r
            .children(CellId(i))
            .unwrap()
            .iter()
            .map(|&c| cell_integral(fine, c, &transferred))
            .sum();
        near(total, cell_integral(&parent, CellId(i), &values), 1e-13);
    }
}
#[test]
fn exact_interpolation_at_midpoints_and_center() {
    let mesh = QuadMesh::new(
        vec![p(0., 0.), p(2., 0.), p(1., 1.), p(0., 1.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    let r = UniformRefinement::new(&mesh).unwrap();
    let v = r.prolongate_q1(&[1., 4., 9., 16.]).unwrap();
    let children = r.children(CellId(0)).unwrap();
    let first = r.mesh().cells()[children[0].0];
    near(v[first[0].0], 1., 1e-14);
    near(v[first[1].0], 2.5, 1e-14);
    near(v[first[2].0], 7.5, 1e-14);
    near(v[first[3].0], 8.5, 1e-14);
    near(r.mesh().vertices()[first[2].0].x, 0.75, 1e-14);
    near(r.mesh().vertices()[first[2].0].y, 0.5, 1e-14);
}
#[test]
fn shared_midpoint_and_transfer_across_common_edge() {
    let mesh = distorted(2);
    let r = UniformRefinement::new(&mesh).unwrap();
    let fine = r.mesh();
    let values: Vec<_> = (0..mesh.vertices().len()).map(|i| (i * i) as f64).collect();
    let v = r.prolongate_q1(&values).unwrap();
    for cell in mesh.cells() {
        for j in 0..4 {
            let a = cell[j];
            let b = cell[(j + 1) % 4];
            let pa = mesh.vertices()[a.0];
            let pb = mesh.vertices()[b.0];
            let x = 0.5 * pa.x + 0.5 * pb.x;
            let y = 0.5 * pa.y + 0.5 * pb.y;
            let found: Vec<_> = fine
                .vertices()
                .iter()
                .enumerate()
                .filter(|(_, p)| p.x == x && p.y == y)
                .collect();
            assert_eq!(found.len(), 1);
            near(v[found[0].0], 0.5 * values[a.0] + 0.5 * values[b.0], 1e-12);
        }
    }
}
#[test]
fn disconnected_components_and_rotated_vertex_numbering() {
    let mesh = QuadMesh::new(
        vec![
            p(0., 1.),
            p(0., 0.),
            p(2., 0.),
            p(1., 1.),
            p(5., 0.),
            p(6., 0.),
            p(6., 1.),
            p(5., 1.),
        ],
        vec![
            [VertexId(1), VertexId(2), VertexId(3), VertexId(0)],
            [VertexId(4), VertexId(5), VertexId(6), VertexId(7)],
        ],
    )
    .unwrap();
    let r = UniformRefinement::new(&mesh).unwrap();
    assert_eq!(r.mesh().vertices().len(), 18);
    assert_eq!(r.mesh().boundary_edges().len(), 16);
    let f = |p: Point2| 2. + 3. * p.x - 4. * p.y;
    let values: Vec<_> = mesh.vertices().iter().copied().map(f).collect();
    for (&u, &point) in r
        .prolongate_q1(&values)
        .unwrap()
        .iter()
        .zip(r.mesh().vertices())
    {
        near(u, f(point), 1e-12);
    }
    for i in 0..2 {
        near(
            r.children(CellId(i))
                .unwrap()
                .iter()
                .map(|&c| area(r.mesh(), c))
                .sum(),
            area(&mesh, CellId(i)),
            1e-13,
        );
    }
}
#[test]
fn rejects_bad_fields_and_unrepresentable_children() {
    let mesh = QuadMesh::unit_square(1).unwrap();
    let r = UniformRefinement::new(&mesh).unwrap();
    assert!(r.prolongate_q1(&[]).is_err());
    assert!(r.prolongate_q1(&[f64::NAN; 4]).is_err());
    assert!(r.prolongate_q1(&[f64::INFINITY; 4]).is_err());
    assert!(r
        .prolongate_q1(&[f64::MAX; 4])
        .unwrap()
        .iter()
        .all(|v| *v == f64::MAX));
    // Width 2 is representable at x=1e16, but its midpoint is not.
    let large = QuadMesh::new(
        vec![p(1e16, 0.), p(1e16 + 2., 0.), p(1e16 + 2., 2.), p(1e16, 2.)],
        vec![[VertexId(0), VertexId(1), VertexId(2), VertexId(3)]],
    )
    .unwrap();
    assert!(UniformRefinement::new(&large).is_err());
}
#[test]
fn repeated_refinement_preserves_transfer_and_lineage() {
    let coarse = distorted(2);
    let r1 = UniformRefinement::new(&coarse).unwrap();
    let r2 = UniformRefinement::new(r1.mesh()).unwrap();
    let values: Vec<_> = (0..coarse.vertices().len())
        .map(|i| (i * i) as f64)
        .collect();
    let fine = r2
        .prolongate_q1(&r1.prolongate_q1(&values).unwrap())
        .unwrap();
    for i in 0..r2.mesh().cells().len() {
        let middle = r2.parent_cell(CellId(i)).unwrap();
        let root = r1.parent_cell(middle).unwrap();
        assert!(root.0 < coarse.cells().len());
    }
    let root_integral: f64 = (0..coarse.cells().len())
        .map(|i| cell_integral(&coarse, CellId(i), &values))
        .sum();
    let fine_integral: f64 = (0..r2.mesh().cells().len())
        .map(|i| cell_integral(r2.mesh(), CellId(i), &fine))
        .sum();
    near(root_integral, fine_integral, 1e-12);
    for i in 0..r2.mesh().boundary_edges().len() {
        let edge = r2.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
        assert!(r1.parent_boundary_edge(edge).unwrap().0 < coarse.boundary_edges().len());
    }
}
fn solution_error(mesh: &QuadMesh) -> (f64, f64) {
    use std::f64::consts::PI;
    let exact = |p: Point2| (PI * p.x).sin() * (PI * p.y).sin();
    let (a, b) = quad::assemble_poisson(mesh, |p| 2. * PI * PI * exact(p)).unwrap();
    let constraints: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
    let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
    let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
    quad::errors(mesh, &u, exact, |p| {
        [
            PI * (PI * p.x).cos() * (PI * p.y).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).cos(),
        ]
    })
    .unwrap()
}
#[test]
fn repeated_refinement_solves_poisson_and_converges() {
    let mut mesh = distorted(4);
    let mut errors = vec![];
    for _ in 0..3 {
        errors.push(solution_error(&mesh));
        mesh = UniformRefinement::new(&mesh).unwrap().mesh().clone();
    }
    for pair in errors.windows(2) {
        let l2 = (pair[0].0 / pair[1].0).log2();
        let h1 = (pair[0].1 / pair[1].1).log2();
        assert!((1.8..2.2).contains(&l2), "L2 {}", l2);
        assert!((0.9..1.1).contains(&h1), "H1 {}", h1);
    }
}

#[test]
fn boundary_labels_inherited_through_multiple_levels() {
    let mut mesh = distorted(3);
    let mut labels: Vec<_> = (0..mesh.boundary_edges().len())
        .map(|i| {
            let e = mesh.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
            if e.points.iter().all(|p| p.x == 0.) {
                1
            } else {
                2
            }
        })
        .collect();
    for _ in 0..3 {
        let r = UniformRefinement::new(&mesh).unwrap();
        let next: Vec<_> = (0..r.mesh().boundary_edges().len())
            .map(|i| labels[r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap().0])
            .collect();
        assert_eq!(
            next.iter().filter(|&&tag| tag == 1).count(),
            2 * labels.iter().filter(|&&tag| tag == 1).count()
        );
        for (i, &tag) in next.iter().enumerate() {
            let e = r.mesh().boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
            assert_eq!(tag == 1, e.points.iter().all(|p| p.x == 0.));
        }
        labels = next;
        mesh = r.mesh().clone();
    }
}
