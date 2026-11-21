use dealii_rust::{
    base::{triangle_determinant, Point2},
    field,
    grid::{BoundaryEdgeId, CellId, Mesh2, NonconformingRefinement, VertexId},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};
use std::collections::BTreeMap;
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() < tol, "{} != {}, tol {}", a, b, tol);
}
fn area(points: [Point2; 3]) -> f64 {
    triangle_determinant(points).unwrap() * 0.5
}
fn key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}
fn opts() -> CgOptions {
    CgOptions {
        relative_tolerance: 1e-12,
        absolute_tolerance: 1e-14,
        ..CgOptions::default()
    }
}
#[test]
fn one_mark_has_exactly_one_hanging_node_and_true_boundary() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let r = NonconformingRefinement::new(&mesh, &[CellId(0)]).unwrap();
    assert_eq!(r.cells().len(), 5);
    assert_eq!(r.vertices().len(), 7);
    assert_eq!(r.hanging_nodes().len(), 1);
    assert_eq!(r.boundary_edges().len(), 6);
    let h = r.hanging_nodes()[0];
    assert_eq!(h.masters, [VertexId(0), VertexId(3)]);
    assert!(!r.boundary_vertices().contains(&h.vertex));
    near(r.vertices()[h.vertex.0].x, 0.5, 1e-15);
    near(r.vertices()[h.vertex.0].y, 0.5, 1e-15);
    assert_eq!(r.children(CellId(0)).unwrap().len(), 4);
    assert_eq!(r.children(CellId(1)).unwrap().len(), 1);
    let c = r.hanging_constraints().unwrap().close().unwrap();
    assert!(!c.free_dofs().contains(&h.vertex.0));
}
#[test]
fn all_markings_have_correct_hanging_edges_parent_areas_and_boundary() {
    let parent = Mesh2::unit_square(2).unwrap();
    let mut incidence: BTreeMap<(VertexId, VertexId), Vec<usize>> = BTreeMap::new();
    for (i, c) in parent.cells().iter().enumerate() {
        for j in 0..3 {
            incidence
                .entry(key(c[j], c[(j + 1) % 3]))
                .or_default()
                .push(i);
        }
    }
    for mask in 0usize..(1 << parent.cells().len()) {
        let marked: Vec<_> = (0..parent.cells().len())
            .filter(|i| mask & (1 << i) != 0)
            .map(CellId)
            .collect();
        let r = NonconformingRefinement::new(&parent, &marked).unwrap();
        assert_eq!(r.cells().len(), parent.cells().len() + 3 * marked.len());
        let expected: Vec<_> = incidence
            .iter()
            .filter_map(|(&edge, c)| {
                if c.len() == 2 && ((mask >> c[0]) & 1) != ((mask >> c[1]) & 1) {
                    Some(edge)
                } else {
                    None
                }
            })
            .collect();
        let actual: Vec<_> = r
            .hanging_nodes()
            .iter()
            .map(|h| key(h.masters[0], h.masters[1]))
            .collect();
        assert_eq!(actual, expected);
        for h in r.hanging_nodes() {
            assert!(!r.boundary_vertices().contains(&h.vertex));
            let p = r.vertices()[h.vertex.0];
            let a = r.vertices()[h.masters[0].0];
            let b = r.vertices()[h.masters[1].0];
            near(p.x, 0.5 * a.x + 0.5 * b.x, 1e-15);
            near(p.y, 0.5 * a.y + 0.5 * b.y, 1e-15);
        }
        for i in 0..parent.cells().len() {
            let children = r.children(CellId(i)).unwrap();
            assert_eq!(children.len(), if mask & (1 << i) != 0 { 4 } else { 1 });
            let sum: f64 = children
                .iter()
                .map(|&child| {
                    assert_eq!(r.parent_cell(child).unwrap(), CellId(i));
                    area(r.cell_points(child).unwrap())
                })
                .sum();
            near(sum, area(parent.cell_points(CellId(i)).unwrap()), 1e-14);
        }
        let mut length = vec![0.; parent.boundary_edges().len()];
        for i in 0..r.boundary_edges().len() {
            let child = r.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
            let parent_id = r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
            let original = parent.boundary_edge_geometry(parent_id).unwrap();
            length[parent_id.0] += child.length;
            near(child.outward_normal[0], original.outward_normal[0], 1e-14);
            near(child.outward_normal[1], original.outward_normal[1], 1e-14);
        }
        for (i, &l) in length.iter().enumerate() {
            near(
                l,
                parent
                    .boundary_edge_geometry(BoundaryEdgeId(i))
                    .unwrap()
                    .length,
                1e-14,
            );
        }
    }
}
#[test]
fn prolongation_satisfies_constraints_and_preserves_piecewise_parent_field() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let values: Vec<_> = (0..mesh.vertices().len())
        .map(|i| ((i * i + 3 * i) % 23) as f64)
        .collect();
    let r = NonconformingRefinement::new(&mesh, &[CellId(0), CellId(7), CellId(12)]).unwrap();
    let fine = r.prolongate_p1(&values).unwrap();
    let c = r.hanging_constraints().unwrap().close().unwrap();
    let reduced: Vec<_> = c.free_dofs().iter().map(|&i| fine[i]).collect();
    let expanded = c.reconstruct(&reduced).unwrap();
    for (a, b) in expanded.iter().zip(&fine) {
        near(*a, *b, 1e-14);
    }
    for (i, cell) in r.cells().iter().enumerate() {
        let parent = r.parent_cell(CellId(i)).unwrap();
        let p = r.cell_points(CellId(i)).unwrap();
        let point = Point2::new(
            (p[0].x + p[1].x + p[2].x) / 3.,
            (p[0].y + p[1].y + p[2].y) / 3.,
        )
        .unwrap();
        let local = cell.iter().map(|v| fine[v.0]).sum::<f64>() / 3.;
        let orig = field::evaluate(&mesh, &values, point, 1e-12)
            .unwrap()
            .unwrap();
        near(local, orig, 1e-10);
        assert!(r.children(parent).unwrap().contains(&CellId(i)));
    }
}
#[test]
fn all_or_no_marks_have_no_hanging_nodes_and_preserve_boundaries() {
    let m = Mesh2::unit_square(2).unwrap();
    for marks in [Vec::new(), (0..m.cells().len()).map(CellId).collect()] {
        let r = NonconformingRefinement::new(&m, &marks).unwrap();
        assert!(r.hanging_nodes().is_empty());
        let c = r.hanging_constraints().unwrap().close().unwrap();
        assert_eq!(c.free_dofs().len(), r.vertices().len());
    }
    let r = NonconformingRefinement::new(&m, &[]).unwrap();
    assert_eq!(r.vertices(), m.vertices());
    assert_eq!(r.cells(), m.cells());
    assert_eq!(r.boundary_edges(), m.boundary_edges());
}
#[test]
fn unrefined_assembly_matches_conforming_implementation() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let r = NonconformingRefinement::new(&mesh, &[]).unwrap();
    let f = |p: Point2| 1. + p.x - 2. * p.y;
    let (a, b) = r.assemble_poisson(f).unwrap();
    let (refa, refb) = numerics::assemble_poisson(&mesh, f).unwrap();
    for i in 0..mesh.vertices().len() {
        near(b[i], refb[i], 1e-13);
        for j in 0..mesh.vertices().len() {
            near(a.get(i, j).unwrap(), refa.get(i, j).unwrap(), 1e-12);
        }
    }
}
#[test]
fn affine_patch_solution_is_continuous_and_matches_manufactured_field() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let r = NonconformingRefinement::new(&mesh, &[CellId(0), CellId(10), CellId(17)]).unwrap();
    assert!(!r.hanging_nodes().is_empty());
    let exact = |p: Point2| 1. + 2. * p.x - 3. * p.y;
    let (a, b) = r.assemble_poisson(|_| 0.).unwrap();
    let mut builder = r.hanging_constraints().unwrap();
    for (i, value) in r.interpolate_boundary(exact).unwrap() {
        builder.add(i, value, &[]).unwrap();
    }
    let c = builder.close().unwrap();
    let (ar, br) = c.condense(&a, &b).unwrap();
    assert!(!c.free_dofs().is_empty());
    let (y, report) = conjugate_gradient(&ar, &br, opts()).unwrap();
    assert!(report.residual_norm < 1e-10);
    let u = c.reconstruct(&y).unwrap();
    for (&value, &p) in u.iter().zip(r.vertices()) {
        near(value, exact(p), 1e-10);
    }
    for h in r.hanging_nodes() {
        near(
            u[h.vertex.0],
            0.5 * (u[h.masters[0].0] + u[h.masters[1].0]),
            1e-14,
        );
        for j in 0..=10 {
            let t = j as f64 / 10.;
            let left = (1. - t) * u[h.masters[0].0] + t * u[h.masters[1].0];
            let right = if t <= 0.5 {
                (1. - 2. * t) * u[h.masters[0].0] + 2. * t * u[h.vertex.0]
            } else {
                (2. - 2. * t) * u[h.vertex.0] + (2. * t - 1.) * u[h.masters[1].0]
            };
            near(left, right, 1e-12);
        }
    }
}
#[test]
fn boundary_and_cell_physical_groups_follow_parent_ids() {
    use dealii_rust::gmsh;
    let imported = gmsh::parse(include_str!("fixtures/square41.msh")).unwrap();
    let r = NonconformingRefinement::new(imported.mesh(), &[CellId(0)]).unwrap();
    for i in 0..r.cells().len() {
        let parent = r.parent_cell(CellId(i)).unwrap();
        assert_eq!(imported.cell_physical_tags(parent).unwrap(), &[100, 101]);
    }
    for i in 0..r.boundary_edges().len() {
        let parent = r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
        assert!(!imported.boundary_physical_tags(parent).unwrap().is_empty());
    }
}
#[test]
fn inputs_and_nonrepresentable_geometry_are_rejected() {
    let mesh = Mesh2::unit_square(1).unwrap();
    assert!(NonconformingRefinement::new(&mesh, &[CellId(0), CellId(0)]).is_err());
    assert!(NonconformingRefinement::new(&mesh, &[CellId(2)]).is_err());
    let r = NonconformingRefinement::new(&mesh, &[CellId(0)]).unwrap();
    assert!(r.children(CellId(2)).is_err());
    assert!(r.parent_cell(CellId(r.cells().len())).is_err());
    assert!(r.cell_points(CellId(r.cells().len())).is_err());
    assert!(r
        .boundary_edge_geometry(BoundaryEdgeId(r.boundary_edges().len()))
        .is_err());
    assert!(r
        .parent_boundary_edge(BoundaryEdgeId(r.boundary_edges().len()))
        .is_err());
    assert!(r.prolongate_p1(&[]).is_err());
    assert!(r
        .prolongate_p1(&vec![f64::NAN; mesh.vertices().len()])
        .is_err());
    assert!(r.interpolate_boundary(|_| f64::NAN).is_err());
    assert!(r.assemble_poisson(|_| f64::INFINITY).is_err());
    let p = |x, y| Point2::new(x, y).unwrap();
    let x = 1e16;
    let tiny = Mesh2::new(
        vec![p(x, x), p(x + 2., x), p(x, x + 2.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    assert!(NonconformingRefinement::new(&tiny, &[CellId(0)]).is_err());
}

#[test]
fn quadratic_solution_converges_on_one_level_hanging_meshes() {
    use dealii_rust::{
        fe::{self, AffineTriangle},
        quadrature::triangle_degree_5,
    };
    let exact = |p: Point2| p.x * p.x + p.y * p.y;
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let mesh = Mesh2::unit_square(n).unwrap();
        let marks: Vec<_> = (0..mesh.cells().len()).step_by(2).map(CellId).collect();
        let r = NonconformingRefinement::new(&mesh, &marks).unwrap();
        assert!(!r.hanging_nodes().is_empty());
        let (a, b) = r.assemble_poisson(|_| -4.).unwrap();
        let mut builder = r.hanging_constraints().unwrap();
        for (i, v) in r.interpolate_boundary(exact).unwrap() {
            builder.add(i, v, &[]).unwrap();
        }
        let c = builder.close().unwrap();
        let (a, b) = c.condense(&a, &b).unwrap();
        let y = if c.free_dofs().is_empty() {
            Vec::new()
        } else {
            conjugate_gradient(&a, &b, opts()).unwrap().0
        };
        let u = c.reconstruct(&y).unwrap();
        let mut l2 = 0.;
        let mut h1 = 0.;
        for (i, cell) in r.cells().iter().enumerate() {
            let map = AffineTriangle::new(r.cell_points(CellId(i)).unwrap()).unwrap();
            let gradients = map.gradients();
            let mut gradient = [0.; 2];
            for j in 0..3 {
                for (d, component) in gradient.iter_mut().enumerate() {
                    *component += u[cell[j].0] * gradients[j][d];
                }
            }
            for q in triangle_degree_5() {
                let p = map.map(q.point);
                let shape = fe::shape_values(q.point);
                let val = (0..3).map(|j| shape[j] * u[cell[j].0]).sum::<f64>();
                let w = map.determinant() * q.weight;
                l2 += (val - exact(p)).powi(2) * w;
                h1 += ((gradient[0] - 2. * p.x).powi(2) + (gradient[1] - 2. * p.y).powi(2)) * w;
            }
        }
        errors.push((l2.sqrt(), h1.sqrt()));
    }
    for pair in errors.windows(2) {
        assert!(pair[0].0 / pair[1].0 > 3.2, "{:?}", errors);
        assert!(pair[0].1 / pair[1].1 > 1.8, "{:?}", errors);
    }
}

#[test]
fn disconnected_components_refine_independently() {
    let p = |x, y| Point2::new(x, y).unwrap();
    let mesh = Mesh2::new(
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
    .unwrap();
    let r = NonconformingRefinement::new(&mesh, &[CellId(0)]).unwrap();
    assert_eq!(r.children(CellId(0)).unwrap().len(), 4);
    assert_eq!(r.children(CellId(1)).unwrap().len(), 1);
    assert!(r.hanging_nodes().is_empty());
    assert_eq!(r.boundary_edges().len(), 9);
    assert_eq!(
        r.cells()[r.children(CellId(1)).unwrap()[0].0],
        mesh.cells()[1]
    );
}

#[test]
fn specialized_boundary_does_not_treat_interface_as_exterior() {
    let parent = Mesh2::unit_square(1).unwrap();
    let r = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let hanging = r.hanging_nodes()[0].vertex;
    // Deliberately demonstrate why Mesh2's boundary algorithm cannot be
    // used here: it sees unmatched coarse and fine edges as boundary edges.
    let unsuitable = Mesh2::new(r.vertices().to_vec(), r.cells().to_vec()).unwrap();
    assert!(unsuitable.boundary_vertices().contains(&hanging));
    assert!(!r.boundary_vertices().contains(&hanging));
    assert!(unsuitable.boundary_edges().len() > r.boundary_edges().len());
}
