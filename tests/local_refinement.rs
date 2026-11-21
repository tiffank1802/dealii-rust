use dealii_rust::{
    base::{triangle_determinant, Point2},
    field, gmsh,
    grid::{BoundaryEdgeId, CellId, LocalRefinement, Mesh2},
};
fn area(mesh: &Mesh2, id: CellId) -> f64 {
    triangle_determinant(mesh.cell_points(id).unwrap()).unwrap() * 0.5
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-11, "{} != {}", a, b);
}
#[test]
fn all_marking_subsets_preserve_area_and_have_valid_boundary_lineage() {
    let parent = Mesh2::unit_square(2).unwrap();
    for mask in 0..(1usize << parent.cells().len()) {
        let marked: Vec<_> = (0..parent.cells().len())
            .filter(|i| mask & (1 << i) != 0)
            .map(CellId)
            .collect();
        let r = LocalRefinement::new(&parent, &marked).unwrap();
        let mesh = r.mesh();
        assert_eq!(
            &mesh.vertices()[..parent.vertices().len()],
            parent.vertices()
        );
        for i in 0..parent.cells().len() {
            let children = r.children(CellId(i)).unwrap();
            assert!((1..=4).contains(&children.len()));
            if marked.contains(&CellId(i)) {
                assert_eq!(children.len(), 4);
            }
            let sum: f64 = children
                .iter()
                .map(|&id| {
                    assert_eq!(r.parent_cell(id).unwrap(), CellId(i));
                    area(mesh, id)
                })
                .sum();
            near(sum, area(&parent, CellId(i)));
        }
        let mut lengths = vec![0.; parent.boundary_edges().len()];
        for i in 0..mesh.boundary_edges().len() {
            let geom = mesh.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
            let old = r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
            let original = parent.boundary_edge_geometry(old).unwrap();
            lengths[old.0] += geom.length;
            near(geom.outward_normal[0], original.outward_normal[0]);
            near(geom.outward_normal[1], original.outward_normal[1]);
        }
        for (i, length) in lengths.iter().enumerate() {
            near(
                *length,
                parent
                    .boundary_edge_geometry(BoundaryEdgeId(i))
                    .unwrap()
                    .length,
            );
        }
    }
}
#[test]
fn one_mark_is_local_and_green_neighbors_are_conforming() {
    let mesh = Mesh2::unit_square(4).unwrap();
    let refined = LocalRefinement::new(&mesh, &[CellId(0)]).unwrap();
    assert_eq!(refined.children(CellId(0)).unwrap().len(), 4);
    assert!(refined.mesh().cells().len() < mesh.cells().len() * 4);
    assert!(
        refined
            .children(CellId(mesh.cells().len() - 1))
            .unwrap()
            .len()
            == 1
    );
    assert!(refined.mesh().cells().len() > mesh.cells().len() + 3);
}
#[test]
fn prolongation_preserves_piecewise_field_on_all_children() {
    let mesh = Mesh2::unit_square(3).unwrap();
    let values: Vec<_> = (0..mesh.vertices().len())
        .map(|i| ((i * i + 3 * i) % 17) as f64)
        .collect();
    let r = LocalRefinement::new(&mesh, &[CellId(0), CellId(4), CellId(10)]).unwrap();
    let fine = r.prolongate_p1(&values).unwrap();
    for (i, cell) in r.mesh().cells().iter().enumerate() {
        let parent = r.parent_cell(CellId(i)).unwrap();
        let original = mesh.cells()[parent.0];
        let parent_integral =
            area(&mesh, parent) * original.iter().map(|v| values[v.0]).sum::<f64>() / 3.;
        let children_integral: f64 = r
            .children(parent)
            .unwrap()
            .iter()
            .map(|id| {
                area(r.mesh(), *id)
                    * r.mesh().cells()[id.0]
                        .iter()
                        .map(|v| fine[v.0])
                        .sum::<f64>()
                    / 3.
            })
            .sum();
        near(children_integral, parent_integral);
        let points = cell.map(|v| r.mesh().vertices()[v.0]);
        let p = Point2::new(
            (points[0].x + points[1].x + points[2].x) / 3.,
            (points[0].y + points[1].y + points[2].y) / 3.,
        )
        .unwrap();
        near(
            field::evaluate(&mesh, &values, p, 1e-12).unwrap().unwrap(),
            field::evaluate(r.mesh(), &fine, p, 1e-12).unwrap().unwrap(),
        );
    }
}
#[test]
fn no_marks_is_identity_and_invalid_inputs_are_rejected() {
    let mesh = Mesh2::unit_square(2).unwrap();
    let r = LocalRefinement::new(&mesh, &[]).unwrap();
    assert_eq!(r.mesh().vertices(), mesh.vertices());
    assert_eq!(r.mesh().cells(), mesh.cells());
    for i in 0..mesh.cells().len() {
        assert_eq!(r.children(CellId(i)).unwrap(), &[CellId(i)]);
    }
    assert!(LocalRefinement::new(&mesh, &[CellId(0), CellId(0)]).is_err());
    assert!(LocalRefinement::new(&mesh, &[CellId(mesh.cells().len())]).is_err());
    assert!(r.children(CellId(mesh.cells().len())).is_err());
    assert!(r.parent_cell(CellId(mesh.cells().len())).is_err());
    assert!(r
        .parent_boundary_edge(BoundaryEdgeId(mesh.boundary_edges().len()))
        .is_err());
    assert!(r.prolongate_p1(&[]).is_err());
    assert!(r
        .prolongate_p1(&vec![f64::NAN; mesh.vertices().len()])
        .is_err());
}
#[test]
fn boundary_physical_tags_can_be_transferred_and_repeated_refinement_works() {
    let imported = gmsh::parse(include_str!("fixtures/square41.msh")).unwrap();
    let r = LocalRefinement::new(imported.mesh(), &[CellId(0)]).unwrap();
    for i in 0..r.mesh().boundary_edges().len() {
        let parent = r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
        assert!(!imported.boundary_physical_tags(parent).unwrap().is_empty());
    }
    let s = LocalRefinement::new(r.mesh(), &[CellId(0)]).unwrap();
    let vals: Vec<_> = imported
        .mesh()
        .vertices()
        .iter()
        .map(|p| 3. + 2. * p.x - 4. * p.y)
        .collect();
    let fine = s.prolongate_p1(&r.prolongate_p1(&vals).unwrap()).unwrap();
    for (v, p) in fine.iter().zip(s.mesh().vertices()) {
        near(*v, 3. + 2. * p.x - 4. * p.y);
    }
}

#[test]
fn indicator_driven_poisson_refinement_solves_and_reduces_error() {
    use dealii_rust::{
        estimator,
        lac::{conjugate_gradient, CgOptions},
        numerics,
    };
    use std::f64::consts::PI;
    let exact = |p: Point2| (PI * p.x).sin() * (PI * p.y).sin();
    let source = |p: Point2| 2. * PI * PI * exact(p);
    let grad = |p: Point2| {
        [
            PI * (PI * p.x).cos() * (PI * p.y).sin(),
            PI * (PI * p.x).sin() * (PI * p.y).cos(),
        ]
    };
    let mut mesh = Mesh2::unit_square(4).unwrap();
    let mut errors = Vec::new();
    for level in 0..3 {
        let (a, b) = numerics::assemble_poisson(&mesh, source).unwrap();
        let constraints: Vec<_> = mesh.boundary_vertices().iter().map(|v| (v.0, 0.)).collect();
        let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
        let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
        errors.push(numerics::errors(&mesh, &u, exact, grad).unwrap());
        let e = estimator::poisson_dirichlet(&mesh, &u, source).unwrap();
        if level < 2 {
            let marked = estimator::mark_bulk(&e.cell_indicators, 0.5).unwrap();
            let refined = LocalRefinement::new(&mesh, &marked).unwrap();
            assert!(refined.mesh().cells().len() < mesh.cells().len() * 4);
            mesh = refined.mesh().clone();
        }
    }
    for pair in errors.windows(2) {
        assert!(pair[1].0 < pair[0].0);
        assert!(pair[1].1 < pair[0].1);
    }
}
