use dealii_rust::{
    base::{triangle_determinant, Point2},
    estimator, gmsh,
    grid::{BoundaryEdgeId, CellId, Mesh2, UniformRefinement},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};
use std::f64::consts::PI;

fn area(mesh: &Mesh2, id: CellId) -> f64 {
    triangle_determinant(mesh.cell_points(id).unwrap()).unwrap() * 0.5
}
fn close(x: f64, y: f64) {
    assert!((x - y).abs() < 1e-12, "{} != {}", x, y);
}

#[test]
fn uniform_refinement_preserves_topology_area_and_boundary_lineage() {
    for n in [1, 2, 4] {
        let coarse = Mesh2::unit_square(n).unwrap();
        let refined = UniformRefinement::new(&coarse).unwrap();
        let fine = refined.mesh();
        assert_eq!(fine.vertices().len(), (2 * n + 1) * (2 * n + 1));
        assert_eq!(fine.cells().len(), coarse.cells().len() * 4);
        assert_eq!(
            fine.boundary_edges().len(),
            coarse.boundary_edges().len() * 2
        );
        assert_eq!(
            &fine.vertices()[..coarse.vertices().len()],
            coarse.vertices()
        );
        for i in 0..coarse.cells().len() {
            let children = refined.children(CellId(i)).unwrap();
            let child_area: f64 = children
                .iter()
                .map(|&c| {
                    assert_eq!(refined.parent_cell(c).unwrap(), CellId(i));
                    area(fine, c)
                })
                .sum();
            close(child_area, area(&coarse, CellId(i)));
        }
        let mut counts = vec![0; coarse.boundary_edges().len()];
        for i in 0..fine.boundary_edges().len() {
            let child = fine.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
            let parent = refined.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
            let old = coarse.boundary_edge_geometry(parent).unwrap();
            counts[parent.0] += 1;
            close(child.length * 2.0, old.length);
            close(child.outward_normal[0], old.outward_normal[0]);
            close(child.outward_normal[1], old.outward_normal[1]);
        }
        assert!(counts.iter().all(|&c| c == 2));
        assert!(refined.children(CellId(coarse.cells().len())).is_err());
        assert!(refined.parent_cell(CellId(fine.cells().len())).is_err());
        assert!(refined
            .parent_boundary_edge(BoundaryEdgeId(fine.boundary_edges().len()))
            .is_err());
    }
}

#[test]
fn p1_prolongation_preserves_piecewise_field_and_integral() {
    let mesh = Mesh2::unit_square(2).unwrap();
    let values: Vec<_> = mesh.vertices().iter().map(|p| p.x * p.y + p.x).collect();
    let refinement = UniformRefinement::new(&mesh).unwrap();
    let fine = refinement.mesh();
    let result = refinement.prolongate_p1(&values).unwrap();
    let constant = refinement.prolongate_p1(&vec![3.0; values.len()]).unwrap();
    assert!(constant.iter().all(|v| *v == 3.0));
    assert_eq!(&result[..values.len()], values.as_slice());
    for (i, cell) in mesh.cells().iter().enumerate() {
        let parent_integral =
            area(&mesh, CellId(i)) * cell.iter().map(|v| values[v.0]).sum::<f64>() / 3.0;
        let mut child_integral = 0.0;
        for child in refinement.children(CellId(i)).unwrap() {
            let ids = fine.cells()[child.0];
            child_integral +=
                area(fine, child) * ids.iter().map(|v| result[v.0]).sum::<f64>() / 3.0;
        }
        close(child_integral, parent_integral);
    }
    assert!(refinement.prolongate_p1(&[]).is_err());
    assert!(refinement
        .prolongate_p1(&vec![f64::NAN; values.len()])
        .is_err());
}

#[test]
fn imported_mesh_tags_can_be_transferred_explicitly() {
    let imported = gmsh::parse(include_str!("fixtures/square41.msh")).unwrap();
    let refined = UniformRefinement::new(imported.mesh()).unwrap();
    for i in 0..refined.mesh().cells().len() {
        let parent = refined.parent_cell(CellId(i)).unwrap();
        assert_eq!(imported.cell_physical_tags(parent).unwrap(), &[100, 101]);
    }
    for tag in [10, 20, 30, 40] {
        let selected = (0..refined.mesh().boundary_edges().len())
            .filter(|&i| {
                let parent = refined.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
                imported
                    .boundary_physical_tags(parent)
                    .unwrap()
                    .contains(&tag)
            })
            .count();
        assert_eq!(selected, 2);
    }
}

#[test]
fn estimator_analytical_source_and_gradient_jump() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let source_only = estimator::poisson_dirichlet(&mesh, &[0.0; 4], |_| 1.0).unwrap();
    for eta in source_only.cell_indicators {
        close(eta, 1.0);
    }
    close(source_only.global_indicator, 2.0_f64.sqrt());
    // The diagonal jump is sqrt(2), and h_e=sqrt(2). Each cell gets
    // (1/2) h_e² jump² = 2. This tests the local formula algebraically;
    // this synthetic field does not satisfy homogeneous Dirichlet data.
    let jump = estimator::poisson_dirichlet(&mesh, &[0.0, 1.0, 0.0, 0.0], |_| 0.0).unwrap();
    for eta in jump.cell_indicators {
        close(eta, 2.0_f64.sqrt());
    }
    close(jump.global_indicator, 2.0);
    let affine: Vec<_> = mesh
        .vertices()
        .iter()
        .map(|p| 1.0 + 2.0 * p.x - 3.0 * p.y)
        .collect();
    close(
        estimator::poisson_dirichlet(&mesh, &affine, |_| 0.0)
            .unwrap()
            .global_indicator,
        0.0,
    );
}

fn source(p: Point2) -> f64 {
    2.0 * PI * PI * (PI * p.x).sin() * (PI * p.y).sin()
}
#[test]
fn estimator_and_error_decrease_on_uniform_refinement() {
    let mut mesh = Mesh2::unit_square(4).unwrap();
    let mut previous = (f64::INFINITY, f64::INFINITY);
    for _ in 0..3 {
        let (matrix, rhs) = numerics::assemble_poisson(&mesh, source).unwrap();
        let constraints: Vec<_> = mesh
            .boundary_vertices()
            .iter()
            .map(|id| (id.0, 0.0))
            .collect();
        let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints).unwrap();
        let (u, _) = conjugate_gradient(&matrix, &rhs, CgOptions::default()).unwrap();
        let eta = estimator::poisson_dirichlet(&mesh, &u, source).unwrap();
        let (_, h1) = numerics::errors(
            &mesh,
            &u,
            |p| (PI * p.x).sin() * (PI * p.y).sin(),
            |p| {
                [
                    PI * (PI * p.x).cos() * (PI * p.y).sin(),
                    PI * (PI * p.x).sin() * (PI * p.y).cos(),
                ]
            },
        )
        .unwrap();
        assert!(
            eta.global_indicator < previous.0 * 0.65,
            "eta = {}",
            eta.global_indicator
        );
        assert!(h1 < previous.1 * 0.65, "h1 = {}", h1);
        previous = (eta.global_indicator, h1);
        mesh = UniformRefinement::new(&mesh).unwrap().mesh().clone();
    }
}

#[test]
fn estimator_rejects_invalid_numerical_data() {
    let mesh = Mesh2::unit_square(1).unwrap();
    assert!(estimator::poisson_dirichlet(&mesh, &[0.0], |_| 0.0).is_err());
    assert!(estimator::poisson_dirichlet(&mesh, &[f64::NAN; 4], |_| 0.0).is_err());
    assert!(estimator::poisson_dirichlet(&mesh, &[0.0; 4], |_| f64::NAN).is_err());
    assert!(estimator::poisson_dirichlet(&mesh, &[0.0; 4], |_| f64::MAX).is_err());
    assert!(
        estimator::poisson_dirichlet(&mesh, &[f64::MAX, -f64::MAX, 0.0, 0.0], |_| 0.0).is_err()
    );
}

#[test]
fn affine_prolongation_and_skew_disconnected_mesh() {
    use dealii_rust::grid::VertexId;
    let p = |x, y| Point2::new(x, y).unwrap();
    let mesh = Mesh2::new(
        vec![
            p(1., 2.),
            p(4., 2.),
            p(2., 4.),
            p(7., 2.),
            p(9., 2.),
            p(8., 5.),
        ],
        vec![
            [VertexId(0), VertexId(1), VertexId(2)],
            [VertexId(3), VertexId(4), VertexId(5)],
        ],
    )
    .unwrap();
    let field = |p: Point2| 2.0 + 3.0 * p.x - 4.0 * p.y;
    let values: Vec<_> = mesh.vertices().iter().copied().map(field).collect();
    let refined = UniformRefinement::new(&mesh).unwrap();
    let values = refined.prolongate_p1(&values).unwrap();
    assert_eq!(refined.mesh().vertices().len(), 12);
    assert_eq!(refined.mesh().boundary_edges().len(), 12);
    for (value, &point) in values.iter().zip(refined.mesh().vertices()) {
        close(*value, field(point));
    }
    for i in 0..2 {
        let sum: f64 = refined
            .children(CellId(i))
            .unwrap()
            .iter()
            .map(|&c| area(refined.mesh(), c))
            .sum();
        close(sum, area(&mesh, CellId(i)));
    }
}

#[test]
fn prolongation_preserves_every_nonaffine_parent_edge_midpoint() {
    let mesh = Mesh2::unit_square(2).unwrap();
    let values: Vec<_> = (0..mesh.vertices().len()).map(|i| (i * i) as f64).collect();
    let refinement = UniformRefinement::new(&mesh).unwrap();
    let fine = refinement.mesh();
    let prolongated = refinement.prolongate_p1(&values).unwrap();
    for cell in mesh.cells() {
        for &(a, b) in &[(cell[0], cell[1]), (cell[1], cell[2]), (cell[2], cell[0])] {
            let x = (mesh.vertices()[a.0].x + mesh.vertices()[b.0].x) * 0.5;
            let y = (mesh.vertices()[a.0].y + mesh.vertices()[b.0].y) * 0.5;
            let matching: Vec<_> = fine
                .vertices()
                .iter()
                .enumerate()
                .filter(|(_, p)| p.x == x && p.y == y)
                .collect();
            assert_eq!(matching.len(), 1, "midpoint must be shared");
            close(
                prolongated[matching[0].0],
                (values[a.0] + values[b.0]) * 0.5,
            );
        }
    }
}
