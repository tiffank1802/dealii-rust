use dealii_rust::{
    base::Point2,
    gmsh,
    grid::{BoundaryEdgeId, CellId, LocalRefinement, VertexId},
    lac::{conjugate_gradient, CgOptions},
    numerics,
};
const MESH: &str = include_str!("fixtures/square22.msh");
const MESHIO: &str = include_str!("fixtures/meshio_square22.msh");
fn reject(text: &str, contains: &str) {
    let error = gmsh::parse_msh22(text).unwrap_err().to_string();
    assert!(
        error.contains(contains),
        "{} does not contain {}",
        error,
        contains
    );
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{} != {}", a, b);
}
fn exact(p: Point2) -> f64 {
    1. + p.x + 2. * p.y
}
#[test]
fn reads_noncontiguous_node_and_element_tags_and_physical_groups() {
    let imported = gmsh::parse_msh22(MESH).unwrap();
    let mesh = imported.mesh();
    assert_eq!(mesh.vertices().len(), 5);
    assert_eq!(mesh.cells().len(), 4);
    assert_eq!(mesh.boundary_edges().len(), 4);
    assert_eq!(
        (0..5)
            .map(|i| imported.node_tag(VertexId(i)).unwrap())
            .collect::<Vec<_>>(),
        vec![50, 10, 90, 30, 70]
    );
    assert_eq!(imported.vertex_physical_tags(VertexId(1)).unwrap(), &[5]);
    assert_eq!(imported.physical_name(0, 5), Some("corner"));
    assert_eq!(imported.cell_element_tag(CellId(0)).unwrap(), 600);
    assert_eq!(imported.cells_with_tag(100).len(), 4);
    assert_eq!(imported.physical_name(2, 100), Some("domain"));
    for tag in [10, 20, 30, 40] {
        let edges = imported.boundary_edges_with_tag(tag);
        assert_eq!(edges.len(), 1);
        assert_eq!(imported.boundary_physical_tags(edges[0]).unwrap(), &[tag]);
        assert!(imported.boundary_element_tag(edges[0]).unwrap().is_some());
    }
    assert!(imported.node_tag(VertexId(5)).is_err());
    assert!(imported.cell_physical_tags(CellId(4)).is_err());
    assert!(imported.boundary_physical_tags(BoundaryEdgeId(4)).is_err());
}
#[test]
fn meshio_generated_fixture_preserves_independently_written_groups() {
    let imported = gmsh::read_msh22(MESHIO.as_bytes()).unwrap();
    assert_eq!(imported.mesh().vertices().len(), 5);
    assert_eq!(imported.mesh().cells().len(), 4);
    for (tag, name) in [(10, "bottom"), (20, "right"), (30, "top"), (40, "left")] {
        assert_eq!(imported.physical_name(1, tag), Some(name));
        assert_eq!(imported.boundary_edges_with_tag(tag).len(), 1);
    }
    assert_eq!(imported.physical_name(2, 100), Some("domain"));
    assert_eq!(imported.cells_with_tag(100).len(), 4);
}
#[test]
fn mixed_problem_uses_physical_groups_and_solves_affine_solution() {
    for text in [MESH, MESHIO] {
        let imported = gmsh::parse_msh22(text).unwrap();
        let mesh = imported.mesh();
        let mut neumann = Vec::new();
        for tag in [10, 20, 30] {
            neumann.extend(imported.boundary_edges_with_tag(tag));
        }
        let (a, mut b) = numerics::assemble_poisson(mesh, |_| 0.).unwrap();
        let rhs = numerics::assemble_neumann(mesh, &neumann, |_, n| n[0] + 2. * n[1]).unwrap();
        for (x, y) in b.iter_mut().zip(rhs) {
            *x += y;
        }
        let constraints: Vec<_> = imported
            .boundary_vertices_with_tag(40)
            .iter()
            .map(|v| (v.0, exact(mesh.vertices()[v.0])))
            .collect();
        let (a, b) = a.with_dirichlet(&b, &constraints).unwrap();
        let (u, _) = conjugate_gradient(&a, &b, CgOptions::default()).unwrap();
        for (u, p) in u.iter().zip(mesh.vertices()) {
            close(*u, exact(*p));
        }
        let r = LocalRefinement::new(mesh, &[CellId(0)]).unwrap();
        let fine = r.prolongate_p1(&u).unwrap();
        for (x, p) in fine.iter().zip(r.mesh().vertices()) {
            close(*x, exact(*p));
        }
        for i in 0..r.mesh().boundary_edges().len() {
            let parent = r.parent_boundary_edge(BoundaryEdgeId(i)).unwrap();
            assert!(!imported.boundary_physical_tags(parent).unwrap().is_empty());
        }
    }
}
#[test]
fn optional_names_zero_one_or_two_tags_and_missing_boundary_lines() {
    let names = MESH
        .split("$PhysicalNames\n")
        .nth(1)
        .unwrap()
        .split("$EndPhysicalNames\n")
        .next()
        .unwrap();
    let text = MESH
        .replace(&format!("$PhysicalNames\n{}$EndPhysicalNames\n", names), "")
        .replace("200 1 2 10 1 10 30", "200 1 0 10 30")
        .replace("300 1 2 20 2 30 50", "300 1 1 20 30 50")
        .replace("400 1 2 30 3 50 90", "400 1 2 0 3 50 90")
        .replace("9\n100 15", "8\n100 15")
        .replace("500 1 2 40 4 90 10\n", "");
    let m = gmsh::parse_msh22(&text).unwrap();
    assert_eq!(m.physical_name(1, 20), None);
    assert_eq!(m.boundary_edges_with_tag(20).len(), 1);
    assert_eq!(m.boundary_edges_with_tag(10).len(), 0);
    assert_eq!(m.boundary_edges_with_tag(30).len(), 0);
    assert_eq!(
        (0..m.mesh().boundary_edges().len())
            .filter(|&i| m.boundary_element_tag(BoundaryEdgeId(i)).unwrap().is_none())
            .count(),
        1
    );
    assert_eq!(
        (0..m.mesh().boundary_edges().len())
            .filter(
                |&i| m.boundary_element_tag(BoundaryEdgeId(i)).unwrap().is_some()
                    && m.boundary_physical_tags(BoundaryEdgeId(i))
                        .unwrap()
                        .is_empty()
            )
            .count(),
        2
    );
}
#[test]
fn rejects_wrong_version_binary_additional_sections_and_non_utf8() {
    reject(&MESH.replace("2.2 0 8", "4.1 0 8"), "2.2 ASCII");
    reject(&MESH.replace("2.2 0 8", "2.2 1 8"), "2.2 ASCII");
    reject(&MESH.replace("2.2 0 8", "2.2 0 4"), "2.2 ASCII");
    reject(
        &format!("{}\n$Entities\n0 0 0 0\n$EndEntities\n", MESH),
        "unsupported in MSH 2.2",
    );
    reject(
        &format!("{}\n$Periodic\n0\n$EndPeriodic\n", MESH),
        "unsupported Gmsh section",
    );
    assert!(gmsh::read_msh22(&[0xff, 0xfe][..]).is_err());
    assert!(gmsh::parse(MESH).is_err());
    assert!(gmsh::parse_msh22(include_str!("fixtures/square41.msh")).is_err());
}
#[test]
fn rejects_invalid_nodes_counts_and_coordinates() {
    reject(&MESH.replace("50 1 1 0", "10 1 1 0"), "duplicate node tag");
    reject(&MESH.replace("50 1 1 0", "0 1 1 0"), "must be positive");
    reject(&MESH.replace("50 1 1 0", "50 1 1 1"), "planar z=0");
    reject(&MESH.replace("50 1 1 0", "50 NaN 1 0"), "nonfinite");
    reject(&MESH.replace("$Nodes\n5\n", "$Nodes\n6\n"), "truncated");
    reject(
        &MESH.replace("$Nodes\n5\n", "$Nodes\n4\n"),
        "unexpected extra",
    );
    reject(&MESH.replace("50 1 1 0", "50 1 1 0 99"), "unexpected extra");
}
#[test]
fn rejects_unsupported_bad_elements_and_extra_tag_data() {
    reject(
        &MESH.replace("600 2 2 100 1 10 30 70", "600 3 2 100 1 10 30 70"),
        "unsupported element type",
    );
    reject(
        &MESH.replace("600 2 2 100 1 10 30 70", "600 2 3 100 1 0 10 30 70"),
        "more than two tags",
    );
    reject(
        &MESH.replace("600 2 2 100 1 10 30 70", "600 2 2 100 1 10 30 999"),
        "unknown node",
    );
    reject(
        &MESH.replace("700 2 2 100 1 30 50 70", "600 2 2 100 1 30 50 70"),
        "duplicate element tag",
    );
    reject(
        &MESH.replace("700 2 2 100 1 30 50 70", "700 2 2 100 1 30 70 50"),
        "invalid Gmsh triangle mesh",
    );
    reject(
        &MESH.replace("$Elements\n9\n", "$Elements\n10\n"),
        "truncated",
    );
    reject(
        &MESH.replace("$Elements\n9\n", "$Elements\n8\n"),
        "unexpected extra",
    );
    reject(
        &MESH.replace("500 1 2 40 4 90 10", "500 1 2 40 4 10 30"),
        "duplicate line",
    );
    reject(
        &MESH.replace("500 1 2 40 4 90 10", "500 1 2 40 4 10 70"),
        "not an exterior",
    );
    reject(
        &MESH.replace("500 1 2 40 4 90 10", "500 1 2 40 4 10 10"),
        "degenerate line",
    );
}
#[test]
fn preserves_duplicate_point_detection_and_rejects_nonconforming_records() {
    let duplicate = MESH
        .replace("$Elements\n9\n", "$Elements\n10\n")
        .replace("200 1 2 10 1 10 30", "101 15 1 6 10\n200 1 2 10 1 10 30");
    reject(&duplicate, "duplicate point");
    reject(
        &MESH.replace("100 15 2 5 1 10", "100 15 2 -1 1 10"),
        "invalid physical tag",
    );
    reject(
        &MESH.replace("100 15 2 5 1 10", "100 15 2 5 -1 10"),
        "invalid elementary entity tag",
    );
    reject(
        &MESH.replace("100 15 2 5 1 10", "100 15 2 5 1 10 42"),
        "unexpected extra",
    );
}
