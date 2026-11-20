use dealii_rust::{
    base::triangle_determinant,
    gmsh,
    grid::{BoundaryEdgeId, CellId, VertexId},
};

const SQUARE: &str = include_str!("fixtures/square41.msh");
fn reject(text: &str, message: &str) {
    let err = gmsh::parse(text).unwrap_err();
    assert!(
        err.to_string().contains(message),
        "unexpected error: {}",
        err
    );
}

#[test]
fn reads_noncontiguous_tags_physical_groups_and_geometry() {
    let imported = gmsh::parse(SQUARE).unwrap();
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
    assert!(imported
        .vertex_physical_tags(VertexId(0))
        .unwrap()
        .is_empty());
    for i in 0..4 {
        assert_eq!(imported.cell_physical_tags(CellId(i)).unwrap(), &[100, 101]);
    }
    assert_eq!(imported.cell_element_tag(CellId(0)).unwrap(), 600);
    assert_eq!(imported.cells_with_tag(100).len(), 4);
    assert!(imported.cells_with_tag(999).is_empty());
    assert_eq!(imported.physical_name(2, 101), Some("material A"));
    assert_eq!(imported.physical_name(1, 40), Some("left"));
    assert_eq!(imported.physical_name(2, 40), None);
    for (tag, expected) in [(10, "bottom"), (20, "right"), (30, "top"), (40, "left")] {
        let edges = imported.boundary_edges_with_tag(tag);
        assert_eq!(edges.len(), 1, "{}", expected);
        assert_eq!(imported.boundary_physical_tags(edges[0]).unwrap(), &[tag]);
        assert!(imported.boundary_element_tag(edges[0]).unwrap().is_some());
        assert_eq!(imported.boundary_vertices_with_tag(tag).len(), 2);
    }
    let left = imported.boundary_vertices_with_tag(40);
    assert!(left.iter().all(|v| mesh.vertices()[v.0].x == 0.0));
    let area: f64 = (0..mesh.cells().len())
        .map(|i| triangle_determinant(mesh.cell_points(CellId(i)).unwrap()).unwrap() / 2.0)
        .sum();
    assert!((area - 1.0).abs() < 1e-14);
    assert!(imported.node_tag(VertexId(5)).is_err());
    assert!(imported.cell_physical_tags(CellId(4)).is_err());
    assert!(imported.boundary_physical_tags(BoundaryEdgeId(4)).is_err());
}

#[test]
fn read_accepts_streams_and_optional_physical_names() {
    let without_names = SQUARE.replace("$PhysicalNames\n7\n0 5 \"corner\"\n1 10 \"bottom\"\n1 20 \"right\"\n1 30 \"top\"\n1 40 \"left\"\n2 100 \"domain\"\n2 101 \"material A\"\n$EndPhysicalNames\n", "");
    let imported = gmsh::read(without_names.as_bytes()).unwrap();
    assert_eq!(imported.physical_name(1, 40), None);
    assert_eq!(imported.boundary_edges_with_tag(40).len(), 1);
    assert!(gmsh::read(&[0xff, 0xfe][..]).is_err());
}

#[test]
fn rejects_unsupported_versions_and_sections() {
    reject(&SQUARE.replace("4.1 0 8", "2.2 0 8"), "4.1 ASCII");
    reject(&SQUARE.replace("4.1 0 8", "4.1 1 8"), "4.1 ASCII");
    reject(&SQUARE.replace("4.1 0 8", "4.1 0 4"), "4.1 ASCII");
    reject(
        &format!("{}\n$Periodic\n0\n$EndPeriodic\n", SQUARE),
        "unsupported Gmsh section $Periodic",
    );
    reject(&SQUARE.replace("$EndNodes", "$EndFoo"), "missing $EndNodes");
    reject(
        &SQUARE.replace(
            "$Entities\n",
            "$Entities\n4 4 1 0\n$EndEntities\n$Entities\n",
        ),
        "duplicate Gmsh section $Entities",
    );
}

#[test]
fn rejects_unsupported_or_malformed_nodes() {
    reject(
        &SQUARE.replace("2 1 0 1\n70\n", "2 1 1 1\n70\n"),
        "parametric nodes",
    );
    reject(&SQUARE.replace("0.5 0.5 0", "0.5 0.5 0.01"), "planar z=0");
    reject(&SQUARE.replace("0.5 0.5 0", "nan 0.5 0"), "nonfinite");
    reject(
        &SQUARE.replace("0 2 0 1\n30\n", "0 2 0 1\n10\n"),
        "duplicate node tag",
    );
    reject(&SQUARE.replace("5 5 10 90", "5 5 1 90"), "tag range");
    reject(
        &SQUARE.replace("0 3 0 1\n50\n", "0 99 0 1\n50\n"),
        "unknown entity",
    );
    reject(
        &SQUARE.replace("0 3 0 1\n50\n", "0 3 0 1000000\n50\n"),
        "block exceeds",
    );
}

#[test]
fn rejects_invalid_elements_and_connectivity() {
    reject(
        &SQUARE.replace("2 1 2 4\n600", "2 1 3 4\n600"),
        "unsupported dimension 2 / element type 3",
    );
    reject(
        &SQUARE.replace("600 10 30 70", "600 10 30 999"),
        "unknown node 999",
    );
    reject(
        &SQUARE.replace("600 10 30 70", "600 10 70 30"),
        "triangle must be counterclockwise",
    );
    reject(
        &SQUARE.replace("200 10 30", "200 10 70"),
        "not an exterior mesh edge",
    );
    reject(
        &SQUARE.replace("300 30 50", "300 10 30"),
        "duplicate line element",
    );
    reject(
        &SQUARE.replace("300 30 50", "200 30 50"),
        "duplicate element tag",
    );
    reject(
        &SQUARE.replace("100 10\n1 1 1 1", "100 30\n1 1 1 1"),
        "point node belongs to a different entity",
    );
    reject(
        &SQUARE.replace("6 9 100 900", "6 10 100 900"),
        "declared count",
    );
}

#[test]
fn rejects_invalid_entities_and_physical_names() {
    reject(
        &SQUARE.replace("4 4 1 0\n1 0 0 0 1 5", "4 4 1 1\n1 0 0 0 1 5"),
        "3D volumes",
    );
    reject(
        &SQUARE.replace("2 1 -2\n2 1 0 0", "2 1 -99\n2 1 0 0"),
        "unknown bounding entity",
    );
    reject(
        &SQUARE.replace("1 20 \"right\"", "1 10 \"right\""),
        "duplicate dimension/tag",
    );
    reject(
        &SQUARE.replace("2 101 \"material A\"", "2 101 material A"),
        "expected quoted name",
    );
}

#[test]
fn untagged_unlisted_boundary_and_escaped_names_are_explicit() {
    // The left boundary edge is not listed as a line element: its tags are empty.
    let text = SQUARE
        .replace("6 9 100 900", "5 8 100 900")
        .replace("1 4 1 1\n500 90 10\n", "")
        .replace("\"material A\"", "\"material \\\"A\\\"\"");
    let imported = gmsh::parse(&text).unwrap();
    assert!(imported.boundary_edges_with_tag(40).is_empty());
    assert_eq!(imported.physical_name(2, 101), Some("material \"A\""));
    let untagged = (0..imported.mesh().boundary_edges().len())
        .filter(|&i| {
            imported
                .boundary_physical_tags(BoundaryEdgeId(i))
                .unwrap()
                .is_empty()
        })
        .count();
    assert_eq!(untagged, 1);
    let unlisted = (0..4)
        .filter(|&i| {
            imported
                .boundary_element_tag(BoundaryEdgeId(i))
                .unwrap()
                .is_none()
        })
        .count();
    assert_eq!(unlisted, 1);
}

#[test]
fn tagged_mesh_solves_mixed_poisson_end_to_end() {
    use dealii_rust::{
        lac::{conjugate_gradient, CgOptions},
        numerics,
    };
    let imported = gmsh::parse(SQUARE).unwrap();
    let mesh = imported.mesh();
    let (matrix, mut rhs) = numerics::assemble_poisson(mesh, |_| 0.0).unwrap();
    let mut edges = Vec::new();
    for tag in [10, 20, 30] {
        edges.extend(imported.boundary_edges_with_tag(tag));
    }
    let load = numerics::assemble_neumann(mesh, &edges, |_, n| n[0] + 2.0 * n[1]).unwrap();
    for (b, g) in rhs.iter_mut().zip(load) {
        *b += g;
    }
    let constraints: Vec<_> = imported
        .boundary_vertices_with_tag(40)
        .into_iter()
        .map(|v| {
            let p = mesh.vertices()[v.0];
            (v.0, 1.0 + p.x + 2.0 * p.y)
        })
        .collect();
    let (matrix, rhs) = matrix.with_dirichlet(&rhs, &constraints).unwrap();
    let (solution, _) = conjugate_gradient(&matrix, &rhs, CgOptions::default()).unwrap();
    for (value, p) in solution.iter().zip(mesh.vertices()) {
        assert!((*value - (1.0 + p.x + 2.0 * p.y)).abs() < 1e-12);
    }
}
