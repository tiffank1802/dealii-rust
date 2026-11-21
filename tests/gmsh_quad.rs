use dealii_rust::{
    gmsh,
    grid::{BoundaryEdgeId, CellId, VertexId},
    quad,
};

const QUAD: &str = r#"$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
7
0 5 "corner"
1 10 "bottom"
1 20 "right"
1 30 "top"
1 40 "left"
2 100 "domain"
2 101 "material A"
$EndPhysicalNames
$Entities
4 4 1 0
1 0 0 0 1 5
2 1 0 0 1 5
3 1 1 0 1 5
4 0 1 0 1 5
1 0 0 0 1 0 0 1 10 2 1 -2
2 1 0 0 1 1 0 1 20 2 2 -3
3 0 1 0 1 1 0 1 30 2 3 -4
4 0 0 0 0 1 0 1 40 2 4 -1
1 0 0 0 1 1 0 2 100 101 4 1 2 3 4
$EndEntities
$Nodes
4 4 10 40
0 1 0 1
10
0 0 0
0 2 0 1
20
1 0 0
0 3 0 1
30
1 1 0
0 4 0 1
40
0 1 0
$EndNodes
$Elements
5 5 200 600
1 1 1 1
200 10 20
1 2 1 1
300 20 30
1 3 1 1
400 30 40
1 4 1 1
500 40 10
2 1 3 1
600 10 20 30 40
$EndElements
"#;

fn reject(text: &str, message: &str) {
    let err = gmsh::parse_msh41_quads(text).unwrap_err();
    assert!(
        err.to_string().contains(message),
        "unexpected error: {}",
        err
    );
}

#[test]
fn reads_quad_mesh() {
    let imported = gmsh::parse_msh41_quads(QUAD).unwrap();
    let mesh = imported.mesh();
    assert_eq!(mesh.vertices().len(), 4);
    assert_eq!(mesh.cells().len(), 1);
    assert_eq!(mesh.boundary_edges().len(), 4);
    assert_eq!(
        (0..4)
            .map(|i| imported.node_tag(VertexId(i)).unwrap())
            .collect::<Vec<_>>(),
        vec![10, 20, 30, 40]
    );
    for i in 0..4 {
        assert_eq!(imported.vertex_physical_tags(VertexId(i)).unwrap(), &[5]);
    }
    assert_eq!(imported.cell_physical_tags(CellId(0)).unwrap(), &[100, 101]);
    assert_eq!(imported.cell_element_tag(CellId(0)).unwrap(), 600);
    assert_eq!(imported.cells_with_tag(100).len(), 1);
    assert!(imported.cells_with_tag(999).is_empty());
    assert_eq!(imported.physical_name(2, 101), Some("material A"));
    for (tag, expected) in [(10, "bottom"), (20, "right"), (30, "top"), (40, "left")] {
        let edges = imported.boundary_edges_with_tag(tag);
        assert_eq!(edges.len(), 1, "{}", expected);
        assert_eq!(imported.boundary_physical_tags(edges[0]).unwrap(), &[tag]);
        assert!(imported.boundary_element_tag(edges[0]).unwrap().is_some());
        assert_eq!(imported.boundary_vertices_with_tag(tag).len(), 2);
    }
    let left = imported.boundary_vertices_with_tag(40);
    assert!(left.iter().all(|v| mesh.vertices()[v.0].x == 0.0));
    assert!(imported.node_tag(VertexId(4)).is_err());
    assert!(imported.cell_physical_tags(CellId(1)).is_err());
    assert!(imported.boundary_physical_tags(BoundaryEdgeId(4)).is_err());
}

#[test]
fn read_stream() {
    let imported = gmsh::read_msh41_quads(QUAD.as_bytes()).unwrap();
    assert_eq!(imported.mesh().cells().len(), 1);
    assert!(gmsh::read_msh41_quads(&[0xff, 0xfe][..]).is_err());
}

#[test]
fn rejects_wrong_element_type() {
    // triangle type 2 should be rejected
    let text = QUAD.replace("2 1 3 1\n600 10 20 30 40", "2 1 2 4\n600 10 20 30 40");
    reject(&text, "unsupported dimension 2 / element type 2");
}

#[test]
fn rejects_triangles_points_only() {
    // remove quad cell block -> no cells -> mesh must contain vertices and cells
    let text = QUAD.replace("2 1 3 1\n600 10 20 30 40\n", "");
    // element count becomes 4
    let text = text.replace("5 5 200 600", "4 4 200 500");
    reject(&text, "mesh must contain vertices and cells");
}

#[test]
fn rejects_3d_and_binary() {
    reject(&QUAD.replace("4.1 0 8", "4.1 1 8"), "4.1 ASCII");
    reject(&QUAD.replace("4.1 0 8", "4.1 0 4"), "4.1 ASCII");
}

#[test]
fn rejects_partitions_and_extra_sections() {
    reject(
        &format!(
            "{}\n$PartitionedEntities\n0\n$EndPartitionedEntities\n",
            QUAD
        ),
        "unsupported Gmsh section $PartitionedEntities",
    );
}

#[test]
fn rejects_degenerate_quad() {
    let text = QUAD.replace("600 10 20 30 40", "600 10 10 30 40");
    reject(&text, "repeated cell vertex");
}

#[test]
fn rejects_clockwise_quad() {
    let text = QUAD.replace("600 10 20 30 40", "600 10 40 30 20");
    reject(
        &text,
        "quadrilateral must be convex, counterclockwise and nondegenerate",
    );
}

#[test]
fn rejects_unknown_node() {
    reject(
        &QUAD.replace("600 10 20 30 40", "600 10 20 30 999"),
        "unknown node 999",
    );
}

#[test]
fn rejects_line_not_exterior() {
    // point a line across diagonal that is internal
    let text = QUAD.replace("200 10 20", "200 10 30");
    reject(&text, "not an exterior mesh edge");
}

#[test]
fn rejects_duplicate_line() {
    reject(
        &QUAD.replace("300 20 30", "300 10 20"),
        "duplicate line element",
    );
}

#[test]
fn rejects_parametric_nodes() {
    reject(&QUAD.replace("0 1 0 1", "0 1 1 1"), "parametric nodes");
}

#[test]
fn rejects_nonplanar() {
    reject(&QUAD.replace("0 0 0", "0 0 0.01"), "planar z=0");
}

#[test]
fn roundtrip_vtk_quad() {
    let imported = gmsh::parse_msh41_quads(QUAD).unwrap();
    let mesh = imported.mesh();
    let mut buf = Vec::new();
    quad::write_vtk(&mut buf, mesh, &vec![0.0; mesh.vertices().len()]).unwrap();
    let s = String::from_utf8(buf).unwrap();
    assert!(s.contains("CELL_TYPES 1"));
    assert!(s.contains("9"));
    assert!(s.contains("POINTS 4 double"));
}
