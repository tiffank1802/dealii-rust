use dealii_rust::{grid::Mesh2, vtk};
#[test]
fn writes_point_and_cell_data_in_mesh_order() {
    let mesh = Mesh2::unit_square(1).unwrap();
    let values = [1., 2., 3., 4.];
    let cells = [0.125, 0.25];
    let mut out = Vec::new();
    vtk::write_scalar_with_cell_scalar(&mut out, &mesh, &values, "indicator_1", &cells).unwrap();
    let text = String::from_utf8(out).unwrap();
    assert!(text
        .contains("POINT_DATA 4\nSCALARS solution double 1\nLOOKUP_TABLE default\n1\n2\n3\n4\n"));
    assert!(text.ends_with(
        "CELL_DATA 2\nSCALARS indicator_1 double 1\nLOOKUP_TABLE default\n0.125\n0.25\n"
    ));
}
#[test]
fn rejects_invalid_fields_and_names_before_writing() {
    let mesh = Mesh2::unit_square(1).unwrap();
    for name in ["", "1a", "invalid name", "invalid\nname", "é"] {
        let mut out = Vec::new();
        assert!(
            vtk::write_scalar_with_cell_scalar(&mut out, &mesh, &[1.; 4], name, &[1.; 2]).is_err()
        );
        assert!(out.is_empty());
    }
    for (points, cells) in [
        (vec![1.; 3], vec![1.; 2]),
        (vec![1.; 4], vec![1.]),
        (vec![f64::NAN; 4], vec![1.; 2]),
        (vec![1.; 4], vec![f64::INFINITY; 2]),
    ] {
        let mut out = Vec::new();
        assert!(
            vtk::write_scalar_with_cell_scalar(&mut out, &mesh, &points, "tag", &cells).is_err()
        );
        assert!(out.is_empty());
    }
}
