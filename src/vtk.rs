//! VTK Legacy ASCII export of a scalar nodal field on triangles.
use crate::grid::Mesh2;
use std::io::{self, Write};

/// Write a VTK unstructured grid. The field has one finite value per vertex.
pub fn write_scalar<W: Write>(mut out: W, mesh: &Mesh2, values: &[f64]) -> io::Result<()> {
    if values.len() != mesh.vertices().len() || values.iter().any(|v| !v.is_finite()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid nodal field",
        ));
    }
    writeln!(out, "# vtk DataFile Version 3.0\ndealii-rust Poisson solution\nASCII\nDATASET UNSTRUCTURED_GRID")?;
    writeln!(out, "POINTS {} double", mesh.vertices().len())?;
    for p in mesh.vertices() {
        writeln!(out, "{} {} 0", p.x, p.y)?;
    }
    writeln!(
        out,
        "CELLS {} {}",
        mesh.cells().len(),
        mesh.cells().len() * 4
    )?;
    for c in mesh.cells() {
        writeln!(out, "3 {} {} {}", c[0].0, c[1].0, c[2].0)?;
    }
    writeln!(out, "CELL_TYPES {}", mesh.cells().len())?;
    for _ in mesh.cells() {
        writeln!(out, "5")?;
    }
    writeln!(
        out,
        "POINT_DATA {}\nSCALARS solution double 1\nLOOKUP_TABLE default",
        values.len()
    )?;
    for value in values {
        writeln!(out, "{}", value)?;
    }
    Ok(())
}
