//! Experimental scalar Q1 finite elements on conforming 2D quadrilaterals.
//!
//! Reference vertex order is `(0,0), (1,0), (1,1), (0,1)`; a physical cell
//! uses the corresponding counterclockwise cyclic order. Geometry is an
//! isoparametric bilinear map; quadrature is tensor-product Gauss 3×3.
//! This module is separate from the triangle APIs and accepts only convex,
//! nondegenerate cells. Global geometric conformity is the caller's duty.
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

use crate::base::Point2;
use crate::grid::{CellId, VertexId};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::square_degree_5;
use crate::{error, Result};

/// Q1 shape values in counterclockwise reference vertex order.
pub fn shape_values(p: Point2) -> [f64; 4] {
    let (x, y) = (p.x, p.y);
    [(1.0 - x) * (1.0 - y), x * (1.0 - y), x * y, (1.0 - x) * y]
}

/// Reference-coordinate shape gradients in the same order.
pub fn reference_gradients(p: Point2) -> [[f64; 2]; 4] {
    let (x, y) = (p.x, p.y);
    [[y - 1.0, x - 1.0], [1.0 - y, -x], [y, x], [-y, 1.0 - x]]
}

/// The bilinear map from `[0,1]²` to one convex physical quadrilateral.
/// The Jacobian determinant is affine in the reference coordinates, hence
/// checking its four corner values guarantees positivity throughout the cell
/// in exact arithmetic. Its minimum must exceed 64 epsilon times the largest
/// squared boundary-edge length. No condition-number bound is promised.
#[derive(Debug, Clone, Copy)]
pub struct BilinearMap {
    points: [Point2; 4],
}

/// Mapped point, signed determinant and physical shape gradients.
#[derive(Debug, Clone, Copy)]
pub struct MappedPoint {
    pub point: Point2,
    pub determinant: f64,
    pub gradients: [[f64; 2]; 4],
}

impl BilinearMap {
    pub fn new(points: [Point2; 4]) -> Result<Self> {
        if points.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
            return Err(error("quadrilateral coordinates must be finite"));
        }
        let mut scale: f64 = 0.0;
        for i in 0..4 {
            let a = points[i];
            let b = points[(i + 1) % 4];
            scale = scale.max((b.x - a.x).powi(2) + (b.y - a.y).powi(2));
        }
        if !scale.is_finite() || scale <= 0.0 {
            return Err(error("invalid quadrilateral scale"));
        }
        let map = Self { points };
        for reference in [
            Point2 { x: 0.0, y: 0.0 },
            Point2 { x: 1.0, y: 0.0 },
            Point2 { x: 1.0, y: 1.0 },
            Point2 { x: 0.0, y: 1.0 },
        ] {
            let j = map.jacobian(reference);
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            if !det.is_finite() || det <= 64.0 * f64::EPSILON * scale {
                return Err(error(
                    "quadrilateral must be convex, counterclockwise and nondegenerate",
                ));
            }
        }
        Ok(map)
    }
    // Rows physical (x,y), columns reference (xi,eta).
    fn jacobian(&self, p: Point2) -> [[f64; 2]; 2] {
        let [a, b, c, d] = self.points;
        let dx = |a: f64, b: f64, c: f64, d: f64| {
            [
                (1.0 - p.y) * (b - a) + p.y * (c - d),
                (1.0 - p.x) * (d - a) + p.x * (c - b),
            ]
        };
        [dx(a.x, b.x, c.x, d.x), dx(a.y, b.y, c.y, d.y)]
    }
    /// Evaluate geometry and physical gradients at a point within `[0,1]²`.
    pub fn evaluate(&self, p: Point2) -> Result<MappedPoint> {
        if !p.x.is_finite() || !p.y.is_finite() || p.x < 0.0 || p.x > 1.0 || p.y < 0.0 || p.y > 1.0
        {
            return Err(error(
                "reference point must be finite and inside the unit square",
            ));
        }
        let [a, b, c, d] = self.points;
        let map = |a: f64, b: f64, c: f64, d: f64| {
            a + p.x * (b - a) + p.y * ((d - a) + p.x * ((c - d) - (b - a)))
        };
        let point = Point2::new(map(a.x, b.x, c.x, d.x), map(a.y, b.y, c.y, d.y))?;
        let j = self.jacobian(p);
        let determinant = j[0][0] * j[1][1] - j[0][1] * j[1][0];
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(error("invalid quadrilateral Jacobian"));
        }
        let refs = reference_gradients(p);
        let mut gradients = [[0.0; 2]; 4];
        for i in 0..4 {
            gradients[i] = [
                (j[1][1] * refs[i][0] - j[1][0] * refs[i][1]) / determinant,
                (-j[0][1] * refs[i][0] + j[0][0] * refs[i][1]) / determinant,
            ];
        }
        if gradients.iter().flatten().any(|x| !x.is_finite()) {
            return Err(error("nonfinite shape gradient"));
        }
        Ok(MappedPoint {
            point,
            determinant,
            gradients,
        })
    }
}

/// Immutable quadrilateral mesh. `CellId`/`VertexId` values belong to this
/// particular mesh and must not be mixed with IDs from triangle meshes.
#[derive(Debug, Clone)]
pub struct QuadMesh {
    vertices: Vec<Point2>,
    cells: Vec<[VertexId; 4]>,
    boundary_edges: Vec<[VertexId; 2]>,
    boundary_vertices: Vec<VertexId>,
}
impl QuadMesh {
    /// Validate local convexity/orientation, vertex use, duplicate cells and
    /// shared-edge topology. Intersections, hanging nodes and overlaps between
    /// different cells are not detected; a conforming mesh is required.
    pub fn new(vertices: Vec<Point2>, cells: Vec<[VertexId; 4]>) -> Result<Self> {
        if vertices.is_empty() || cells.is_empty() {
            return Err(error("mesh must contain vertices and cells"));
        }
        if vertices
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite())
        {
            return Err(error("mesh coordinates must be finite"));
        }
        let mut used = vec![false; vertices.len()];
        let mut unique = BTreeSet::new();
        let mut edges = BTreeMap::<(VertexId, VertexId), (usize, [VertexId; 2])>::new();
        for cell in &cells {
            if cell.iter().any(|v| v.0 >= vertices.len()) {
                return Err(error("cell vertex index out of bounds"));
            }
            let mut sorted = *cell;
            sorted.sort();
            if sorted.windows(2).any(|p| p[0] == p[1]) {
                return Err(error("repeated cell vertex"));
            }
            if !unique.insert(sorted) {
                return Err(error("duplicate cell"));
            }
            BilinearMap::new([
                vertices[cell[0].0],
                vertices[cell[1].0],
                vertices[cell[2].0],
                vertices[cell[3].0],
            ])?;
            for &v in cell {
                used[v.0] = true;
            }
            for i in 0..4 {
                let a = cell[i];
                let b = cell[(i + 1) % 4];
                let key = if a < b { (a, b) } else { (b, a) };
                let (count, first) = edges.entry(key).or_insert((0, [a, b]));
                *count += 1;
                if *count > 2 {
                    return Err(error("non-manifold edge"));
                }
                if *count == 2 && *first == [a, b] {
                    return Err(error("inconsistent orientation across shared edge"));
                }
            }
        }
        if used.iter().any(|u| !u) {
            return Err(error("mesh has unused vertices"));
        }
        let boundary_edges: Vec<_> = edges
            .iter()
            .filter(|(_, (count, _))| *count == 1)
            .map(|(&(a, b), _)| [a, b])
            .collect();
        let boundary_vertices = boundary_edges
            .iter()
            .flat_map(|edge| edge.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            vertices,
            cells,
            boundary_edges,
            boundary_vertices,
        })
    }
    /// Uniform square with `n × n` Q1 quadrilateral cells.
    pub fn unit_square(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(error("subdivisions must be positive"));
        }
        let side = n
            .checked_add(1)
            .ok_or_else(|| error("mesh size overflow"))?;
        let nv = side
            .checked_mul(side)
            .ok_or_else(|| error("mesh size overflow"))?;
        let nc = n
            .checked_mul(n)
            .ok_or_else(|| error("mesh size overflow"))?;
        let mut vertices = Vec::with_capacity(nv);
        for j in 0..side {
            for i in 0..side {
                vertices.push(Point2 {
                    x: i as f64 / n as f64,
                    y: j as f64 / n as f64,
                });
            }
        }
        let mut cells = Vec::with_capacity(nc);
        for j in 0..n {
            for i in 0..n {
                let a = j * side + i;
                cells.push([
                    VertexId(a),
                    VertexId(a + 1),
                    VertexId(a + side + 1),
                    VertexId(a + side),
                ]);
            }
        }
        Self::new(vertices, cells)
    }
    pub fn vertices(&self) -> &[Point2] {
        &self.vertices
    }
    pub fn cells(&self) -> &[[VertexId; 4]] {
        &self.cells
    }
    /// Boundary edges in canonical endpoint order, not an orientation.
    pub fn boundary_edges(&self) -> &[[VertexId; 2]] {
        &self.boundary_edges
    }
    pub fn boundary_vertices(&self) -> &[VertexId] {
        &self.boundary_vertices
    }
    pub fn cell_points(&self, id: CellId) -> Result<[Point2; 4]> {
        let cell = self
            .cells
            .get(id.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([
            self.vertices[cell[0].0],
            self.vertices[cell[1].0],
            self.vertices[cell[2].0],
            self.vertices[cell[3].0],
        ])
    }
}

/// Assemble unconstrained Q1 volume terms for `-Δu=f`. Source is evaluated
/// at 3×3 tensor Gauss points. For nonaffine mappings the stiffness integrand
/// is generally rational, so this quadrature is approximate even for f=0.
/// A conforming, non-overlapping mesh is a precondition.
pub fn assemble_poisson<F: Fn(Point2) -> f64>(
    mesh: &QuadMesh,
    source: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let n = mesh.vertices.len();
    let mut builder = CooBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    for (i, cell) in mesh.cells.iter().enumerate() {
        let mapping = BilinearMap::new(mesh.cell_points(CellId(i))?)?;
        let mut matrix = [[0.0; 4]; 4];
        let mut local_rhs = [0.0; 4];
        for q in square_degree_5() {
            let mapped = mapping.evaluate(q.point)?;
            let value = source(mapped.point);
            if !value.is_finite() {
                return Err(error("source function returned a nonfinite value"));
            }
            let shape = shape_values(q.point);
            let weight = mapped.determinant * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            for a in 0..4 {
                local_rhs[a] += value * shape[a] * weight;
                for b in a..4 {
                    matrix[a][b] += (mapped.gradients[a][0] * mapped.gradients[b][0]
                        + mapped.gradients[a][1] * mapped.gradients[b][1])
                        * weight;
                }
            }
        }
        for a in 0..4 {
            rhs[cell[a].0] += local_rhs[a];
            if !rhs[cell[a].0].is_finite() {
                return Err(error("nonfinite assembled right hand side"));
            }
            for b in a..4 {
                builder.add(cell[a].0, cell[b].0, matrix[a][b])?;
                if a != b {
                    builder.add(cell[b].0, cell[a].0, matrix[a][b])?;
                }
            }
        }
    }
    Ok((builder.into_csr()?, rhs))
}

/// Compute integrated L2 and H1-seminorm errors with the 3×3 rule.
pub fn errors<F, G>(
    mesh: &QuadMesh,
    solution: &[f64],
    exact: F,
    exact_gradient: G,
) -> Result<(f64, f64)>
where
    F: Fn(Point2) -> f64,
    G: Fn(Point2) -> [f64; 2],
{
    if solution.len() != mesh.vertices.len() || solution.iter().any(|v| !v.is_finite()) {
        return Err(error("invalid solution vector"));
    }
    let mut l2 = 0.0;
    let mut h1 = 0.0;
    for (i, cell) in mesh.cells.iter().enumerate() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i))?)?;
        for q in square_degree_5() {
            let mapped = map.evaluate(q.point)?;
            let shape = shape_values(q.point);
            let value = exact(mapped.point);
            let gradient = exact_gradient(mapped.point);
            if !value.is_finite() || gradient.iter().any(|v| !v.is_finite()) {
                return Err(error("exact solution returned a nonfinite value"));
            }
            let mut uh = 0.0;
            let mut grad = [0.0; 2];
            for j in 0..4 {
                let u = solution[cell[j].0];
                uh += u * shape[j];
                grad[0] += u * mapped.gradients[j][0];
                grad[1] += u * mapped.gradients[j][1];
            }
            let du = value - uh;
            let dx = gradient[0] - grad[0];
            let dy = gradient[1] - grad[1];
            let weight = mapped.determinant * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            l2 += du * du * weight;
            h1 += (dx * dx + dy * dy) * weight;
        }
    }
    if !l2.is_finite() || !h1.is_finite() {
        return Err(error("nonfinite integrated error"));
    }
    Ok((l2.sqrt(), h1.sqrt()))
}

/// Write VTK Legacy ASCII quadrilaterals (VTK cell type 9) and a scalar
/// per vertex. This is separate from the triangle exporter.
pub fn write_vtk<W: Write>(mut out: W, mesh: &QuadMesh, values: &[f64]) -> io::Result<()> {
    if values.len() != mesh.vertices.len() || values.iter().any(|v| !v.is_finite()) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid nodal field",
        ));
    }
    let entries = mesh
        .cells
        .len()
        .checked_mul(5)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "mesh size overflow"))?;
    writeln!(
        out,
        "# vtk DataFile Version 3.0\ndealii-rust Q1 solution\nASCII\nDATASET UNSTRUCTURED_GRID"
    )?;
    writeln!(out, "POINTS {} double", mesh.vertices.len())?;
    for p in &mesh.vertices {
        writeln!(out, "{} {} 0", p.x, p.y)?;
    }
    writeln!(out, "CELLS {} {}", mesh.cells.len(), entries)?;
    for c in &mesh.cells {
        writeln!(out, "4 {} {} {} {}", c[0].0, c[1].0, c[2].0, c[3].0)?;
    }
    writeln!(out, "CELL_TYPES {}", mesh.cells.len())?;
    for _ in &mesh.cells {
        writeln!(out, "9")?;
    }
    writeln!(
        out,
        "POINT_DATA {}\nSCALARS solution double 1\nLOOKUP_TABLE default",
        values.len()
    )?;
    for v in values {
        writeln!(out, "{}", v)?;
    }
    Ok(())
}
