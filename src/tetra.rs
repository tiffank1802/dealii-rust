//! Scalar P1 Poisson on conforming, positively oriented tetrahedral meshes.
//!
//! Reference vertices are `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`.
//! A positive Jacobian is required. Global geometric conformity is assumed.
use crate::base::Point3;
use crate::grid::{CellId, VertexId};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::{tetra_degree_5, triangle_degree_5};
use crate::{error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

fn sub(a: Point3, b: Point3) -> [f64; 3] {
    [a.x - b.x, a.y - b.y, a.z - b.z]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// P1 reference shape values, in reference vertex order.
pub fn shape_values(p: Point3) -> [f64; 4] {
    [1.0 - p.x - p.y - p.z, p.x, p.y, p.z]
}

/// Affine tetrahedron map, with constant physical P1 gradients.
#[derive(Debug, Clone, Copy)]
pub struct AffineTet {
    points: [Point3; 4],
    determinant: f64,
    gradients: [[f64; 3]; 4],
}
impl AffineTet {
    /// The signed determinant must exceed 64 eps times the largest edge
    /// length cubed. Inverted, degenerate and nonfinite geometry is rejected.
    pub fn new(points: [Point3; 4]) -> Result<Self> {
        if points
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
        {
            return Err(error("tetrahedron coordinates must be finite"));
        }
        let [a, b, c, d] = points;
        let e1 = sub(b, a);
        let e2 = sub(c, a);
        let e3 = sub(d, a);
        let determinant = dot(e1, cross(e2, e3));
        let mut scale: f64 = 0.0;
        for i in 0..4 {
            for j in i + 1..4 {
                let e = sub(points[i], points[j]);
                scale = scale.max(dot(e, e));
            }
        }
        let scale3 = scale * scale.sqrt();
        if !determinant.is_finite()
            || !scale3.is_finite()
            || scale <= 0.0
            || determinant <= 64.0 * f64::EPSILON * scale3
        {
            return Err(error(
                "tetrahedron must be positively oriented and nondegenerate",
            ));
        }
        let mut gradients = [[0.0; 3]; 4];
        for (i, vector) in [cross(e2, e3), cross(e3, e1), cross(e1, e2)]
            .iter()
            .enumerate()
        {
            for (k, value) in vector.iter().enumerate() {
                gradients[i + 1][k] = value / determinant;
            }
        }
        for k in 0..3 {
            gradients[0][k] = -gradients[1][k] - gradients[2][k] - gradients[3][k];
        }
        if gradients.iter().flatten().any(|x| !x.is_finite()) {
            return Err(error("nonfinite tetrahedral gradient"));
        }
        Ok(Self {
            points,
            determinant,
            gradients,
        })
    }
    pub fn determinant(&self) -> f64 {
        self.determinant
    }
    pub fn gradients(&self) -> [[f64; 3]; 4] {
        self.gradients
    }
    /// Map a finite reference point by the affine map, including extrapolation.
    pub fn map(&self, p: Point3) -> Result<Point3> {
        let [a, b, c, d] = self.points;
        let e1 = sub(b, a);
        let e2 = sub(c, a);
        let e3 = sub(d, a);
        Point3::new(
            a.x + p.x * e1[0] + p.y * e2[0] + p.z * e3[0],
            a.y + p.x * e1[1] + p.y * e2[1] + p.z * e3[1],
            a.z + p.x * e1[2] + p.y * e2[2] + p.z * e3[2],
        )
    }
}

/// Index into a tetrahedral mesh's boundary-face list.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BoundaryFaceId(pub usize);
/// Outward-oriented triangle and physical geometry.
#[derive(Debug, Copy, Clone)]
pub struct BoundaryFaceGeometry {
    pub vertices: [VertexId; 3],
    pub points: [Point3; 3],
    pub outward_normal: [f64; 3],
    pub area: f64,
}

/// Immutable mesh with tetrahedral cells. IDs belong only to this mesh.
#[derive(Debug, Clone)]
pub struct TetMesh {
    vertices: Vec<Point3>,
    cells: Vec<[VertexId; 4]>,
    boundary_faces: Vec<[VertexId; 3]>,
    boundary_vertices: Vec<VertexId>,
}
fn odd(face: [VertexId; 3]) -> bool {
    let mut inversions = 0;
    for i in 0..3 {
        for j in i + 1..3 {
            if face[i] > face[j] {
                inversions += 1;
            }
        }
    }
    inversions % 2 != 0
}
impl TetMesh {
    /// Validate indices, orientations, degeneracy, duplicates, unused vertices
    /// and face incidence/orientation. Intersections, overlaps and hanging
    /// nodes are not detected; the caller must supply a conforming mesh.
    pub fn new(vertices: Vec<Point3>, cells: Vec<[VertexId; 4]>) -> Result<Self> {
        if vertices.is_empty() || cells.is_empty() {
            return Err(error("mesh must contain vertices and cells"));
        }
        if vertices
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
        {
            return Err(error("mesh coordinates must be finite"));
        }
        let mut used = vec![false; vertices.len()];
        let mut unique = BTreeSet::new();
        let mut faces = BTreeMap::<[VertexId; 3], (usize, [VertexId; 3])>::new();
        for cell in &cells {
            if cell.iter().any(|v| v.0 >= vertices.len()) {
                return Err(error("cell vertex index out of bounds"));
            }
            let mut sorted = *cell;
            sorted.sort();
            if sorted.windows(2).any(|v| v[0] == v[1]) {
                return Err(error("repeated cell vertex"));
            }
            if !unique.insert(sorted) {
                return Err(error("duplicate cell"));
            }
            AffineTet::new([
                vertices[cell[0].0],
                vertices[cell[1].0],
                vertices[cell[2].0],
                vertices[cell[3].0],
            ])?;
            for v in cell {
                used[v.0] = true;
            }
            let [a, b, c, d] = *cell;
            for oriented in [[b, c, d], [a, d, c], [a, b, d], [a, c, b]] {
                let mut key = oriented;
                key.sort();
                let (count, first) = faces.entry(key).or_insert((0, oriented));
                *count += 1;
                if *count > 2 {
                    return Err(error("non-manifold face"));
                }
                if *count == 2 && odd(*first) == odd(oriented) {
                    return Err(error("inconsistent orientation across shared face"));
                }
            }
        }
        if used.iter().any(|v| !v) {
            return Err(error("mesh has unused vertices"));
        }
        let boundary_faces: Vec<_> = faces
            .values()
            .filter(|(n, _)| *n == 1)
            .map(|(_, f)| *f)
            .collect();
        let boundary_vertices = boundary_faces
            .iter()
            .flat_map(|f| f.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            vertices,
            cells,
            boundary_faces,
            boundary_vertices,
        })
    }
    /// Split each structured cube into six conforming tetrahedra around its
    /// body diagonal. Each tetrahedron is reoriented to positive determinant.
    pub fn unit_cube(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(error("subdivisions must be positive"));
        }
        let side = n
            .checked_add(1)
            .ok_or_else(|| error("mesh size overflow"))?;
        let nv = side
            .checked_mul(side)
            .and_then(|x| x.checked_mul(side))
            .ok_or_else(|| error("mesh size overflow"))?;
        let nc = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(n))
            .and_then(|x| x.checked_mul(6))
            .ok_or_else(|| error("mesh size overflow"))?;
        let mut vertices = Vec::with_capacity(nv);
        for k in 0..side {
            for j in 0..side {
                for i in 0..side {
                    vertices.push(Point3 {
                        x: i as f64 / n as f64,
                        y: j as f64 / n as f64,
                        z: k as f64 / n as f64,
                    });
                }
            }
        }
        let mut cells = Vec::with_capacity(nc);
        let permutations = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    for order in permutations {
                        let mut position = [i, j, k];
                        let index = |p: [usize; 3]| VertexId((p[2] * side + p[1]) * side + p[0]);
                        let a = index(position);
                        position[order[0]] += 1;
                        let b = index(position);
                        position[order[1]] += 1;
                        let c = index(position);
                        position[order[2]] += 1;
                        let d = index(position);
                        let mut inversions = 0;
                        for x in 0..3 {
                            for y in x + 1..3 {
                                if order[x] > order[y] {
                                    inversions += 1;
                                }
                            }
                        }
                        if inversions % 2 == 0 {
                            cells.push([a, b, c, d]);
                        } else {
                            cells.push([a, c, b, d]);
                        }
                    }
                }
            }
        }
        Self::new(vertices, cells)
    }
    pub fn vertices(&self) -> &[Point3] {
        &self.vertices
    }
    pub fn cells(&self) -> &[[VertexId; 4]] {
        &self.cells
    }
    /// Outward-oriented boundary-face vertex IDs; position is BoundaryFaceId.
    pub fn boundary_faces(&self) -> &[[VertexId; 3]] {
        &self.boundary_faces
    }
    pub fn boundary_vertices(&self) -> &[VertexId] {
        &self.boundary_vertices
    }
    pub fn cell_points(&self, id: CellId) -> Result<[Point3; 4]> {
        let c = self
            .cells
            .get(id.0)
            .ok_or_else(|| error("cell index out of bounds"))?;
        Ok([
            self.vertices[c[0].0],
            self.vertices[c[1].0],
            self.vertices[c[2].0],
            self.vertices[c[3].0],
        ])
    }
    pub fn boundary_face_geometry(&self, id: BoundaryFaceId) -> Result<BoundaryFaceGeometry> {
        let v = *self
            .boundary_faces
            .get(id.0)
            .ok_or_else(|| error("boundary face index out of bounds"))?;
        let points = [
            self.vertices[v[0].0],
            self.vertices[v[1].0],
            self.vertices[v[2].0],
        ];
        let normal = cross(sub(points[1], points[0]), sub(points[2], points[0]));
        let length = normal[0].hypot(normal[1]).hypot(normal[2]);
        let area = length * 0.5;
        if !length.is_finite() || !area.is_finite() || area <= 0.0 {
            return Err(error("invalid boundary face geometry"));
        }
        Ok(BoundaryFaceGeometry {
            vertices: v,
            points,
            outward_normal: [normal[0] / length, normal[1] / length, normal[2] / length],
            area,
        })
    }
}

/// Assemble the P1 Laplace matrix and volume load for `-Δu=f`. Uses a
/// degree-five tetrahedron quadrature for the source. No boundary condition
/// is imposed. The mesh must be geometrically conforming and non-overlapping.
pub fn assemble_poisson<F: Fn(Point3) -> f64>(
    mesh: &TetMesh,
    source: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let n = mesh.vertices.len();
    let mut builder = CooBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    let quadrature = tetra_degree_5();
    for (i, cell) in mesh.cells.iter().enumerate() {
        let map = AffineTet::new(mesh.cell_points(CellId(i))?)?;
        let grads = map.gradients();
        let volume = map.determinant() / 6.0;
        if !volume.is_finite() || volume <= 0.0 {
            return Err(error("invalid tetrahedral volume"));
        }
        for a in 0..4 {
            for b in 0..4 {
                builder.add(cell[a].0, cell[b].0, volume * dot(grads[a], grads[b]))?;
            }
        }
        for q in &quadrature {
            let point = map.map(q.point)?;
            let value = source(point);
            if !value.is_finite() {
                return Err(error("source function returned a nonfinite value"));
            }
            let shapes = shape_values(q.point);
            let weight = map.determinant() * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            for a in 0..4 {
                rhs[cell[a].0] += value * shapes[a] * weight;
            }
        }
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        return Err(error("nonfinite assembled right hand side"));
    }
    Ok((builder.into_csr()?, rhs))
}

/// L2 and H1-seminorm errors using degree-five tetrahedron quadrature.
pub fn errors<F, G>(mesh: &TetMesh, solution: &[f64], exact: F, gradient: G) -> Result<(f64, f64)>
where
    F: Fn(Point3) -> f64,
    G: Fn(Point3) -> [f64; 3],
{
    if solution.len() != mesh.vertices.len() || solution.iter().any(|x| !x.is_finite()) {
        return Err(error("invalid solution vector"));
    }
    let mut l2 = 0.0;
    let mut h1 = 0.0;
    let quadrature = tetra_degree_5();
    for (i, cell) in mesh.cells.iter().enumerate() {
        let map = AffineTet::new(mesh.cell_points(CellId(i))?)?;
        let grads = map.gradients();
        let mut uh_gradient = [0.0; 3];
        for j in 0..4 {
            for (k, value) in uh_gradient.iter_mut().enumerate() {
                *value += solution[cell[j].0] * grads[j][k];
            }
        }
        for q in &quadrature {
            let point = map.map(q.point)?;
            let value = exact(point);
            let target_gradient = gradient(point);
            if !value.is_finite() || target_gradient.iter().any(|x| !x.is_finite()) {
                return Err(error("exact solution returned a nonfinite value"));
            }
            let shape = shape_values(q.point);
            let mut uh = 0.0;
            for j in 0..4 {
                uh += solution[cell[j].0] * shape[j];
            }
            let diff = value - uh;
            let mut grad_diff = 0.0;
            for k in 0..3 {
                grad_diff += (target_gradient[k] - uh_gradient[k]).powi(2);
            }
            let weight = map.determinant() * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            l2 += diff * diff * weight;
            h1 += grad_diff * weight;
        }
    }
    if !l2.is_finite() || !h1.is_finite() {
        return Err(error("nonfinite integrated error"));
    }
    Ok((l2.sqrt(), h1.sqrt()))
}

/// Write a scalar P1 field as VTK Legacy ASCII tetrahedra (cell type 10).
pub fn write_vtk<W: Write>(mut out: W, mesh: &TetMesh, values: &[f64]) -> io::Result<()> {
    if values.len() != mesh.vertices.len() || values.iter().any(|x| !x.is_finite()) {
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
        "# vtk DataFile Version 3.0\ndealii-rust 3D solution\nASCII\nDATASET UNSTRUCTURED_GRID"
    )?;
    writeln!(out, "POINTS {} double", mesh.vertices.len())?;
    for p in &mesh.vertices {
        writeln!(out, "{} {} {}", p.x, p.y, p.z)?;
    }
    writeln!(out, "CELLS {} {}", mesh.cells.len(), entries)?;
    for c in &mesh.cells {
        writeln!(out, "4 {} {} {} {}", c[0].0, c[1].0, c[2].0, c[3].0)?;
    }
    writeln!(out, "CELL_TYPES {}", mesh.cells.len())?;
    for _ in &mesh.cells {
        writeln!(out, "10")?;
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

/// Interpolate Dirichlet data on selected triangular boundary faces. Returns
/// sorted `(vertex index, value)` constraints for `CsrMatrix::with_dirichlet`;
/// vertices shared by faces are sampled once. Face IDs must belong to this
/// mesh, be in bounds and be unique. Nonfinite function values are rejected.
pub fn interpolate_boundary<F: Fn(Point3) -> f64>(
    mesh: &TetMesh,
    faces: &[BoundaryFaceId],
    value: F,
) -> Result<Vec<(usize, f64)>> {
    let mut seen = BTreeSet::new();
    let mut vertices = BTreeSet::new();
    for &id in faces {
        if !seen.insert(id) {
            return Err(error("duplicate boundary face"));
        }
        let face = mesh
            .boundary_faces
            .get(id.0)
            .ok_or_else(|| error("boundary face index out of bounds"))?;
        vertices.extend(face.iter().copied());
    }
    vertices
        .into_iter()
        .map(|v| {
            let result = value(mesh.vertices[v.0]);
            if !result.is_finite() {
                Err(error("boundary function returned a nonfinite value"))
            } else {
                Ok((v.0, result))
            }
        })
        .collect()
}

/// Assemble `∫ g φ_i dS` over selected tetrahedral boundary faces.
/// For `-Δu=f`, the callback receives the physical point and the unit
/// **outward** normal and must return `g = ∂u/∂n`. Its contribution is added
/// to the volume right-hand side before Dirichlet elimination. Face IDs must
/// belong to this mesh, be in bounds and be unique. Nonfinite values and
/// nonrepresentable results are rejected.
///
/// The P1 face basis and flux use a degree-five triangle rule. The integral is
/// exact for polynomial fluxes of degree ≤ 4, apart from floating-point error.
/// This is only a load assembly routine: pure Neumann problems remain singular
/// and require compatibility and normalization. Use [`crate::neumann3d::solve`]
/// for full-boundary Neumann Poisson, or CG only once the system is positive definite.
pub fn assemble_neumann<F>(mesh: &TetMesh, faces: &[BoundaryFaceId], flux: F) -> Result<Vec<f64>>
where
    F: Fn(Point3, [f64; 3]) -> f64,
{
    let mut seen = BTreeSet::new();
    for &id in faces {
        if !seen.insert(id) {
            return Err(error("duplicate Neumann boundary face"));
        }
        mesh.boundary_face_geometry(id)?;
    }
    let quadrature = triangle_degree_5();
    let mut load = vec![0.0; mesh.vertices.len()];
    for &id in faces {
        let face = mesh.boundary_face_geometry(id)?;
        let [a, b, c] = face.points;
        let e1 = sub(b, a);
        let e2 = sub(c, a);
        let jacobian = 2.0 * face.area;
        for q in &quadrature {
            let p = Point3::new(
                a.x + q.point.x * e1[0] + q.point.y * e2[0],
                a.y + q.point.x * e1[1] + q.point.y * e2[1],
                a.z + q.point.x * e1[2] + q.point.y * e2[2],
            )?;
            let g = flux(p, face.outward_normal);
            if !g.is_finite() {
                return Err(error("Neumann flux returned a nonfinite value"));
            }
            let weight = jacobian * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid boundary quadrature weight"));
            }
            let shape = crate::fe::shape_values(q.point);
            for i in 0..3 {
                load[face.vertices[i].0] += g * shape[i] * weight;
            }
        }
    }
    if load.iter().any(|x| !x.is_finite()) {
        return Err(error("nonfinite Neumann load"));
    }
    Ok(load)
}
