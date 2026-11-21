//! Scalar P1 Poisson on conforming, positively oriented tetrahedral meshes.
//!
//! Reference vertices are `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`.
//! A positive Jacobian is required. Global geometric conformity is assumed.
use crate::base::{Point3, UpdateFlags};
use crate::grid::{CellId, VertexId};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::{tetra_degree_5, triangle_degree_5};
pub use crate::tetra_refinement::UniformRefinement;
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

/// Number of shape functions of the P1 tetrahedron.
const P1_TETRA_N_SHAPES: usize = 4;

/// Constant reference gradients of the P1 basis functions
/// `1 - x - y - z`, `x`, `y`, `z` on the reference tetrahedron.
const P1_TETRA_REFERENCE_GRADIENTS: [[f64; 3]; P1_TETRA_N_SHAPES] = [
    [-1.0, -1.0, -1.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// FEValues-style evaluation of the P1 tetrahedron on an affine cell.
///
/// The evaluation is built once for a cell, a quadrature rule and a set of
/// [`UpdateFlags`], and then exposes the requested quantities. Quantities whose
/// flag is false are not computed at all: the matching accessor returns `None`.
///
/// Arrays are stored shape-major, so `phi[i][q]` lives at index
/// `i * n_quadrature_points() + q`; the per-entry accessors ([`value`],
/// [`grad_ref_at`], [`grad_phys_at`], [`jxw_at`]) hide that layout.
///
/// Physical quadrature points are always computed and are not flag-gated.
///
/// [`value`]: FeValuesP1Tetrahedron::value
/// [`grad_ref_at`]: FeValuesP1Tetrahedron::grad_ref_at
/// [`grad_phys_at`]: FeValuesP1Tetrahedron::grad_phys_at
/// [`jxw_at`]: FeValuesP1Tetrahedron::jxw_at
#[derive(Debug, Clone)]
pub struct FeValuesP1Tetrahedron {
    n_quadrature_points: usize,
    flags: UpdateFlags,
    quadrature_points: Vec<Point3>,
    values: Option<Vec<f64>>,
    grad_ref: Option<Vec<[f64; 3]>>,
    grad_phys: Option<Vec<[f64; 3]>>,
    jxw: Option<Vec<f64>>,
}

impl FeValuesP1Tetrahedron {
    /// Evaluate the P1 tetrahedron on `cell` at the given reference points.
    ///
    /// `reference_points` and `weights` must have the same non-zero length. All
    /// reference coordinates and all weights must be finite, and every weight
    /// must be strictly positive. `cell` must be a positively oriented,
    /// non-degenerate affine tetrahedron, as enforced by
    /// [`AffineTet::new`].
    ///
    /// Input validation always runs, including when `flags` is empty.
    pub fn new(
        cell: &AffineTet,
        reference_points: &[Point3],
        weights: &[f64],
        flags: UpdateFlags,
    ) -> Result<Self> {
        if reference_points.len() != weights.len() {
            return Err(error("quadrature point and weight counts must match"));
        }
        if reference_points.is_empty() {
            return Err(error("quadrature rule must contain at least one point"));
        }
        for point in reference_points {
            if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                return Err(error("quadrature point coordinates must be finite"));
            }
        }
        for weight in weights {
            if !weight.is_finite() || *weight <= 0.0 {
                return Err(error("quadrature weights must be finite and positive"));
            }
        }
        let determinant = cell.determinant();
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(error(
                "cell must be a positively oriented nondegenerate tetrahedron",
            ));
        }

        let n = reference_points.len();
        let quadrature_points = reference_points
            .iter()
            .map(|point| cell.map(*point))
            .collect::<Result<Vec<_>>>()?;

        let values = if flags.values {
            let mut data = vec![0.0; P1_TETRA_N_SHAPES * n];
            for (q, point) in reference_points.iter().enumerate() {
                for (i, value) in shape_values(*point).iter().enumerate() {
                    data[i * n + q] = *value;
                }
            }
            Some(data)
        } else {
            None
        };

        let grad_ref = if flags.gradients {
            let mut data = vec![[0.0; 3]; P1_TETRA_N_SHAPES * n];
            for (i, gradient) in P1_TETRA_REFERENCE_GRADIENTS.iter().enumerate() {
                data[i * n..(i + 1) * n].fill(*gradient);
            }
            Some(data)
        } else {
            None
        };

        let grad_phys = if flags.gradients {
            let physical = cell.gradients();
            let mut data = vec![[0.0; 3]; P1_TETRA_N_SHAPES * n];
            for (i, gradient) in physical.iter().enumerate() {
                data[i * n..(i + 1) * n].fill(*gradient);
            }
            Some(data)
        } else {
            None
        };

        let jxw = if flags.jxw {
            let abs_determinant = determinant.abs();
            Some(
                weights
                    .iter()
                    .map(|weight| abs_determinant * weight)
                    .collect(),
            )
        } else {
            None
        };

        Ok(Self {
            n_quadrature_points: n,
            flags,
            quadrature_points,
            values,
            grad_ref,
            grad_phys,
            jxw,
        })
    }

    /// Number of quadrature points in this evaluation.
    pub fn n_quadrature_points(&self) -> usize {
        self.n_quadrature_points
    }

    /// The flags this evaluation was built with.
    pub fn update_flags(&self) -> UpdateFlags {
        self.flags
    }

    /// Physical quadrature points `x_q = map(p_q)`, always computed.
    pub fn quadrature_points(&self) -> &[Point3] {
        &self.quadrature_points
    }

    /// Shape values `phi[i][q]`, shape-major; `None` unless `values` was requested.
    pub fn values(&self) -> Option<&[f64]> {
        self.values.as_deref()
    }

    /// Shape value `phi[shape][q]`; `None` if values were not requested or out of range.
    pub fn value(&self, shape: usize, q: usize) -> Option<f64> {
        let index = self.flat_index(shape, q)?;
        self.values.as_ref()?.get(index).copied()
    }

    /// Reference gradients `dphi_i/d(ref coords)`, shape-major;
    /// `None` unless `gradients` was requested.
    pub fn grad_ref(&self) -> Option<&[[f64; 3]]> {
        self.grad_ref.as_deref()
    }

    /// Reference gradient of `shape` at point `q`; `None` if gradients were not
    /// requested or out of range.
    pub fn grad_ref_at(&self, shape: usize, q: usize) -> Option<[f64; 3]> {
        let index = self.flat_index(shape, q)?;
        self.grad_ref.as_ref()?.get(index).copied()
    }

    /// Physical gradients `J^{-T} grad_ref`, shape-major;
    /// `None` unless `gradients` was requested.
    pub fn grad_phys(&self) -> Option<&[[f64; 3]]> {
        self.grad_phys.as_deref()
    }

    /// Physical gradient of `shape` at point `q`; `None` if gradients were not
    /// requested or out of range.
    pub fn grad_phys_at(&self, shape: usize, q: usize) -> Option<[f64; 3]> {
        let index = self.flat_index(shape, q)?;
        self.grad_phys.as_ref()?.get(index).copied()
    }

    /// Jacobian-times-weight factors `|det J| * w_q`; `None` unless `jxw` was requested.
    pub fn jxw(&self) -> Option<&[f64]> {
        self.jxw.as_deref()
    }

    /// `JxW` at point `q`; `None` if `jxw` was not requested or out of range.
    pub fn jxw_at(&self, q: usize) -> Option<f64> {
        self.jxw.as_ref()?.get(q).copied()
    }

    /// Flattened shape-major index of `(shape, q)`, or `None` when out of range.
    fn flat_index(&self, shape: usize, q: usize) -> Option<usize> {
        if shape >= P1_TETRA_N_SHAPES || q >= self.n_quadrature_points {
            return None;
        }
        Some(shape * self.n_quadrature_points + q)
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
/// **outward** normal and must return `g = ∂u/∂n` (or the conormal flux
/// `g = n · A∇u` for anisotropic diffusion). Its contribution is added
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::Point3;
    use crate::quadrature::tetra_degree_5;

    fn p(x: f64, y: f64, z: f64) -> Point3 {
        Point3::new(x, y, z).unwrap()
    }

    fn near(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{} != {} (tolerance {})", a, b, tol);
    }

    fn unit_cell() -> AffineTet {
        AffineTet::new([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ])
        .unwrap()
    }

    fn scaled_cell() -> AffineTet {
        AffineTet::new([
            p(1.0, 2.0, 3.0),
            p(4.0, 2.0, 3.0),
            p(1.0, 5.0, 3.0),
            p(1.0, 2.0, 7.0),
        ])
        .unwrap()
    }

    fn degree_five_rule() -> (Vec<Point3>, Vec<f64>) {
        let rule = tetra_degree_5();
        (
            rule.iter().map(|q| q.point).collect(),
            rule.iter().map(|q| q.weight).collect(),
        )
    }

    fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn sub3(a: Point3, b: Point3) -> [f64; 3] {
        [a.x - b.x, a.y - b.y, a.z - b.z]
    }

    /// J^{-T} applied to a reference gradient, computed independently from edges.
    fn apply_jtinv(g: [f64; 3], e1: [f64; 3], e2: [f64; 3], e3: [f64; 3], det: f64) -> [f64; 3] {
        let c0 = cross3(e2, e3);
        let c1 = cross3(e3, e1);
        let c2 = cross3(e1, e2);
        [
            (c0[0] * g[0] + c1[0] * g[1] + c2[0] * g[2]) / det,
            (c0[1] * g[0] + c1[1] * g[1] + c2[1] * g[2]) / det,
            (c0[2] * g[0] + c1[2] * g[1] + c2[2] * g[2]) / det,
        ]
    }

    /// ADR-0003 invariant 1.
    #[test]
    fn values_satisfy_partition_of_unity() {
        let (points, weights) = degree_five_rule();
        for cell in [unit_cell(), scaled_cell()].iter() {
            let fe = FeValuesP1Tetrahedron::new(cell, &points, &weights, UpdateFlags::ALL).unwrap();
            assert_eq!(fe.n_quadrature_points(), points.len());
            assert_eq!(fe.values().unwrap().len(), 4 * points.len());
            for q in 0..fe.n_quadrature_points() {
                let sum: f64 = (0..4).map(|i| fe.value(i, q).unwrap()).sum();
                near(sum, 1.0, 1e-12);
            }
        }
    }

    /// ADR-0003 invariant 2.
    #[test]
    fn values_are_nodal_at_reference_vertices() {
        let vertices = [
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ];
        let weights = [1.0; 4];
        let flags = UpdateFlags {
            values: true,
            gradients: false,
            jxw: false,
        };
        let fe = FeValuesP1Tetrahedron::new(&unit_cell(), &vertices, &weights, flags).unwrap();
        for (q, vertex) in vertices.iter().enumerate() {
            let expected = shape_values(*vertex);
            for i in 0..4 {
                let delta = if i == q { 1.0 } else { 0.0 };
                near(fe.value(i, q).unwrap(), delta, 1e-14);
                near(fe.value(i, q).unwrap(), expected[i], 1e-14);
            }
        }
    }

    /// ADR-0003 invariant 3.
    #[test]
    fn reference_gradients_sum_to_zero() {
        let (points, weights) = degree_five_rule();
        let flags = UpdateFlags {
            values: false,
            gradients: true,
            jxw: false,
        };
        let fe = FeValuesP1Tetrahedron::new(&scaled_cell(), &points, &weights, flags).unwrap();
        assert!(fe.values().is_none());
        for q in 0..fe.n_quadrature_points() {
            let mut sum = [0.0, 0.0, 0.0];
            for i in 0..4 {
                let gradient = fe.grad_ref_at(i, q).unwrap();
                let expected = P1_TETRA_REFERENCE_GRADIENTS[i];
                for k in 0..3 {
                    near(gradient[k], expected[k], 1e-15);
                    sum[k] += gradient[k];
                }
            }
            for k in 0..3 {
                near(sum[k], 0.0, 1e-15);
            }
        }
    }

    /// ADR-0003 invariants 4, 5 and 6.
    #[test]
    fn jxw_and_physical_gradients_match_the_affine_mapping() {
        let (points, weights) = degree_five_rule();
        let cell = scaled_cell();
        let fe = FeValuesP1Tetrahedron::new(&cell, &points, &weights, UpdateFlags::ALL).unwrap();

        let determinant = cell.determinant();
        assert!(determinant > 0.0);

        // reconstruct edges from scaled_cell vertices
        let a = p(1.0, 2.0, 3.0);
        let b = p(4.0, 2.0, 3.0);
        let c = p(1.0, 5.0, 3.0);
        let d = p(1.0, 2.0, 7.0);
        let e1 = sub3(b, a);
        let e2 = sub3(c, a);
        let e3 = sub3(d, a);

        for (q, weight) in weights.iter().enumerate() {
            let jxw = fe.jxw_at(q).unwrap();
            assert!(jxw > 0.0, "JxW must be positive, got {}", jxw);
            near(jxw, determinant.abs() * weight, 1e-12);
            near(fe.jxw().unwrap()[q], jxw, 1e-15);

            for i in 0..4 {
                let reference = fe.grad_ref_at(i, q).unwrap();
                let expected = apply_jtinv(reference, e1, e2, e3, determinant);
                let physical = fe.grad_phys_at(i, q).unwrap();
                for k in 0..3 {
                    near(physical[k], expected[k], 1e-12);
                }
            }
        }

        for (q, point) in points.iter().enumerate() {
            assert_eq!(fe.quadrature_points()[q], cell.map(*point).unwrap());
        }
    }

    /// ADR-0003 invariant 7.
    #[test]
    fn flags_gate_computation() {
        let (points, weights) = degree_five_rule();
        let cell = unit_cell();

        let none =
            FeValuesP1Tetrahedron::new(&cell, &points, &weights, UpdateFlags::EMPTY).unwrap();
        assert_eq!(none.update_flags(), UpdateFlags::EMPTY);
        assert!(none.values().is_none());
        assert!(none.grad_ref().is_none());
        assert!(none.grad_phys().is_none());
        assert!(none.jxw().is_none());
        assert!(none.value(0, 0).is_none());
        assert!(none.grad_ref_at(0, 0).is_none());
        assert!(none.grad_phys_at(0, 0).is_none());
        assert!(none.jxw_at(0).is_none());
        assert_eq!(none.quadrature_points().len(), points.len());

        let only_values = FeValuesP1Tetrahedron::new(
            &cell,
            &points,
            &weights,
            UpdateFlags {
                values: true,
                gradients: false,
                jxw: false,
            },
        )
        .unwrap();
        assert!(only_values.values().is_some());
        assert!(only_values.value(2, 0).is_some());
        assert!(only_values.grad_ref().is_none());
        assert!(only_values.grad_phys().is_none());
        assert!(only_values.jxw().is_none());

        let only_gradients = FeValuesP1Tetrahedron::new(
            &cell,
            &points,
            &weights,
            UpdateFlags {
                values: false,
                gradients: true,
                jxw: false,
            },
        )
        .unwrap();
        assert!(only_gradients.values().is_none());
        assert!(only_gradients.grad_ref().is_some());
        assert!(only_gradients.grad_phys().is_some());
        assert!(only_gradients.jxw().is_none());

        let only_jxw = FeValuesP1Tetrahedron::new(
            &cell,
            &points,
            &weights,
            UpdateFlags {
                values: false,
                gradients: false,
                jxw: true,
            },
        )
        .unwrap();
        assert!(only_jxw.values().is_none());
        assert!(only_jxw.grad_ref().is_none());
        assert!(only_jxw.grad_phys().is_none());
        assert!(only_jxw.jxw().is_some());

        let all = FeValuesP1Tetrahedron::new(&cell, &points, &weights, UpdateFlags::all()).unwrap();
        assert!(all.values().is_some());
        assert!(all.grad_ref().is_some());
        assert!(all.grad_phys().is_some());
        assert!(all.jxw().is_some());

        assert!(all.value(4, 0).is_none());
        assert!(all.value(0, all.n_quadrature_points()).is_none());
        assert!(all.grad_ref_at(0, usize::MAX).is_none());
        assert!(all.grad_phys_at(usize::MAX, 0).is_none());
        assert!(all.jxw_at(usize::MAX).is_none());
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let cell = unit_cell();
        let points = [p(0.25, 0.25, 0.25)];
        let weights = [0.5];

        // A reversed or degenerate cell never reaches the evaluator.
        assert!(AffineTet::new([
            p(0.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 0.0, 1.0)
        ])
        .is_err());
        assert!(AffineTet::new([
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(2.0, 0.0, 0.0),
            p(0.0, 0.0, 1.0)
        ])
        .is_err());
        assert!(AffineTet::new([
            p(0.0, 0.0, 0.0),
            Point3 {
                x: f64::NAN,
                y: 0.0,
                z: 0.0
            },
            p(0.0, 1.0, 0.0),
            p(0.0, 0.0, 1.0),
        ])
        .is_err());

        // Count mismatches and empty rules.
        assert!(FeValuesP1Tetrahedron::new(&cell, &points, &[], UpdateFlags::ALL).is_err());
        assert!(FeValuesP1Tetrahedron::new(
            &cell,
            &[p(0.2, 0.2, 0.2), p(0.3, 0.3, 0.3)],
            &weights,
            UpdateFlags::ALL
        )
        .is_err());
        assert!(FeValuesP1Tetrahedron::new(&cell, &[], &[], UpdateFlags::ALL).is_err());

        // Invalid weights.
        for weight in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY].iter() {
            assert!(
                FeValuesP1Tetrahedron::new(&cell, &points, &[*weight], UpdateFlags::ALL).is_err(),
                "weight {} must be rejected",
                weight
            );
            assert!(
                FeValuesP1Tetrahedron::new(&cell, &points, &[*weight], UpdateFlags::EMPTY).is_err(),
                "weight {} must be rejected even with empty flags",
                weight
            );
        }

        // Non-finite reference coordinates.
        let bad_points = [
            Point3 {
                x: f64::NAN,
                y: 0.0,
                z: 0.0,
            },
            Point3 {
                x: 0.0,
                y: f64::INFINITY,
                z: 0.0,
            },
            Point3 {
                x: 0.0,
                y: 0.0,
                z: f64::NEG_INFINITY,
            },
        ];
        for point in bad_points.iter() {
            assert!(
                FeValuesP1Tetrahedron::new(&cell, &[*point], &weights, UpdateFlags::ALL).is_err()
            );
            assert!(
                FeValuesP1Tetrahedron::new(&cell, &[*point], &weights, UpdateFlags::EMPTY).is_err()
            );
        }

        // determinant check in FeValues (cannot construct bad cell, but ensure
        // a valid cell passes with empty flags).
        assert!(FeValuesP1Tetrahedron::new(&cell, &points, &weights, UpdateFlags::EMPTY).is_ok());
    }
}
