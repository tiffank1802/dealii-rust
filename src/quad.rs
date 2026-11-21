//! Experimental scalar Q1 finite elements on conforming 2D quadrilaterals.
//!
//! Reference vertex order is `(0,0), (1,0), (1,1), (0,1)`; a physical cell
//! uses the corresponding counterclockwise cyclic order. Geometry is an
//! isoparametric bilinear map; quadrature is tensor-product Gauss 3×3.
//! This module is separate from the triangle APIs and accepts only convex,
//! nondegenerate cells. Global geometric conformity is the caller's duty.
pub use crate::quad_refinement::UniformRefinement;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

use crate::base::{Point2, UpdateFlags};
use crate::grid::{BoundaryEdgeGeometry, BoundaryEdgeId, CellId, VertexId};
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::{line_degree_5, square_degree_5};
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
    oriented_boundary_edges: Vec<[VertexId; 2]>,
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
        let oriented_boundary_edges = edges
            .values()
            .filter(|(count, _)| *count == 1)
            .map(|(_, vertices)| *vertices)
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
            oriented_boundary_edges,
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
    /// Directed edge as in its counterclockwise cell, with right-hand
    /// outward normal. The ID is valid only for this quadrilateral mesh.
    pub fn boundary_edge_geometry(&self, id: BoundaryEdgeId) -> Result<BoundaryEdgeGeometry> {
        let vertices = *self
            .oriented_boundary_edges
            .get(id.0)
            .ok_or_else(|| error("boundary edge index out of bounds"))?;
        let points = [self.vertices[vertices[0].0], self.vertices[vertices[1].0]];
        let dx = points[1].x - points[0].x;
        let dy = points[1].y - points[0].y;
        let length = dx.hypot(dy);
        if !length.is_finite() || length <= 0.0 {
            return Err(error("invalid boundary edge length"));
        }
        Ok(BoundaryEdgeGeometry {
            vertices,
            points,
            outward_normal: [dy / length, -dx / length],
            length,
        })
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

/// Interpolate Dirichlet data on selected boundary edges. Returns sorted,
/// deduplicated `(vertex index, value)` constraints for `with_dirichlet`.
/// Edge IDs must belong to this mesh, be unique and in bounds. Nonfinite
/// function values are rejected; a shared vertex is sampled once.
pub fn interpolate_boundary<F: Fn(Point2) -> f64>(
    mesh: &QuadMesh,
    edges: &[BoundaryEdgeId],
    value: F,
) -> Result<Vec<(usize, f64)>> {
    let mut seen = BTreeSet::new();
    let mut vertices = BTreeSet::new();
    for &id in edges {
        if !seen.insert(id) {
            return Err(error("duplicate boundary edge"));
        }
        let edge = mesh
            .boundary_edges
            .get(id.0)
            .ok_or_else(|| error("boundary edge index out of bounds"))?;
        vertices.extend(edge.iter().copied());
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

/// Assemble `∫g φ_i ds` on selected Q1 boundary edges for `-Δu=f`.
/// The callback receives the physical point and the unit **outward** normal;
/// it supplies `g=∂u/∂n`, or `g=n·A∇u` for anisotropic diffusion.
/// Add this load before Dirichlet elimination.
/// Edge IDs must belong to this mesh, be unique and in bounds. Boundary
/// edges of a bilinear quadrilateral are straight and the Q1 trace is linear;
/// the degree-five Gauss rule is exact for polynomial fluxes through degree 4
/// (up to rounding). Nonfinite flux or accumulated loads are rejected.
/// Pure Neumann systems require compatibility and normalization; this routine
/// assembles only the load and does not solve the singular system.
pub fn assemble_neumann<F>(mesh: &QuadMesh, edges: &[BoundaryEdgeId], flux: F) -> Result<Vec<f64>>
where
    F: Fn(Point2, [f64; 2]) -> f64,
{
    let mut seen = BTreeSet::new();
    for &id in edges {
        if !seen.insert(id) {
            return Err(error("duplicate Neumann boundary edge"));
        }
        mesh.boundary_edge_geometry(id)?;
    }
    let mut load = vec![0.0; mesh.vertices.len()];
    for &id in edges {
        let edge = mesh.boundary_edge_geometry(id)?;
        for q in line_degree_5() {
            let t = q.point;
            let p = Point2::new(
                (1.0 - t) * edge.points[0].x + t * edge.points[1].x,
                (1.0 - t) * edge.points[0].y + t * edge.points[1].y,
            )?;
            let g = flux(p, edge.outward_normal);
            if !g.is_finite() {
                return Err(error("Neumann flux returned a nonfinite value"));
            }
            let weight = edge.length * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid boundary quadrature weight"));
            }
            load[edge.vertices[0].0] += g * (1.0 - t) * weight;
            load[edge.vertices[1].0] += g * t * weight;
        }
    }
    if load.iter().any(|x| !x.is_finite()) {
        return Err(error("nonfinite Neumann load"));
    }
    Ok(load)
}

/// Number of shape functions of the Q1 quadrilateral.
const Q1_QUAD_N_SHAPES: usize = 4;

/// FEValues-style evaluation of the Q1 quadrilateral on a bilinear cell.
///
/// The evaluation is built once for a cell (via [`BilinearMap`]), a quadrature
/// rule and a set of [`UpdateFlags`], and then exposes the requested
/// quantities. Quantities whose flag is false are not computed at all: the
/// matching accessor returns `None`.
///
/// Arrays are stored shape-major, so `phi[i][q]` lives at index
/// `i * n_quadrature_points() + q`; the per-entry accessors ([`value`],
/// [`grad_ref_at`], [`grad_phys_at`], [`jxw_at`]) hide that layout.
///
/// Physical quadrature points are always computed and are not flag-gated.
/// The Jacobian varies per quadrature point, so physical gradients and `JxW`
/// are evaluated pointwise.
///
/// [`value`]: FeValuesQ1Quad::value
/// [`grad_ref_at`]: FeValuesQ1Quad::grad_ref_at
/// [`grad_phys_at`]: FeValuesQ1Quad::grad_phys_at
/// [`jxw_at`]: FeValuesQ1Quad::jxw_at
#[derive(Debug, Clone)]
pub struct FeValuesQ1Quad {
    n_quadrature_points: usize,
    flags: UpdateFlags,
    quadrature_points: Vec<Point2>,
    values: Option<Vec<f64>>,
    grad_ref: Option<Vec<[f64; 2]>>,
    grad_phys: Option<Vec<[f64; 2]>>,
    jxw: Option<Vec<f64>>,
}

impl FeValuesQ1Quad {
    /// Evaluate the Q1 quadrilateral on `map` at the given reference points.
    ///
    /// `reference_points` and `weights` must have the same non-zero length. All
    /// reference coordinates must be finite and inside the unit square
    /// `[0, 1]^2`, and all weights must be finite and strictly positive. `map`
    /// must be a valid convex, counterclockwise, non-degenerate bilinear map,
    /// as enforced by [`BilinearMap::new`].
    ///
    /// Input validation always runs, including when `flags` is empty.
    pub fn new(
        map: &BilinearMap,
        reference_points: &[Point2],
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
            if !point.x.is_finite() || !point.y.is_finite() {
                return Err(error("quadrature point coordinates must be finite"));
            }
        }
        for weight in weights {
            if !weight.is_finite() || *weight <= 0.0 {
                return Err(error("quadrature weights must be finite and positive"));
            }
        }

        let n = reference_points.len();

        // Validate every reference point against the mapping (finiteness, unit
        // square, positive Jacobian) and cache the per-point geometry.
        let mut quadrature_points = Vec::with_capacity(n);
        let mut determinants = Vec::with_capacity(n);
        let mut mapped_gradients = Vec::with_capacity(n);
        for point in reference_points {
            let mapped = map.evaluate(*point)?;
            quadrature_points.push(mapped.point);
            determinants.push(mapped.determinant);
            mapped_gradients.push(mapped.gradients);
        }

        let values = if flags.values {
            let mut data = vec![0.0; Q1_QUAD_N_SHAPES * n];
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
            let mut data = vec![[0.0; 2]; Q1_QUAD_N_SHAPES * n];
            for (q, point) in reference_points.iter().enumerate() {
                let refs = reference_gradients(*point);
                for i in 0..Q1_QUAD_N_SHAPES {
                    data[i * n + q] = refs[i];
                }
            }
            Some(data)
        } else {
            None
        };

        let grad_phys = if flags.gradients {
            let mut data = vec![[0.0; 2]; Q1_QUAD_N_SHAPES * n];
            for q in 0..n {
                for i in 0..Q1_QUAD_N_SHAPES {
                    data[i * n + q] = mapped_gradients[q][i];
                }
            }
            Some(data)
        } else {
            None
        };

        let jxw = if flags.jxw {
            Some(
                determinants
                    .iter()
                    .zip(weights.iter())
                    .map(|(det, w)| det.abs() * w)
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
    pub fn quadrature_points(&self) -> &[Point2] {
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
    pub fn grad_ref(&self) -> Option<&[[f64; 2]]> {
        self.grad_ref.as_deref()
    }

    /// Reference gradient of `shape` at point `q`; `None` if gradients were not
    /// requested or out of range.
    pub fn grad_ref_at(&self, shape: usize, q: usize) -> Option<[f64; 2]> {
        let index = self.flat_index(shape, q)?;
        self.grad_ref.as_ref()?.get(index).copied()
    }

    /// Physical gradients `J^{-T} grad_ref`, shape-major;
    /// `None` unless `gradients` was requested.
    pub fn grad_phys(&self) -> Option<&[[f64; 2]]> {
        self.grad_phys.as_deref()
    }

    /// Physical gradient of `shape` at point `q`; `None` if gradients were not
    /// requested or out of range.
    pub fn grad_phys_at(&self, shape: usize, q: usize) -> Option<[f64; 2]> {
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
        if shape >= Q1_QUAD_N_SHAPES || q >= self.n_quadrature_points {
            return None;
        }
        Some(shape * self.n_quadrature_points + q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::UpdateFlags;
    use crate::quadrature::square_degree_5;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y).unwrap()
    }

    fn near(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{} != {} (tolerance {})", a, b, tol);
    }

    /// Reference vertices of the Q1 quadrilateral, in vertex order.
    fn reference_vertices() -> [Point2; 4] {
        [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]
    }

    /// Counterclockwise, convex and non-affine: the Jacobian varies per point.
    fn nonaffine_points() -> [Point2; 4] {
        [p(0.0, 0.0), p(2.0, 0.0), p(1.8, 1.5), p(0.2, 1.0)]
    }

    fn unit_square_points() -> [Point2; 4] {
        [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)]
    }

    fn degree_five_rule() -> (Vec<Point2>, Vec<f64>) {
        let rule = square_degree_5();
        (
            rule.iter().map(|q| q.point).collect(),
            rule.iter().map(|q| q.weight).collect(),
        )
    }

    /// Jacobian derived independently from the bilinear map
    /// `x(xi,eta) = a + xi (b-a) + eta (d-a) + xi eta ((c-d) - (b-a))`,
    /// written differently from the production formula on purpose.
    /// Rows are physical `(x, y)`, columns reference `(xi, eta)`.
    fn independent_jacobian(points: [Point2; 4], q: Point2) -> [[f64; 2]; 2] {
        let [a, b, c, d] = points;
        let derivatives = |a: f64, b: f64, c: f64, d: f64| {
            let cross = (c - d) - (b - a);
            [(b - a) + q.y * cross, (d - a) + q.x * cross]
        };
        [
            derivatives(a.x, b.x, c.x, d.x),
            derivatives(a.y, b.y, c.y, d.y),
        ]
    }

    /// Independent `J^{-T} r`: invert, transpose, then multiply.
    fn apply_inverse_transpose(jacobian: [[f64; 2]; 2], r: [f64; 2]) -> [f64; 2] {
        let det = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
        let inverse = [
            [jacobian[1][1] / det, -jacobian[0][1] / det],
            [-jacobian[1][0] / det, jacobian[0][0] / det],
        ];
        let inverse_transpose = [
            [inverse[0][0], inverse[1][0]],
            [inverse[0][1], inverse[1][1]],
        ];
        [
            inverse_transpose[0][0] * r[0] + inverse_transpose[0][1] * r[1],
            inverse_transpose[1][0] * r[0] + inverse_transpose[1][1] * r[1],
        ]
    }

    /// ADR-0003 invariant 1, generalised to the Q1 quadrilateral.
    #[test]
    fn values_satisfy_partition_of_unity() {
        let (points, weights) = degree_five_rule();
        for cell in [unit_square_points(), nonaffine_points()].iter() {
            let map = BilinearMap::new(*cell).unwrap();
            let fe = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::ALL).unwrap();
            assert_eq!(fe.n_quadrature_points(), points.len());
            assert_eq!(fe.values().unwrap().len(), 4 * points.len());
            for q in 0..fe.n_quadrature_points() {
                let sum: f64 = (0..4).map(|i| fe.value(i, q).unwrap()).sum();
                near(sum, 1.0, 1e-12);
            }
        }
    }

    /// ADR-0003 invariant 2 at the four reference vertices.
    #[test]
    fn values_are_nodal_at_reference_vertices() {
        let vertices = reference_vertices();
        let weights = [0.25; 4];
        let flags = UpdateFlags {
            values: true,
            gradients: false,
            jxw: false,
        };
        let map = BilinearMap::new(unit_square_points()).unwrap();
        let fe = FeValuesQ1Quad::new(&map, &vertices, &weights, flags).unwrap();
        for (q, vertex) in vertices.iter().enumerate() {
            let expected = shape_values(*vertex);
            for (i, value) in expected.iter().enumerate() {
                let delta = if i == q { 1.0 } else { 0.0 };
                near(fe.value(i, q).unwrap(), delta, 1e-14);
                near(fe.value(i, q).unwrap(), *value, 1e-14);
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
        let map = BilinearMap::new(nonaffine_points()).unwrap();
        let fe = FeValuesQ1Quad::new(&map, &points, &weights, flags).unwrap();
        assert!(fe.values().is_none());
        for q in 0..fe.n_quadrature_points() {
            let mut sum = [0.0, 0.0];
            for i in 0..4 {
                let gradient = fe.grad_ref_at(i, q).unwrap();
                sum[0] += gradient[0];
                sum[1] += gradient[1];
            }
            near(sum[0], 0.0, 1e-15);
            near(sum[1], 0.0, 1e-15);
        }
    }

    /// ADR-0003 invariant 4 on a valid convex, non-degenerate quadrilateral.
    #[test]
    fn jxw_is_positive_on_a_valid_cell() {
        let (points, weights) = degree_five_rule();
        for cell in [unit_square_points(), nonaffine_points()].iter() {
            let map = BilinearMap::new(*cell).unwrap();
            let fe = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::all()).unwrap();
            for q in 0..fe.n_quadrature_points() {
                let jxw = fe.jxw_at(q).unwrap();
                assert!(jxw > 0.0, "JxW must be positive, got {}", jxw);
            }
        }
    }

    /// ADR-0003 invariant 6, with the Jacobian recomputed independently.
    #[test]
    fn jxw_equals_absolute_determinant_times_weight() {
        let (points, weights) = degree_five_rule();
        let cell = nonaffine_points();
        let map = BilinearMap::new(cell).unwrap();
        let fe = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::all()).unwrap();
        for (q, point) in points.iter().enumerate() {
            let j = independent_jacobian(cell, *point);
            let det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
            let expected = det.abs() * weights[q];
            let jxw = fe.jxw_at(q).unwrap();
            near(jxw, expected, 1e-13);
            near(fe.jxw().unwrap()[q], jxw, 1e-15);
        }
    }

    /// ADR-0003 invariant 5, with `J^{-T}` rebuilt independently.
    #[test]
    fn physical_gradients_apply_inverse_transpose_jacobian() {
        let (points, weights) = degree_five_rule();
        let cell = nonaffine_points();
        let map = BilinearMap::new(cell).unwrap();
        let fe = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::all()).unwrap();
        for (q, point) in points.iter().enumerate() {
            let jacobian = independent_jacobian(cell, *point);
            for i in 0..4 {
                let reference = fe.grad_ref_at(i, q).unwrap();
                let expected = apply_inverse_transpose(jacobian, reference);
                let physical = fe.grad_phys_at(i, q).unwrap();
                near(physical[0], expected[0], 1e-13);
                near(physical[1], expected[1], 1e-13);
            }
        }
        // Physical quadrature points are always available.
        for (q, point) in points.iter().enumerate() {
            assert_eq!(
                fe.quadrature_points()[q],
                map.evaluate(*point).unwrap().point
            );
        }
    }

    /// ADR-0003 invariant 7.
    #[test]
    fn flags_gate_computation() {
        let (points, weights) = degree_five_rule();
        let map = BilinearMap::new(unit_square_points()).unwrap();

        let none = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::EMPTY).unwrap();
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

        let only_values = FeValuesQ1Quad::new(
            &map,
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
        assert!(only_values.value(3, 0).is_some());
        assert!(only_values.grad_ref().is_none());
        assert!(only_values.grad_phys().is_none());
        assert!(only_values.jxw().is_none());

        let only_gradients = FeValuesQ1Quad::new(
            &map,
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

        let only_jxw = FeValuesQ1Quad::new(
            &map,
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

        let all = FeValuesQ1Quad::new(&map, &points, &weights, UpdateFlags::all()).unwrap();
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

    /// Invalid cells never produce a map; invalid rules never produce an
    /// evaluation, and rule validation runs even with empty flags.
    #[test]
    fn invalid_inputs_are_rejected() {
        // Reversed (clockwise), degenerate and non-finite cells.
        assert!(BilinearMap::new([p(0.0, 0.0), p(0.0, 1.0), p(1.0, 1.0), p(1.0, 0.0)]).is_err());
        assert!(BilinearMap::new([p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0), p(1.0, 0.0)]).is_err());
        assert!(BilinearMap::new([
            p(0.0, 0.0),
            Point2 {
                x: f64::NAN,
                y: 0.0,
            },
            p(1.0, 1.0),
            p(0.0, 1.0),
        ])
        .is_err());

        let map = BilinearMap::new(unit_square_points()).unwrap();
        let points = [p(0.25, 0.25)];
        let weights = [0.5];

        // Count mismatches and empty rules.
        assert!(FeValuesQ1Quad::new(&map, &points, &[], UpdateFlags::ALL).is_err());
        assert!(FeValuesQ1Quad::new(&map, &points, &[], UpdateFlags::EMPTY).is_err());
        assert!(FeValuesQ1Quad::new(
            &map,
            &[p(0.2, 0.2), p(0.3, 0.3)],
            &weights,
            UpdateFlags::ALL
        )
        .is_err());
        assert!(FeValuesQ1Quad::new(&map, &[], &[], UpdateFlags::ALL).is_err());

        // Invalid weights.
        for weight in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY].iter() {
            assert!(
                FeValuesQ1Quad::new(&map, &points, &[*weight], UpdateFlags::ALL).is_err(),
                "weight {} must be rejected",
                weight
            );
            assert!(
                FeValuesQ1Quad::new(&map, &points, &[*weight], UpdateFlags::EMPTY).is_err(),
                "weight {} must be rejected even with empty flags",
                weight
            );
        }

        // Non-finite reference coordinates and points outside the unit square.
        let bad_points = [
            Point2 {
                x: f64::NAN,
                y: 0.0,
            },
            Point2 {
                x: 0.0,
                y: f64::INFINITY,
            },
            p(-0.1, 0.5),
            p(0.5, 1.1),
        ];
        for point in bad_points.iter() {
            assert!(FeValuesQ1Quad::new(&map, &[*point], &weights, UpdateFlags::ALL).is_err());
            assert!(FeValuesQ1Quad::new(&map, &[*point], &weights, UpdateFlags::EMPTY).is_err());
        }
    }
}
