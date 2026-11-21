//! Interpolation, L2 projection and physical-point evaluation of scalar Q1 fields.
use crate::base::Point2;
use crate::grid::CellId;
use crate::lac::{conjugate_gradient, CgOptions, CgReport, CooBuilder, CsrMatrix};
use crate::quad::{self, BilinearMap, QuadMesh};
use crate::quadrature::square_degree_5;
use crate::{error, Result};

/// Sample a scalar function at the vertices, in the mesh's Q1 DoF order.
/// Nonfinite function values are rejected. Boundary interpolation on selected
/// edges is available separately as [`quad::interpolate_boundary`].
pub fn interpolate<F: Fn(Point2) -> f64>(mesh: &QuadMesh, f: F) -> Result<Vec<f64>> {
    mesh.vertices()
        .iter()
        .map(|&p| {
            let value = f(p);
            if !value.is_finite() {
                Err(error("interpolation function returned a nonfinite value"))
            } else {
                Ok(value)
            }
        })
        .collect()
}

/// Assemble the consistent Q1 mass matrix `M_ij = ∫φ_i φ_j` with no
/// constraints. The 3×3 tensor Gauss rule exactly integrates these entries
/// for bilinear quadrilateral mappings in exact arithmetic. Positive
/// quadrature and use of every vertex yield an SPD matrix in exact arithmetic
/// for conforming nonoverlapping meshes; floating-point conditioning is not
/// guaranteed. Unrepresentable or nonfinite local mass entries are rejected.
pub fn mass_matrix(mesh: &QuadMesh) -> Result<CsrMatrix> {
    let n = mesh.vertices().len();
    let mut builder = CooBuilder::new(n, n);
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i))?)?;
        let mut local = [[0.0; 4]; 4];
        for q in square_degree_5() {
            let mapped = map.evaluate(q.point)?;
            let weight = mapped.determinant * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            let shape = quad::shape_values(q.point);
            for a in 0..4 {
                for b in a..4 {
                    local[a][b] += shape[a] * shape[b] * weight;
                }
            }
        }
        for a in 0..4 {
            for b in a..4 {
                let value = local[a][b];
                if !value.is_finite() || value <= 0.0 {
                    return Err(error("mass matrix entry cannot be represented"));
                }
                builder.add(cell[a].0, cell[b].0, value)?;
                if a != b {
                    builder.add(cell[b].0, cell[a].0, value)?;
                }
            }
        }
    }
    builder.into_csr()
}

/// Assemble the consistent system `M u=b`, with `b_i=∫ f φ_i`. Coefficients
/// are sampled at the 3×3 tensor Gauss points; no boundary values are imposed.
/// The rule is exact up to roundoff for physical polynomials of total degree
/// ≤3 on bilinear quadrilaterals (≤4 on affine cells); other functions incur
/// quadrature error. The mesh must be conforming and nonoverlapping. Nonfinite
/// callback values, quadrature weights or accumulated loads are errors.
pub fn projection_system<F: Fn(Point2) -> f64>(
    mesh: &QuadMesh,
    f: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let matrix = mass_matrix(mesh)?;
    let mut rhs = vec![0.0; mesh.vertices().len()];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = BilinearMap::new(mesh.cell_points(CellId(i))?)?;
        for q in square_degree_5() {
            let mapped = map.evaluate(q.point)?;
            let value = f(mapped.point);
            if !value.is_finite() {
                return Err(error("projection function returned a nonfinite value"));
            }
            let weight = mapped.determinant * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            let shape = quad::shape_values(q.point);
            for j in 0..4 {
                rhs[cell[j].0] += value * shape[j] * weight;
            }
        }
    }
    if rhs.iter().any(|x| !x.is_finite()) {
        return Err(error("nonfinite projection right hand side"));
    }
    Ok((matrix, rhs))
}

/// Project a function onto the unconstrained continuous Q1 space by CG.
/// The quadrature limitations of [`projection_system`] apply. Absolute CG
/// tolerance can dominate for small data; select options for the problem's
/// scale. For a constrained projection, use [`projection_system`], impose
/// constraints via `CsrMatrix::with_dirichlet`, then solve explicitly.
pub fn project_l2<F: Fn(Point2) -> f64>(
    mesh: &QuadMesh,
    f: F,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    let (matrix, rhs) = projection_system(mesh, f)?;
    conjugate_gradient(&matrix, &rhs, options)
}

/// Location of a physical point in a quadrilateral.
#[derive(Debug, Clone, Copy)]
pub struct PointLocation {
    /// Cell ID in the mesh passed to [`locate`].
    pub cell: CellId,
    /// Coordinates on `[0,1]²`, in the Q1 reference vertex convention.
    pub reference: Point2,
    /// Physical point used for evaluation; projected onto an edge when
    /// [`projected`](Self::projected) is true.
    pub point: Point2,
    /// True for a tolerance-only match outside the selected cell.
    pub projected: bool,
    /// Physical distance from the query to `point`; zero for an inside match.
    pub distance: f64,
}

// Edge-relative calculations avoid products of large physical coordinates.
// None of these tests proves global mesh conformity or correct classification
// at every floating-point boundary, especially on very thin/translated cells.
fn edge_geometry(a: Point2, b: Point2, point: Point2) -> Result<(f64, f64, f64)> {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let length = dx.hypot(dy);
    let qx = point.x - a.x;
    let qy = point.y - a.y;
    if !length.is_finite() || length <= 0.0 || !qx.is_finite() || !qy.is_finite() {
        return Err(error("nonfinite point-location geometry"));
    }
    let tangent = (dx / length, dy / length);
    let along = tangent.0 * qx + tangent.1 * qy;
    let inward = -tangent.1 * qx + tangent.0 * qy;
    if !along.is_finite() || !inward.is_finite() {
        return Err(error("nonfinite point-location arithmetic"));
    }
    Ok((length, along, inward))
}

// Invert the bilinear map for a point already classified in the convex cell.
// Horizontal reference lines do not cross inside the cell. The signed
// distance to the line at eta changes sign at the inverse eta, so bisection
// avoids an unconstrained Newton iteration and remains inside [0,1].
fn inverse_inside(points: [Point2; 4], point: Point2) -> Result<Point2> {
    let [a, b, c, d] = points;
    let q = [point.x - a.x, point.y - a.y];
    let bottom = [b.x - a.x, b.y - a.y];
    let left = [d.x - a.x, d.y - a.y];
    let top = [c.x - d.x, c.y - d.y];
    let line = |eta: f64| -> Result<([f64; 2], [f64; 2], f64)> {
        let direction = [
            (1.0 - eta) * bottom[0] + eta * top[0],
            (1.0 - eta) * bottom[1] + eta * top[1],
        ];
        let from_left = [q[0] - eta * left[0], q[1] - eta * left[1]];
        let length = direction[0].hypot(direction[1]);
        if !length.is_finite() || length <= 0.0 || from_left.iter().any(|x| !x.is_finite()) {
            return Err(error("nonfinite inverse mapping geometry"));
        }
        let side = -(direction[1] / length) * from_left[0] + (direction[0] / length) * from_left[1];
        if !side.is_finite() {
            return Err(error("nonfinite inverse mapping arithmetic"));
        }
        Ok((direction, from_left, side))
    };
    let f0 = line(0.0)?.2;
    let f1 = line(1.0)?.2;
    let eta = if f0 <= 0.0 {
        0.0
    } else if f1 >= 0.0 {
        1.0
    } else {
        let mut lo = 0.0;
        let mut hi = 1.0;
        for _ in 0..64 {
            let mid = (lo + hi) * 0.5;
            if mid == lo || mid == hi {
                break;
            }
            if line(mid)?.2 > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        (lo + hi) * 0.5
    };
    let (direction, diff, _) = line(eta)?;
    let length = direction[0].hypot(direction[1]);
    let xi = ((direction[0] / length) * diff[0] + (direction[1] / length) * diff[1]) / length;
    if !xi.is_finite() {
        return Err(error("nonfinite inverse mapping coordinates"));
    }
    let reference = Point2 {
        x: xi.clamp(0.0, 1.0),
        y: eta,
    };
    let mapped = BilinearMap::new(points)?.evaluate(reference)?.point;
    let residual = (mapped.x - point.x).hypot(mapped.y - point.y);
    let scale = points
        .iter()
        .fold(point.x.abs().max(point.y.abs()), |m, p| {
            m.max(p.x.abs()).max(p.y.abs())
        });
    // Includes translation-dependent representability in the comparison.
    let tolerance = 256.0 * f64::EPSILON * scale.max(length);
    if !residual.is_finite() || !tolerance.is_finite() || residual > tolerance {
        return Err(error("inverse mapping failed to reproduce physical point"));
    }
    Ok(reference)
}

/// Locate a physical point by a linear scan of convex Q1 quadrilaterals.
/// The inverse of the bilinear mapping is computed inside the cell by bounded
/// bisection. If several cells contain the point, the first cell in mesh order
/// is returned. Exact numerical containment is preferred over tolerance-only
/// matches, even if an earlier cell matched the tolerance.
///
/// `tolerance >= 0` is dimensionless: a point outside a cell is accepted
/// within `tolerance * max_edge_length` of that cell. A tolerance-only match
/// is **projected to the nearest point on its boundary**, not extrapolated;
/// the first eligible cell is selected, not necessarily the nearest cell.
/// This differs from the triangle API `field::locate`. `Ok(None)` indicates
/// no match. Invalid inputs and nonrepresentable/inconsistent arithmetic are
/// errors. Floating-point classifications on edges and very thin or highly
/// translated cells are not guaranteed. The caller supplies a conforming,
/// nonoverlapping mesh; no spatial index is used.
pub fn locate(mesh: &QuadMesh, point: Point2, tolerance: f64) -> Result<Option<PointLocation>> {
    if !point.x.is_finite() || !point.y.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
        return Err(error("invalid point or Q1 location tolerance"));
    }
    let mut approximate = None;
    for (i, _) in mesh.cells().iter().enumerate() {
        let id = CellId(i);
        let points = mesh.cell_points(id)?;
        let mut edges = [(0.0, 0.0, 0.0); 4];
        let mut scale: f64 = 0.0;
        for j in 0..4 {
            edges[j] = edge_geometry(points[j], points[(j + 1) % 4], point)?;
            scale = scale.max(edges[j].0);
        }
        if edges.iter().all(|e| e.2 >= 0.0) {
            return Ok(Some(PointLocation {
                cell: id,
                reference: inverse_inside(points, point)?,
                point,
                projected: false,
                distance: 0.0,
            }));
        }
        if approximate.is_some() || tolerance == 0.0 {
            continue;
        }
        let limit = tolerance * scale;
        if !limit.is_finite() {
            return Err(error("nonfinite Q1 location tolerance scale"));
        }
        let mut nearest: Option<PointLocation> = None;
        for j in 0..4 {
            let a = points[j];
            let b = points[(j + 1) % 4];
            let t = (edges[j].1 / edges[j].0).clamp(0.0, 1.0);
            let projected = Point2::new((1.0 - t) * a.x + t * b.x, (1.0 - t) * a.y + t * b.y)?;
            let distance = (point.x - projected.x).hypot(point.y - projected.y);
            if !distance.is_finite() {
                return Err(error("nonfinite point projection distance"));
            }
            let reference = match j {
                0 => Point2 { x: t, y: 0.0 },
                1 => Point2 { x: 1.0, y: t },
                2 => Point2 { x: 1.0 - t, y: 1.0 },
                _ => Point2 { x: 0.0, y: 1.0 - t },
            };
            if nearest.as_ref().map_or(true, |old| distance < old.distance) {
                nearest = Some(PointLocation {
                    cell: id,
                    reference,
                    point: projected,
                    projected: true,
                    distance,
                });
            }
        }
        if let Some(location) = nearest {
            if location.distance <= limit {
                approximate = Some(location);
            }
        }
    }
    Ok(approximate)
}

/// Evaluate a scalar Q1 nodal field at a physical point. The search, tolerance
/// and boundary projection behavior are those of [`locate`]. On a conforming
/// shared edge, either neighboring cell gives the same field value in exact
/// arithmetic. Invalid field dimensions or nonfinite values are rejected.
pub fn evaluate(
    mesh: &QuadMesh,
    values: &[f64],
    point: Point2,
    tolerance: f64,
) -> Result<Option<f64>> {
    if values.len() != mesh.vertices().len() || values.iter().any(|x| !x.is_finite()) {
        return Err(error("invalid Q1 field"));
    }
    let location = match locate(mesh, point, tolerance)? {
        Some(location) => location,
        None => return Ok(None),
    };
    let cell = mesh.cells()[location.cell.0];
    let shape = quad::shape_values(location.reference);
    let value: f64 = (0..4).map(|i| shape[i] * values[cell[i].0]).sum();
    if !value.is_finite() {
        return Err(error("nonfinite Q1 field evaluation"));
    }
    Ok(Some(value))
}
