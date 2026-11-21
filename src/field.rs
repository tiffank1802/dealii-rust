//! Interpolation, L2 projection and point evaluation of scalar continuous P1 fields.
use std::collections::BTreeSet;

use crate::base::Point2;
use crate::fe::{shape_values, AffineTriangle};
use crate::grid::{BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::lac::{conjugate_gradient, CgOptions, CgReport, CooBuilder, CsrMatrix};
use crate::quadrature::triangle_degree_5;
use crate::{error, Result};

/// Interpolate a function at every mesh vertex, in P1 degree-of-freedom order.
/// Nonfinite function values are rejected.
pub fn interpolate<F: Fn(Point2) -> f64>(mesh: &Mesh2, f: F) -> Result<Vec<f64>> {
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

/// Interpolate Dirichlet data on a selection of boundary edges. The result is
/// a sorted, deduplicated list of `(vertex index, value)` constraints usable by
/// `CsrMatrix::with_dirichlet`. Repeated/invalid edge IDs are errors; shared
/// corners are sampled only once. The closure may not return nonfinite values.
pub fn interpolate_boundary<F: Fn(Point2) -> f64>(
    mesh: &Mesh2,
    edges: &[BoundaryEdgeId],
    f: F,
) -> Result<Vec<(usize, f64)>> {
    let mut seen = BTreeSet::new();
    let mut vertices = BTreeSet::new();
    for &id in edges {
        if !seen.insert(id) {
            return Err(error("duplicate boundary edge"));
        }
        let edge = mesh.boundary_edge_geometry(id)?;
        vertices.extend(edge.vertices);
    }
    vertices
        .into_iter()
        .map(|id| {
            let value = f(mesh.vertices()[id.0]);
            if !value.is_finite() {
                Err(error("boundary function returned a nonfinite value"))
            } else {
                Ok((id.0, value))
            }
        })
        .collect()
}

/// Assemble the consistent P1 mass matrix `M_ij = ∫ φ_i φ_j`.
/// No boundary constraints are imposed. Mesh vertices are required to belong
/// to a cell by `Mesh2`, making the mass matrix SPD in exact arithmetic.
pub fn mass_matrix(mesh: &Mesh2) -> Result<CsrMatrix> {
    let n = mesh.vertices().len();
    let mut coo = CooBuilder::new(n, n);
    for (i, cell) in mesh.cells().iter().enumerate() {
        let det = AffineTriangle::new(mesh.cell_points(CellId(i))?)?.determinant();
        let off_diagonal = det / 24.0;
        if off_diagonal == 0.0 || !off_diagonal.is_finite() {
            return Err(error("mass matrix entry cannot be represented"));
        }
        for (local_i, row) in cell.iter().enumerate() {
            for (local_j, col) in cell.iter().enumerate() {
                coo.add(
                    row.0,
                    col.0,
                    off_diagonal * if local_i == local_j { 2.0 } else { 1.0 },
                )?;
            }
        }
    }
    coo.into_csr()
}

/// Assemble the consistent mass system `M u = b`, where `b_i = ∫ f φ_i`.
/// The right-hand side uses degree-five triangle quadrature (exact for
/// polynomials f of degree at most four up to floating-point roundoff).
/// No boundary values are imposed; finite values at quadrature points are required.
pub fn projection_system<F: Fn(Point2) -> f64>(
    mesh: &Mesh2,
    f: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let matrix = mass_matrix(mesh)?;
    let mut rhs = vec![0.0; mesh.vertices().len()];
    let quadrature = triangle_degree_5();
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = AffineTriangle::new(mesh.cell_points(CellId(i))?)?;
        for q in &quadrature {
            let value = f(map.map(q.point));
            if !value.is_finite() {
                return Err(error("projection function returned a nonfinite value"));
            }
            let shape = shape_values(q.point);
            let weight = map.determinant() * q.weight;
            for j in 0..3 {
                rhs[cell[j].0] += value * shape[j] * weight;
            }
        }
    }
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(error("nonfinite projection right hand side"));
    }
    Ok((matrix, rhs))
}

/// L2-project `f` onto the unconstrained continuous P1 space. Uses the
/// consistent mass matrix and CG; errors and convergence options are those of
/// [`conjugate_gradient`]. Absolute tolerance can dominate for small data;
/// choose options appropriate to the scale. Projection is quadrature-based and
/// imposes no Dirichlet constraints. For constrained projection, obtain the
/// system with [`projection_system`] and apply constraints before solving.
pub fn project_l2<F: Fn(Point2) -> f64>(
    mesh: &Mesh2,
    f: F,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    let (matrix, rhs) = projection_system(mesh, f)?;
    conjugate_gradient(&matrix, &rhs, options)
}

/// Location of a point in a triangle; barycentric coordinates follow the
/// cell's vertex order. IDs are valid only for the mesh passed to [`locate`].
#[derive(Debug, Clone, Copy)]
pub struct PointLocation {
    pub cell: CellId,
    pub barycentric: [f64; 3],
}

/// Locate a point using a linear scan of the triangles, with dimensionless
/// barycentric tolerance `tolerance >= 0`. An exact numerical containment
/// (`λ_i >= 0`) is preferred to a tolerance-only match, then the first cell
/// in mesh order is chosen. A tolerance-only match may lie outside the mesh;
/// its barycentric weights are **not** clamped, so evaluation extrapolates.
/// Floating-point arithmetic prevents guaranteed boundary classification,
/// especially on very thin/large translated cells; select tolerance according
/// to the application. A point outside returns `Ok(None)`. Nonfinite inputs,
/// tolerance, or intermediate arithmetic return an error. The mesh is assumed
/// conforming and nonoverlapping, as required by `Mesh2`.
pub fn locate(mesh: &Mesh2, point: Point2, tolerance: f64) -> Result<Option<PointLocation>> {
    if !point.x.is_finite() || !point.y.is_finite() || !tolerance.is_finite() || tolerance < 0.0 {
        return Err(error("invalid point or barycentric tolerance"));
    }
    let mut approximate = None;
    for (i, _) in mesh.cells().iter().enumerate() {
        let [a, b, c] = mesh.cell_points(CellId(i))?;
        let determinant = AffineTriangle::new([a, b, c])?.determinant();
        let dx = point.x - a.x;
        let dy = point.y - a.y;
        let l1 = (dx * (c.y - a.y) - dy * (c.x - a.x)) / determinant;
        let l2 = ((b.x - a.x) * dy - (b.y - a.y) * dx) / determinant;
        let barycentric = [1.0 - l1 - l2, l1, l2];
        if barycentric.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite point location arithmetic"));
        }
        let location = PointLocation {
            cell: CellId(i),
            barycentric,
        };
        if barycentric.iter().all(|&v| v >= 0.0) {
            return Ok(Some(location));
        }
        if approximate.is_none() && barycentric.iter().all(|&v| v >= -tolerance) {
            approximate = Some(location);
        }
    }
    Ok(approximate)
}

/// Evaluate a continuous scalar P1 field at a point, with the same search and
/// tolerance behavior as [`locate`]. Returns `None` if no triangle matches.
/// On shared boundaries, either adjacent cell represents the same value for a
/// conforming field. Vector length and finiteness are checked before searching.
pub fn evaluate(
    mesh: &Mesh2,
    values: &[f64],
    point: Point2,
    tolerance: f64,
) -> Result<Option<f64>> {
    if values.len() != mesh.vertices().len() || values.iter().any(|v| !v.is_finite()) {
        return Err(error("invalid P1 field"));
    }
    let location = match locate(mesh, point, tolerance)? {
        Some(p) => p,
        None => return Ok(None),
    };
    let vertices: [VertexId; 3] = mesh.cells()[location.cell.0];
    let value = (0..3)
        .map(|i| location.barycentric[i] * values[vertices[i].0])
        .sum::<f64>();
    if !value.is_finite() {
        return Err(error("nonfinite point evaluation"));
    }
    Ok(Some(value))
}
