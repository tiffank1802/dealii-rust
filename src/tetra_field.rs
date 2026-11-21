//! Interpolation, L2 projection and point evaluation of scalar P1 tetrahedral fields.
use crate::base::Point3;
use crate::grid::CellId;
use crate::lac::{conjugate_gradient, CgOptions, CgReport, CooBuilder, CsrMatrix};
use crate::quadrature::tetra_degree_5;
use crate::tetra::{self, AffineTet, TetMesh};
use crate::{error, Result};

/// Sample a scalar function at the mesh vertices in P1 degree-of-freedom order.
/// Nonfinite values are rejected. For boundary data on selected faces, use
/// [`tetra::interpolate_boundary`].
pub fn interpolate<F: Fn(Point3) -> f64>(mesh: &TetMesh, f: F) -> Result<Vec<f64>> {
    mesh.vertices()
        .iter()
        .map(|&p| {
            let value = f(p);
            if value.is_finite() {
                Ok(value)
            } else {
                Err(error("interpolation function returned a nonfinite value"))
            }
        })
        .collect()
}

/// Assemble the unconstrained consistent P1 mass matrix `M_ij = ∫ φ_i φ_j`.
/// The local entries are `volume/20` off diagonal and `volume/10` on the
/// diagonal. For a conforming, nonoverlapping mesh where all vertices belong
/// to a cell, it is SPD in exact arithmetic; conditioning is not guaranteed.
/// Unrepresentable entries and nonfinite sums are rejected.
pub fn mass_matrix(mesh: &TetMesh) -> Result<CsrMatrix> {
    let mut builder = CooBuilder::new(mesh.vertices().len(), mesh.vertices().len());
    for (i, cell) in mesh.cells().iter().enumerate() {
        let determinant = AffineTet::new(mesh.cell_points(CellId(i))?)?.determinant();
        let off = determinant / 120.0;
        let diagonal = off * 2.0;
        if off <= 0.0 || !off.is_finite() || !diagonal.is_finite() {
            return Err(error("mass matrix entry cannot be represented"));
        }
        for a in 0..4 {
            for b in 0..4 {
                builder.add(cell[a].0, cell[b].0, if a == b { diagonal } else { off })?;
            }
        }
    }
    builder.into_csr()
}

/// Assemble the unconstrained mass system `M u = b`, where `b_i = ∫ f φ_i`.
/// The degree-five tetrahedral quadrature is exact (up to roundoff) for
/// physical polynomial `f` of total degree at most four; other functions incur
/// quadrature error. Callback values, weights and accumulated loads must be
/// finite and representable.
pub fn projection_system<F: Fn(Point3) -> f64>(
    mesh: &TetMesh,
    f: F,
) -> Result<(CsrMatrix, Vec<f64>)> {
    let matrix = mass_matrix(mesh)?;
    let mut rhs = vec![0.0; mesh.vertices().len()];
    let quadrature = tetra_degree_5();
    for (i, cell) in mesh.cells().iter().enumerate() {
        let map = AffineTet::new(mesh.cell_points(CellId(i))?)?;
        for q in &quadrature {
            let value = f(map.map(q.point)?);
            if !value.is_finite() {
                return Err(error("projection function returned a nonfinite value"));
            }
            let weight = map.determinant() * q.weight;
            if !weight.is_finite() || weight <= 0.0 {
                return Err(error("invalid quadrature weight"));
            }
            let shape = tetra::shape_values(q.point);
            for j in 0..4 {
                rhs[cell[j].0] += value * shape[j] * weight;
            }
        }
    }
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(error("nonfinite projection right hand side"));
    }
    Ok((matrix, rhs))
}

/// L2-project onto the unconstrained continuous P1 space using CG and the
/// quadrature of [`projection_system`]. Adjust CG tolerances to the scale of
/// the data. For a constrained projection, apply `CsrMatrix::with_dirichlet`
/// to the assembled system before solving.
pub fn project_l2<F: Fn(Point3) -> f64>(
    mesh: &TetMesh,
    f: F,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    let (matrix, rhs) = projection_system(mesh, f)?;
    conjugate_gradient(&matrix, &rhs, options)
}

/// Location of a physical point in a tetrahedron. Barycentric coordinates
/// follow the local vertex order; the cell ID belongs to the input mesh.
#[derive(Debug, Clone, Copy)]
pub struct PointLocation {
    /// Cell ID in the mesh passed to [`locate`].
    pub cell: CellId,
    /// Barycentric weights in the local vertex order (not clamped).
    pub barycentric: [f64; 4],
}

/// Locate a physical point by a linear scan of the tetrahedra. The tolerance
/// is dimensionless and applies to each barycentric weight: `λ_i >= -tolerance`.
/// The first numerically containing cell (`λ_i >= 0`) is preferred even over
/// an earlier tolerance-only match; otherwise the first compatible cell wins.
/// Tolerance-only matches **extrapolate**, with no clamping or projection,
/// as in [`crate::field::locate`] and unlike [`crate::quad_field::locate`].
/// `Ok(None)` means no cell matched. Nonfinite inputs and arithmetic are
/// errors. Boundary classification can be sensitive to rounding, particularly
/// for thin or strongly translated cells. A conforming, nonoverlapping mesh is
/// assumed; no spatial index or global intersection check is provided.
pub fn locate(mesh: &TetMesh, point: Point3, tolerance: f64) -> Result<Option<PointLocation>> {
    if !point.x.is_finite()
        || !point.y.is_finite()
        || !point.z.is_finite()
        || !tolerance.is_finite()
        || tolerance < 0.0
    {
        return Err(error("invalid point or barycentric tolerance"));
    }
    let mut approximate = None;
    for i in 0..mesh.cells().len() {
        let id = CellId(i);
        let points = mesh.cell_points(id)?;
        let gradients = AffineTet::new(points)?.gradients();
        let diff = [
            point.x - points[0].x,
            point.y - points[0].y,
            point.z - points[0].z,
        ];
        if diff.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite point location arithmetic"));
        }
        let lambda = |j: usize| -> f64 {
            gradients[j][0] * diff[0] + gradients[j][1] * diff[1] + gradients[j][2] * diff[2]
        };
        let l1 = lambda(1);
        let l2 = lambda(2);
        let l3 = lambda(3);
        let barycentric = [1.0 - l1 - l2 - l3, l1, l2, l3];
        if barycentric.iter().any(|v| !v.is_finite()) {
            return Err(error("nonfinite point location arithmetic"));
        }
        let location = PointLocation {
            cell: id,
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

/// Evaluate a scalar P1 nodal field at a physical point. Search, tolerance and
/// extrapolation follow [`locate`]; on a shared face, either cell agrees in
/// exact arithmetic for a conforming field. Invalid dimensions, nonfinite
/// values and nonfinite computed results are rejected.
pub fn evaluate(
    mesh: &TetMesh,
    values: &[f64],
    point: Point3,
    tolerance: f64,
) -> Result<Option<f64>> {
    if values.len() != mesh.vertices().len() || values.iter().any(|v| !v.is_finite()) {
        return Err(error("invalid P1 tetrahedral field"));
    }
    let location = match locate(mesh, point, tolerance)? {
        Some(location) => location,
        None => return Ok(None),
    };
    let cell = mesh.cells()[location.cell.0];
    let value = (0..4)
        .map(|i| location.barycentric[i] * values[cell[i].0])
        .sum::<f64>();
    if !value.is_finite() {
        return Err(error("nonfinite point evaluation"));
    }
    Ok(Some(value))
}
