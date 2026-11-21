//! Scalar P1 triangle and affine mapping.
use crate::base::{triangle_determinant, Point2, UpdateFlags};
use crate::{error, Result};

/// Reference shape values in vertex order (0,0), (1,0), (0,1).
pub fn shape_values(p: Point2) -> [f64; 3] {
    [1.0 - p.x - p.y, p.x, p.y]
}

/// Affine map from the reference triangle to a counterclockwise physical triangle.
#[derive(Debug, Clone, Copy)]
pub struct AffineTriangle {
    points: [Point2; 3],
    det: f64,
}
impl AffineTriangle {
    pub fn new(points: [Point2; 3]) -> Result<Self> {
        let det = triangle_determinant(points)?;
        Ok(Self { points, det })
    }
    pub fn determinant(&self) -> f64 {
        self.det
    }
    pub fn map(&self, p: Point2) -> Point2 {
        let [a, b, c] = self.points;
        Point2 {
            x: a.x + (b.x - a.x) * p.x + (c.x - a.x) * p.y,
            y: a.y + (b.y - a.y) * p.x + (c.y - a.y) * p.y,
        }
    }
    /// Gradients transformed by the inverse transpose of the Jacobian.
    pub fn gradients(&self) -> [[f64; 2]; 3] {
        let [a, b, c] = self.points;
        let j00 = b.x - a.x;
        let j10 = b.y - a.y;
        let j01 = c.x - a.x;
        let j11 = c.y - a.y;
        let refs = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];
        let mut result = [[0.0; 2]; 3];
        for i in 0..3 {
            result[i] = [
                (j11 * refs[i][0] - j10 * refs[i][1]) / self.det,
                (-j01 * refs[i][0] + j00 * refs[i][1]) / self.det,
            ];
        }
        result
    }
}

/// Number of shape functions of the P1 triangle.
const P1_TRIANGLE_N_SHAPES: usize = 3;

/// Constant reference gradients of the P1 basis functions `1 - x - y`, `x`, `y`.
const P1_REFERENCE_GRADIENTS: [[f64; 2]; P1_TRIANGLE_N_SHAPES] =
    [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]];

/// FEValues-style evaluation of the P1 triangle on an affine cell.
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
/// [`value`]: FeValuesP1Triangle::value
/// [`grad_ref_at`]: FeValuesP1Triangle::grad_ref_at
/// [`grad_phys_at`]: FeValuesP1Triangle::grad_phys_at
/// [`jxw_at`]: FeValuesP1Triangle::jxw_at
#[derive(Debug, Clone)]
pub struct FeValuesP1Triangle {
    n_quadrature_points: usize,
    flags: UpdateFlags,
    quadrature_points: Vec<Point2>,
    values: Option<Vec<f64>>,
    grad_ref: Option<Vec<[f64; 2]>>,
    grad_phys: Option<Vec<[f64; 2]>>,
    jxw: Option<Vec<f64>>,
}

impl FeValuesP1Triangle {
    /// Evaluate the P1 triangle on `cell` at the given reference points.
    ///
    /// `reference_points` and `weights` must have the same non-zero length. All
    /// reference coordinates and all weights must be finite, and every weight
    /// must be strictly positive. `cell` must be a counterclockwise,
    /// non-degenerate affine triangle, as enforced by
    /// [`AffineTriangle::new`].
    ///
    /// Input validation always runs, including when `flags` is empty.
    pub fn new(
        cell: &AffineTriangle,
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
        let determinant = cell.determinant();
        if !determinant.is_finite() || determinant <= 0.0 {
            return Err(error(
                "cell must be a counterclockwise nondegenerate triangle",
            ));
        }

        let n = reference_points.len();
        let quadrature_points = reference_points
            .iter()
            .map(|point| cell.map(*point))
            .collect();

        let values = if flags.values {
            let mut data = vec![0.0; P1_TRIANGLE_N_SHAPES * n];
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
            let mut data = vec![[0.0; 2]; P1_TRIANGLE_N_SHAPES * n];
            for (i, gradient) in P1_REFERENCE_GRADIENTS.iter().enumerate() {
                data[i * n..(i + 1) * n].fill(*gradient);
            }
            Some(data)
        } else {
            None
        };

        let grad_phys = if flags.gradients {
            let physical = cell.gradients();
            let mut data = vec![[0.0; 2]; P1_TRIANGLE_N_SHAPES * n];
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
        if shape >= P1_TRIANGLE_N_SHAPES || q >= self.n_quadrature_points {
            return None;
        }
        Some(shape * self.n_quadrature_points + q)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::quadrature::triangle_degree_5;

    fn p(x: f64, y: f64) -> Point2 {
        Point2::new(x, y).unwrap()
    }

    fn near(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{} != {} (tolerance {})", a, b, tol);
    }

    fn unit_cell() -> AffineTriangle {
        AffineTriangle::new([p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0)]).unwrap()
    }

    fn scaled_cell() -> AffineTriangle {
        AffineTriangle::new([p(1.0, 2.0), p(3.0, 2.0), p(1.0, 5.0)]).unwrap()
    }

    fn degree_five_rule() -> (Vec<Point2>, Vec<f64>) {
        let rule = triangle_degree_5();
        (
            rule.iter().map(|q| q.point).collect(),
            rule.iter().map(|q| q.weight).collect(),
        )
    }

    /// ADR-0003 invariant 1.
    #[test]
    fn values_satisfy_partition_of_unity() {
        let (points, weights) = degree_five_rule();
        for cell in [unit_cell(), scaled_cell()].iter() {
            let fe = FeValuesP1Triangle::new(cell, &points, &weights, UpdateFlags::ALL).unwrap();
            assert_eq!(fe.n_quadrature_points(), points.len());
            assert_eq!(fe.values().unwrap().len(), 3 * points.len());
            for q in 0..fe.n_quadrature_points() {
                let sum: f64 = (0..3).map(|i| fe.value(i, q).unwrap()).sum();
                near(sum, 1.0, 1e-12);
            }
        }
    }

    /// ADR-0003 invariant 2.
    #[test]
    fn values_are_nodal_at_reference_vertices() {
        let vertices = [p(0.0, 0.0), p(1.0, 0.0), p(0.0, 1.0)];
        let weights = [1.0 / 6.0; 3];
        let flags = UpdateFlags {
            values: true,
            gradients: false,
            jxw: false,
        };
        let fe = FeValuesP1Triangle::new(&unit_cell(), &vertices, &weights, flags).unwrap();
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
    fn reference_gradients_are_constant_and_sum_to_zero() {
        let (points, weights) = degree_five_rule();
        let flags = UpdateFlags {
            values: false,
            gradients: true,
            jxw: false,
        };
        let fe = FeValuesP1Triangle::new(&scaled_cell(), &points, &weights, flags).unwrap();
        assert!(fe.values().is_none());
        for q in 0..fe.n_quadrature_points() {
            let mut sum = [0.0, 0.0];
            for i in 0..3 {
                let gradient = fe.grad_ref_at(i, q).unwrap();
                let expected = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]][i];
                near(gradient[0], expected[0], 1e-15);
                near(gradient[1], expected[1], 1e-15);
                sum[0] += gradient[0];
                sum[1] += gradient[1];
            }
            near(sum[0], 0.0, 1e-15);
            near(sum[1], 0.0, 1e-15);
        }
    }

    /// ADR-0003 invariants 4, 5 and 6, plus the non-gated physical points.
    #[test]
    fn jxw_and_physical_gradients_match_the_affine_mapping() {
        let (points, weights) = degree_five_rule();
        let cell = scaled_cell();
        let fe = FeValuesP1Triangle::new(&cell, &points, &weights, UpdateFlags::ALL).unwrap();

        let determinant = cell.determinant();
        near(determinant, 6.0, 1e-14);

        // Constant Jacobian columns and the resulting J^{-T} = 1/det * [[j11, -j10], [-j01, j00]].
        let a = p(1.0, 2.0);
        let b = p(3.0, 2.0);
        let c = p(1.0, 5.0);
        let (j00, j10) = (b.x - a.x, b.y - a.y);
        let (j01, j11) = (c.x - a.x, c.y - a.y);

        for (q, weight) in weights.iter().enumerate() {
            let jxw = fe.jxw_at(q).unwrap();
            assert!(jxw > 0.0, "JxW must be positive, got {}", jxw);
            near(jxw, determinant.abs() * weight, 1e-13);
            near(fe.jxw().unwrap()[q], jxw, 1e-15);

            for i in 0..3 {
                let reference = fe.grad_ref_at(i, q).unwrap();
                let expected = [
                    (j11 * reference[0] - j10 * reference[1]) / determinant,
                    (-j01 * reference[0] + j00 * reference[1]) / determinant,
                ];
                let physical = fe.grad_phys_at(i, q).unwrap();
                near(physical[0], expected[0], 1e-13);
                near(physical[1], expected[1], 1e-13);
            }
        }

        for (q, point) in points.iter().enumerate() {
            assert_eq!(fe.quadrature_points()[q], cell.map(*point));
        }
    }

    /// ADR-0003 invariant 7.
    #[test]
    fn flags_gate_computation() {
        let (points, weights) = degree_five_rule();
        let cell = unit_cell();

        let none = FeValuesP1Triangle::new(&cell, &points, &weights, UpdateFlags::EMPTY).unwrap();
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

        let only_values = FeValuesP1Triangle::new(
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

        let only_gradients = FeValuesP1Triangle::new(
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

        let only_jxw = FeValuesP1Triangle::new(
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

        let all = FeValuesP1Triangle::new(&cell, &points, &weights, UpdateFlags::all()).unwrap();
        assert!(all.values().is_some());
        assert!(all.grad_ref().is_some());
        assert!(all.grad_phys().is_some());
        assert!(all.jxw().is_some());

        assert!(all.value(3, 0).is_none());
        assert!(all.value(0, all.n_quadrature_points()).is_none());
        assert!(all.grad_ref_at(0, usize::MAX).is_none());
        assert!(all.grad_phys_at(usize::MAX, 0).is_none());
        assert!(all.jxw_at(usize::MAX).is_none());
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let cell = unit_cell();
        let points = [p(0.25, 0.25)];
        let weights = [0.5];

        // A reversed or degenerate cell never reaches the evaluator.
        assert!(AffineTriangle::new([p(0.0, 0.0), p(0.0, 1.0), p(1.0, 0.0)]).is_err());
        assert!(AffineTriangle::new([p(0.0, 0.0), p(1.0, 1.0), p(2.0, 2.0)]).is_err());
        assert!(AffineTriangle::new([
            p(0.0, 0.0),
            Point2 {
                x: f64::NAN,
                y: 0.0
            },
            p(0.0, 1.0)
        ])
        .is_err());

        // Count mismatches and empty rules.
        assert!(FeValuesP1Triangle::new(&cell, &points, &[], UpdateFlags::ALL).is_err());
        assert!(FeValuesP1Triangle::new(
            &cell,
            &[p(0.2, 0.2), p(0.3, 0.3)],
            &weights,
            UpdateFlags::ALL
        )
        .is_err());
        assert!(FeValuesP1Triangle::new(&cell, &[], &[], UpdateFlags::ALL).is_err());

        // Invalid weights.
        for weight in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY].iter() {
            assert!(
                FeValuesP1Triangle::new(&cell, &points, &[*weight], UpdateFlags::ALL).is_err(),
                "weight {} must be rejected",
                weight
            );
            assert!(
                FeValuesP1Triangle::new(&cell, &points, &[*weight], UpdateFlags::EMPTY).is_err(),
                "weight {} must be rejected even with empty flags",
                weight
            );
        }

        // Non-finite reference coordinates.
        let bad_points = [
            Point2 {
                x: f64::NAN,
                y: 0.0,
            },
            Point2 {
                x: 0.0,
                y: f64::INFINITY,
            },
        ];
        for point in bad_points.iter() {
            assert!(FeValuesP1Triangle::new(&cell, &[*point], &weights, UpdateFlags::ALL).is_err());
            assert!(
                FeValuesP1Triangle::new(&cell, &[*point], &weights, UpdateFlags::EMPTY).is_err()
            );
        }
    }
}
