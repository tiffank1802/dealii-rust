//! Local P1 triangle integrators for a symmetric diffusion-reaction equation.
use crate::base::Point2;
use crate::fe::{shape_values, AffineTriangle};
use crate::quadrature::triangle_degree_5;
use crate::{error, Result};

/// Symmetric diffusion tensor `[[xx, xy], [xy, yy]]` in physical coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiffusionTensor {
    pub xx: f64,
    pub xy: f64,
    pub yy: f64,
}

impl DiffusionTensor {
    /// Identity tensor, recovering the Laplace operator.
    pub fn identity() -> Self {
        Self {
            xx: 1.0,
            xy: 0.0,
            yy: 1.0,
        }
    }

    pub(crate) fn validate(self) -> Result<()> {
        if !self.xx.is_finite()
            || !self.xy.is_finite()
            || !self.yy.is_finite()
            || self.xx <= 0.0
            || self.yy <= 0.0
        {
            return Err(error(
                "diffusion tensor must be finite and positive definite",
            ));
        }
        // Equivalent to xy² < xx*yy without forming products that can
        // overflow. In ill-conditioned cases this check is subject to rounding.
        let normalized_off_diagonal = (self.xy / self.xx.sqrt()) / self.yy.sqrt();
        if !normalized_off_diagonal.is_finite() || normalized_off_diagonal.abs() >= 1.0 {
            return Err(error(
                "diffusion tensor must be finite and positive definite",
            ));
        }
        Ok(())
    }
}

/// A local bilinear form and right-hand side in the input vertex order.
#[derive(Debug, Clone)]
pub struct LocalSystem {
    pub matrix: [[f64; 3]; 3],
    pub rhs: [f64; 3],
}

/// Integrate `∫ (∇φ_i · A ∇φ_j + c φ_i φ_j)` and `∫ f φ_i` on one
/// positively oriented P1 triangle. All three functions are sampled at the
/// seven degree-five quadrature points. `A` must be symmetric positive
/// definite and `c` nonnegative at those points; source values must be finite.
/// These sample checks do not establish positivity between quadrature points,
/// uniform ellipticity, or adequate conditioning. For polynomial coefficients,
/// the rule is exact up to total degree 5 for A, 3 for c, and 4 for f (apart
/// from floating-point error). Other inputs introduce quadrature error.
///
/// The result contains volume terms only; boundary terms are assembled
/// separately. Unsupported overflow/nonfinite local contributions are errors.
pub fn diffusion_reaction_p1<D, C, F>(
    points: [Point2; 3],
    diffusion: D,
    reaction: C,
    source: F,
) -> Result<LocalSystem>
where
    D: Fn(Point2) -> DiffusionTensor,
    C: Fn(Point2) -> f64,
    F: Fn(Point2) -> f64,
{
    let mapping = AffineTriangle::new(points)?;
    let gradients = mapping.gradients();
    if gradients.iter().flatten().any(|x| !x.is_finite()) {
        return Err(error("nonfinite basis gradient"));
    }
    let mut matrix = [[0.0; 3]; 3];
    let mut rhs = [0.0; 3];
    for q in triangle_degree_5() {
        let point = mapping.map(q.point);
        if !point.x.is_finite() || !point.y.is_finite() {
            return Err(error("nonfinite quadrature point"));
        }
        let tensor = diffusion(point);
        tensor.validate()?;
        let c = reaction(point);
        let f = source(point);
        if !c.is_finite() || c < 0.0 {
            return Err(error("reaction must be finite and nonnegative"));
        }
        if !f.is_finite() {
            return Err(error("source function returned a nonfinite value"));
        }
        let shape = shape_values(q.point);
        let weight = mapping.determinant() * q.weight;
        if !weight.is_finite() || weight <= 0.0 {
            return Err(error("invalid quadrature weight"));
        }
        for i in 0..3 {
            rhs[i] += f * shape[i] * weight;
            for j in i..3 {
                let [gx, gy] = gradients[j];
                let ag = [
                    tensor.xx * gx + tensor.xy * gy,
                    tensor.xy * gx + tensor.yy * gy,
                ];
                let diffusion_term = gradients[i][0] * ag[0] + gradients[i][1] * ag[1];
                matrix[i][j] += (diffusion_term + c * shape[i] * shape[j]) * weight;
            }
        }
    }
    for i in 0..3 {
        if !rhs[i].is_finite() {
            return Err(error("nonfinite local right hand side"));
        }
        for j in i..3 {
            if !matrix[i][j].is_finite() {
                return Err(error("nonfinite local matrix entry"));
            }
            matrix[j][i] = matrix[i][j];
        }
    }
    Ok(LocalSystem { matrix, rhs })
}
