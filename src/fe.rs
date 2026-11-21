//! Scalar P1 triangle and affine mapping.
use crate::base::{triangle_determinant, Point2};
use crate::Result;

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
