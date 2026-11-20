//! Two-dimensional geometry.

use crate::{error, Result};

/// Point in physical or reference coordinates, depending on the caller.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    /// Construct a point with finite coordinates.
    pub fn new(x: f64, y: f64) -> Result<Self> {
        if !x.is_finite() || !y.is_finite() {
            return Err(error("point coordinates must be finite"));
        }
        Ok(Self { x, y })
    }
}

/// Signed twice-area; positive means counterclockwise orientation.
pub fn signed_double_area(a: Point2, b: Point2, c: Point2) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

/// Validate a triangle and return its positive determinant.
///
/// Rejects clockwise and nearly degenerate triangles: the determinant must
/// exceed 64 machine epsilons times the largest squared edge length.
pub fn triangle_determinant(points: [Point2; 3]) -> Result<f64> {
    if points.iter().any(|p| !p.x.is_finite() || !p.y.is_finite()) {
        return Err(error("triangle coordinates must be finite"));
    }
    let [a, b, c] = points;
    let det = signed_double_area(a, b, c);
    let edge2 = |p: Point2, q: Point2| (p.x - q.x).powi(2) + (p.y - q.y).powi(2);
    let scale = edge2(a, b).max(edge2(b, c)).max(edge2(c, a));
    if !det.is_finite() || !scale.is_finite() || det <= 64.0 * f64::EPSILON * scale {
        return Err(error("triangle must be counterclockwise and nondegenerate"));
    }
    Ok(det)
}
