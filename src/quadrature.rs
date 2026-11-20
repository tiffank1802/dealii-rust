//! Quadrature on the reference triangle and unit interval.
use crate::base::Point2;

/// A point and its weight; weights integrate over the reference triangle of area 1/2.
#[derive(Clone, Copy, Debug)]
pub struct QuadraturePoint {
    pub point: Point2,
    pub weight: f64,
}

/// Seven-point symmetric rule, exact for polynomials through degree five.
pub fn triangle_degree_5() -> Vec<QuadraturePoint> {
    let a = 0.059_715_871_789_770;
    let b = 0.470_142_064_105_115;
    let c = 0.797_426_985_353_087;
    let d = 0.101_286_507_323_456;
    let w1 = 0.066_197_076_394_253;
    let w2 = 0.062_969_590_272_414;
    vec![
        QuadraturePoint {
            point: Point2 {
                x: 1.0 / 3.0,
                y: 1.0 / 3.0,
            },
            weight: 0.1125,
        },
        QuadraturePoint {
            point: Point2 { x: a, y: b },
            weight: w1,
        },
        QuadraturePoint {
            point: Point2 { x: b, y: a },
            weight: w1,
        },
        QuadraturePoint {
            point: Point2 { x: b, y: b },
            weight: w1,
        },
        QuadraturePoint {
            point: Point2 { x: c, y: d },
            weight: w2,
        },
        QuadraturePoint {
            point: Point2 { x: d, y: c },
            weight: w2,
        },
        QuadraturePoint {
            point: Point2 { x: d, y: d },
            weight: w2,
        },
    ]
}

/// A point and weight on the unit interval [0, 1].
#[derive(Clone, Copy, Debug)]
pub struct LineQuadraturePoint {
    pub point: f64,
    pub weight: f64,
}

/// Three-point Gauss rule, exact for polynomials through degree five.
pub fn line_degree_5() -> [LineQuadraturePoint; 3] {
    let offset = (3.0_f64 / 5.0).sqrt() / 2.0;
    [
        LineQuadraturePoint {
            point: 0.5 - offset,
            weight: 5.0 / 18.0,
        },
        LineQuadraturePoint {
            point: 0.5,
            weight: 4.0 / 9.0,
        },
        LineQuadraturePoint {
            point: 0.5 + offset,
            weight: 5.0 / 18.0,
        },
    ]
}
