//! Quadrature on the reference triangle, square, tetrahedron and unit interval.
use crate::base::Point2;

/// A reference point and its quadrature weight.
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

/// Tensor product of the three-point Gauss rule on `[0,1]²`.
/// Exact for polynomials of degree at most five in each coordinate.
pub fn square_degree_5() -> [QuadraturePoint; 9] {
    let q = line_degree_5();
    let mut result = [QuadraturePoint {
        point: Point2 { x: 0.0, y: 0.0 },
        weight: 0.0,
    }; 9];
    for i in 0..3 {
        for j in 0..3 {
            result[3 * i + j] = QuadraturePoint {
                point: Point2 {
                    x: q[i].point,
                    y: q[j].point,
                },
                weight: q[i].weight * q[j].weight,
            };
        }
    }
    result
}

/// Quadrature point on the reference tetrahedron.
#[derive(Clone, Copy, Debug)]
pub struct QuadraturePoint3 {
    pub point: crate::base::Point3,
    pub weight: f64,
}

/// A 4×4×4 Gauss-Duffy rule on `x,y,z >= 0`, `x+y+z <= 1`.
/// The weights sum to 1/6 and polynomials of total degree ≤ 5 are
/// integrated exactly in exact arithmetic. This rule favors simplicity over
/// minimal point count.
pub fn tetra_degree_5() -> Vec<QuadraturePoint3> {
    let a = ((3.0 + 2.0 * (6.0_f64 / 5.0).sqrt()) / 7.0).sqrt() / 2.0;
    let b = ((3.0 - 2.0 * (6.0_f64 / 5.0).sqrt()) / 7.0).sqrt() / 2.0;
    let outer = (18.0 - 30.0_f64.sqrt()) / 72.0;
    let inner = (18.0 + 30.0_f64.sqrt()) / 72.0;
    let q = [
        (0.5 - a, outer),
        (0.5 - b, inner),
        (0.5 + b, inner),
        (0.5 + a, outer),
    ];
    let mut result = Vec::with_capacity(64);
    for &(r, wr) in &q {
        for &(s, ws) in &q {
            for &(t, wt) in &q {
                result.push(QuadraturePoint3 {
                    point: crate::base::Point3 {
                        x: r,
                        y: (1.0 - r) * s,
                        z: (1.0 - r) * (1.0 - s) * t,
                    },
                    weight: wr * ws * wt * (1.0 - r) * (1.0 - r) * (1.0 - s),
                });
            }
        }
    }
    result
}
