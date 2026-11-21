//! Two- and three-dimensional geometry.

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

/// A point in three-dimensional physical or reference coordinates.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Point3 {
    /// Construct a point with finite coordinates.
    pub fn new(x: f64, y: f64, z: f64) -> Result<Self> {
        if !x.is_finite() || !y.is_finite() || !z.is_finite() {
            return Err(error("point coordinates must be finite"));
        }
        Ok(Self { x, y, z })
    }
}

/// Reference vertices of the reference triangle, in vertex order.
static TRIANGLE_VERTICES: [Point2; 3] = [
    Point2 { x: 0.0, y: 0.0 },
    Point2 { x: 1.0, y: 0.0 },
    Point2 { x: 0.0, y: 1.0 },
];

/// Reference vertices of the reference quadrilateral, in vertex order.
static QUADRILATERAL_VERTICES: [Point2; 4] = [
    Point2 { x: 0.0, y: 0.0 },
    Point2 { x: 1.0, y: 0.0 },
    Point2 { x: 1.0, y: 1.0 },
    Point2 { x: 0.0, y: 1.0 },
];

/// Reference vertices of the reference tetrahedron, in vertex order.
static TETRAHEDRON_VERTICES: [Point3; 4] = [
    Point3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    Point3 {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    },
    Point3 {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    },
    Point3 {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    },
];

/// Local vertex indices of each triangle edge, as `[[0,1],[1,2],[2,0]]`.
static TRIANGLE_FACES: [[usize; 2]; 3] = [[0, 1], [1, 2], [2, 0]];

/// Local vertex indices of each quadrilateral edge, as `[[0,1],[1,2],[2,3],[3,0]]`.
static QUADRILATERAL_FACES: [[usize; 2]; 4] = [[0, 1], [1, 2], [2, 3], [3, 0]];

/// Local vertex indices of each tetrahedron face; face `i` is opposite vertex `i`.
static TETRAHEDRON_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

/// The topological shape of a reference cell.
///
/// The variants describe only the reference topology and its fixed vertex
/// order. Orientation and validity of a *physical* cell are the concern of the
/// grid seam; this type makes no statement about them.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum ReferenceCell {
    /// Reference triangle with vertices `(0,0)`, `(1,0)`, `(0,1)`.
    Triangle,
    /// Reference quadrilateral with vertices `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`.
    Quadrilateral,
    /// Reference tetrahedron with vertices `(0,0,0)`, `(1,0,0)`, `(0,1,0)`, `(0,0,1)`.
    Tetrahedron,
}

impl ReferenceCell {
    /// Spatial dimension: 2 for `Triangle` and `Quadrilateral`, 3 for `Tetrahedron`.
    pub fn dim(&self) -> usize {
        match self {
            ReferenceCell::Triangle | ReferenceCell::Quadrilateral => 2,
            ReferenceCell::Tetrahedron => 3,
        }
    }

    /// Number of vertices: 3, 4, 4.
    pub fn n_vertices(&self) -> usize {
        match self {
            ReferenceCell::Triangle => 3,
            ReferenceCell::Quadrilateral | ReferenceCell::Tetrahedron => 4,
        }
    }

    /// Number of faces (codimension-1 facets): 3 edges, 4 edges, 4 triangles.
    pub fn n_faces(&self) -> usize {
        match self {
            ReferenceCell::Triangle => 3,
            ReferenceCell::Quadrilateral => 4,
            ReferenceCell::Tetrahedron => 4,
        }
    }

    /// Local vertex indices forming the face `face_index`.
    ///
    /// Returns `None` when `face_index` is out of range, that is when it is not
    /// smaller than [`n_faces`](Self::n_faces). The returned slice has length 2
    /// for the 2D cells and length 3 for the tetrahedron.
    pub fn face_vertex_indices(&self, face_index: usize) -> Option<&'static [usize]> {
        match self {
            ReferenceCell::Triangle => TRIANGLE_FACES.get(face_index).map(|face| &face[..]),
            ReferenceCell::Quadrilateral => {
                QUADRILATERAL_FACES.get(face_index).map(|face| &face[..])
            }
            ReferenceCell::Tetrahedron => TETRAHEDRON_FACES.get(face_index).map(|face| &face[..]),
        }
    }

    /// Reference vertex coordinates for the 2D cells; `None` for `Tetrahedron`.
    pub fn reference_vertices_2d(&self) -> Option<&'static [Point2]> {
        match self {
            ReferenceCell::Triangle => Some(&TRIANGLE_VERTICES[..]),
            ReferenceCell::Quadrilateral => Some(&QUADRILATERAL_VERTICES[..]),
            ReferenceCell::Tetrahedron => None,
        }
    }

    /// Reference vertex coordinates for the 3D cell; `None` for the 2D cells.
    pub fn reference_vertices_3d(&self) -> Option<&'static [Point3]> {
        match self {
            ReferenceCell::Tetrahedron => Some(&TETRAHEDRON_VERTICES[..]),
            ReferenceCell::Triangle | ReferenceCell::Quadrilateral => None,
        }
    }
}

/// Quantities an FEValues evaluation should compute.
///
/// An evaluation computes only the requested quantities; when a flag is false
/// the corresponding output is `None`. Normal vectors are intentionally absent:
/// they belong to face evaluation in the grid seam.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct UpdateFlags {
    /// Compute shape function values.
    pub values: bool,
    /// Compute reference and physical shape function gradients.
    pub gradients: bool,
    /// Compute the Jacobian-times-weight factors.
    pub jxw: bool,
}

impl UpdateFlags {
    /// No quantity requested.
    pub const EMPTY: Self = Self {
        values: false,
        gradients: false,
        jxw: false,
    };

    /// Every quantity requested.
    pub const ALL: Self = Self {
        values: true,
        gradients: true,
        jxw: true,
    };

    /// No quantity requested; equivalent to [`UpdateFlags::EMPTY`].
    pub fn empty() -> Self {
        Self::EMPTY
    }

    /// Every quantity requested; equivalent to [`UpdateFlags::ALL`].
    pub fn all() -> Self {
        Self::ALL
    }

    /// Whether no quantity is requested.
    pub fn is_empty(&self) -> bool {
        !self.values && !self.gradients && !self.jxw
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-15, "{} != {}", a, b);
    }

    #[test]
    fn reference_cell_dimensions_and_counts() {
        assert_eq!(ReferenceCell::Triangle.dim(), 2);
        assert_eq!(ReferenceCell::Quadrilateral.dim(), 2);
        assert_eq!(ReferenceCell::Tetrahedron.dim(), 3);

        assert_eq!(ReferenceCell::Triangle.n_vertices(), 3);
        assert_eq!(ReferenceCell::Quadrilateral.n_vertices(), 4);
        assert_eq!(ReferenceCell::Tetrahedron.n_vertices(), 4);

        assert_eq!(ReferenceCell::Triangle.n_faces(), 3);
        assert_eq!(ReferenceCell::Quadrilateral.n_faces(), 4);
        assert_eq!(ReferenceCell::Tetrahedron.n_faces(), 4);
    }

    #[test]
    fn reference_cell_face_vertex_indices() {
        let triangle = [[0, 1], [1, 2], [2, 0]];
        for (face, expected) in triangle.iter().enumerate() {
            assert_eq!(
                ReferenceCell::Triangle.face_vertex_indices(face),
                Some(&expected[..])
            );
        }
        let quadrilateral = [[0, 1], [1, 2], [2, 3], [3, 0]];
        for (face, expected) in quadrilateral.iter().enumerate() {
            assert_eq!(
                ReferenceCell::Quadrilateral.face_vertex_indices(face),
                Some(&expected[..])
            );
        }
        let tetrahedron = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
        for (face, expected) in tetrahedron.iter().enumerate() {
            assert_eq!(
                ReferenceCell::Tetrahedron.face_vertex_indices(face),
                Some(&expected[..])
            );
        }

        for cell in [
            ReferenceCell::Triangle,
            ReferenceCell::Quadrilateral,
            ReferenceCell::Tetrahedron,
        ]
        .iter()
        {
            let expected_len = if cell.dim() == 3 { 3 } else { 2 };
            for face in 0..cell.n_faces() {
                let indices = cell.face_vertex_indices(face).unwrap();
                assert_eq!(indices.len(), expected_len);
                assert!(indices.iter().all(|&v| v < cell.n_vertices()));
            }
            assert_eq!(cell.face_vertex_indices(cell.n_faces()), None);
            assert_eq!(cell.face_vertex_indices(usize::MAX), None);
        }
    }

    #[test]
    fn reference_cell_vertex_coordinates() {
        let triangle = ReferenceCell::Triangle.reference_vertices_2d().unwrap();
        assert_eq!(triangle.len(), 3);
        for (vertex, expected) in triangle.iter().zip([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]) {
            near(vertex.x, expected[0]);
            near(vertex.y, expected[1]);
        }

        let quadrilateral = ReferenceCell::Quadrilateral
            .reference_vertices_2d()
            .unwrap();
        assert_eq!(quadrilateral.len(), 4);
        for (vertex, expected) in
            quadrilateral
                .iter()
                .zip([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
        {
            near(vertex.x, expected[0]);
            near(vertex.y, expected[1]);
        }

        let tetrahedron = ReferenceCell::Tetrahedron.reference_vertices_3d().unwrap();
        assert_eq!(tetrahedron.len(), 4);
        for (vertex, expected) in tetrahedron.iter().zip([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]) {
            near(vertex.x, expected[0]);
            near(vertex.y, expected[1]);
            near(vertex.z, expected[2]);
        }
    }

    #[test]
    fn reference_cell_dimension_mismatch_returns_none() {
        assert!(ReferenceCell::Triangle.reference_vertices_3d().is_none());
        assert!(ReferenceCell::Quadrilateral
            .reference_vertices_3d()
            .is_none());
        assert!(ReferenceCell::Tetrahedron.reference_vertices_2d().is_none());
    }

    #[test]
    fn update_flags_constructors_and_emptiness() {
        assert_eq!(UpdateFlags::empty(), UpdateFlags::EMPTY);
        assert_eq!(UpdateFlags::all(), UpdateFlags::ALL);
        assert!(UpdateFlags::EMPTY.is_empty());
        assert!(UpdateFlags::empty().is_empty());
        assert!(!UpdateFlags::ALL.is_empty());
        assert!(!UpdateFlags::all().is_empty());

        let empty = UpdateFlags::empty();
        assert!(!empty.values);
        assert!(!empty.gradients);
        assert!(!empty.jxw);
        let all = UpdateFlags::all();
        assert!(all.values);
        assert!(all.gradients);
        assert!(all.jxw);

        for flags in [
            UpdateFlags {
                values: true,
                gradients: false,
                jxw: false,
            },
            UpdateFlags {
                values: false,
                gradients: true,
                jxw: false,
            },
            UpdateFlags {
                values: false,
                gradients: false,
                jxw: true,
            },
        ]
        .iter()
        {
            assert!(!flags.is_empty());
        }
    }
}
