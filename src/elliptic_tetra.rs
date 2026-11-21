//! Scalar anisotropic diffusion-reaction on conforming P1 tetrahedra.
use crate::base::Point3;
use crate::grid::CellId;
use crate::lac::{CooBuilder, CsrMatrix};
use crate::quadrature::tetra_degree_5;
use crate::tetra::{self, AffineTet, TetMesh};
use crate::{error, Error, Result};

/// Symmetric 3D diffusion tensor in physical coordinates.
/// Its entries represent `[[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiffusionTensor3 {
    pub xx: f64,
    pub xy: f64,
    pub xz: f64,
    pub yy: f64,
    pub yz: f64,
    pub zz: f64,
}
impl DiffusionTensor3 {
    pub fn identity() -> Self {
        Self {
            xx: 1.0,
            xy: 0.0,
            xz: 0.0,
            yy: 1.0,
            yz: 0.0,
            zz: 1.0,
        }
    }
    fn validate(self) -> Result<()> {
        let entries = [self.xx, self.xy, self.xz, self.yy, self.yz, self.zz];
        let invalid = || error("diffusion tensor must be finite and positive definite");
        if entries.iter().any(|x| !x.is_finite())
            || self.xx <= 0.0
            || self.yy <= 0.0
            || self.zz <= 0.0
        {
            return Err(invalid());
        }
        // Positive definiteness is invariant under diagonal scaling. Work on
        // the correlation matrix to avoid overflowing diagonal products.
        // Divide by the larger square root first to avoid intermediate overflow.
        let normalized = |v: f64, a: f64, b: f64| (v / a.max(b)) / a.min(b);
        let x = self.xx.sqrt();
        let y = self.yy.sqrt();
        let z = self.zz.sqrt();
        let a = normalized(self.xy, x, y);
        let b = normalized(self.xz, x, z);
        let c = normalized(self.yz, y, z);
        if !a.is_finite()
            || !b.is_finite()
            || !c.is_finite()
            || a.abs() >= 1.0
            || b.abs() >= 1.0
            || c.abs() >= 1.0
        {
            return Err(invalid());
        }
        // Cholesky pivots of [[1,a,b],[a,1,c],[b,c,1]]. Near-singular
        // tensors are subject to rounding; no condition number is guaranteed.
        let pivot = 1.0 - a * a;
        let off = (c - a * b) / pivot.sqrt();
        let last = 1.0 - b * b - off * off;
        if !last.is_finite() || last <= 0.0 {
            return Err(invalid());
        }
        Ok(())
    }
}

/// Local P1 volume matrix and load in the input vertex order.
#[derive(Debug, Clone)]
pub struct LocalSystem {
    pub matrix: [[f64; 4]; 4],
    pub rhs: [f64; 4],
}

/// Integrate `∫(∇φ_i·A∇φ_j + cφ_iφ_j)` and `∫fφ_i` on a positively
/// oriented affine tetrahedron using a degree-five Gauss-Duffy rule.
/// At quadrature points `A` must be finite, symmetric positive definite,
/// `c` finite and nonnegative, and `f` finite. Sampling does not establish
/// positivity between points, uniform ellipticity, coercivity or conditioning.
/// Polynomial coefficients of total degree ≤5 for A, ≤3 for c, and ≤4 for f
/// are integrated exactly up to rounding. Other inputs incur quadrature error.
/// No boundary conditions are included; nonfinite arithmetic is rejected.
pub fn diffusion_reaction_p1<D, C, F>(
    points: [Point3; 4],
    diffusion: D,
    reaction: C,
    source: F,
) -> Result<LocalSystem>
where
    D: Fn(Point3) -> DiffusionTensor3,
    C: Fn(Point3) -> f64,
    F: Fn(Point3) -> f64,
{
    let map = AffineTet::new(points)?;
    let gradients = map.gradients();
    let mut matrix = [[0.0; 4]; 4];
    let mut rhs = [0.0; 4];
    for q in tetra_degree_5() {
        let point = map.map(q.point)?;
        let a = diffusion(point);
        a.validate()?;
        let c = reaction(point);
        let f = source(point);
        if !c.is_finite() || c < 0.0 {
            return Err(error("reaction must be finite and nonnegative"));
        }
        if !f.is_finite() {
            return Err(error("source function returned a nonfinite value"));
        }
        let weight = map.determinant() * q.weight;
        if !weight.is_finite() || weight <= 0.0 {
            return Err(error("invalid quadrature weight"));
        }
        let shape = tetra::shape_values(q.point);
        for i in 0..4 {
            rhs[i] += f * shape[i] * weight;
            for j in i..4 {
                let [x, y, z] = gradients[j];
                let ag = [
                    a.xx * x + a.xy * y + a.xz * z,
                    a.xy * x + a.yy * y + a.yz * z,
                    a.xz * x + a.yz * y + a.zz * z,
                ];
                let term =
                    gradients[i][0] * ag[0] + gradients[i][1] * ag[1] + gradients[i][2] * ag[2];
                matrix[i][j] += (term + c * shape[i] * shape[j]) * weight;
            }
        }
    }
    for i in 0..4 {
        if !rhs[i].is_finite() {
            return Err(error("nonfinite local right hand side"));
        }
        for j in i..4 {
            if !matrix[i][j].is_finite() {
                return Err(error("nonfinite local matrix entry"));
            }
            matrix[j][i] = matrix[i][j];
        }
    }
    Ok(LocalSystem { matrix, rhs })
}

/// Assemble `-div(A∇u)+cu=f` volume terms on a conforming tetrahedral mesh.
/// Callbacks receive the mesh-local cell ID and the physical quadrature point.
/// No boundary conditions are imposed. For mixed boundaries add the load from
/// `tetra::assemble_neumann` using the outward **conormal** flux `n·A∇u`,
/// before applying Dirichlet constraints. See [`diffusion_reaction_p1`] for
/// sample validation and quadrature limitations.
///
/// The matrix is not automatically positive definite. Pure diffusion with
/// Neumann everywhere has a constant nullspace per connected component; use
/// [`crate::neumann3d::solve_diffusion`] with compatible data for zero reaction.
/// Positive reaction can remove the nullspace, but no global coercivity or
/// condition-number guarantee is made here.
pub fn assemble<D, C, F>(
    mesh: &TetMesh,
    diffusion: D,
    reaction: C,
    source: F,
) -> Result<(CsrMatrix, Vec<f64>)>
where
    D: Fn(CellId, Point3) -> DiffusionTensor3,
    C: Fn(CellId, Point3) -> f64,
    F: Fn(CellId, Point3) -> f64,
{
    let n = mesh.vertices().len();
    let mut matrix = CooBuilder::new(n, n);
    let mut rhs = vec![0.0; n];
    for (i, cell) in mesh.cells().iter().enumerate() {
        let id = CellId(i);
        let local = diffusion_reaction_p1(
            mesh.cell_points(id)?,
            |p| diffusion(id, p),
            |p| reaction(id, p),
            |p| source(id, p),
        )
        .map_err(|e| Error(format!("cell {}: {}", i, e)))?;
        for r in 0..4 {
            rhs[cell[r].0] += local.rhs[r];
            if !rhs[cell[r].0].is_finite() {
                return Err(error("nonfinite assembled right hand side"));
            }
            for c in 0..4 {
                matrix.add(cell[r].0, cell[c].0, local.matrix[r][c])?;
            }
        }
    }
    Ok((matrix.into_csr()?, rhs))
}
