//! Minimal serial COO/CSR algebra and conjugate gradients.
use crate::{error, Result};

/// Incremental sparse matrix builder. Duplicate entries are summed on compression.
pub struct CooBuilder {
    rows: usize,
    cols: usize,
    entries: Vec<(usize, usize, f64)>,
}
impl CooBuilder {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            entries: Vec::new(),
        }
    }
    pub fn add(&mut self, row: usize, col: usize, value: f64) -> Result<()> {
        if row >= self.rows || col >= self.cols {
            return Err(error("matrix index out of bounds"));
        }
        if !value.is_finite() {
            return Err(error("matrix entry must be finite"));
        }
        self.entries.push((row, col, value));
        Ok(())
    }
    pub fn into_csr(mut self) -> Result<CsrMatrix> {
        self.entries.sort_by_key(|e| (e.0, e.1));
        let mut merged: Vec<(usize, usize, f64)> = Vec::new();
        for (row, col, value) in self.entries {
            if let Some(last) = merged.last_mut() {
                if last.0 == row && last.1 == col {
                    last.2 += value;
                    if !last.2.is_finite() {
                        return Err(error("nonfinite sum of matrix entries"));
                    }
                    continue;
                }
            }
            merged.push((row, col, value));
        }
        let row_ptr_len = self
            .rows
            .checked_add(1)
            .ok_or_else(|| error("matrix size overflow"))?;
        let mut row_ptr = vec![0; row_ptr_len];
        let mut col_idx = Vec::with_capacity(merged.len());
        let mut values = Vec::with_capacity(merged.len());
        for (row, col, value) in merged {
            row_ptr[row + 1] += 1;
            col_idx.push(col);
            values.push(value);
        }
        for i in 0..self.rows {
            row_ptr[i + 1] += row_ptr[i];
        }
        Ok(CsrMatrix {
            rows: self.rows,
            cols: self.cols,
            row_ptr,
            col_idx,
            values,
        })
    }
}

/// Compressed sparse row matrix with sorted column indices.
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    rows: usize,
    cols: usize,
    row_ptr: Vec<usize>,
    col_idx: Vec<usize>,
    values: Vec<f64>,
}
impl CsrMatrix {
    pub fn rows(&self) -> usize {
        self.rows
    }
    pub fn cols(&self) -> usize {
        self.cols
    }
    /// Number of explicitly stored entries, including stored numerical zeros.
    pub fn stored_entries(&self) -> usize {
        self.values.len()
    }
    /// Logical bytes occupied by the CSR row pointers, column indices and
    /// values at their current lengths. Excludes vector capacity, headers,
    /// allocator overhead, matrix assembly buffers and process peak memory.
    pub fn logical_storage_bytes(&self) -> u128 {
        let index_bytes = std::mem::size_of::<usize>() as u128;
        let value_bytes = std::mem::size_of::<f64>() as u128;
        (self.row_ptr.len() as u128 + self.col_idx.len() as u128) * index_bytes
            + self.values.len() as u128 * value_bytes
    }
    pub fn get(&self, row: usize, col: usize) -> Result<f64> {
        if row >= self.rows || col >= self.cols {
            return Err(error("matrix index out of bounds"));
        }
        for k in self.row_ptr[row]..self.row_ptr[row + 1] {
            if self.col_idx[k] == col {
                return Ok(self.values[k]);
            }
        }
        Ok(0.0)
    }
    pub fn mul(&self, x: &[f64]) -> Result<Vec<f64>> {
        if x.len() != self.cols {
            return Err(error("matrix/vector dimension mismatch"));
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(error("vector must be finite"));
        }
        let mut result = vec![0.0; self.rows];
        for (row, output) in result.iter_mut().enumerate() {
            for k in self.row_ptr[row]..self.row_ptr[row + 1] {
                *output += self.values[k] * x[self.col_idx[k]];
            }
            if !output.is_finite() {
                return Err(error("nonfinite matrix/vector product"));
            }
        }
        Ok(result)
    }
    /// Symmetric elimination for Dirichlet values. Multiple constraints are
    /// applied simultaneously; duplicate indices are rejected. Neither input is
    /// modified on failure.
    pub fn with_dirichlet(
        &self,
        rhs: &[f64],
        constraints: &[(usize, f64)],
    ) -> Result<(Self, Vec<f64>)> {
        if self.rows != self.cols || rhs.len() != self.rows {
            return Err(error("Dirichlet dimension mismatch"));
        }
        if rhs.iter().any(|x| !x.is_finite()) {
            return Err(error("right hand side must be finite"));
        }
        let mut prescribed = vec![None; self.rows];
        for &(i, value) in constraints {
            if i >= self.rows || !value.is_finite() {
                return Err(error("invalid Dirichlet constraint"));
            }
            if prescribed[i].is_some() {
                return Err(error("duplicate Dirichlet constraint"));
            }
            prescribed[i] = Some(value);
        }
        let mut b = rhs.to_vec();
        let mut builder = CooBuilder::new(self.rows, self.cols);
        for row in 0..self.rows {
            if let Some(value) = prescribed[row] {
                builder.add(row, row, 1.0)?;
                b[row] = value;
            } else {
                for k in self.row_ptr[row]..self.row_ptr[row + 1] {
                    let col = self.col_idx[k];
                    let value = self.values[k];
                    if let Some(g) = prescribed[col] {
                        b[row] -= value * g;
                    } else {
                        builder.add(row, col, value)?;
                    }
                }
                if !b[row].is_finite() {
                    return Err(error("nonfinite constrained right hand side"));
                }
            }
        }
        Ok((builder.into_csr()?, b))
    }
}

/// Conjugate-gradient stopping options: ||r|| <= max(abs_tol, rel_tol ||b||).
#[derive(Debug, Clone, Copy)]
pub struct CgOptions {
    pub max_iterations: usize,
    pub relative_tolerance: f64,
    pub absolute_tolerance: f64,
}
impl Default for CgOptions {
    fn default() -> Self {
        Self {
            max_iterations: 10000,
            relative_tolerance: 1e-10,
            absolute_tolerance: 1e-12,
        }
    }
}
/// Successful solver diagnostics; residual is recalculated from b - Ax.
#[derive(Debug, Clone, Copy)]
pub struct CgReport {
    pub iterations: usize,
    pub residual_norm: f64,
}

/// Diagonal Jacobi preconditioner `M = diag(A)` for an SPD system.
///
/// The diagonal is copied at construction, so rebuild it after changing the
/// matrix or its Dirichlet constraints. `new` rejects rectangular matrices,
/// missing, zero or negative diagonal entries. A positive diagonal alone is
/// not proof that the matrix is symmetric positive definite. A preconditioner
/// can be reused with a same-sized system, but it will retain its old diagonal.
#[derive(Debug, Clone)]
pub struct JacobiPreconditioner {
    diagonal: Vec<f64>,
}

impl JacobiPreconditioner {
    /// Extract the diagonal of a square CSR matrix and require finite positive
    /// entries. This does not check symmetry or positive definiteness of `a`.
    pub fn new(a: &CsrMatrix) -> Result<Self> {
        if a.rows != a.cols {
            return Err(error("Jacobi matrix must be square"));
        }
        let mut diagonal = Vec::with_capacity(a.rows);
        for row in 0..a.rows {
            let value = a.get(row, row)?;
            if !value.is_finite() || value <= 0.0 {
                return Err(error("Jacobi diagonal must be finite and positive"));
            }
            diagonal.push(value);
        }
        Ok(Self { diagonal })
    }
    /// Apply `M⁻¹r`, validating dimensions and finite input/output. Dividing
    /// by the stored diagonal avoids explicitly forming a possibly overflowing
    /// reciprocal; unrepresentable results are errors.
    pub fn apply(&self, r: &[f64]) -> Result<Vec<f64>> {
        if r.len() != self.diagonal.len() {
            return Err(error("Jacobi dimension mismatch"));
        }
        if r.iter().any(|v| !v.is_finite()) {
            return Err(error("Jacobi input must be finite"));
        }
        let values: Vec<_> = r.iter().zip(&self.diagonal).map(|(v, d)| v / d).collect();
        if values.iter().any(|x| !x.is_finite()) {
            return Err(error("nonfinite Jacobi result"));
        }
        Ok(values)
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |acc, x| acc.hypot(*x))
}

/// Solve Ax=b from zero. A must be symmetric positive definite; breakdown is
/// detected but does not constitute a complete test of this precondition.
pub fn conjugate_gradient(
    a: &CsrMatrix,
    b: &[f64],
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    cg_impl(a, b, None, None, options)
}

/// Solve `Ax=b` starting from an explicit initial guess in the matrix's
/// vector ordering. Neither `initial` nor the other inputs are modified.
/// The stopping threshold remains `max(abs_tol, rel_tol * ||b||)`, not a
/// fraction of the initial residual. An already converged initial guess
/// returns zero iterations; even for `b=0`, a nonzero initial guess must be
/// checked against the actual residual. The reported residual is recomputed
/// from `b-Ax`. As for [`conjugate_gradient`], `A` must be SPD and the
/// implementation is unpreconditioned. Large or tiny values can make inner
/// products unrepresentable and result in a numerical breakdown error.
///
/// When using a Dirichlet-eliminated system, callers should initialize its
/// constrained entries to the prescribed values; no constraints are imposed
/// automatically. A good initial guess is not guaranteed to reduce iterations.
pub fn conjugate_gradient_with_initial_guess(
    a: &CsrMatrix,
    b: &[f64],
    initial: &[f64],
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    cg_impl(a, b, Some(initial), None, options)
}

/// Solve an SPD system from zero with a Jacobi preconditioner. The stopping
/// test and reported residual use the *unpreconditioned* norm `||b-Ax||` and
/// `max(abs_tol, rel_tol * ||b||)`. Jacobi must have the same dimension as the
/// system. Positive diagonal does not establish SPD; breakdown is detected but
/// is not a complete SPD test. Jacobi does not guarantee faster convergence.
pub fn conjugate_gradient_jacobi(
    a: &CsrMatrix,
    b: &[f64],
    jacobi: &JacobiPreconditioner,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    cg_impl(a, b, None, Some(jacobi), options)
}

/// Jacobi-preconditioned CG from a caller-supplied initial guess. Neither the
/// guess nor the preconditioner is modified. See [`conjugate_gradient_jacobi`]
/// and [`conjugate_gradient_with_initial_guess`] for the stopping criteria,
/// SPD assumptions, and Dirichlet constraints on the initial guess.
pub fn conjugate_gradient_jacobi_with_initial_guess(
    a: &CsrMatrix,
    b: &[f64],
    initial: &[f64],
    jacobi: &JacobiPreconditioner,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    cg_impl(a, b, Some(initial), Some(jacobi), options)
}

fn precondition(r: &[f64], jacobi: Option<&JacobiPreconditioner>) -> Result<(Vec<f64>, f64)> {
    let z = match jacobi {
        Some(m) => m.apply(r)?,
        None => r.to_vec(),
    };
    let rho = dot(r, &z);
    if !rho.is_finite() || rho <= 0.0 {
        return Err(error(
            "CG breakdown: unrepresentable or nonpositive preconditioned residual",
        ));
    }
    Ok((z, rho))
}

fn cg_impl(
    a: &CsrMatrix,
    b: &[f64],
    initial: Option<&[f64]>,
    jacobi: Option<&JacobiPreconditioner>,
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    if a.rows != a.cols
        || a.rows != b.len()
        || initial.map_or(false, |x| x.len() != b.len())
        || jacobi.map_or(false, |m| m.diagonal.len() != b.len())
    {
        return Err(error("CG dimension mismatch"));
    }
    if b.iter().any(|v| !v.is_finite())
        || initial.map_or(false, |x| x.iter().any(|v| !v.is_finite()))
        || !options.relative_tolerance.is_finite()
        || !options.absolute_tolerance.is_finite()
        || options.relative_tolerance < 0.0
        || options.absolute_tolerance < 0.0
    {
        return Err(error("invalid CG input or tolerance"));
    }
    let b_norm = norm(b);
    if !b_norm.is_finite() {
        return Err(error("nonfinite right hand side norm"));
    }
    let threshold = options
        .absolute_tolerance
        .max(options.relative_tolerance * b_norm);
    if !threshold.is_finite() {
        return Err(error("nonfinite CG tolerance"));
    }
    let mut x = initial.map_or_else(|| vec![0.0; b.len()], |x| x.to_vec());
    let mut r = match initial {
        Some(_) => b.iter().zip(a.mul(&x)?).map(|(bi, ai)| bi - ai).collect(),
        None => b.to_vec(),
    };
    let initial_norm = norm(&r);
    if !initial_norm.is_finite() {
        return Err(error("nonfinite initial CG residual"));
    }
    if initial_norm <= threshold {
        return Ok((
            x,
            CgReport {
                iterations: 0,
                residual_norm: initial_norm,
            },
        ));
    }
    let (mut p, mut rho) = precondition(&r, jacobi)?;
    for iteration in 1..=options.max_iterations {
        let ap = a.mul(&p)?;
        let curvature = dot(&p, &ap);
        if !curvature.is_finite() || curvature <= 0.0 {
            return Err(error("CG breakdown: nonpositive or nonfinite curvature"));
        }
        let alpha = rho / curvature;
        if !alpha.is_finite() {
            return Err(error("CG breakdown: nonfinite step"));
        }
        for i in 0..x.len() {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let next = dot(&r, &r);
        if !next.is_finite() {
            return Err(error("CG breakdown: nonfinite residual"));
        }
        if next.sqrt() <= threshold {
            let ax = a.mul(&x)?;
            let true_r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            let actual = norm(&true_r);
            if !actual.is_finite() {
                return Err(error("CG breakdown: nonfinite true residual"));
            }
            if actual <= threshold {
                return Ok((
                    x,
                    CgReport {
                        iterations: iteration,
                        residual_norm: actual,
                    },
                ));
            }
            // Restart if roundoff made the recurrence optimistic.
            r = true_r;
            (p, rho) = precondition(&r, jacobi)?;
            continue;
        }
        let (z, next_rho) = precondition(&r, jacobi)?;
        let beta = next_rho / rho;
        if !beta.is_finite() {
            return Err(error("CG breakdown: nonfinite direction coefficient"));
        }
        for i in 0..p.len() {
            p[i] = z[i] + beta * p[i];
        }
        rho = next_rho;
    }
    let ax = a.mul(&x)?;
    let residual = norm(
        &b.iter()
            .zip(ax.iter())
            .map(|(bi, ai)| bi - ai)
            .collect::<Vec<_>>(),
    );
    Err(crate::Error(format!(
        "CG did not converge after {} iterations (residual {})",
        options.max_iterations, residual
    )))
}
