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

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solve Ax=b from zero. A must be symmetric positive definite; breakdown is
/// detected but does not constitute a complete test of this precondition.
pub fn conjugate_gradient(
    a: &CsrMatrix,
    b: &[f64],
    options: CgOptions,
) -> Result<(Vec<f64>, CgReport)> {
    if a.rows != a.cols || a.rows != b.len() {
        return Err(error("CG dimension mismatch"));
    }
    if b.iter().any(|v| !v.is_finite())
        || !options.relative_tolerance.is_finite()
        || !options.absolute_tolerance.is_finite()
        || options.relative_tolerance < 0.0
        || options.absolute_tolerance < 0.0
    {
        return Err(error("invalid CG input or tolerance"));
    }
    let b_norm = dot(b, b).sqrt();
    if !b_norm.is_finite() {
        return Err(error("nonfinite right hand side norm"));
    }
    let threshold = options
        .absolute_tolerance
        .max(options.relative_tolerance * b_norm);
    if !threshold.is_finite() {
        return Err(error("nonfinite CG tolerance"));
    }
    let mut x = vec![0.0; b.len()];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    if rr.sqrt() <= threshold {
        return Ok((
            x,
            CgReport {
                iterations: 0,
                residual_norm: rr.sqrt(),
            },
        ));
    }
    for iteration in 1..=options.max_iterations {
        let ap = a.mul(&p)?;
        let curvature = dot(&p, &ap);
        if !curvature.is_finite() || curvature <= 0.0 {
            return Err(error("CG breakdown: nonpositive or nonfinite curvature"));
        }
        let alpha = rr / curvature;
        if !alpha.is_finite() {
            return Err(error("CG breakdown: nonfinite step"));
        }
        for i in 0..x.len() {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let mut next = dot(&r, &r);
        if !next.is_finite() {
            return Err(error("CG breakdown: nonfinite residual"));
        }
        if next.sqrt() <= threshold {
            let ax = a.mul(&x)?;
            let true_r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
            let actual = dot(&true_r, &true_r).sqrt();
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
            next = dot(&r, &r);
            p = r.clone();
            rr = next;
            continue;
        }
        let beta = next / rr;
        for i in 0..p.len() {
            p[i] = r[i] + beta * p[i];
        }
        rr = next;
    }
    let ax = a.mul(&x)?;
    let residual = b
        .iter()
        .zip(ax.iter())
        .map(|(bi, ai)| (bi - ai) * (bi - ai))
        .sum::<f64>()
        .sqrt();
    Err(crate::Error(format!(
        "CG did not converge after {} iterations (residual {})",
        options.max_iterations, residual
    )))
}
