//! Algebraic affine constraints and reduced sparse systems.
//!
//! This module implements `x = C y + g` and the reduced equations
//! `Cᵀ A C y = Cᵀ(b − A g)`. It does not discover hanging nodes or change
//! mesh topology; the caller supplies the complete set of relations.
use crate::lac::{CooBuilder, CsrMatrix};
use crate::{error, Result};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone)]
struct Rule {
    offset: f64,
    terms: Vec<(usize, f64)>,
}
#[derive(Debug, Clone)]
struct Row {
    offset: f64,
    terms: Vec<(usize, f64)>,
}

/// Builder for rules `x[slave] = offset + Σ coefficient * x[master]`.
/// Each slave may have one rule. Masters may themselves be constrained; call
/// [`close`](Self::close) to expand dependencies and detect cycles. Indices
/// refer to the original vector, independent of the finite element mesh.
#[derive(Debug, Clone)]
pub struct AffineConstraints {
    rules: Vec<Option<Rule>>,
}
impl AffineConstraints {
    /// Create a builder for `size` original degrees of freedom.
    pub fn new(size: usize) -> Self {
        Self {
            rules: vec![None; size],
        }
    }
    /// Add a rule, or a fixed value when `terms` is empty. Rejects duplicate
    /// slaves or masters, out-of-range indices and nonfinite numbers. Zero
    /// coefficients are ignored after duplicate validation. Errors do not
    /// mutate the builder; nonzero self-dependencies are rejected at `close`.
    pub fn add(&mut self, slave: usize, offset: f64, terms: &[(usize, f64)]) -> Result<()> {
        if slave >= self.rules.len() || !offset.is_finite() {
            return Err(error("invalid constraint slave or offset"));
        }
        if self.rules[slave].is_some() {
            return Err(error("duplicate constraint slave"));
        }
        let mut seen = BTreeSet::new();
        for &(master, coefficient) in terms {
            if master >= self.rules.len() || !coefficient.is_finite() {
                return Err(error("invalid constraint master or coefficient"));
            }
            if !seen.insert(master) {
                return Err(error("duplicate constraint master"));
            }
        }
        self.rules[slave] = Some(Rule {
            offset,
            terms: terms.iter().copied().filter(|&(_, w)| w != 0.0).collect(),
        });
        Ok(())
    }
    /// Resolve all dependencies without recursion. Cycles and nonfinite
    /// arithmetic are errors. Free DoFs are ordered by increasing original
    /// index. Expressions can grow substantially when many rules are chained.
    pub fn close(self) -> Result<ClosedConstraints> {
        let n = self.rules.len();
        let free_dofs: Vec<_> = (0..n).filter(|&i| self.rules[i].is_none()).collect();
        let mut rows: Vec<Option<Row>> = vec![None; n];
        for (reduced, &original) in free_dofs.iter().enumerate() {
            rows[original] = Some(Row {
                offset: 0.0,
                terms: vec![(reduced, 1.0)],
            });
        }
        let mut dependents = vec![Vec::new(); n];
        let mut indegree = vec![0usize; n];
        for (slave, rule) in self.rules.iter().enumerate() {
            if let Some(rule) = rule {
                for &(master, _) in &rule.terms {
                    if self.rules[master].is_some() {
                        indegree[slave] += 1;
                        dependents[master].push(slave);
                    }
                }
            }
        }
        let mut ready = BTreeSet::new();
        for (i, &count) in indegree.iter().enumerate() {
            if self.rules[i].is_some() && count == 0 {
                ready.insert(i);
            }
        }
        let mut completed = 0;
        while let Some(&slave) = ready.iter().next() {
            ready.remove(&slave);
            let rule = self.rules[slave]
                .as_ref()
                .ok_or_else(|| error("invalid constraint state"))?;
            let mut offset = rule.offset;
            let mut combined = BTreeMap::<usize, f64>::new();
            for &(master, weight) in &rule.terms {
                let row = rows[master]
                    .as_ref()
                    .ok_or_else(|| error("invalid constraint dependency"))?;
                offset += weight * row.offset;
                if !offset.is_finite() {
                    return Err(error("nonfinite expanded constraint offset"));
                }
                for &(column, value) in &row.terms {
                    let entry = combined.entry(column).or_insert(0.0);
                    *entry += weight * value;
                    if !entry.is_finite() {
                        return Err(error("nonfinite expanded constraint coefficient"));
                    }
                }
            }
            rows[slave] = Some(Row {
                offset,
                terms: combined.into_iter().filter(|&(_, v)| v != 0.0).collect(),
            });
            completed += 1;
            for &next in &dependents[slave] {
                indegree[next] -= 1;
                if indegree[next] == 0 {
                    ready.insert(next);
                }
            }
        }
        if completed != n - free_dofs.len() {
            return Err(error("cyclic affine constraints"));
        }
        let rows = rows
            .into_iter()
            .map(|row| row.ok_or_else(|| error("incomplete affine constraints")))
            .collect::<Result<Vec<_>>>()?;
        Ok(ClosedConstraints { rows, free_dofs })
    }
}

/// Expanded, immutable transformation `x = C y + g`. Reduced coordinates are
/// the free original DoFs in the order exposed by [`free_dofs`](Self::free_dofs).
#[derive(Debug, Clone)]
pub struct ClosedConstraints {
    rows: Vec<Row>,
    free_dofs: Vec<usize>,
}
impl ClosedConstraints {
    /// Original vector dimension.
    pub fn size(&self) -> usize {
        self.rows.len()
    }
    /// Original indices corresponding to the reduced unknowns, in order.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }
    /// Expand reduced values into the original vector. Length, finiteness and
    /// computed values are checked. A fully constrained system takes `[]`.
    pub fn reconstruct(&self, reduced: &[f64]) -> Result<Vec<f64>> {
        if reduced.len() != self.free_dofs.len() || reduced.iter().any(|v| !v.is_finite()) {
            return Err(error("invalid reduced vector"));
        }
        let mut result = Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let mut value = row.offset;
            for &(j, w) in &row.terms {
                value += w * reduced[j];
            }
            if !value.is_finite() {
                return Err(error("nonfinite constraint reconstruction"));
            }
            result.push(value);
        }
        Ok(result)
    }
    /// Assemble `(Cᵀ A C, Cᵀ(b − A g))` without modifying the inputs. If all
    /// DoFs are constrained, returns an empty 0×0 system; reconstruct with an
    /// empty reduced vector rather than invoking a solver. Nonsymmetric input
    /// is allowed. For SPD `A`, full-column-rank `C` gives SPD reduced `A` in
    /// exact arithmetic; this routine does not verify SPD or conditioning.
    /// Products may cause substantial fill-in or overflow, and are accumulated
    /// in COO before compression. Detected nonfinite arithmetic is an error.
    pub fn condense(&self, matrix: &CsrMatrix, rhs: &[f64]) -> Result<(CsrMatrix, Vec<f64>)> {
        let n = self.rows.len();
        let m = self.free_dofs.len();
        if matrix.rows() != n
            || matrix.cols() != n
            || rhs.len() != n
            || rhs.iter().any(|v| !v.is_finite())
        {
            return Err(error(
                "invalid constraint system dimensions or right hand side",
            ));
        }
        let offsets: Vec<_> = self.rows.iter().map(|row| row.offset).collect();
        let ag = matrix.mul(&offsets)?;
        let mut reduced_rhs = vec![0.0; m];
        let mut builder = CooBuilder::new(m, m);
        for (i, row) in self.rows.iter().enumerate() {
            let residual = rhs[i] - ag[i];
            if !residual.is_finite() {
                return Err(error("nonfinite shifted right hand side"));
            }
            for &(k, weight) in &row.terms {
                reduced_rhs[k] += weight * residual;
                if !reduced_rhs[k].is_finite() {
                    return Err(error("nonfinite reduced right hand side"));
                }
            }
            for (j, a) in matrix.row_entries(i) {
                if a == 0.0 {
                    continue;
                }
                for &(k, left) in &row.terms {
                    for &(l, right) in &self.rows[j].terms {
                        builder.add(k, l, left * a * right)?;
                    }
                }
            }
        }
        Ok((builder.into_csr()?, reduced_rhs))
    }
}
