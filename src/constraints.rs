//! Algebraic affine constraints and reduced sparse systems.
//!
//! This module implements `x = C y + g` and the reduced equations
//! `Cᵀ A C y = Cᵀ(b − A g)`. It does not discover hanging nodes or change
//! mesh topology; the caller supplies the complete set of relations.
//!
//! Frozen contract: constraints are expressed as `x = C y + g`; [`AffineConstraints::close`]
//! develops dependency chains iteratively (without recursion) and detects cycles;
//! [`ClosedConstraints::condense`] assembles the reduced system and
//! [`ClosedConstraints::reconstruct`] recovers `x`. A fully constrained system is
//! 0×0 and reconstructs from an empty reduced vector.
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
///
/// Frozen contract: each slave may have one rule; masters may themselves be
/// constrained. [`close`](Self::close) expands dependency chains iteratively
/// (without recursion) and detects cycles. Indices refer to the original
/// vector, independent of the finite element mesh.
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
    ///
    /// Frozen policy: dependency chains are developed iteratively (no
    /// recursion) using topological order; any cycle, including a nonzero
    /// self-dependency, is reported as an error. A zero-coefficient
    /// self-dependency is ignored.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lac::{conjugate_gradient, CgOptions, CooBuilder, CsrMatrix};

    fn near(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-9, "{} != {}", a, b);
    }

    fn identity(size: usize) -> CsrMatrix {
        let mut builder = CooBuilder::new(size, size);
        for i in 0..size {
            builder.add(i, i, 1.0).unwrap();
        }
        builder.into_csr().unwrap()
    }

    #[test]
    fn close_develops_deep_dependency_chains_without_recursion() {
        let n = 5_000;
        let mut constraints = AffineConstraints::new(n);
        for i in 1..n {
            constraints.add(i, 1.0, &[(i - 1, 1.0)]).unwrap();
        }
        let closed = constraints.close().unwrap();
        assert_eq!(closed.free_dofs(), &[0]);
        let values = closed.reconstruct(&[2.0]).unwrap();
        near(values[n - 1], n as f64 + 1.0);
    }

    #[test]
    fn close_detects_direct_and_indirect_cycles() {
        let mut constraints = AffineConstraints::new(3);
        constraints.add(0, 0.0, &[(1, 1.0)]).unwrap();
        constraints.add(1, 0.0, &[(2, 1.0)]).unwrap();
        constraints.add(2, 0.0, &[(0, 1.0)]).unwrap();
        assert!(constraints.close().is_err());

        let mut self_loop = AffineConstraints::new(1);
        self_loop.add(0, 0.0, &[(0, 1.0)]).unwrap();
        assert!(self_loop.close().is_err());
    }

    #[test]
    fn condense_and_reconstruct_round_trip_satisfies_affine_relation() {
        let mut constraints = AffineConstraints::new(5);
        constraints.add(2, 3.0, &[(1, 0.5)]).unwrap();
        constraints.add(1, 2.0, &[(0, 2.0)]).unwrap();
        constraints
            .add(4, -1.0, &[(2, 1.0), (1, -1.0), (3, 3.0)])
            .unwrap();
        let closed = constraints.close().unwrap();
        assert_eq!(closed.free_dofs(), &[0, 3]);

        let reduced = [4.0, 5.0];
        let original = closed.reconstruct(&reduced).unwrap();
        assert_eq!(original, vec![4.0, 10.0, 8.0, 5.0, 12.0]);

        let matrix = identity(5);
        let (reduced_matrix, reduced_rhs) = closed.condense(&matrix, &original).unwrap();
        let (solution, _) =
            conjugate_gradient(&reduced_matrix, &reduced_rhs, CgOptions::default()).unwrap();
        let reconstructed = closed.reconstruct(&solution).unwrap();
        for (value, expected) in reconstructed.iter().zip(&original) {
            near(*value, *expected);
        }
    }

    #[test]
    fn fully_constrained_system_reconstructs_from_empty_reduced_vector() {
        let mut constraints = AffineConstraints::new(3);
        constraints.add(0, 3.0, &[]).unwrap();
        constraints.add(1, 1.0, &[(0, 2.0)]).unwrap();
        constraints.add(2, -1.0, &[(1, 0.5)]).unwrap();
        let closed = constraints.close().unwrap();
        assert!(closed.free_dofs().is_empty());
        assert_eq!(closed.reconstruct(&[]).unwrap(), vec![3.0, 7.0, 2.5]);

        let matrix = identity(3);
        let (reduced_matrix, reduced_rhs) = closed.condense(&matrix, &[1.0, 2.0, 3.0]).unwrap();
        assert_eq!(reduced_matrix.rows(), 0);
        assert_eq!(reduced_matrix.cols(), 0);
        assert!(reduced_rhs.is_empty());
    }

    #[test]
    fn zero_dimensional_constraints_reconstruct_empty() {
        let closed = AffineConstraints::new(0).close().unwrap();
        assert!(closed.reconstruct(&[]).unwrap().is_empty());
        let matrix = CooBuilder::new(0, 0).into_csr().unwrap();
        assert!(closed.condense(&matrix, &[]).unwrap().1.is_empty());
    }
}
