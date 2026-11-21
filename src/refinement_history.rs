//! Persistent one-level transition history for conforming P1 triangle meshes.
use crate::grid::{CellId, LocalRefinement, Mesh2};
use crate::{error, Result};

/// A stack of conforming red-green refinements, retaining every mesh level.
///
/// `refine` appends one [`LocalRefinement`] after it succeeds. `coarsen_last`
/// removes the *entire most recent refinement level*, returning to its exact
/// stored parent mesh; it is not selective element coarsening. Mesh IDs and
/// vertex fields belong to their level and must not be reused at another
/// level, even for unchanged cells. After popping, later IDs can be reused for
/// unrelated cells. No metadata, solution vectors, hanging-node constraints,
/// mesh-quality guarantee, or coarsening error estimator is managed here.
#[derive(Debug, Clone)]
pub struct RefinementHistory {
    base: Mesh2,
    transitions: Vec<LocalRefinement>,
}
impl RefinementHistory {
    /// Start with a validated mesh. Geometric conformity and nonoverlap remain
    /// caller obligations, just as for `Mesh2` and [`LocalRefinement`].
    pub fn new(mesh: Mesh2) -> Self {
        Self {
            base: mesh,
            transitions: Vec::new(),
        }
    }
    /// Index of the current level; the original mesh is level zero.
    pub fn level(&self) -> usize {
        self.transitions.len()
    }
    /// Current mesh. IDs are meaningful only at this level.
    pub fn mesh(&self) -> &Mesh2 {
        self.transitions
            .last()
            .map_or(&self.base, LocalRefinement::mesh)
    }
    /// Mesh at a retained level, from zero through [`level`](Self::level).
    pub fn mesh_at(&self, level: usize) -> Result<&Mesh2> {
        if level == 0 {
            return Ok(&self.base);
        }
        self.transitions
            .get(level - 1)
            .map(LocalRefinement::mesh)
            .ok_or_else(|| error("refinement level out of bounds"))
    }
    /// Refinement from `parent_level` to `parent_level + 1`.
    pub fn transition(&self, parent_level: usize) -> Result<&LocalRefinement> {
        self.transitions
            .get(parent_level)
            .ok_or_else(|| error("refinement transition out of bounds"))
    }
    /// Refine the current level. Errors leave the history unchanged.
    /// Empty selections are rejected to avoid indistinguishable levels.
    pub fn refine(&mut self, marked: &[CellId]) -> Result<()> {
        if marked.is_empty() {
            return Err(error("cannot append an empty refinement"));
        }
        let refined = LocalRefinement::new(self.mesh(), marked)?;
        self.transitions.push(refined);
        Ok(())
    }
    /// Remove the entire latest refinement level. Errors at level zero leave
    /// the history unchanged. Retained old fields are caller-managed.
    pub fn coarsen_last(&mut self) -> Result<()> {
        if self.transitions.pop().is_none() {
            return Err(error("cannot coarsen the base level"));
        }
        Ok(())
    }
    /// Restrict a child-level P1 nodal vector to its parent by taking the
    /// values at retained parent vertices. This is nodal injection, **not** an
    /// L2 projection: detail in the child-only vertices is discarded. It
    /// exactly reverses prolongation for an unchanged P1 parent field in real
    /// arithmetic. Call before popping the relevant child level.
    pub fn restrict_injection(
        &self,
        parent_level: usize,
        child_values: &[f64],
    ) -> Result<Vec<f64>> {
        let parent = self.mesh_at(parent_level)?;
        let transition = self.transition(parent_level)?;
        if child_values.len() != transition.mesh().vertices().len()
            || child_values.iter().any(|v| !v.is_finite())
        {
            return Err(error("invalid child P1 field"));
        }
        Ok(child_values[..parent.vertices().len()].to_vec())
    }
}
