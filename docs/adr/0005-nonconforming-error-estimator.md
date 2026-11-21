---
id: ADR-0005
status: accepted
created: 2026-11-21
decided: 2026-11-21
supersedes: []
superseded_by: null
rfc: RFC-0002
ai-usage: ai-assisted
tags: [adaptivity, error-estimator, nonconforming, poisson]
---

# ADR-0005: Residual error estimator for nonconforming P1 triangles

**Status:** Accepted
**Date:** 2026-11-21
**Deciders:** project maintainers

## Context

RFC-0001 / ADR-0004 delivered multi-level nonconforming P1 refinement with an
active-only `NonconformingMesh`, vertex-indexed DoFs and hanging constraints,
but no a posteriori estimator. The conforming estimator
`estimator::poisson_dirichlet` (`src/estimator.rs:30-102`) accepts only `&Mesh2`,
documents that nonconforming grids are unsupported (`:28-29`), and counts jumps
only for edges with exactly two incident cells (`src/estimator.rs:39,77-89`).
On a hanging interface the coarse edge and both fine half-edges each have one
incident cell, so all interface terms vanish silently. Driving Dörfler marking
(`src/estimator.rs:113-148`) with such values would under-refine exactly the
adapted region. The assembled load vector cannot supply the volume term: it
holds `∫f φᵢ`, not `∫f²` (`src/nonconforming_refinement.rs:1018-1044`).

Constraints: ADR-0001 P0-P2 parity is frozen; ADR-0003 §3–§4 require additive
changes and index-based, level-scoped IDs; ADR-0002 forbids verbatim upstream
reuse; `docs/contracts/seams.md:32,35-43` assigns error estimation to the
`assembly` seam with one owner per seam; MSRV 1.70, no `unsafe`, no new deps.

## Decision

Accept [RFC-0002](../rfcs/0002-nonconforming-error-estimator/README.md) with two
interdependent sub-decisions.

### D1 — Additive API in the assembly seam

One new public function in `src/estimator.rs`:

```rust
pub fn poisson_dirichlet_nonconforming<F>(
    mesh: &NonconformingMesh,
    solution: &[f64],
    source: F,
) -> Result<ResidualEstimate>
where
    F: Fn(Point2) -> f64;
```

It reuses `ResidualEstimate` and `mark_bulk` unchanged, with per-cell
indicators aligned to `mesh.cells` (active order). It calls `mesh.validate()?`,
requires the reconstructed solution (each hanging rule satisfied to
`1e-12 · (1 + max|u|)` over the slave and both masters), finite
values/source/gradients/lengths, and resamples `f` with the existing
degree-five quadrature. That tolerance scales with the largest magnitude in the
rule; it must not be tightened to an unscaled equality test.
`η_K² = h_K²∫_K f² + ½Σ_e h_e∫_e[∂ₙu]²` is a residual **indicator**, not a
certified bound. Private helpers handle the interface walk and jump
accumulation; no method is added to a grid type.

### D2 — Interface jumps use two-sided edges paired by `edge_splits`

For `η_K²`, interior edges with two incident cells keep the existing formula,
except that a key in `edge_splits` with two incident active cells returns `Err`
as a defensive guard: the `new` / `refine` / `coarsen` API cannot produce that
state, and rejecting it prevents a malformed mesh from skipping the interface
term silently. A key in `edge_splits` with one incident cell is a coarse
interface root: walk its sub-edges recursively with a visited set (the
`collect_subedges` pattern, `src/nonconforming_refinement.rs:22-55`) to leaf
segments, require exactly one
fine incident cell per leaf, and charge, for `jump = (∇u_c − ∇u_f)·n`:

- coarse cell: `0.5 · L_master · l_leaf · jump²`;
- fine cell: `0.5 · l_leaf² · jump²`.

An incidence mismatch is `Err`, never a skipped term. Split exterior
half-edges never become roots because their master key is not incident to an
active cell. Tangential traces of a constrained P1 solution agree at the
interface, but normal derivatives do not, so the jump is real and must be
counted. This works for unbalanced, multi-level chains because the root walk
reaches sub-midpoints that have no direct active-cell incidence.

## Options considered

- **A. Two-sided edges via `edge_splits` (chosen).** Counts the real interface
  jump; handles deep chains; stays inside the existing seams.
- **B. Skip interface jumps with a documented bias.** Smallest change, but
  silently under-estimates on the adapted patch and corrupts Dörfler marking.
  Rejected as the primary failure mode.
- **C. Project via `AffineConstraints`.** `close()`/`reconstruct()` only expose
  `x = Cy + g` (`src/constraints.rs:78-150,171-187`), not an edge projection;
  unwrapping deep chains would add a seam surface and amplify row-fill risk.
  Rejected.
- **D. Method on `NonconformingMesh`.** Convenient but creates a grid→assembly
  dependency and violates one-owner-per-seam. Rejected.

## Rationale

Option A is the only one that preserves the meaning of the Dörfler marker:
marking depends on *relative* cell indicators, so dropping the largest
interface contribution systematically biases refinement away from hanging
interfaces. Pairing sides through `edge_splits` reuses the single source of
truth already used by `discover_hanging` (`src/nonconforming_refinement.rs:62-76`)
and needs no new ID, state or geometric search. The per-side `h_e` scaling is
the direct generalization of `½ h_e ∫_e jump²` when the two sides have
different resolution. The implementation cost is two private helpers and a
differential test, which is proportionate to the correctness risk.

## Consequences

### Positive

- Interface contributions are counted; Dörfler marking acts on the refined
  patch instead of avoiding it.
- No P0-P2 numerical behaviour changes; the frozen conforming estimator stays
  untouched and a test asserts equality on uniformly refined meshes.
- No seam, ID, dependency or quadrature-table changes; ADR-0002 provenance is
  preserved (formula derived from first principles, no upstream constants).

### Negative

- The source closure must be supplied again; the assembled load vector is not
  reusable for `∫f²`.
- The caller must pass a reconstructed solution; the tolerance check turns a
  bad solve into `Err` rather than a plausible wrong number.
- Scope is narrow: homogeneous Dirichlet, scalar unit diffusion, P1 triangles,
  indicator only.

### Risks

- Silent under-estimation is mitigated by mandatory one-fine-owner incidence,
  `Err` on mismatch, and `hanging_interface_jump_is_counted`.
- Solver-dependence in tests is avoided with crafted constraint-satisfying
  vectors; the rate test uses `x²+y²` with a `0.25–0.65` ratio band.
- Constraint row-fill is not introduced: the estimator reads `mesh.hanging`
  only; the iterative leaf walk has a visited set.

## Scope

**IN:** additive `poisson_dirichlet_nonconforming` in `src/estimator.rs`;
interface jump policy; reuse of `ResidualEstimate`/`mark_bulk`; example update;
`tests/nonconforming_estimator.rs` with nine contract tests including a
manufactured-solution rate test and a hanging-interface differential case.

**OUT:** global/hierarchical estimators, variable diffusion/reaction, Neumann
residuals, ZZ/goal-oriented estimators, coarsening-error estimators,
anisotropy, higher order, quads/tetra/3D, distributed/MPI, matrix-free, new
quadrature, persistent IDs, ADR-0003 seam rework, changes to RFC-0001 or
ADR-0004.

## References

- [RFC-0002: Residual a posteriori error estimator for nonconforming P1 triangles](../rfcs/0002-nonconforming-error-estimator/README.md)
- [RFC-0001: Multi-level non-conforming refinement](../rfcs/0001-multi-level-nonconforming/README.md)
- [ADR-0004: Multi-level non-conforming refinement and selective derefinement](0004-multi-level-nonconforming.md)
- [ADR-0003: Seam contracts](0003-seam-contracts.md) — §3 additive, §4 index-based IDs
- [ADR-0001: Parity scope](0001-parity-scope.md), [ADR-0002: Licensing and provenance](0002-licensing-provenance.md)
- [Seam contracts](../../docs/contracts/seams.md) — assembly owns error estimation
