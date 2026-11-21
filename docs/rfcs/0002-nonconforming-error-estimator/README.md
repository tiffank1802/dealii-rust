---
id: RFC-0002
status: accepted
revision: 1
created: 2026-11-21
decided: 2026-11-21
author: carter-architect-local
adr: ADR-0005
plan: null
supersedes: []
related: [RFC-0001, ADR-0004]
ai-usage: ai-assisted
tags: [adaptivity, error-estimator, nonconforming, hanging-nodes, poisson]
---

# RFC-0002: Residual a posteriori error estimator for nonconforming P1 triangles

**Status:** Accepted
**Revision:** 1
**Author:** carter-architect-local

> This RFC records the design accepted on 2026-11-21 for the workstream following
> [RFC-0001](../0001-multi-level-nonconforming/README.md) / [ADR-0004](../../adr/0004-multi-level-nonconforming.md).
> It is a design decision only: no `src/`, `tests/`, `examples/` or plan file is
> changed here. The durable decision is [ADR-0005](../../adr/0005-nonconforming-error-estimator.md).

## Summary

Add one additive function to the existing `estimator` seam that computes a
cellwise residual indicator for the active cells of a `NonconformingMesh`, with
an explicit jump policy for coarse/fine interfaces derived from `edge_splits`.
Reuse the existing `ResidualEstimate` result type and `mark_bulk` Dörfler marker
unchanged. The estimator is an *indicator*, not a guaranteed bound.

## 1. Context and motivation

RFC-0001 revision 5 / ADR-0004 delivered multi-level nonconforming P1
refinement and selective derefinement: `NonconformingMesh` is active-only with
vertex-indexed DoFs, and hanging relations live in `HangingNode` /
`AffineConstraints` (`src/nonconforming_refinement.rs:104-133`, `:961-971`).
The adaptive loop (estimate → mark → refine → solve) still lacks an estimator.

The conforming estimator cannot be reused:

- `estimator::poisson_dirichlet` takes `&Mesh2` and states that nonconforming
  grids are unsupported (`src/estimator.rs:28-30`).
- It builds edge→cell incidence keyed by an exact canonical vertex pair and
  contributes a jump only when `adjacent.len() == 2`
  (`src/estimator.rs:39`, `src/estimator.rs:77-89`). On a hanging interface the
  coarse edge `(a,b)` is split at midpoint `m`: `(a,b)` has one incident active
  cell and the fine edges `(a,m)`, `(m,b)` also have one each. **Every interface
  contribution is silently skipped**, under-estimating precisely on the refined
  patch that adaptivity is meant to control. Passing such values to
  `mark_bulk` (`src/estimator.rs:113-148`) would under-mark that patch.
- The load vector from `NonconformingRefinement::assemble_poisson`
  (`src/nonconforming_refinement.rs:1018-1044`) contains `∫f φᵢ`, not `∫f²`, so
  the volume term cannot be reconstructed from it. The source closure must be
  supplied to the estimator again.

Finally, `hanging_constraints()` returns an *unclosed* builder
(`src/nonconforming_refinement.rs:961-971`); only `close()` / `condense()` /
`reconstruct()` (`src/constraints.rs:78-150`, `:171-187`, `:195-234`) yield the
continuous trace. The estimator must therefore consume the reconstructed
solution vector, not a raw assembled one.

## 2. Goals / non-goals

**IN:** P1 triangles only; active-only `NonconformingMesh`; homogeneous
Dirichlet and scalar unit diffusion (the same assumptions as
`src/estimator.rs:19-29`); residual *indicator* (not a guaranteed bound);
interface jumps paired through `edge_splits`; reuse of `ResidualEstimate` and
`mark_bulk`; one example update and one new test file in a later plan.

**OUT:** global/hierarchical estimators, variable diffusion or reaction,
Neumann boundary residuals, recovery-based (ZZ) and goal-oriented estimators,
coarsening-error estimators, anisotropic meshes, higher order / P2+, quads,
tetrahedra, 3D, distributed/MPI, matrix-free execution, new quadrature rules,
persistent cell identity, and any rework of ADR-0003 seams.

## 3. Estimator definition

For each active cell `K`:

```
η_K² = h_K² ∫_K f² + ½ Σ_{e ⊂ ∂K} h_e ∫_e [∂ₙ u_h]²
```

with `h_K` the longest edge of `K`, `h_e` the edge length as seen from the
charging cell, `Δu_h = 0` inside a P1 triangle, and a constant jump per edge
segment. `∫_K f²` is approximated by the crate's existing degree-five triangle
quadrature (`src/quadrature.rs:12`); gradients come from `AffineTriangle`.
This reproduces the frozen conforming formula at `src/estimator.rs:23-29`.

### Interface policy (two-sided edges + `edge_splits` pairing)

1. Build active-cell edge incidence exactly as `discover_hanging` does
   (`src/nonconforming_refinement.rs:62-76`).
2. `incidence.len() == 2`: ordinary interior edge — apply the existing formula
   unchanged. A key in `edge_splits` with two incident active cells is instead
   an `Err`, a defensive guard for malformed input that the `new` / `refine` /
   `coarsen` API never produces: splitting a shared edge always leaves a single
   coarse owner. The guard keeps such a mesh from silently skipping the
   interface term.
3. A key in `edge_splits` with `incidence.len() == 1` is an interface root; its
   single incident cell is the **coarse owner**. Recursively split the master
   edge through `edge_splits` to leaf sub-segments, using a visited set, the
   same iterative pattern as `collect_subedges`
   (`src/nonconforming_refinement.rs:22-55`). This also reaches sub-midpoints of
   unbalanced, multi-level chains, which are emitted by the root walk rather
   than by their own incidence key.
4. Each leaf segment must have exactly one incident active **fine** cell. A
   second incident cell, or none, is an error (`Err`), never a silent skip.
   Split exterior half-edges are excluded because their master key is not
   incident to any active cell, so they never become interface roots.
5. Per leaf segment with `jump = (∇u_coarse − ∇u_fine) · n`:

   - coarse cell receives `0.5 · L_master · l_leaf · jump²`;
   - fine cell receives `0.5 · l_leaf² · jump²`.

   When `L_master == l_leaf` (conforming edge) both reduce to the existing
   `0.5 · h_e² · jump²`. The per-side scaling is the direct generalization of
   `½ h_e ∫_e jump²` to an edge whose two sides have different resolutions.

A constrained P1 trace is continuous across an interface: coarse and fine
traces agree at `a`, `m`, `b`, so their **tangential** derivatives match. Their
**normal** derivatives need not (e.g. coarse `(1,1)` vs fine `(1,-0.5)` across a
horizontal interface, jump `1.5`). Interface jumps are therefore real and must
be charged.

## 4. API

```rust
// src/estimator.rs (assembly seam) — one new public function
pub fn poisson_dirichlet_nonconforming<F>(
    mesh: &NonconformingMesh,
    solution: &[f64],
    source: F,
) -> Result<ResidualEstimate>
where
    F: Fn(Point2) -> f64;
```

Private helpers in the same file: an interface leaf-segment walker and a jump
accumulator. Reuse the existing `ResidualEstimate`
(`src/estimator.rs:10-17`): `cell_indicators` is aligned to `mesh.cells`
(active-only order), and `global_indicator = sqrt(Σ η_K²)`.

Input contract and validation:

- `mesh.validate()?` (`src/nonconforming_refinement.rs:143-195`);
- `solution.len() == mesh.vertices.len()` and all values finite;
- the solution must be the **reconstructed** vector: for every
  `mesh.hanging` entry,
  `|u[slave] − ½(u[m0]+u[m1])| ≤ 1e-12 · (1 + max(|u[slave]|, |u[m0]|, |u[m1]|))`,
  otherwise `Err`. The scale is the largest magnitude among the three values in
  that rule, not a bare `1 + |u|`; do not tighten it to an unscaled equality
  test. Do not silently reconstruct — that would hide a bad solve and would
  interact with Dirichlet rules;
- source values, gradients, edge lengths and all accumulators must be finite;
- inconsistent interface incidence is `Err`.

Marking is unchanged: `mark_bulk(&est.cell_indicators, theta)` returns
`Vec<CellId>` in the current active index space. Do not add a `mark_dorfler`
alias (YAGNI).

## 5. Seam discipline

- The estimator belongs to the `assembly` seam, which already owns error
  estimation (`docs/contracts/seams.md:32`); `estimator.rs` may consume
  `NonconformingMesh` from `grid` without changing seam direction.
- No method is added to `NonconformingMesh` or `NonconformingRefinement`: a
  grid-owned estimator would introduce a grid→assembly dependency and violate
  one-owner-per-seam (`docs/contracts/seams.md:35-43`). The existing
  `assemble_poisson` wrapper on the transition is a pre-existing exception,
  not a precedent to extend.
- Additive only: no existing public signature changes; P0-P2 conforming tests
  and `poisson_dirichlet` stay frozen; `RefinementHistory` untouched.
- No ADR-0003 seam rework: no new typed IDs, no persistent IDs, no new
  quadrature tables (ADR-0002 provenance is respected — the formula is derived
  from first principles, no upstream code or constants are reused).
- MSRV 1.70 (no let-chains), `#![forbid(unsafe_code)]`, zero new dependencies,
  serial. No new module or `lib.rs` change is required; `estimator` is already
  public.

## 6. Adaptive loop

```text
transition = NonconformingRefinement::new(&parent, &[])?;   // or a seeded patch
loop ≤ max_cycles:
    (A, b) = transition.assemble_poisson(source)?;
    c = transition.hanging_constraints()?;                  // unclosed builder
    c.add(id, value, &[]) for (id, value) in transition.interpolate_boundary(exact)?;
    closed = c.close()?;
    (A_r, b_r) = closed.condense(&A, &b)?;
    u = closed.reconstruct(&cg(&A_r, &b_r)?)?;              // slaves correct here
    est = poisson_dirichlet_nonconforming(transition.mesh(), &u, source)?;
    marked = mark_bulk(&est.cell_indicators, theta)?;
    if marked.is_empty() || est.global_indicator <= tol { break }
    transition = transition.refine(&marked)?;               // all ids now stale
```

Notes:

- Only `refine` is driven by η. Naive coarsening is **not**: `coarsen` requires
  all four children of a sibling group, no descendants, and no surviving
  reference to removed midpoints (`src/nonconforming_refinement.rs:598-616`);
  an indicator alone cannot supply that. `NonconformingRefinementHistory`
  offers only level-pop `coarsen_last` (`:1170-1175`). Keep the existing
  deterministic coarsen demonstration in the example separate from the
  estimator-driven loop (`examples/poisson_nonconforming_adaptive.rs:101-151`).
- Marked `CellId`s belong to the current active index space and must be
  discarded after every structural change (ADR-0003 §4,
  `docs/adr/0003-seam-contracts.md:172-180`; coarsen invalidation at
  `src/nonconforming_refinement.rs:537-546`). Re-solve after each change; no
  solution transfer is defined for nonconforming refinement.
- `DofHandler::from_vertex_count` is count-only (`src/dofs.rs:58-60,78-85,
  107-114`): the estimator indexes vertices directly and needs no DofHandler.

## 7. Validation plan

New file `tests/nonconforming_estimator.rs`; names are illustrative, assertions
are the contract.

| Test | Assertion |
|---|---|
| `nonconforming_estimator_matches_conforming_on_uniform_refinement` | Mark all cells via `NonconformingRefinement::new`; per-cell and global indicators equal `estimator::poisson_dirichlet(UniformRefinement::new(&parent).mesh(), …)` within `1e-10`. No interfaces → identical arithmetic. |
| `uniform_refinement_halves_global_indicator` | Three levels with `u = x²+y²`, `f = -4`, Dirichlet from `exact`: `0.25 < η_{k+1}/η_k < 0.65` (theory ≈ 0.5), all η finite and positive. Upper bound mirrors `tests/refinement.rs:170-174`. |
| `hanging_interface_jump_is_counted` | Solve a one-level hanging mesh; clone the mesh and clear `edge_splits`, `hanging`, `midpoint_parents` (public fields) to make it interface-blind; assert `η_aware > 1.01 · η_blind`. Directly detects the primary failure mode. |
| `hanging_interface_vanishes_after_both_sides_refined` | Refine cell 0, then the hanging neighbour; build `Mesh2::new(mesh.vertices.clone(), mesh.cells.clone())`; new estimator equals conforming `poisson_dirichlet` on that mesh. |
| `unbalanced_chain_estimates_finite` | Refine a child on an existing interface to create a sub-midpoint; `Ok`, finite positive η, hanging rules satisfied. |
| `estimator_rejects_unreconstructed_hanging_values` | Perturb one slave value by `0.1` → `Err`. |
| `estimator_rejects_invalid_numerical_data` | Short vector, NaN vector, NaN source, `\|f\| = f64::MAX`, invalid mesh → `Err` (mirrors `tests/refinement.rs:181-191`). |
| `mark_bulk_on_nonconforming_indicators_returns_active_cell_ids` | Every returned `CellId < cells.len()`; `refine(&marked)` succeeds; `theta = 0` returns empty. |
| `affine_solution_has_zero_indicator` | `u = 1 + 2x − 3y`, `f = 0` → all indicators `0`. |

Where possible use crafted vectors satisfying the hanging rules (exact
quadratic values with each slave overwritten by its masters' average) instead
of solver output, so assertions are deterministic.

## 8. Risks

- **Silent under-estimation (primary).** Killed by the two-sided interface
  policy, the "one fine owner per leaf" incidence check, and the differential
  test `hanging_interface_jump_is_counted`. Any incidence mismatch is `Err`,
  never a skip.
- **Hanging-value mismatch.** Consumed solution must be reconstructed; the
  tolerance check turns a silent gradient error into an `Err`.
- **Solver dependence.** Avoided by crafted constraint-satisfying vectors in
  the hanging-specific tests; the rate test uses `x²+y²` with a documented
  ratio band.
- **Constraint row-fill is a non-issue here.** The estimator does not call
  `close()` — it only reads `mesh.hanging`; the documented deep-chain fill risk
  (`docs/adr/0004-multi-level-nonconforming.md:117-120`) stays confined to the
  solve. The iterative `edge_splits` walk has a visited set, so deep interfaces
  cannot loop.
- **`h_e` convention.** Coarse/fine edge lengths differ; the per-side scaling
  above is documented explicitly so future maintainers do not "simplify" it
  to one length and reintroduce a bias.
- **Cost.** The source is resampled by degree-five quadrature in both the
  assembly and the estimator; this duplicates work but keeps the two seams
  decoupled and the frozen conforming path unchanged.

## Alternatives considered

| Option | Pros | Cons |
|---|---|---|
| Two-sided edges via `edge_splits` (chosen) | Counts the real interface jump; works for unbalanced/multi-level chains; no seam change | Requires a leaf-segment walk and per-side scaling |
| Skip interface jumps with a documented bias | Trivial implementation | Silent under-estimation on the refined patch; breaks Dörfler marking — rejected |
| Algebraic projection through `AffineConstraints` | Reuses the constraint machinery | Only exposes `x = Cy + g`, not an edge projection; would unwrap deep chains and add a seam surface — rejected |
| Method on `NonconformingMesh` | Convenient call site | Grid→assembly dependency; violates one-owner-per-seam — rejected |
| Refactor the conforming estimator to accept both meshes | One code path | Changes frozen P0-P2 numerical behaviour; `&Mesh2` contract — rejected |

## Trade-offs / Impact

Touches `src/estimator.rs` (one public function plus private helpers),
`examples/poisson_nonconforming_adaptive.rs`, and one new test file. No
dependency, seam, ID or quadrature changes. The design is deliberately narrow:
an indicator for homogeneous-Dirichlet scalar Poisson, not a certified bound.

## Open questions

None blocking. The per-side `h_e` scaling and the hanging-rule tolerance are
documented conventions; revisit only with evidence from the planned rate and
hanging-interface tests.

## Decision

**Accepted** on 2026-11-21 as the basis for
[ADR-0005](../../adr/0005-nonconforming-error-estimator.md). RFC-0001 is not
reopened; ADR-0003 seams remain frozen. Implementation is sequenced by a
separate plan (not created by this RFC).
