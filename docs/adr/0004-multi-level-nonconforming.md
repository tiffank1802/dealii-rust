# ADR-0004: Multi-level non-conforming refinement and selective derefinement

- **Status:** Accepted
- **Date:** 2026-11-21
- **Amended:** 2026-11-21 (revision 5 implementation-review amendment)
- **Deciders:** project maintainers

## Context

Adaptive loops (estimate → mark → refine → solve → coarsen) run for more than
one cycle only if refinement is multi-level and coarsening is selective; the
current `NonconformingRefinement` is a one-level transition over a conforming
`Mesh2` (`src/nonconforming_refinement.rs`). The P0-P2 scope of
[ADR-0001](0001-parity-scope.md) deliberately includes only "one-level
non-conforming (hanging-node) refinement for triangles", and the project
commits to deferred capabilities under their own seams and ADRs/RFCs. This
work expands that scope via
[RFC-0001](../rfcs/0001-multi-level-nonconforming/README.md).

The seam contracts of [ADR-0003](0003-seam-contracts.md) constrain the design:
typed IDs stay index-based, valid only for the mesh instance they came from
(ADR-0003 §4); refinements stay additive (ADR-0003 §3); no external crates,
no unsafe, `rust-version = 1.70`.

## Decision

Accept **RFC-0001 (revision 5), Option C** — multi-level non-conforming
meshes with hanging nodes plus selective derefinement — recorded as
**ADR-0004** on 2026-11-21 and amended 2026-11-21
([README.md](../rfcs/0001-multi-level-nonconforming/README.md)).

Key decisions (preserved from the RFC for decision memory):

- **Multi-level, active-only cell model.** New grid-seam type
  `NonconformingMesh` carries the flattened active (leaf) cell list,
  level-scoped `parent_cells` provenance for the frozen accessor (source-cell
  index in the mesh that produced this level; **not** a self-parent pointer or
  a cell level), the cross-level midpoint source of truth `edge_splits`,
  `hanging`, exterior boundary edges and original-vertex provenance.
  `NonconformingRefinement` wraps it; `cells()` returns active cells only and
  is the single source of truth for assembly and `hanging_constraints`.
- **Vertex reuse across levels.** Splitting an edge reuses the stored
  midpoint via the canonical vertex-pair key; a midpoint created at level *k*
  keeps its index at level *k+1*. The split level is deliberately not stored.
- **Interface detection by `edge_splits`, no cell level.** An edge is a
  coarse/fine interface iff its canonical key is in `edge_splits` and the
  unsplit edge is still incident to an active cell (RFC-0001 §3). No cell
  level is stored or derived; `parent_cells` cannot yield one. A transition's
  ordinal level is its history step index. The `parent_cells` chain is
  ambiguous in the active-only model (a carried root points at an earlier
  index, not at itself), so deriving a level from it produced wrong interface
  classifications and an unchecked-index/cycle-prone `depth()` helper.
- **Coarsening midpoint provenance.** Each transition privately records
  `introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>>`: the midpoints
  *newly created* for each refined sibling group, excluding midpoints reused
  from an earlier level via `edge_splits`. `coarsen` removes only those
  midpoints (plus an orphan cascade); it never infers them from
  `midpoint_parents` master membership within the group. This is level-free
  and index-based; `edge_splits` still holds no level.
- **Unbalanced refinement allowed.** No 2:1 neighbour expansion; continuity is
  enforced through constraints, not mesh quality.
- **Dedicated history.** `NonconformingRefinementHistory` (additive) stores
  transitions (`NonconformingRefinement` per level) with `refine`,
  `coarsen_last`, `level`, `mesh()`, `mesh_at()`; `RefinementHistory` stays
  frozen to `Mesh2`/`LocalRefinement`.
- **No persistent IDs.** Parent/child and boundary-edge provenance is
  level-scoped and index-based; meaningful only within one transition.
  `edge_splits` and `introduced_midpoints` both store vertex indices (and
  level-scoped parent `CellId`s for the latter), never cell paths or persistent
  identifiers; both are invalidated/renumbered on coarsening.
- **Transition-level boundary provenance.** The transition-level parent-edge
  map is named `parent_boundary_edges: Vec<BoundaryEdgeId>` (accessor
  `parent_boundary_edge`); the mesh-level owning-cell field of the earlier
  draft is not materialized — `boundary_edge_geometry` needs the oriented
  exterior edge, not an owner.
- **Additive entry points.** `NonconformingRefinement::new(parent: &Mesh2,
  marked)` and all frozen seam signatures unchanged; new: `refine`,
  `coarsen`, `NonconformingMesh`, `NonconformingRefinementHistory`,
  `DofHandler::from_vertex_count` (scope: vertex-indexed DoFs; no
  `cell_dofs`/`boundary_dofs` via it).
- **Empty markings rejected** by `refine`/`coarsen` (avoid indistinguishable
  levels).
- **Environment unchanged:** serial, triangles only, no new dependencies,
  `#![forbid(unsafe_code)]`, MSRV 1.70. Scope is additive: P0-P2 examples and
  tests and `RefinementHistory` untouched.

## Amendment (2026-11-21)

The accepted decision above is unchanged; two implementation-review defects
narrowed two of its clauses (RFC-0001 revision 5):

1. **Drop "cell level derived from `parent_cells`".** The P0-P2 convention
   stores a previous-level index in `parent_cells`, not a self-parent (a carried
   root computes an arbitrary "depth", and the helper had no bounds or cycle
   guard). Interface detection is now defined by `edge_splits` plus active-cell
   incidence (RFC-0001 §3); no depth function is part of the design. Frozen P0-P2
   accessor behaviour is preserved (ADR-0001).
2. **Track midpoint introduction for coarsening.** `introduced_midpoints` on the
   transition distinguishes midpoints created by the sibling group being
   coarsened from midpoints reused at an earlier level, so `coarsen` can remove
   exactly its own midpoints without storing a level in `edge_splits`. Hanging
   nodes are recomputed from `edge_splits` after coarsening so pre-existing
   interfaces reappear.

Both changes are additive and stay inside ADR-0003 §3–§4. Rationale, exact API
changes and document line references are in
[decision-2026-11-21.md](../rfcs/0001-multi-level-nonconforming/decision-2026-11-21.md).

## Consequences

**Positive**

- Adaptive loops over multiple refine/coarsen cycles run without rebuilding
  or losing provenance; `mesh_at(level)` keeps per-level views.
- `edge_splits` vertex reuse composes multi-level constraints through the
  iterative `AffineConstraints::close`; no recursion needed.
- Deep levels on existing meshes are possible; the deep-level validation
  tests (`deep_refinement_row_width_bounded`,
  `deep_refinement_completes_in_budget`, RFC-0001 §13) measure constraint-fill
  risk.

**Negative / accepted risks**

- ID invalidation on `coarsen` is legal per ADR-0003 §4: all
  `VertexId`/`CellId`/`DofId`/`BoundaryEdgeId` become invalid when indices
  compact.
- `coarsen` restores a parent only when all children are marked, none has
  descendants, and no surviving reference (cell or `HangingNode`) uses the
  excised midpoints; the operation is immutable (returns a new mesh).
- Deep levels can widen constraint rows; documented in RFC-0001 §13.

**Re-evaluation trigger**

Re-open this decision if the adaptivity story later needs persistent
cross-level cell identity (deal.II-style triangulation paths) — per
ADR-0001's re-evaluation triggers.

## Alternatives considered

- **A. Keep one-level only, rebuild each cycle** (RFC-0001 option A): simple,
  matches the current immutable types, but loses provenance and is expensive.
- **B. Full deal.II-style triangulation with persistent cell IDs** (option B):
  matches deal.II but violates ADR-0003 §4 and was rejected at decision time.
- **C. Multi-level with index-based IDs + dedicated history** (chosen):
  builds on the existing one-level code, stays inside the ADR-0003 seam
  contracts, and enables multi-cycle adaptivity without persistent IDs.
