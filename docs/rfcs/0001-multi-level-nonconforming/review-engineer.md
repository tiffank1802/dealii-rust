# Review: RFC-0001 — engineer

## Round 4 — 2026-11-21 — reviewed revision 4 — verdict: approve

### Summary

Revision 4 adds two non-blocking clarifications on top of the approved revision 3:
(1) `edge_splits` renumbering on coarsening is now stated explicitly — surviving
entries are retained and their vertex keys/values are remapped with the same
compaction map as the vertices, alongside the existing cascade drop of entries
whose midpoint or key vertex was removed; (2) the refinement-level definition is
now pinned down — a cell's level is its depth from the root mesh (root = 0),
a transition's level is parent level + 1, and the level metadata lives in
`NonconformingRefinementHistory`, not in `NonconformingMesh`. Both are consistent
with the substrate already reviewed (index-based IDs, active-only cells,
`AffineConstraints::close`) and introduce no new blocking issues. A handful of
Round 3 residual notes (children in the §9 retained-field list, the
`DofHandler::from_vertex_count` storage representation, the `boundary_parents`
naming collision, an explicit "split level intentionally omitted" sentence, and
the illustrative parenthetical in §3) are still not folded in, but they remain
cosmetic/clarificatory and do not block implementation or test writing.

### Verification of revision 4 clarifications

| Clarification | Location | Status |
|---|---|---|
| `edge_splits` renumbering on coarsen | §8 (lines 166-174) | Consistent and complete: survivors keep original vertex order, removed midpoints excised/compacted, all IDs invalidated; `edge_splits` drops entries whose midpoint or key vertex was removed (cascade) and remaps surviving keys/values with the same index map; pre-existing interface edges persist while their midpoints survive. No self-inconsistency found. |
| Refinement-level definition | §1 (lines 73-75) | Consistent with §3 interface discovery: level is per-cell depth from the root, derivable by walking `parent_cells` across the retained transitions in history; transition level is the history step counter. No contradiction with the active-only model or with `edge_splits` (which stores vertex indices, not levels). Minor note: "a transition's level = parent level + 1" reads as a single counter even under selective refinement where marked cells may have differing depths; recommend one sentence clarifying that the transition level is the history step index and that per-cell depth is derived, not the transition level. Non-blocking. |

### Carry-over residual notes (non-blocking, do not block implementation)

1. **§9 retained-field list vs. `children`.** §1 says the transition retains
   parent→child links (see §9), but §9's enumeration
   (`midpoint_parents`, `oriented_boundary_edges`, `boundary_parents`,
   `boundary_vertices`, `original_vertices`) still omits `children`/parent→child
   links. Recommend adding it for implementer clarity. Cosmetic.
2. **`edge_splits` split level intentionally omitted.** The type signature
   (`BTreeMap<(VertexId, VertexId), VertexId>`) makes the omission visible, but
   an explicit "level is deliberately not stored; reuse and coarsening operate on
   midpoint identity" sentence (per engineer clarification D) would prevent
   re-asking. Cosmetic.
3. **`DofHandler::from_vertex_count` representation.** §6 fixes the contract
   (vertex-indexed DoFs; `cell_dofs`/`boundary_dofs` unavailable) but not how a
   handler stores a bare count given `DofHandler<'a>` borrows `&Mesh2`
   (`src/dofs.rs:18-25`). Implementer will need an owned/count variant or
   `Option<&Mesh2>`; core contract unaffected. Non-blocking.
4. **`boundary_parents` naming collision.** `NonconformingMesh.boundary_parents`
   is `Vec<Option<CellId>>` (owning active cell), while the existing one-level
   type stores `boundary_parents: Vec<BoundaryEdgeId>` for `parent_boundary_edge`
   (`src/nonconforming_refinement.rs:51`, `src/nonconforming_refinement.rs:204-205`).
   §1's note distinguishes the two meanings, but both fields share a name; an
   implementer (and the transition, which §9 says retains `boundary_parents`)
   cannot tell which is meant. Recommend renaming the mesh-level field
   (e.g. `boundary_owners`) or the transition-level one. Non-blocking.
5. **§3 interface-discovery parenthetical.** "(one refined this step, one not)"
   is illustrative; the actual condition is differing cell levels, which also
   covers pre-existing interface edges whose cells are not marked this step.
   Recommend rewording to "at different levels (typically one refined this step,
   one not)" so pre-existing interfaces are not misread as excluded. Non-blocking.

### MSRV 1.70 / edge cases

Unchanged from Round 3: triangles only; serial only; no `hp`; no new crates; no
unsafe; `rust-version = 1.70` (`Cargo.toml`); `#![forbid(unsafe_code)]`
(`src/lib.rs:5`); no let-chains/`try` blocks; `BTreeMap`/`BTreeSet` ordering;
guarded indexing; `checked_mul`/`checked_add`. Validation plan (§13) and §8/§11
edge cases (already-split edge reuse, empty marking, non-finite midpoint,
malformed rules, cumulative overflow, partial sibling group, descendants,
neighbour dependency, surviving-vertex reference, ID invalidation) remain covered.

### Verdict

**approve.** Revision 4's two clarifications are consistent and complete enough to
implement against (with the minor wording note on the transition-level counter),
and no new blocking or seam issues were introduced. The carry-over notes above
are cosmetic/clarificatory and can be folded in during implementation; P0-P2
examples/tests remain untouched per ADR-0001.


## Round 3 — 2026-11-21 — reviewed revision 3 — verdict: approve

### Summary

Revision 3 closes the single residual blocking item from Round 2 (the `children`
field that was dead under the active-only cell model) by removing
`children: Vec<Vec<CellId>>` from `NonconformingMesh` and stating that
parent→child links are retained only inside the transition stored by the
history (`README.md` §1, §9). The four non-blocking clarifications from Round 2
are all addressed at least in substance: two explicitly (B, C) and two via the
type contract and global renumbering statement (A, D). The MSRV 1.70 guardrails
and the edge-case list remain present and consistent with the codebase. The
engineering substrate is unchanged and still supports the design
(`AffineConstraints::close` iterative Kahn sort, `NonconformingRefinement`,
`RefinementHistory`). Verdict is **approve**; implementation and test writing
can start. The notes below are cosmetic and do not block.

### Resolution of the Round 2 residual blocking item

| Round 2 item | Resolved in rev 3 | Status |
|---|---|---|
| `children` dead under active-only model (`README.md` §1, §4) | §1 now states "There is **no** `children: Vec<Vec<CellId>>` field: under the active-only model children are not retained in the mesh (they are either active or removed). The transition type retains parent→child links during a single refinement (see §9)." | Resolved. Provenance children are derivable by inverting `parent_cells` within a mesh, and the transition/history retains links across a refinement. No field-level contradiction remains. |

### Resolution of the Round 2 non-blocking clarifications (A–D)

| # | Clarification | Resolved in rev 3 | Status |
|---|---|---|---|
| A | `edge_splits` renumbering on coarsen | §8 states "Renumbering: survivors keep their original vertex order; removed midpoints are excised and indices compacted." and "`edge_splits`: drop entries whose midpoint was removed and any entry whose canonical key contains a removed vertex (cascade). Surviving entries retained." The global compaction statement covers retained entries; their keys/values are implicitly remapped with the same index map because a retained `VertexId` must refer to the compacted index space to remain valid. | Addressed in substance; recommend one explicit sentence ("surviving `edge_splits` entries are renumbered with the same vertex index map") to remove the last ambiguity. Non-blocking. |
| B | `DofHandler::from_vertex_count` cell-connectivity gap | §6 states "`cell_dofs`/`boundary_dofs` are **not** available for nonconforming meshes via this constructor (callers use the vertex-indexed numbering directly)." | Resolved explicitly. |
| C | Hanging discovery restated in edge terms | §3 now defines an edge→cell-incidence map and declares an edge interface iff its two incident active cells are at different refinement levels, with boundary (0/1 incident cells) excluded and each sub-edge key in `edge_splits` tested independently for level difference > 1. | Resolved explicitly. |
| D | Intentional omission of split level | §1's `edge_splits: BTreeMap<(VertexId, VertexId), VertexId>` carries only a midpoint; no level field exists, and reuse (§2) and coarsening (§8) operate on midpoint identity rather than level. The omission is visible in the type signature. | Addressed in substance via the type contract; an explicit "level is intentionally omitted" sentence would prevent re-asking. Non-blocking. |

### MSRV 1.70 constraints (still covered)

§12 guardrails carry: triangles only; serial only (no MPI); no `hp`; no order-mismatch constraints; no new external crates; no unsafe; `rust-version = 1.70` (verified `Cargo.toml:5` `rust-version = "1.70"` and `src/lib.rs:5` `#![forbid(unsafe_code)]`); no let-chains; no `try` blocks; `?`/explicit `return Err`; `checked_mul`/`checked_add` for cumulative counts; no `get_unchecked`. Consistent with the codebase and Round 1/2 lists (BTreeMap/BTreeSet ordering, guarded indexing via `.get`). No further constraints to add. Verified against `src/constraints.rs` (iterative `close`, no recursion) and `src/nonconforming_refinement.rs` (`checked_mul` at `src/nonconforming_refinement.rs:75-79`).

### Edge cases (still covered)

§13 validation plan and §11/§8 cover: already-split edge reuse (§2), empty marking (`new` allows, `refine`/`coarsen` reject `Err`, §11), non-finite midpoint (`0.5*a + 0.5*b` form per `src/local_refinement.rs:64`), malformed-rule rejection (`AffineConstraints::add` + `close` cycle detection), cumulative overflow (`checked_mul`/`checked_add`), partial sibling group, descendants present, neighbour dependency, surviving-vertex reference, and ID invalidation. Deep-level stress tests (`deep_refinement_row_width_bounded`, `deep_refinement_completes_in_budget`) are retained. No new edge cases identified.

### Minor residual notes (non-blocking, do not block implementation)

1. **Transition retained-field list vs. `children()`.** §1 removes `children` from `NonconformingMesh` and says the transition retains parent→child links (§1, §9), but §9's explicit retained-field list (`midpoint_parents`, `oriented_boundary_edges`, `boundary_parents`, `boundary_vertices`, `original_vertices`) omits `children`/parent→child links. Since children are derivable by inverting `parent_cells`, and `NonconformingRefinement::children()` (used by `tests/nonconforming_refinement.rs`) must be served from the transition, recommend adding `children` (or "parent→child links") to §9's retained list for implementer clarity. Non-blocking.
2. **`DofHandler::from_vertex_count` representation.** §6 specifies the additive constructor but not how a handler stores a bare count given the current `DofHandler<'a>` borrows `&Mesh2` (`src/dofs.rs:18-25`). Implementer will need an owned/count variant or `Option<&Mesh2>`; the core contract (vertex-indexed DoFs) is unaffected. Non-blocking.
3. **`mesh()` accessor.** §1 mentions `mesh() -> &NonconformingMesh` on the wrapper and §9 `mesh()`/`mesh_at()` on the history; keep these as the single view entry points so `children()` is clearly transition-served. Non-blocking.

### Verdict

**approve.** The residual blocking item (dead `children` field under the active-only model) is resolved, the four Round 2 clarifications are addressed in substance (B and C explicitly; A and D via the type/global-renumbering contract, with one recommended sentence each), and the MSRV 1.70 constraints and edge-case coverage are intact. Implementation, tests, and the new adaptive Poisson example can proceed; P0-P2 examples/tests remain untouched per ADR-0001.


## Round 2 — 2026-11-21 — reviewed revision 2 — verdict: request-changes

### Summary

Revision 2 resolves every blocking item from Round 1 (interface-edge metadata,
vertex reuse, multi-node/2:1 policy, DoF numbering, active-cell model, coarsening
renumbering, dedicated history type, cross-level provenance scoping) and
addresses the edge-case and MSRV lists. The engineering substrate is unchanged and
still supports the design (`AffineConstraints::close` iterative Kahn sort,
one-level `NonconformingRefinement`, `RefinementHistory`). Composition without
recursion remains feasible on the stated contract.

One residual data-model contradiction remains (§1 `children` vs. §4 active-only
cells), plus three clarifications needed before implementation. Verdict is
therefore **request-changes**, but narrowly: the Round 1 blockers are closed.

### Resolution of Round 1 blocking items

| # | Round 1 item | Resolved in | Status |
|---|---|---|---|
| 1 | Interface-edge metadata (`src/nonconforming_refinement.rs:44-56`) | §1 `edge_splits: BTreeMap<(VertexId, VertexId), VertexId>` canonical vertex-pair → midpoint; single source of truth; accumulates on refine | Resolved |
| 2 | Vertex reuse across levels (`src/nonconforming_refinement.rs:60-91`) | §2 explicit reuse via canonical key, reusing `edge_key` (`src/nonconforming_refinement.rs:10-16`); create only when absent | Resolved |
| 3 | Multiple hanging nodes per edge / 2:1 balance (`src/nonconforming_refinement.rs:22-25, 242-251`) | §5 commits to **unbalanced allowed**; no neighbour expansion; one `HangingNode` per midpoint; `hanging_constraints` unchanged; module-doc warning | Resolved |
| 4 | DoF numbering for nonconforming meshes (`src/dofs.rs:18-25`) | §6 vertex-indexed (incl. hanging), `AffineConstraints` unchanged; additive `DofHandler::from_vertex_count(n)`; `DofHandler::new(&Mesh2)` unchanged | Resolved |
| 5 | Active-cell model (`src/nonconforming_refinement.rs:44-56, 290-316`) | §4 `cells()` returns active (leaf) only; no active flag; provenance in parallel vectors | Resolved |
| 6 | Coarsening renumbering/validity (`src/nonconforming_refinement.rs:62-67`) | §8 immutable `Result<Self>`, one level per call; survivor ordering (original order, removed midpoints excised, compacted); rejection conditions (partial siblings, descendants, neighbour dependency); ID invalidation; `edge_splits` drops removed midpoints | Resolved (see clarification A) |
| 7 | Dedicated history type (`src/refinement_history.rs:15-18, 56-71`) | §9 additive `NonconformingRefinementHistory` mirroring `RefinementHistory`; `RefinementHistory` untouched | Resolved |
| 8 | Cross-level provenance scoping (ADR-0003 §4) | §10 level-scoped IDs in history; `edge_splits` only cross-level structure storing vertex indices | Resolved |

### Residual blocking items

1. **`[README.md §1, §4]` — blocking — `children` is dead under the active-only model.**
   §1 defines `children: Vec<Vec<CellId>>` "per active cell, its children (empty if
   leaf)", and §4 defines active cells as leaves. The conjunction makes `children`
   empty for every entry in every mesh, so the field can never carry non-trivial
   data. Provenance children are derivable by inverting `parent_cells` within the
   new mesh (or from the history's transition record), so the field is either
   redundant or its intended non-empty scenario is unstated. An implementer cannot
   tell whether to store it, drop it, or populate it during `refine`.
   Required: either remove `children` from `NonconformingMesh` and document the
   inversion, or define the one case in which an active cell has children (and how
   that coexists with "active = leaf"). Until then the data model is internally
   contradictory and tests cannot assert on `children`.

### Non-blocking clarifications

- **A. `edge_splits` renumbering on coarsening.** §8 states survivors keep order
  and removed midpoints are excised/compacted, and that `edge_splits` drops
  entries whose midpoint was removed. For retained entries to stay valid, their
  vertex keys/values must be compacted with the same index map. State explicitly
  that `edge_splits` is renumbered alongside vertices (or that it stores
  post-compaction indices), so `coarsen`'s output is self-consistent and
  `coarsen_invalidates_ids_documented` is testable.
- **B. `DofHandler::from_vertex_count`.** §6's additive constructor yields a count
  but `DofHandler::cell_dofs` needs cell connectivity, which a vertex count cannot
  provide. Either widen the constructor to accept a `&NonconformingMesh`-like
  view (or a cell slice), or document that nonconforming callers map
  `VertexId.0 → DofId.0` directly and use the constructor only for `count()`.
  The core contract (vertex-indexed DoFs) is unaffected.
- **C. Hanging discovery precision.** §3's rule ("an edge is interface if its two
  endpoint vertices belong to cells at different levels") is vertex-based and may
  misclassify edges whose endpoints are shared by multiple cells at mixed levels.
  Recommend restating it in edge terms: an edge is interface iff the (up to two)
  active cells sharing it are at different refinement depths, with boundary edges
  excluded. Implementable from `parent_cells` plus history depths; clarifies the
  "computable without a marking" claim.
- **D. `edge_splits` level.** The map stores only a midpoint, not the split level.
  Acceptable because reuse is key-based and coarsening drops by midpoint, but
  state that level is intentionally omitted so reviewers do not re-ask for it.

### Edge cases

All Round 1 edge cases are now covered by §12's validation plan:
already-split edge reuse, empty marking (§3 rejects `Err`, matching
`RefinementHistory::refine`), non-finite midpoint (keep `0.5*a + 0.5*b` form per
`src/local_refinement.rs:64`), malformed-rule rejection (`AffineConstraints::add`
duplicate/slave checks + `close` cycle detection), cumulative overflow
(`checked_mul`/`checked_add`), partial sibling group, descendants present, neighbour
dependency, and ID invalidation. No new edge cases identified.

### MSRV 1.70 constraints

§11 guardrails carry the Round 1 list verbatim: triangles only, serial only, no
`hp`, no new crates, no unsafe, `rust-version = 1.70` (`Cargo.toml`),
`#![forbid(unsafe_code)]` (`src/lib.rs:5`), no let-chains, no `try` blocks, use
`?`/explicit `return Err`, `BTreeMap`/`BTreeSet` ordering, guarded indexing
(`.get`), `checked_mul`/`checked_add` for cumulative counts, no `get_unchecked`.
Consistent with the codebase. No further constraints to add.

### Verdict

**request-changes.** Round 1 blocking items are all resolved and the design is
implementable. One residual data-model contradiction (`children` under an
active-only cell model) must be settled, and clarifications A–D should be folded
into the text, before implementation and test writing start. Once those are
written in, a third round can approve.

## Round 1 — 2026-11-21 — reviewed revision 1 — verdict: request-changes

### Summary

The chosen option (C) is implementable on top of the existing code, and the
frozen `AffineConstraints::close` (`src/constraints.rs:78-150`) already does the
hard part iteratively and cycle-safely. The one-level machinery
(`NonconformingRefinement`, `src/nonconforming_refinement.rs`) is a good
substrate. However, revision 1 is not yet implementable or testable: it does
not pin down the per-mesh metadata needed to *discover* hanging nodes on an
already nonconforming mesh, the vertex-reuse rule across levels, DoF numbering
for nonconforming meshes, or coarsening's renumbering. Until those are specified,
an implementer must invent APIs and tests cannot be written against a stable
contract. The direction is sound; the text needs a second pass.

### Feasibility

**Constraint composition without recursion: yes, with conditions.**

- `AffineConstraints::close` (`src/constraints.rs:78-150`) resolves dependency
  chains with an iterative Kahn topological sort over `dependents`/`indegree`;
  it already passes a 5 000-deep chain test
  (`src/constraints.rs:254-265`). Multi-level composition does not require
  recursion in the caller: one `AffineConstraints::new(n_vertices)` plus one
  `add(slave, 0.0, &[(m0, 0.5), (m1, 0.5)])` per hanging midpoint, then one
  `close()`, is sufficient. `close` substitutes each constrained master into its
  dependents and deduplicates by reduced index through a `BTreeMap`, so chains
  collapse automatically.
- Feasible only if every hanging midpoint is a *distinct vertex index* that is
  a slave in exactly one rule. That holds today for one level
  (`src/nonconforming_refinement.rs:83-91` creates one midpoint per split edge).
  For multi-level it requires **vertex reuse**: a midpoint created at level *k*
  must keep its index at level *k+1* and must not be recreated when its edge is
  split again. The current `NonconformingRefinement::new`
  (`src/nonconforming_refinement.rs:60-167`) always appends fresh midpoints and has
  no deduplication or "edge already split" check. This is the single biggest
  implementation gap; without it, rules refer to different vertices and `close`
  cannot compose across levels.
- `HangingNode` (`src/nonconforming_refinement.rs:22-25`) is one midpoint with
  two masters. That is fine as a *per-midpoint* record; an edge may carry
  several of them. The RFC must state whether 2:1 balance is enforced. If it is
  not, `hanging_constraints` (`src/nonconforming_refinement.rs:242-251`) still
  works unchanged because it iterates `self.hanging` and emits one rule per
  entry — but discovery must produce all of them. If balance *is* enforced,
  `refine` must expand the marked set (deal.II-style), which is a much larger
  change and is not in the current code at all.
- `close` is `self`-consuming. Rebuilding constraints after a structural change
  (refinement or coarsening) means rebuilding from the new mesh's vertex set.
  That matches the immutable refinement types (`LocalRefinement`,
  `UniformRefinement`) and is acceptable; it just must be documented so callers
  do not cache a `ClosedConstraints` across a refinement.

**Memory / performance risks.**

- Row width. `Row.terms` is a `Vec<(usize, f64)>` built via a `BTreeMap`
  (`src/constraints.rs:113,132-133`). For midpoint rules the support of a slave
  is the set of free vertices reachable through masters. The `BTreeMap` dedups
  by reduced index, so a midpoint on a repeatedly subdivided edge stays at 2
  masters in the common case; blow-up occurs when a midpoint's two masters
  expand to *disjoint* free-vertex sets. On a 2D triangular mesh this can grow
  with depth but is bounded by the number of original vertices in the support
  patch, not exponentially, provided vertices are deduped across levels. This
  must be measured, not assumed.
- `condense` (`src/constraints.rs:195-234`) is O(n · w²) in the worst case
  because it nests `row.terms` × `self.rows[j].terms`. For deep levels with
  wide rows this dominates; deep refinement is usually local, but a fully
  refined small mesh at depth 6 is the stress case. Needs a benchmark.
- History storage. Retaining every level (as `RefinementHistory` does,
  `src/refinement_history.rs:15-18`) costs O(total cells across levels). For
  selective refinement this is fine; for full refinement to depth 6 on a small
  mesh it is ~4⁶ cells and is the main memory risk. The RFC asks for a
  deep-level memory test; see below.
- Overflow. Existing code guards `cells.len() * 4` with `checked_mul`
  (`src/nonconforming_refinement.rs:75-79`). Multi-level code must guard the
  *cumulative* count before allocating vertices/cells, because a depth-6 full
  refinement of even a modest mesh overflows `usize` on 64-bit only at ~2⁶³
  cells, but the allocation itself OOMs long before that. Guard early and
  return `Err`.

### Testability

Tests that must exist before this can be merged, grouped by concern. Names are
suggested; the point is the contract each one pins down.

**Continuity (multi-level).**
- `multi_level_hanging_node_is_average_of_masters`: build a 2- and 3-level mesh
  with mixed marking, interpolate a smooth function at every vertex, and assert
  for each `HangingNode` that `value[vertex] == 0.5*(value[m0] + value[m1])`
  within 1e-12. This directly checks the P1 continuity condition.
- `multi_level_no_duplicate_hanging_slaves`: assert every hanging vertex appears
  as a slave in at most one rule and that no vertex is both a slave and absent
  from `hanging_nodes()`.
- `multi_level_boundary_edges_are_parent_boundary`: extend the existing
  one-level boundary assertions to two levels; interface edges must never appear
  in `boundary_edges()`.

**Composition.**
- `constraints_close_resolves_to_free_dofs`: for a 3-level mesh, call `close()`,
  then for every slave index assert it is *not* in `free_dofs()` and that its
  row's masters are all free DoF indices (no slave leaks into a final row).
- `constraints_reconstruct_satisfies_original_rules`: choose arbitrary finite
  values for `free_dofs()`, `reconstruct`, and re-check every midpoint rule
  `x[slave] == 0.5*(x[m0]+x[m1])`. This is the end-to-end composition test.
- `constraints_deep_chain_still_terminates`: reuse the existing 5 000-deep test
  pattern (`src/constraints.rs:254-265`) but seeded from a multi-level mesh's
  actual rules, to guard against accidental recursion or stack growth.
- `fully_constrained_submesh_reconstructs_empty`: a mesh where every vertex is
  hanging (degenerate but possible in tests) must produce a 0×0 system and
  reconstruct from `[]`, as `src/constraints.rs:305-320` requires.

**Derefinement round-trip.**
- `coarsen_full_sibling_group_restores_parent`: refine then coarsen the same
  cells, assert `vertices()` and `cells()` equal the parent's (exact, including
  order) and `hanging_nodes()` is empty.
- `coarsen_rejects_partial_sibling_group`: mark only 3 of 4 siblings, expect
  `Err`.
- `coarsen_rejects_cell_with_descendants`: refine a child, then attempt to
  coarsen its parent, expect `Err`.
- `coarsen_rejects_neighbor_dependency`: a parent whose midpoint is a hanging
  master for an adjacent unrefined cell must not be restorable; expect `Err`.
- `prolongate_then_restrict_recovers_parent_field`: multi-level analogue of
  `tests/refinement_history.rs:44-60`, including a non-constant field.
- `coarsen_invalidates_ids_documented`: after coarsening, old `VertexId`/
  `CellId`/`DofId` must not be assumed valid; tests should use only fresh IDs
  from the returned mesh.

**Deep-level memory.**
- `deep_refinement_row_width_bounded`: refine one cell to depth 6 on a tiny
  mesh, build and close constraints, assert `max(row.terms.len())` is below a
  documented threshold (recommend recording the observed value and setting the
  threshold at 2× it).
- `deep_refinement_completes_in_budget`: same setup, assert wall time and
  allocation stay within a generous bound (e.g. < 1 s, < 64 MiB) so regressions
  in fill are caught. Use `#[ignore]` if it is too slow for default CI, but
  keep it in the suite.

### Blocking issues

1. **`[src/nonconforming_refinement.rs:44-56]` — blocking — no multi-level
   input representation.** The type stores only active cells and hanging nodes;
   there is no interface-edge metadata, split level, or "which edge is split"
   record. The RFC's "tracking which edges are interface edges across levels"
   (`README.md` §Proposal 2) has no data structure, invariant, or construction
   rule. Without it `refine(&self, marked)` cannot be implemented or tested.
   Required: specify a per-mesh edge-split record (e.g. a map from canonical
   vertex-pair to split level / midpoint vertex) and state how it is derived
   when an edge carries more than one hanging node.

2. **`[src/nonconforming_refinement.rs:60-91]` — blocking — vertex reuse across
   levels.** `new` appends a fresh midpoint for every split edge and never
   checks for an existing vertex at that position. Multi-level refinement needs
   the opposite: reuse the stored midpoint when an edge is already split, and
   only create one when newly splitting. The RFC must state the reuse rule and
   the equality test used to detect "same edge" (recommend: canonical
   `(VertexId, VertexId)` key, exactly as `edge_key`, `src/nonconforming_refinement.rs:10-16`).

3. **`[src/nonconforming_refinement.rs:22-25, 242-251]` — blocking — multiple
   hanging nodes per edge.** `HangingNode` is one midpoint and
   `hanging_constraints` emits one rule per entry. That is fine *if* the RFC
   commits to "one `HangingNode` per midpoint, edges may have several". It must
   also state the 2:1 balance policy: either enforce it (and specify the
   neighbour-splitting expansion, which is new code with no precedent here) or
   explicitly allow unbalanced meshes and rely on constraints for continuity.
   Recommendation: allow unbalanced, document it, and test all-midpoint
   constraint satisfaction. Do not silently assume balance.

4. **`[src/dofs.rs:18-25]` — blocking — DoF numbering for nonconforming
   meshes.** `DofHandler` borrows `&Mesh2` only; a nonconforming mesh is not a
   `Mesh2` and conversion is forbidden
   (`src/nonconforming_refinement.rs:36-37`). The RFC touches `dofs.rs` but
   never defines numbering. Recommendation (engineer): keep DoFs vertex-indexed
   (one per vertex, including hanging nodes) and express hanging relations via
   `AffineConstraints`, exactly as `hanging_constraints` already does. If a
   handler is wanted, add an additive constructor taking a vertex count or a
   `&NonconformingRefinement`-like view; do not replace `DofHandler::new`.

5. **`[src/nonconforming_refinement.rs:44-56, 290-316]` — blocking — active-cell
   model.** If inactive parents are retained in `cells`, `assemble_poisson`
   (`src/nonconforming_refinement.rs:297-314`) and every `cells()` caller
   double-counts. The architect recommends the flattened active list; I agree.
   The RFC must commit to it (or to a documented filter) so that assembly and
   tests have one source of truth. Recommend: `cells()` returns active cells only,
   provenance in parallel vectors, no active flag in cell records.

6. **`[src/nonconforming_refinement.rs:62-67]` (proposed `coarsen`) — blocking —
   renumbering and validity.** Removing midpoint vertices renumbers every later
   vertex and invalidates all IDs (legal per ADR-0003 §4, but must be stated).
   There is no compacting helper in the codebase. The RFC must specify: (a)
   immutable operation returning a new mesh (recommended, matching
   `LocalRefinement`/`UniformRefinement`), (b) survivor ordering (recommend:
   original vertex order with removed midpoints excised, indices compacted), and
   (c) rejection conditions — partial sibling groups, descendants present, or a
   neighbour depending on the midpoint. Without (b) round-trip tests cannot be
   written.

7. **`[src/refinement_history.rs:15-18, 56-71]` — blocking — history reuse.**
   `RefinementHistory` stores `Mesh2` levels and calls `LocalRefinement::new`.
   It cannot hold nonconforming levels. The RFC says "reuse the existing
   pattern"; the existing type cannot. Required: a dedicated
   `NonconformingRefinementHistory` (additive, recommended) rather than
   generalising `RefinementHistory`, which would change a frozen public type and
   needs an ADR entry.

### Missing edge cases

- **Edge already split then split again.** Must reuse the existing midpoint and
  not create a duplicate vertex; otherwise constraints do not chain and the mesh
  has coincident vertices that `Mesh2::new`-style validation would reject.
- **Coarsening a parent whose midpoint is a master for a neighbour.** Must be
  rejected; otherwise the neighbour's hanging node dangles.
- **Coarsening with descendants.** A parent whose children were themselves
  refined cannot be restored; reject.
- **Partial sibling group.** Reject unless all siblings are marked, per the RFC.
- **Empty marking.** `RefinementHistory::refine` rejects empty selections
  (`src/refinement_history.rs:57-58`); decide the same policy for the
  nonconforming history to avoid indistinguishable levels.
- **Cumulative overflow.** Guard `vertices.len()` / `cells.len()` growth before
  allocation, not only `cells * 4`.
- **Non-finite midpoints.** `Point2::new` already rejects non-finite; keep the
  `0.5*a + 0.5*b` form (not `0.5*(a+b)`) to avoid overflow, as
  `src/local_refinement.rs:64` does.
- **Cycle safety.** `close` detects cycles (`src/constraints.rs:142-144`), but
  a multi-level builder could create a self-dependency if a midpoint's master list
  accidentally includes itself. The builder's `add` rejects nonzero
  self-dependency only at `close` (`src/constraints.rs:47`); tests should
  include a malformed-rule rejection case.

### MSRV 1.70 constraints

- Target is `rust-version = "1.70"` (`Cargo.toml`) with `#![forbid(unsafe_code)]`
  (`src/lib.rs:5`). No unsafe, no new crates.
- Do not use let-chains / `if let` chains in conditions (stabilised later than
  1.70); use nested `if`/`match`.
- No `try` blocks; use `?` and explicit `return Err(...)`.
- Sorting/dedup must use `BTreeMap`/`BTreeSet` (already the convention,
  `src/nonconforming_refinement.rs:8`); no `slice::sort` reliance on unstable
  APIs.
- Avoid indexing that assumes bounds; use `.get` or guarded loops, consistent
  with `cell_points`/`children` (`src/nonconforming_refinement.rs:190-202`).
- `checked_mul`/`checked_add` for all size growth (already used); extend to
  cumulative counts.
- `forbid(unsafe_code)` also rules out `get_unchecked` for performance in deep
  loops; if profiling shows a bottleneck, document and escalate rather than
  bypass.

### Risks

- **Fill-in / memory growth** (`README.md` §Trade-offs) is real but unquantified.
  Multi-level midpoint rules can widen rows; `condense` is quadratic in row
  width. Needs the deep-level benchmark proposed above before acceptance.
- **Seam leak growth.** `assemble_poisson` and `hanging_constraints` already
  place assembly/DoFs behaviour on a grid type. Multi-level extends that
  surface. Either keep them as documented grid-owned wrappers or move them under
  the assembly seam; decide before the API freezes.
- **ID invalidation.** Coarsening renumbers vertices; all `VertexId`/`CellId`/
  `DofId`/`BoundaryEdgeId` become invalid. Legal under ADR-0003 §4, but must be
  stated in `coarsen`'s docs exactly as `UniformRefinement` documents it
  (`src/grid.rs:225-230`).

### Non-blocking suggestions

- Keep `NonconformingRefinement::new(parent: &Mesh2, marked)` unchanged; add
  `refine(&self, marked) -> Result<Self>` and `coarsen(&self, marked) ->
  Result<Self>` additively, so existing callers and tests are unaffected.
- Name the history type `NonconformingRefinementHistory`; do not overload
  `RefinementHistory`.
- State explicitly that the adaptive Poisson example is a new file and that
  P0-P2 examples/tests are untouched, per ADR-0001.
- Record observed max row width and memory in the deep-level test so future
  regressions are visible.

### Verdict

**request-changes.** The engineering substrate (`AffineConstraints::close`,
one-level `NonconformingRefinement`, `RefinementHistory`) can carry this work,
and composition without recursion is feasible. But revision 1 must resolve the
blocking items — interface-edge metadata, vertex reuse across levels, the
multi-node/2:1-balance policy, DoF numbering for nonconforming meshes, the
active-cell model, coarsening renumbering, and a dedicated history type — and
specify the edge cases above, before implementation and test writing can start.
Once those decisions are written into the RFC, a second review round can
approve.