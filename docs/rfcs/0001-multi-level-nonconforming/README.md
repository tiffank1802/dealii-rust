---
id: RFC-0001
status: accepted
revision: 5
created: 2026-11-21
decided: 2026-11-21
amended: 2026-11-21
author: carter-maestro
adr: ADR-0004
plan: 001-multi-level-nonconforming
supersedes: []
ai-usage: ai-assisted
tags: [adaptivity, refinement, nonconforming, hanging-nodes]
---

# RFC-0001: Multi-level non-conforming refinement and selective derefinement

**Status:** Accepted
**Revision:** 5
**Amended:** 2026-11-21 (implementation-review amendment, see
[decision-2026-11-21.md](decision-2026-11-21.md))
**Author:** carter-maestro

## Summary

Extend the existing one-level triangular non-conforming refinement
(`NonconformingRefinement`, ADR-0003/P0-P2 scope) to **multi-level**
non-conforming meshes with hanging nodes, plus **selective derefinement**, while
keeping the crate serial, zero-dependency and `forbid(unsafe_code)`. This moves
the project beyond the P0-P2 scope defined in ADR-0001 and therefore requires
this RFC. Revisions 1–4 resolved the Round 1–4 design reviews; revision 5
amends the parent-cell/level and coarsening-midpoint contracts found during
implementation review (2026-11-21). The accepted Option C direction, the data
model and all public signatures are unchanged.

## Problem

The current implementation (`src/nonconforming_refinement.rs`) supports a single
refinement level from one conforming parent mesh: it detects hanging nodes from
parent connectivity, generates P1 midpoint constraints, and exposes volume
assembly + prolongation. It cannot:

- refine an already non-conforming mesh (multi-level),
- coarsen selected cells back to a parent level (selective derefinement),
- keep constraints consistent across more than one hanging level,
- track provenance so a cell knows its parent/children across levels.

Without these, adaptive loops (estimate → mark → refine → solve → coarsen) cannot
run for more than one cycle, and the library cannot reproduce deal.II's
`Triangulation`/`DoFHandler` adaptivity story beyond a single step.

## Proposal

### 1. Mesh representation: a dedicated `NonconformingMesh`

Introduce a new grid-seam type `NonconformingMesh` (additive, does not replace
`Mesh2`) carrying the active mesh plus provenance needed for multi-level work:

- `vertices: Vec<Point2>` — all vertices (original + midpoints), compacted.
- `cells: Vec<[VertexId; 3]>` — **active (leaf) cells only** (flattened list).
- `parent_cells: Vec<CellId>` — per active cell, the `CellId` of the **source
  cell in the mesh that produced this level**: for `new`, the index in the
  input `Mesh2`; for `refine`, the index in `self.mesh.cells`; for `coarsen`,
  the parent-level group index recorded in the transition's `children`. This
  keeps the existing `Vec<CellId>` convention so the frozen
  `parent_cell(child) -> Result<CellId>` accessor, its P0-P2 physical-tag
  lookups and its tests are unchanged (`src/nonconforming_refinement.rs:96`,
  `tests/nonconforming_refinement.rs:96,225`). `None` is not used. **These ids
  are opaque, level-scoped provenance tokens.** They are *not* self-parent
  pointers (a carried root keeps its index in the input mesh, which generally
  differs from its index in the active list), *not* a cell level, and must not
  be used for level comparisons or to index another level's cell list. See
  §3 for the actual interface rule.
- `edge_splits: BTreeMap<(VertexId, VertexId), VertexId>` — canonical (sorted)
  vertex-pair → midpoint vertex, for every split edge at every level. Single
  source of truth for "is this edge split, and by which vertex".
- `hanging: Vec<HangingNode>` — one entry per midpoint vertex (slave), with its
  two masters. An edge may carry several midpoints (level difference > 1).
- `boundary_edges: Vec<[VertexId; 2]>` — exterior edge vertex pairs only.
  Interface edges are **never** boundary edges.
- The transition-level parent-edge map is named `parent_boundary_edges:
  Vec<BoundaryEdgeId>` — transition-only, retained in `NonconformingRefinement`
  and mapped by the accessor `parent_boundary_edge(child) ->
  Result<BoundaryEdgeId>`. The mesh-level owning-cell field of the earlier
  one-level draft is not part of this design: `boundary_edge_geometry` needs
  the oriented exterior edge, which the transition retains
  (`oriented_boundary_edges`), so the owning-cell field is dropped as unused
  (YAGNI).
- **No per-cell level is stored or derived.** Interface detection (§3) uses
  `edge_splits` plus active-cell edge incidence, not a cell level; the
  `parent_cells` chain is provenance only and cannot yield a level (a carried
  root may point at an earlier cell that is not its own active-list index).
  A transition's ordinal level is its history step index in
  `NonconformingRefinementHistory`; this is an operation counter, not a
  property of any cell. Cells of one mesh may have different effective depths
  under unbalanced refinement (§5), and no depth value is computed anywhere.
- `boundary_vertices: Vec<VertexId>` — exterior vertices (for `boundary_vertices()`).
- `original_vertices: usize` — count of level-0 vertices, for provenance/tests.
- `midpoint_parents: BTreeMap<VertexId, (VertexId, VertexId)>` — per midpoint
  vertex, its two master vertices. Needed by `prolongate_p1` and
  `interpolate_boundary`.

There is **no** `children: Vec<Vec<CellId>>` field: under the active-only model
children are not retained in the mesh (they are either active or removed). The
transition type retains parent→child links during a single refinement (see §9).

`NonconformingRefinement` is refactored to **wrap** a `NonconformingMesh`
(additive: `mesh() -> &NonconformingMesh` accessor; existing `vertices()`/
`cells()`/`hanging_nodes()`/etc. delegate to the inner mesh). The existing
constructor `NonconformingRefinement::new(parent: &Mesh2, marked)` stays
unchanged. New entry points `refine(&self, marked) -> Result<Self>` and
`coarsen(&self, marked) -> Result<Self>` are additive.

### 2. Vertex reuse across levels (blocking item 2)

When splitting an edge, look up the canonical key `(min(v0,v1), max(v0,v1))` in
`edge_splits`. If present, reuse the stored midpoint vertex; only create a new
vertex when the key is absent. This is what makes constraint composition across
levels possible: a midpoint created at level *k* keeps its index at level *k+1*.
Reuse the existing `edge_key` helper (`src/nonconforming_refinement.rs:10-16`).

The split level is deliberately not stored: vertex reuse is key-based and
coarsening drops by midpoint identity, so a stored level could go stale.
Midpoint *provenance for coarsening* is recorded separately on the transition
in `introduced_midpoints` (§8); `edge_splits` itself holds only the current
midpoint for each key.

### 3. Interface / hanging-edge discovery (blocking item 1, 5)

`edge_splits` is the single source of truth for interface detection. An edge
(canonical vertex pair) is **interface** iff its key is present in `edge_splits`
**and** the unsplit key still occurs among the edges of at least one active
cell. No cell level, `parent_cells` chain, or marking comparison is consulted.
Concretely, during `refine` (and again after `coarsen`, §8):

1. Build an edge→cell-incidence map: for each active cell, for each of its 3
   edges (canonical key), append the cell's id to that edge's list.
2. An edge is an interface root iff `edge_splits.contains_key(edge)` and the
   edge occurs in that incidence map. This one rule covers both newly created
   interfaces and **pre-existing interfaces** whose cells were not marked this
   step: their midpoint remains in `edge_splits`, so they are re-tested and
   their `HangingNode`s re-emitted automatically. An interface disappears once
   both sides are refined and the unsplit key is no longer an active-cell edge,
   even though the midpoint stays in `edge_splits` for reuse (§2).
3. For each interface edge, walk its split sub-edge keys in `edge_splits`
   (iteratively, with a visited set): emit one `HangingNode` per split key
   found. An edge may carry several midpoints (level difference > 1); each
   sub-edge key is tested independently.

Edges with 0 or 1 incident active cells are boundary edges (never interface;
exterior splitting is handled by the retained oriented boundary edges).
Edges with 2 incident active cells whose key is absent from `edge_splits` are
interior (no hanging node). The discovery procedure is shared by `refine` and
`coarsen` so a pre-existing interface reappears when both sides become coarse
again.

### 4. Active-cell model (blocking item 5)

`cells()` returns **only active (leaf) cells**. There is no active flag in cell
records and no inactive-inclusive list. Provenance (parents) lives in the
parallel `parent_cells` vector. `assemble_poisson` and `hanging_constraints`
iterate `cells()` and therefore see only active cells. This is the single
source of truth for assembly.

### 5. 2:1 balance policy (blocking item 3)

**Unbalanced refinement is allowed.** No neighbour-splitting expansion is
performed. Continuity is enforced by constraints, not by mesh quality. The RFC
documents this explicitly: meshes may have level differences > 1 across an edge,
and `hanging_constraints` emits one midpoint rule per `HangingNode`. Users are
warned (in module docs) that unbalanced meshes can be ill-conditioned.

### 6. DoF numbering for nonconforming meshes (blocking item 4)

DoFs remain **vertex-indexed** (one per vertex, including hanging nodes), exactly
as the one-level code does. Hanging relations are expressed through
`AffineConstraints` via `hanging_constraints`. No change to `DofHandler`'s
signature. Additively, a constructor
`DofHandler::from_vertex_count(n: usize) -> Self` is added for callers that
only have a vertex count; `cell_dofs`/`boundary_dofs` are **not** available for
nonconforming meshes via this constructor (callers use the vertex-indexed
numbering directly). Existing `DofHandler::new(&Mesh2)` is unchanged. Storage
shape is the implementer's choice: the current `DofHandler<'a>` borrows
`&Mesh2` (`src/dofs.rs:18-20`), so a count-only constructor implies an
owned-count or `Option<&Mesh2>` variant; the contract — no
`cell_dofs`/`boundary_dofs` via this constructor — is fixed.

### 7. Constraint composition (no recursion)

Multi-level constraints are built by a single `AffineConstraints::new(n_vertices)`
plus one `add(slave, 0.0, &[(m0, 0.5), (m1, 0.5)])` per `HangingNode`, then
one `close()`. `close()` (iterative Kahn sort, `src/constraints.rs:78-150`)
resolves chains automatically and detects cycles. Constraints are rebuilt from
scratch after every refinement/coarsening (they are not cached across
structural changes, matching the immutable refinement types).

### 8. Selective derefinement (`coarsen`)

- Input: `&self` + `marked: &[CellId]` (active cells to coarsen).
- Each transition records midpoint provenance in a private field
  `introduced_midpoints: BTreeMap<CellId, BTreeSet<VertexId>>` on
  `NonconformingRefinement` (never in `edge_splits` or in `NonconformingMesh`).
  For a refined sibling group `p`, the set contains the `VertexId`s of the
  midpoints **newly created** while refining `p`: exactly the midpoints of the
  group's three canonical edge keys that were **absent** from `edge_splits` at
  the start of that transition. A midpoint reused from an earlier level (§2)
  is deliberately **excluded**, even when both of its masters happen to be
  vertices of the sibling group. A midpoint introduced for an edge shared by
  two marked groups belongs to both groups' sets. This is a level-free,
  index-based provenance record; no split level is stored anywhere.
- The candidate removal set for a restored group is its `introduced_midpoints`
  entry. A midpoint is **never** selected for removal by inspecting
  `midpoint_parents` master membership within the group's vertices: that
  heuristic misclassifies a reused interface midpoint and would reject a legal
  derefinement (the level-1 interface midpoint reused when the coarse neighbour
  is refined in the next level). After all selected groups are collapsed,
  midpoints that have become geometrically unused are additionally removed by
  the orphan cascade below (this also removes nested sub-edge midpoints).
- A parent is restored **only if** all of the following hold:
  1. all its children are in `marked`;
  2. none of those children has descendants;
  3. no surviving active cell references any midpoint in the candidate removal
     set as an ordinary vertex;
  4. no remaining `HangingNode` (as slave or master) references any midpoint in
     the candidate removal set.
- Operation is **immutable**: returns a new `NonconformingMesh`.
- **One level per call** (callers loop), matching `RefinementHistory::coarsen_last`.
- Renumbering: survivors keep their original vertex order; removed midpoints are
  excised and indices compacted. All `VertexId`/`CellId`/`DofId`/
  `BoundaryEdgeId` are invalidated (legal per ADR-0003 §4); documented in
  `coarsen`'s rustdoc.
- `edge_splits`: drop entries whose midpoint was removed **and** any entry whose
  canonical key contains a removed vertex (cascade). Surviving entries retained;
  their vertex indices are remapped to the new compacted indices alongside the
  vertices. Pre-existing interface edges (already split before this coarsening)
  persist as long as their midpoints survive; a reused midpoint is not in any
  later transition's `introduced_midpoints`, so it survives until the
  introducing transition is undone or the orphan cascade removes it.
- The returned transition carries the surviving `introduced_midpoints` entries,
  renumbered with the same vertex compaction map (entries whose midpoint was
  removed are dropped), so a subsequent coarsen of the coarser mesh keeps the
  same provenance rule.
- Hanging nodes of the result are recomputed with the §3 discovery routine on
  the coarsened active mesh, rather than carried forward by filtering the
  input's `hanging` list. This restores a pre-existing interface whose midpoint
  had become inactive while both neighbours were refined at the deeper level.

### 9. History storage (blocking item 7)

Introduce a **dedicated** `NonconformingRefinementHistory` (additive), mirroring
`RefinementHistory` (`src/refinement_history.rs`): stores one
`NonconformingRefinement` (transition) per level, with `refine(marked)`,
`coarsen_last()`, `level()`, `mesh()` (current level) and `mesh_at(level)`.
The history stores **transitions** (which retain `children` / parent→child
links, `midpoint_parents`, `introduced_midpoints`, `oriented_boundary_edges`,
`parent_boundary_edges`, `boundary_vertices`, `original_vertices`), not bare
meshes; `mesh()`/`mesh_at()` expose the `NonconformingMesh` view. Mesh-only
fields (`edge_splits`, `hanging`, `boundary_edges`, `parent_cells`) are not
transition-retained; `parent_cells` remains a field of the exposed mesh for the
frozen accessor, not a level source.
Do **not** overload `RefinementHistory`, which is frozen to
`Mesh2`/`LocalRefinement`.

### 10. Cross-level provenance scoping (blocking item 8)

Parent/child IDs are **level-scoped**: they are meaningful only within one
transition and are stored in the level-retaining history, not as long-lived
cross-level identifiers. This respects ADR-0003 §4 (no persistent/path-based IDs
in P0-P2). `edge_splits` maps canonical vertex pairs to vertex indices, and
`introduced_midpoints` maps level-scoped parent `CellId`s to vertex indices;
both are level-scoped index-based records, never cell paths or persistent IDs.

### 11. Empty-marking policy (clarification 7)

- `NonconformingRefinement::new(parent, marked)` keeps accepting empty markings
  (existing tests rely on this).
- `refine` and `coarsen` **reject** empty markings (`Err`), to avoid
  indistinguishable levels in the history. `NonconformingRefinementHistory::refine`
  delegates to `refine` and thus inherits the rejection.

### 12. Guardrails

- Triangles only; quads/tetra adaptivity deferred.
- Serial only; no MPI.
- No `hp`, no order-mismatch constraints.
- No new external crates; no unsafe; MSRV 1.70 (no let-chains, no try blocks,
  no `get_unchecked`; use `checked_mul`/`checked_add` for cumulative counts).
- `NonconformingRefinement::new(parent, marked)` signature unchanged; new
  methods are additive.
- Assembly/dofs convenience methods on the grid type stay as documented
  grid-owned wrappers (no seam change this RFC).
- New example file for adaptive Poisson; P0-P2 examples/tests untouched.

### 13. Validation plan

Tests (names are illustrative; the contract is the assertion):
- Continuity: `multi_level_hanging_node_is_average_of_masters`,
  `multi_level_no_duplicate_hanging_slaves`,
  `multi_level_boundary_edges_are_parent_boundary`.
- Composition: `constraints_close_resolves_to_free_dofs`,
  `constraints_reconstruct_satisfies_original_rules`,
  `constraints_deep_chain_still_terminates`,
  `fully_constrained_submesh_reconstructs_empty`.
- Derefinement: `coarsen_full_sibling_group_restores_parent`,
  `coarsen_rejects_partial_sibling_group`,
  `coarsen_rejects_cell_with_descendants`,
  `coarsen_rejects_neighbor_dependency`,
  `coarsen_reused_interface_midpoint_restores_parent` (replaces the earlier
  `coarsen_rejects_surviving_vertex_reference`: a midpoint reused from an
  earlier level is not removed, so the legal derefinement must succeed and
  restore the parent mesh, including its hanging node),
  `prolongate_then_restrict_recovers_parent_field`,
  `coarsen_invalidates_ids_documented`,
  `coarsen_children_and_parent_cells_share_index_space`.
- Deep-level: `deep_refinement_row_width_bounded`,
  `deep_refinement_completes_in_budget` (ignored by default if slow).
- Edge cases: already-split edge reuse, empty marking, non-finite midpoint,
  malformed rule rejection, cumulative overflow,
  `mesh_validate_rejects_out_of_range_ids_and_nonfinite_coordinates`.

## Alternatives considered

| Option | Pros | Cons |
|---|---|---|
| A. Keep one-level only, rebuild mesh each cycle | Simple, matches current immutable types | Loses provenance, expensive, error-prone |
| B. Full deal.II-style triangulation with persistent cell IDs | Matches deal.II exactly | Violates ADR-0003 §4; too large |
| C. Multi-level with index-based IDs + dedicated history (chosen) | Builds on existing code, stays in ADR-0003 | Less faithful to deal.II; manual provenance |

## Trade-offs / Impact

- Touches `src/nonconforming_refinement.rs`, new `NonconformingRefinementHistory`,
  `src/constraints.rs` (docs only), `src/dofs.rs` (additive constructor),
  `src/grid.rs` (NonconformingMesh type), and new tests/examples.
- Risk: constraint fill/memory growth on deep levels — must be measured (deep-level tests).
- Risk: derefinement invalidates all IDs — documented per ADR-0003 §4.
- No new external crates; no unsafe; MSRV 1.70.

## Open questions

1. Maximum supported depth before constraint composition becomes impractical?
   → Documented, not API-capped; test up to 4 levels (unit) / 6 levels (bench).
2. Derefinement across a level boundary in one call?
   → One level per call; callers loop.
3. Active cells view?
   → Flattened active list; `cells()` returns active only.

## Decision

**Accepted.** RFC-0001 revision 5, **Option C** — multi-level with
index-based IDs + a dedicated `NonconformingRefinementHistory` — is accepted
as **ADR-0004** on 2026-11-21 (amended 2026-11-21)
([ADR-0004](../adr/0004-multi-level-nonconforming.md)).

Brief rationale: the adaptive loop (estimate → mark → refine → solve →
coarsen) requires a multi-level mesh and selective derefinement; the work is
additive per ADR-0003 §3 and expands the P0-P2 scope per ADR-0001. Vertex
reuse via the canonical `edge_splits` key map composes constraints across
levels through the existing iterative `AffineConstraints::close`; level-scoped
IDs keep cross-level provenance inside the history, per ADR-0003 §4.

Round 4 residual clarifications were folded into revision 4 (split level not
stored, transition level is the history step index, `parent_boundary_edges`
naming, parent→child links in §9, `DofHandler::from_vertex_count` storage
note). Revision 5 amends two implementation-review defects without changing
this decision or any signature: (1) interface detection is now defined solely
by `edge_splits` and active-cell incidence — there is no derived cell level,
and `parent_cells` is documented as opaque P0-P2 provenance (§1, §3); (2)
coarsening selects removed midpoints via the private per-transition
`introduced_midpoints` map, so reused interface midpoints survive (§2, §8, §9).
See [decision-2026-11-21.md](decision-2026-11-21.md) for rationale and exact
implementation instructions.