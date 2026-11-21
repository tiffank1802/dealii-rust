---
id: RFC-0001
status: draft
revision: 4
created: 2026-11-21
decided: null
author: carter-maestro
adr: null
plan: null
supersedes: []
ai-usage: ai-assisted
tags: [adaptivity, refinement, nonconforming, hanging-nodes]
---

# RFC-0001: Multi-level non-conforming refinement and selective derefinement

**Status:** Draft
**Revision:** 4
**Author:** carter-maestro

## Summary

Extend the existing one-level triangular non-conforming refinement
(`NonconformingRefinement`, ADR-0003/P0-P2 scope) to **multi-level**
non-conforming meshes with hanging nodes, plus **selective derefinement**, while
keeping the crate serial, zero-dependency and `forbid(unsafe_code)`. This moves
the project beyond the P0-P2 scope defined in ADR-0001 and therefore requires
this RFC. Revision 3 resolves all blocking items from the Round 1 and Round 2
reviews.

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
- `parent_cells: Vec<CellId>` — per active cell, its parent cell in the previous
  level; self (the cell's own id) for original/root cells. **Keeps the existing
  `Vec<CellId>` convention** (`src/nonconforming_refinement.rs:47`) so the
  `parent_cell(child) -> Result<CellId>` accessor and its tests are unchanged.
  `None` is **not** used; roots are self-parents.
- `edge_splits: BTreeMap<(VertexId, VertexId), VertexId>` — canonical (sorted)
  vertex-pair → midpoint vertex, for every split edge at every level. Single
  source of truth for "is this edge split, and by which vertex".
- `hanging: Vec<HangingNode>` — one entry per midpoint vertex (slave), with its
  two masters. An edge may carry several midpoints (level difference > 1).
- `boundary_edges: Vec<[VertexId; 2]>` — exterior edge vertex pairs only.
  Interface edges are **never** boundary edges.
- `boundary_parents: Vec<Option<CellId>>` — per boundary edge, the active cell
  that owns it (or `None` if orphaned). Needed by `boundary_edge_geometry`.
  Note: this is distinct from `parent_boundary_edge` (a transition-level field
  retained in `NonconformingRefinement` that maps a boundary edge to its parent
  level's boundary edge id); `boundary_parents` is mesh-level ownership only.
- A cell's **refinement level** is its depth from the root mesh (root cells =
  level 0). A transition's level = parent level + 1. Stored per transition in
  `NonconformingRefinementHistory`, not in `NonconformingMesh`.
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

### 3. Interface / hanging-edge discovery (blocking item 1, 5)

An edge (canonical vertex pair) is **interface** iff the two active cells that
share it are at different refinement levels. Concretely, during `refine`:

1. Build an edge→cell-incidence map: for each active cell, for each of its 3
   edges (canonical key), append the cell's id to that edge's list.
2. An edge with exactly 2 incident active cells at different levels (one
   refined this step, one not) is an interface edge.
3. For each interface edge, if its key is in `edge_splits`, its midpoint is a
   hanging node: emit one `HangingNode { slave: midpoint, masters: (v0, v1) }`.
   If an edge carries several midpoints (level difference > 1), each sub-edge
   key in `edge_splits` is tested independently and emits its own `HangingNode`.

Edges with 0 or 1 incident cells are boundary edges (never interface).
Edges with 2 incident cells at the **same** level are interior (no hanging node).

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
numbering directly). Existing `DofHandler::new(&Mesh2)` is unchanged.

### 7. Constraint composition (no recursion)

Multi-level constraints are built by a single `AffineConstraints::new(n_vertices)`
plus one `add(slave, 0.0, &[(m0, 0.5), (m1, 0.5)])` per `HangingNode`, then
one `close()`. `close()` (iterative Kahn sort, `src/constraints.rs:78-150`)
resolves chains automatically and detects cycles. Constraints are rebuilt from
scratch after every refinement/coarsening (they are not cached across
structural changes, matching the immutable refinement types).

### 8. Selective derefinement (`coarsen`)

- Input: `&self` + `marked: &[CellId]` (active cells to coarsen).
- A parent is restored **only if** all of the following hold:
  1. all its children are in `marked`;
  2. none of those children has descendants;
  3. no surviving active cell references any of the children's midpoints as an
     ordinary vertex;
  4. no remaining `HangingNode` (as slave or master) references any of those
     midpoints.
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
  persist as long as their midpoints survive.

### 9. History storage (blocking item 7)

Introduce a **dedicated** `NonconformingRefinementHistory` (additive), mirroring
`RefinementHistory` (`src/refinement_history.rs`): stores one
`NonconformingRefinement` (transition) per level, with `refine(marked)`,
`coarsen_last()`, `level()`, `mesh()` (current level) and `mesh_at(level)`.
The history stores **transitions** (which retain `midpoint_parents`,
`oriented_boundary_edges`, `boundary_parents`, `boundary_vertices`,
`original_vertices`), not bare meshes; `mesh()`/`mesh_at()` expose the
`NonconformingMesh` view. Do **not** overload `RefinementHistory`, which is
frozen to `Mesh2`/`LocalRefinement`.

### 10. Cross-level provenance scoping (blocking item 8)

Parent/child IDs are **level-scoped**: they are meaningful only within one
transition and are stored in the level-retaining history, not as long-lived
cross-level identifiers. This respects ADR-0003 §4 (no persistent/path-based IDs
in P0-P2). The `edge_splits` map is the only cross-level structure, and it
stores vertex indices, not cell paths.

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
  `coarsen_rejects_surviving_vertex_reference`,
  `prolongate_then_restrict_recovers_parent_field`,
  `coarsen_invalidates_ids_documented`.
- Deep-level: `deep_refinement_row_width_bounded`,
  `deep_refinement_completes_in_budget` (ignored by default if slow).
- Edge cases: already-split edge reuse, empty marking, non-finite midpoint,
  malformed rule rejection, cumulative overflow.

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

(Filled when accepted/rejected.)
EOF
echo done