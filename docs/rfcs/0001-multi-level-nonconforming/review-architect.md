# Review: RFC-0001 — architect

## Round 4 — 2026-11-21 — reviewed revision 4 — verdict: approve

### Summary

Revision 4 folds in the non-blocking clarifications called out in Round 3
(architect residuals 1, 3, 5 and engineer residuals 1–3) without changing any
design decision. The data model, signatures, ADR-0003 alignment, and
implementability conclusions from Round 3 still hold. No new seam violations,
ADR-0003 contradictions, or implementation blockers were introduced. The items
below that remain open are clarifications only; they can be folded in during
implementation without touching the RFC's decisions.

### Resolution check of Round 3 residual clarifications

| # | Round 3 item (source) | Status in rev 4 | Notes |
|---|---|---|---|
| A1 | `edge_splits` renumbering on coarsen (architect #1) | Resolved | §8 now states surviving entries are "renumbered with the same vertex index map" alongside vertices, in addition to the drop/cascade rules. |
| A2 | `edge_splits` level intentionally omitted (architect #2) | Not addressed | Still only implicit in the type signature (`BTreeMap<(VertexId, VertexId), VertexId>`). Recommend one sentence in §2 or §8; non-blocking. |
| A3 | Cell level definition (architect #3) | Resolved (minor) | §1 now defines a cell's level as its depth from the root mesh (root = 0), with a transition's level = parent + 1. Recommend also stating it is derived by counting hops to a self-parent in `parent_cells` (no stored field) and that "stored per transition" refers to the operation level. |
| A4 | Pre-existing interface persistence during `refine` (architect #4) | Partially | §3's first sentence defines interface by level difference, which covers pre-existing interfaces, but step 2's parenthetical "(one refined this step, one not)" reads as the only case and would let an implementer drop pre-existing hanging nodes. Reword step 2 to apply the level-difference test to every edge (see residual 1). |
| A5 | `boundary_parents` naming/type ambiguity (architect #5) | Partially | §1's note now distinguishes mesh-level ownership from the transition-level `parent_boundary_edge` map. But §9 still lists `boundary_parents` without qualification, and the transition-level field is not named, leaving a collision risk with the mesh-level `Vec<Option<CellId>>`. Also the stated consumer (`boundary_edge_geometry`) does not use an owning cell in the current code; justify or drop the field. |
| A6 | §9 retained-field list includes `children` (architect #6) | Not addressed | §9 still omits `children`/parent→child links even though §1 says the transition retains them. Add for implementer clarity. |
| A7 | Empty-marking deviation note (architect #7) | Resolved/accepted | §11 unchanged and consistent with `RefinementHistory::refine`. |
| E1 | `children()` served from transition (engineer #1) | Not addressed | Same as A6; §9 should list parent→child links. |
| E2 | `DofHandler::from_vertex_count` representation (engineer #2) | Not addressed | §6 fixes scope but not storage (owned count vs. `Option<&Mesh2>`). Minor; implementer decision. |
| E3 | `mesh()` as single view entry point (engineer #3) | Resolved | §1 (`mesh() -> &NonconformingMesh`) and §9 (`mesh()`/`mesh_at()`) consistently present. |

### New issues introduced by revision 4

None. The revision-4 wording is consistent with the Round 3 decisions:
- `boundary_parents` note keeps the mesh-level/transition-level distinction explicit;
- the cell-level definition is consistent with `parent_cells` (self for roots) and §10 (no per-cell stored level);
- the §8 renumbering/cascade clause is consistent with survivor ordering and the
  "pre-existing interface edges persist as long as their midpoints survive" clause;
- `from_vertex_count` scope, `mesh()` accessor, and empty-marking policy are
  unchanged and still consistent.

The open points below are pre-existing residuals, not regressions introduced by
revision 4.

### Residual clarifications (non-blocking)

1. **Pre-existing interface persistence in §3 (carried, A4).** Step 2 should read
   "an edge with exactly 2 incident active cells at *different levels*" and drop
   the "(one refined this step, one not)" parenthetical, or add a sentence that
   pre-existing interfaces (levels already differ, neither cell marked this step)
   are re-tested against `edge_splits` and their `HangingNode`s are carried
   forward. As written, a strict reading of step 2 would silently drop pre-existing
   hanging nodes. Non-blocking; the first sentence already gives the correct rule.
2. **`boundary_parents` field naming and necessity (carried, A5).** Name the
   transition-level parent-edge map explicitly (e.g. `parent_boundary_edges:
   Vec<BoundaryEdgeId>`) so it cannot be confused with the mesh-level
   `boundary_parents: Vec<Option<CellId>>`, and list only the transition-only
   fields in §9. Also either identify a concrete consumer for the owning-cell
   field or remove it (YAGNI); `boundary_edge_geometry` as it exists needs
   orientation, not an owning cell.
3. **§9 retained-field list (carried, A6/E1).** Add `children` / parent→child
   links to the list; they are referenced by §1 but absent from the enumeration.
4. **Cell-level derivability (carried, A3).** State that the level is computed by
   counting hops to a self-parent in `parent_cells`, with no stored field, and that
   "stored per transition" refers to the operation level in
   `NonconformingRefinementHistory`.
5. **`edge_splits` level omission (carried, A2).** Add one sentence in §2/§8
   that the split level is deliberately not stored because reuse is key-based and
   coarsening drops by midpoint.
6. **`from_vertex_count` storage (carried, E2).** Note whether the handler stores
   an owned count or `Option<&Mesh2>`; the scope contract (no `cell_dofs`/
   `boundary_dofs`) is already fixed.

### Verdict

**approve.** Revision 4's clarifications are consistent with the Round 3 design
and introduce no new issues or blockers. The six residual items above are
non-blocking clarifications (four carried from Round 3, two wording/precision
notes) that can be folded into the text during implementation; none requires a
change to the RFC's decisions, signatures, or ADR-0003 alignment. Implementation,
tests, and the adaptive Poisson example may proceed; P0–P2 examples/tests remain
untouched per ADR-0001.

## Round 3 — 2026-11-21 — reviewed revision 3 — verdict: approve

### Summary

Revision 3 resolves all four Round 2 blocking items and addresses (or explicitly
defers) the Round 2 non-blocking clarifications. The data model is now internally
consistent: `NonconformingMesh` is the single source of truth for the active mesh,
`parent_cells` keeps the existing `Vec<CellId>` convention (self for roots),
`edge_splits` is the cross-level midpoint source of truth with vertex reuse stated,
and derefinement now guards geometric vertex references. The direction (Option C,
additive, index-based IDs) remains sound and no new seam violations or ADR-0003
contradictions were introduced. Residual items below are clarifications, not
blockers; they can be folded in during implementation.

### Resolution check of Round 2 blocking items

| # | Round 2 item | Status in rev 3 | Notes |
|---|---|---|---|
| 1 | Discovery rule for interface/hanging edges | Resolved | §3 now defines edge→cell incidence map; interface iff the two sharing active cells are at different levels; one `HangingNode` per midpoint; multi-midpoint edges tested per sub-edge key. |
| 2 | Container split / transition metadata | Resolved | §1 `NonconformingMesh` carries `midpoint_parents`, `boundary_parents`, `boundary_vertices`, `original_vertices`; §9 transitions retain `oriented_boundary_edges` (and `children`); `mesh()`/`mesh_at()` expose the mesh view. No `children` field in the mesh (active-only). |
| 3 | `parent_cells` `Vec<CellId>` vs `Option` | Resolved | §1 keeps `Vec<CellId>`, self for roots, matching `parent_cell(child) -> Result<CellId>` and existing tests. |
| 4 | Derefinement geometric vertex dependency | Resolved | §8 now requires (3) no surviving active cell references excised midpoints as ordinary vertices, and (4) no remaining `HangingNode` references them (slave or master). |

### Resolution check of Round 2 non-blocking clarifications

| # | Round 2 item | Status in rev 3 | Notes |
|---|---|---|---|
| 5 | `boundary_edges` field type | Resolved | §1: `Vec<[VertexId; 2]>` (geometry); interface edges never boundary. |
| 6 | `edge_splits` cascade on coarsen | Resolved | §8: drop entries whose midpoint removed **and** any entry whose key contains a removed vertex. |
| 7 | Empty-marking policy | Addressed (deviation) | §11: `new` accepts empty; `refine`/`coarsen` reject empty. Consistent with `RefinementHistory::refine`; differs from Round 2's softer recommendation (refine accepts empty) but is documented and defensible. |
| 8 | `DofHandler::from_vertex_count` scope | Resolved | §6: `cell_dofs`/`boundary_dofs` unavailable for nonconforming meshes via this constructor; callers use vertex-indexed numbering directly. |
| 9 | Wrapper accessor | Resolved | §1: `mesh() -> &NonconformingMesh` plus delegation of existing accessors. |
| 10 | Seam leak (carried over) | Resolved/accepted | §12 keeps assembly/dofs convenience methods as documented grid-owned wrappers. |

Engineer Round 2 items: `children` contradiction resolved (field removed from mesh, §1); clarifications A–D: B and C resolved; A (`edge_splits` renumbering) and D (level intentionally omitted) partially addressed (see residual).

### Residual clarifications (non-blocking)

1. **`edge_splits` renumbering on coarsen (engineer A).** §8 compacts vertex indices and says surviving entries are retained, but does not explicitly state that the vertex keys/values of retained `edge_splits` entries are remapped by the same compaction map. Add one sentence so `coarsen`'s output is self-consistent and `coarsen_invalidates_ids_documented` is testable.
2. **`edge_splits` level intentionally omitted (engineer D).** State explicitly (as §2/§8 note) that the split level is deliberately not stored, since reuse is key-based and coarsening drops by midpoint. Avoids re-asking.
3. **Cell level definition.** §3 compares cells "at different refinement levels" but `NonconformingMesh` stores no per-cell level. Define level as parent-chain depth (hops to a self-parent in `parent_cells`); the "one refined this step, one not" operational rule then follows. Needed for implementability of pre-existing interface persistence.
4. **Pre-existing interface persistence during `refine`.** §3's operational criterion ("one refined this step, one not") catches newly created interfaces; ensure `refine` also preserves `HangingNode`s for interface edges whose cells are not marked this step (their midpoints remain in `edge_splits`). Otherwise a refine pass could silently drop pre-existing hanging nodes.
5. **`boundary_parents` naming/type ambiguity (new).** §1 types `boundary_parents: Vec<Option<CellId>>` ("owning active cell, or None if orphaned") and says it is needed by `boundary_edge_geometry` and `parent_boundary_edge`. But `parent_boundary_edge` currently returns a `BoundaryEdgeId` (parent edge), and `boundary_edge_geometry` does not use an owning cell. The transition also retains a `boundary_parents` (§9) needed for `parent_boundary_edge`. Reconcile: keep `boundary_parents` as the parent-edge map (`Vec<BoundaryEdgeId>`, possibly `Option` for orphaned) for `parent_boundary_edge`, and name the owning-cell field separately (e.g. `boundary_owners: Vec<Option<CellId>>`), or drop the owning-cell field if unused. As written, an implementer cannot tell which meaning applies.
6. **§9 retained-field list.** Include `children` (parent→child links) explicitly in the transition-retained list; §1 references it ("transition type retains parent→child links") but §9's enumeration omits it. Cosmetic.
7. **Empty-marking deviation note.** §11's choice (refine/coarsen reject empty) is consistent with `RefinementHistory::refine` and is accepted; no change required, just documenting the resolution of Round 2 item 7.

### Verdict

**approve.** All four Round 2 blocking items are resolved and the design is
implementable on the existing substrate (`AffineConstraints::close`,
`NonconformingRefinement`, `RefinementHistory`). The residual items above are
clarifications to fold in during implementation; none requires a change to the
RFC's decisions or signatures. Implementation may proceed.

## Round 2 — 2026-11-21 — reviewed revision 2 — verdict: request-changes

### Summary

Revision 2 is a substantial improvement and directly addresses all eight blocking
items from Round 1. The direction (Option C, additive, index-based IDs) remains
sound and no new seam violations or ADR-0003 contradictions were introduced
beyond the ones noted below. However, four items are still not implementable as
written, and several field-level inconsistencies would force an engineer to guess
at the data model. Verdict remains **request-changes**; the items in §Blocking
must be resolved before implementation, the §Non-blocking items can be folded in
opportunistically.

### Resolution check of Round 1 blocking items

| # | Round 1 item | Status in rev 2 | Notes |
|---|---|---|---|
| 1 | Input representation / interface-edge metadata | Resolved | `NonconformingMesh` + `edge_splits: BTreeMap<(VertexId,VertexId), VertexId>` is the single source of truth; `edge_key` reuse stated. |
| 2 | History storage | Resolved | Dedicated `NonconformingRefinementHistory`, additive, mirrors `RefinementHistory`; no overload. |
| 3 | DoF contract | Resolved | Vertex-indexed DoFs retained; `AffineConstraints` composition; additive `DofHandler::from_vertex_count`. |
| 4 | Active-cell model | Resolved | `cells()` active only; no active flag; parallel provenance vectors. |
| 5 | Hanging-node detection on nonconforming input | Partially | `edge_splits` + rebuild rule given, but the discovery rule is not well-defined (see Blocking 1). |
| 6 | Multi hanging nodes per edge / 2:1 | Resolved | Unbalanced allowed; one `HangingNode` rule per midpoint; `Vec` permits multiple midpoints per geometric edge. |
| 7 | Derefinement validity / API shape | Partially | Rejection conditions, immutability, one-level-per-call, renumbering, `edge_splits` cleanup stated, but geometric vertex dependency missing (see Blocking 4). |
| 8 | Cross-level provenance scoping | Resolved | Level-scoped IDs, history retains levels, `edge_splits` is the only cross-level structure (vertex indices only). |

### Blocking issues (must resolve before implementation)

1. **Discovery rule for interface/hanging edges is not well-defined.**
   `README.md` §3 states an edge is "interface" if "its two endpoint vertices
   belong to cells at different levels". A vertex can belong to several cells at
   several levels, so "belongs to a level" is ambiguous, and the criterion does not
   match the actual condition (the two *cells sharing the edge* are at different
   refinement levels). Rewrite the rule in terms of the two cells adjacent to the
   edge, and state how edge→cell incidence is built during `refine` (it is not
   derivable from `parent_cells`/`children` alone). Also specify detection for
   edges carrying more than one hanging node (level difference > 1): each sub-edge
   key in `edge_splits` is tested independently.

2. **Container split loses transition metadata.**
   `NonconformingMesh` (§1) carries `vertices, cells, parent_cells, children,
   edge_splits, hanging, boundary_edges` but omits `midpoint_parents`,
   `oriented_boundary_edges`, `boundary_parents`, `boundary_vertices` and
   `original_vertices`, all of which the current `NonconformingRefinement` needs
   for `prolongate_p1`, `boundary_edge_geometry`, `parent_boundary_edge`,
   `boundary_vertices` and `interpolate_boundary`. Either fold these into
   `NonconformingMesh`, or state that `NonconformingRefinement` (the transition)
   retains them and that `NonconformingRefinementHistory` stores transitions
   (mirroring `RefinementHistory` storing `LocalRefinement`, not bare meshes) with
   `mesh()`/`mesh_at()` views. Reconcile §1 vs §9 so a retained level is
   self-sufficient.

3. **`parent_cells: Vec<Option<CellId>>` conflicts with the existing accessor and
   tests.**
   The current type uses `Vec<CellId>` with a self-parent for unrefined cells, and
   `parent_cell(child) -> Result<CellId>` is exercised by existing tests for
   every active cell. Switching to `Option` changes a frozen-ish accessor and
   breaks additive-only / untouched P0-P2 tests. Either keep `Vec<CellId>`
   (self for roots, matching current behaviour) or explicitly specify the new
   `parent_cell() -> Result<Option<CellId>>` signature and the migration. Recommend
   keeping `Vec<CellId>`.

4. **Derefinement validity omits geometric vertex dependency.**
   §8 allows restoring a parent only if all children are marked, have no
   descendants, and no neighbour depends on their midpoints as a hanging master.
   It does not forbid a surviving active cell from still referencing one of those
   midpoints as an ordinary vertex; removing it would then leave a dangling
   `VertexId`. Add the condition that every midpoint being excised must not be
   referenced by any surviving active cell or by any remaining `HangingNode`
   (as slave or master). Otherwise the coarsened mesh is invalid.

### Non-blocking issues (clarifications / suggestions)

5. **`boundary_edges: Vec<BoundaryEdgeId>` in §1 reads as a type error.** Boundary
   edges must carry geometry (`Vec<[VertexId; 2]>`) plus the parent mapping; a
   `Vec<BoundaryEdgeId>` is self-referential. Clarify the field type.

6. **`edge_splits` cleanup on coarsen is incomplete.** §8 drops entries whose
   midpoint was removed; entries whose key contains a removed vertex must also be
   dropped (cascade), otherwise stale sub-edges survive.

7. **Empty-marking policy is inconsistent.** `NonconformingRefinement::new`
   allows an empty marking (used by existing tests), but `refine` rejects empty.
   Reconcile: recommend the history layer rejects empty (to avoid indistinguishable
   levels) while `refine` itself accepts empty for API consistency with `new`.

8. **`DofHandler::from_vertex_count` scope.** As specified it only yields a count;
   `cell_dofs`/`boundary_dofs` need cells/boundary data. State whether those are
   unavailable for nonconforming meshes, or add a nonconforming-aware constructor
   taking `&NonconformingMesh`.

9. **Wrapper accessor.** Since `NonconformingRefinement` wraps `NonconformingMesh`,
   specify the accessor (`mesh() -> &NonconformingMesh`) and delegation of the
   existing `vertices`/`cells`/... methods, and that `NonconformingRefinementHistory`
   stores transitions with `mesh()`/`mesh_at()` views.

10. **Seam leak (carried over).** `assemble_poisson`/`hanging_constraints` on the
    grid type is pre-existing; §11 keeps it as a grid-owned wrapper. Acceptable,
    but keep the convenience methods thin and documented as such.

### Verdict

**request-changes.** The eight Round 1 blocking decisions are all present in rev 2
and the chosen option (C) remains the right one, but Blocking items 1–4 must be
resolved in the RFC text before an engineer can implement without guessing at the
data model or producing an invalid coarsened mesh. Once those are written in,
this can move to approve.

## Round 1 — 2026-11-21 — reviewed revision 1 — verdict: request-changes

### Summary

The direction (Option C: multi-level, index-based IDs, reuse existing patterns) is
acceptable and is the only option of the three that stays inside ADR-0003. Option B
would violate the deferred persistent-ID decision and Option A discards the
provenance the RFC itself asks for. However, revision 1 leaves the core data model
unspecified in ways that cannot be implemented against the frozen seam contracts
without guessing. The issues below are blocking: they each require a decision in
the RFC before an engineer can write code without inventing an API.

### Design coherence with the seam model

- **Grid seam:** `NonconformingRefinement` is correctly a grid-seam type and owns
  mesh topology, boundary edges, and hanging-node discovery
  (`src/nonconforming_refinement.rs`). Keeping multi-level refinement there is
  consistent with `docs/contracts/seams.md` (grid owns refinement). The proposal
  to compose constraints via `AffineConstraints::close` is also consistent with
  the dofs seam (`src/constraints.rs`).
- **Dofs seam:** the RFC says it "touches `src/dofs.rs`" but never states how
  `DofHandler` participates. `DofHandler` is hard-bound to `&Mesh2`
  (`src/dofs.rs:18-25`) and returns `DofId(vertex index)`. Nonconforming meshes
  are *not* `Mesh2` and there is deliberately no conversion
  (`src/nonconforming_refinement.rs:36-37`). The RFC must define DoF numbering for
  multi-level nonconforming meshes (recommendation: keep vertex-indexed DoFs and
  build constraints over them, as the one-level code already does) and state the
  invalidation rule on derefinement.
- **Assembly seam:** `NonconformingRefinement::assemble_poisson`
  (`src/nonconforming_refinement.rs:290-316`) already lets a grid type perform
  assembly, which is a pre-existing seam leak. The RFC expands the surface that
  assembly depends on (active-cell view, hanging-node sets). It should not make
  this worse: either keep the convenience method as a thin wrapper documented as
  grid-seam-owned, or move assembly ownership fully into the assembly seam. This
  is not blocking, but must be decided before the multi-level API freezes.
- **Additive-only (ADR-0003 §3):** the proposal is implementable additively only
  if existing signatures survive. In particular `NonconformingRefinement::new(parent: &Mesh2, marked)`
  (`src/nonconforming_refinement.rs:60`) must stay, with new entry points added
  (`refine`, `coarsen`) rather than changing the constructor to accept an enum or
  a nonconforming mesh.

### Does it respect ADR-0003 (index-based IDs, no persistent IDs, additive-only)?

- **Index-based IDs:** yes in intent. But the proposal's "parent cell id, child cell
  ids, level, and active flag" per cell (`README.md` §Proposal 1) risks becoming
  persistent/path-based identifiers if those IDs are retained across refinement
  levels. ADR-0003 §4 forbids that in P0-P2 and defers it to the adaptivity seam.
  The IDs in `parent_cells`/`children` (`src/nonconforming_refinement.rs:47-48`)
  are valid only within one transition, as `LocalRefinement` already documents
  (`src/grid.rs:225-230`). The RFC must state that provenance is **level-scoped**
  and stored in a level-retaining history, not as long-lived cross-level IDs.
- **No persistent IDs:** the recommended fix (below) is to mirror
  `RefinementHistory` (`src/refinement_history.rs`), which retains every level and
  keeps IDs meaningful only at their level. That is the existing, accepted pattern.
- **Additive-only:** achievable, but only if the RFC explicitly commits to new
  methods/constructors rather than generalising `NonconformingRefinement::new`.

### Alternatives

- **Option A (rebuild each cycle):** simpler and actually matches the current
  immutable refinement types, but loses provenance and is expensive. Acceptable as
  a fallback if the RFC's data model cannot be pinned down; the RFC's rejection of
  A is reasonable for adaptive loops.
- **Option B (persistent cell IDs):** rejected correctly; violates ADR-0003 §4 and
  would force a breaking change to typed IDs across grid/dofs/assembly.
- **Option C (chosen):** right choice, but revision 1 under-specifies the mesh
  representation, history storage, and derefinement validity. Approve the
  direction, not the current text.

### Blocking issues

1. **`[src/nonconforming_refinement.rs:60]` — blocking — input representation.**
   `NonconformingRefinement::new` takes `&Mesh2`, but multi-level refinement needs
   to refine an already nonconforming mesh. `NonconformingRefinement` is not a
   `Mesh2`, and converting one to `Mesh2` is explicitly forbidden because
   `Mesh2::new` misclassifies interface edges as boundary edges
   (`src/nonconforming_refinement.rs:36-37`, `tests/nonconforming_refinement.rs:347-358`).
   The RFC must specify either a new `NonconformingMesh` type carrying interface-edge
   metadata, or a `refine(&self, marked)` method that operates on stored
   provenance. "Generalise `NonconformingRefinement` to accept a mesh that may
   already contain hanging nodes" is not implementable as stated.

2. **`[src/refinement_history.rs:15-18,56-62]` — blocking — history reuse.**
   `RefinementHistory` stores `Mesh2` levels and calls `LocalRefinement::new(&Mesh2)`.
   It cannot store `NonconformingRefinement` levels or produce nonconforming
   transitions. The RFC says "reuse the existing `RefinementHistory` pattern" but
   the existing type cannot hold this data. Required decision: introduce a
   dedicated `NonconformingRefinementHistory` (mirroring `RefinementHistory`) or
   generalise `RefinementHistory` over a mesh enum. The latter changes a frozen
   public type and needs an ADR entry; the former is additive and preferred.

3. **`[src/dofs.rs:18-25]` — blocking — DoF numbering contract.**
   `DofHandler` borrows `&Mesh2` only. The RFC touches `dofs.rs` but does not
   define DoF numbering for nonconforming meshes, nor invalidation on
   derefinement. Recommendation: keep DoFs vertex-indexed (one per vertex,
   including hanging nodes) and express hanging relations through
   `AffineConstraints`, exactly as `hanging_constraints` already does
   (`src/nonconforming_refinement.rs:242-251`). Any change to `DofHandler`'s
   signature must be additive (new constructor/method), not a replacement.

4. **`[src/nonconforming_refinement.rs:44-56]` — blocking — active vs. stored cells.**
   The current type stores only active cells in `cells` and derives parents via
   `parent_cells`. The RFC proposes a per-cell "active" flag and an
   "active cells" view (open question 3) but does not settle it. If inactive
   parents are retained in `cells`, then `assemble_poisson`
   (`src/nonconforming_refinement.rs:297-314`) and every caller iterating
   `cells()` will double-count. Recommendation: keep the flattened active-cell
   list as the single source of truth (current model), store provenance in
   parallel vectors, and expose `cells()` as active cells only. Drop the active
   flag from the cell record.

5. **`[src/nonconforming_refinement.rs:70-91,111-122]` — blocking — hanging-node
   detection on nonconforming input.** Discovery depends on a parent split map
   built from the marking. For an already nonconforming mesh there is no marking
   and no parent; interface edges must be classified from stored metadata. The RFC
   mentions "tracking interface edges" but defines no data structure, invariant,
   or how an edge's split level is determined. This must be specified, including
   the case of an edge carrying more than one hanging node (level difference > 1).

6. **`[src/nonconforming_refinement.rs:22-25,242-251]` — blocking — multiple
   hanging nodes per edge / 2:1 balance.** `HangingNode` represents a single
   midpoint with two masters, and `hanging_constraints` emits one midpoint rule.
   Pure red refinement applied repeatedly can create edges with several hanging
   nodes (quarter points, etc.) and the RFC does not state whether 2:1 balance is
   enforced. Without that decision, constraint generation and mesh validity are
   undefined. Recommendation: either (a) enforce 2:1 balance and document it as a
   precondition of `refine`, or (b) generalise `HangingNode` to a list of slaves
   along an edge and generate all midpoint rules. The RFC must pick one.

7. **`[src/nonconforming_refinement.rs:62-67]` (proposed `coarsen`) — blocking —
   derefinement validity and API shape.** "Remove children of a parent when all
   siblings are marked" is insufficient: a parent cannot be restored if a
   neighbouring cell depends on one of its midpoints as a hanging node, or if any
   child has descendants. The RFC must specify rejection conditions, whether the
   operation is immutable (returns a new mesh, recommended, matching
   `LocalRefinement`/`UniformRefinement`) or mutating, how vertex indices are
   renumbered, and whether one call may cross multiple levels (recommendation:
   one level per call, loop in the caller, matching `RefinementHistory::coarsen_last`).

8. **`[src/nonconforming_refinement.rs:47-48]` + ADR-0003 §4 — blocking —
   cross-level provenance.** Storing parent/child IDs that outlive a transition
   effectively creates persistent IDs, which ADR-0003 defers. The RFC's "per cell"
   record must be scoped to a single transition (as it already is via
   `parent_cells`/`children`) and cross-level provenance must live in the
   level-retaining history, not in the mesh. The RFC text should be corrected to
   say this explicitly.

### Risks

- **Constraint fill / memory growth** (`README.md` §Trade-offs): acknowledged but
  unquantified. Multi-level midpoint composition can produce rows with more
  nonzeros than the current two-master rules; deep levels may explode. Add a
  benchmark/test that records row width and peak memory for, e.g., 4 and 6 levels.
- **DoF invalidation on derefinement:** removing midpoint vertices renumbers all
  later vertices; all `VertexId`/`CellId`/`DofId`/`BoundaryEdgeId` become invalid.
  This is legal under ADR-0003 §4, but must be stated in the API docs of
  `coarsen`, exactly as `UniformRefinement` documents invalidation.
- **Seam leak growth:** `assemble_poisson` and `hanging_constraints` already put
  assembly/dofs behaviour on a grid type. Extending them to multi-level increases
  the leak. Consider keeping the convenience methods but documenting them as
  grid-owned wrappers, or moving them under the assembly seam.
- **No 2:1 balance:** unbalanced refinement can produce invalid or very poor
  meshes; if the RFC keeps "no quality guarantees" it must at least state the
  balance policy (or absence of one) so users know what they get.

### Open questions — recommendations

1. **Maximum depth.** Do not impose an API-level cap initially. Instead, test
   correctness up to a documented depth (recommend 4 levels in unit tests, 6 in
   benchmarks) and document observed constraint fill. Revisit if memory growth
   exceeds an agreed threshold.
2. **Derefinement across a level boundary in one call.** Recommend **one level per
   call**, matching `RefinementHistory::coarsen_last`
   (`src/refinement_history.rs:66-71`). Callers can loop. This keeps provenance and
   invalidation understandable.
3. **Active cells view.** Recommend the **flattened active-cell list** (current
   model): `cells()` returns only active cells; provenance lives in parallel
   vectors. Do not add an active flag to cell records and do not expose a filter
   over an inactive-inclusive list.

### Non-blocking suggestions

- Keep `NonconformingRefinement::new(parent: &Mesh2, marked)` unchanged and add
  `refine(&self, marked) -> Result<Self>` and `coarsen(&self, marked) -> Result<Self>`
  as additive methods, so existing callers and tests are unaffected.
- Name the new history type `NonconformingRefinementHistory` (or similar) rather
  than overloading `RefinementHistory`, to avoid changing a frozen public type.
- Add round-trip tests: refine → coarsen → mesh equality with the parent, and
  multi-level constraint composition resolving to level-0 vertices.
- State explicitly that the adaptive Poisson example is additive (new example
  file) and that existing P0-P2 examples/tests are not modified, per ADR-0001.

### Verdict

**request-changes.** The chosen option is sound, but revision 1 must resolve the
blocking items (input representation, history storage, DoF contract, active-cell
model, hanging-node detection, multi-node edges, derefinement validity, and
cross-level provenance scoping) before implementation. Once those decisions are
written into the RFC, a second review round can approve.