# ADR-0003: Seam contracts for parallel implementation

- **Status:** Accepted
- **Date:** 2026-11-21
- **Deciders:** project maintainers

## Context

The reimplementation is split into seven acyclic seams (see
`docs/contracts/seams.md`) so that multiple workers can implement in parallel.
Parallel work requires frozen, concrete contracts: ambiguous wording here is a
defect. This ADR freezes the contracts for the first implementation slice
(`ReferenceCell` + `UpdateFlags` in `base`, then FEValues P1 triangle) and the
cross-cutting typed-ID and documentation-language decisions.

All contracts below target Rust 1.70, zero external crates, and
`#![forbid(unsafe_code)]`. They are intentionally concrete rather than generic.

## Decision

### 1. `ReferenceCell`

An enum in the `base` seam. No generics, no unsafe.

```rust
// src/base.rs (base seam)
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum ReferenceCell {
    Triangle,       // 2D
    Quadrilateral,  // 2D
    Tetrahedron,    // 3D
}

impl ReferenceCell {
    /// Spatial dimension: 2 for Triangle/Quadrilateral, 3 for Tetrahedron.
    pub fn dim(&self) -> usize;

    /// Number of vertices: 3, 4, 4.
    pub fn n_vertices(&self) -> usize;

    /// Number of faces (codimension-1 facets): 3 edges, 4 edges, 4 triangles.
    pub fn n_faces(&self) -> usize;

    /// Local vertex indices forming face `face_index` (0..n_faces()).
    /// Returns `None` when `face_index` is out of range.
    /// Slice length is 2 for 2D cells and 3 for Tetrahedron.
    pub fn face_vertex_indices(&self, face_index: usize) -> Option<&'static [usize]>;

    /// Reference vertex coordinates for 2D cells; `None` for Tetrahedron.
    pub fn reference_vertices_2d(&self) -> Option<&'static [Point2]>;

    /// Reference vertex coordinates for the 3D cell; `None` for 2D cells.
    pub fn reference_vertices_3d(&self) -> Option<&'static [Point3]>;
}
```

Fixed reference vertex order (must match existing code and tests):

- `Triangle`: `(0,0), (1,0), (0,1)`.
- `Quadrilateral`: `(0,0), (1,0), (1,1), (0,1)`.
- `Tetrahedron`: `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`.

Fixed face-to-vertex indices:

- `Triangle`: `[[0,1],[1,2],[2,0]]`.
- `Quadrilateral`: `[[0,1],[1,2],[2,3],[3,0]]`.
- `Tetrahedron` (face `i` is opposite vertex `i`): `[[1,2,3],[0,2,3],[0,1,3],[0,1,2]]`.

Notes for the implementer:

- Return `&'static` slices into `const` arrays; no allocation, no unsafe.
- `Point2`/`Point3` are the existing types in `src/base.rs`.
- No trait abstraction is introduced; `ReferenceCell` is a plain enum.
- Orientation and validity of a *physical* cell are the grid seam's concern;
  `ReferenceCell` only describes the reference topology.

### 2. `UpdateFlags`

A simple concrete type in the `base` seam. No external crate, no bitflags
dependency. Normals are **not** in scope now (see below).

```rust
// src/base.rs (base seam)
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct UpdateFlags {
    pub values: bool,
    pub gradients: bool,
    pub jxw: bool,
}

impl UpdateFlags {
    pub const EMPTY: Self = Self { values: false, gradients: false, jxw: false };
    pub const ALL:   Self = Self { values: true,  gradients: true,  jxw: true };

    pub fn empty() -> Self { Self::EMPTY }
    pub fn all()   -> Self { Self::ALL }

    pub fn is_empty(&self) -> bool {
        !self.values && !self.gradients && !self.jxw
    }
}
```

Normals: **out of scope for this slice.** Normal vectors are needed for
Neumann/Robin boundary assembly, but in the current architecture they come from
`BoundaryEdgeGeometry` (grid seam), not from cell FEValues. A future additive
flag `normals` may be added when a face-evaluation seam is introduced; adding it
is a backwards-compatible change. Until then, no worker may assume
`UpdateFlags` carries normals.

Rule: an evaluation computes **only** requested quantities. If a flag is false,
the corresponding output field is `None` and the computation is skipped. Skipping
must be observable in tests (assert the field is `None`).

### 3. FEValues evaluation contract (P1 triangle first slice)

The first implementation slice is FEValues for the P1 triangle. The contract
below is written for that slice; generalisation to other cells follows the same
shape and is additive.

Mapping: affine, as in the existing `AffineTriangle` (`src/fe.rs`). Constant
Jacobian `J` (2x2), `det J` constant, `J^{-T}` constant.

For each quadrature point `q` on the reference triangle, with reference point
`p_q`, reference weight `w_q`, physical point `x_q = map(p_q)`:

- Shape values (flag `values`): `phi[i][q]` for `i in 0..3`, where `phi` are
  the P1 basis functions `1 - x - y, x, y` evaluated at `p_q`.
- Reference gradients (flag `gradients`): `grad_ref[i][q] = dphi_i / d(ref coords)`.
- Physical gradients (flag `gradients`): `grad_phys[i][q] = J^{-T} grad_ref[i][q]`.
- JxW (flag `jxw`): `JxW[q] = |det J| * w_q`.

Mandatory vs flag-gated:

| Quantity | When computed |
|---|---|
| Input validation (finite coords, valid cell, finite positive weights) | always, including empty flags |
| Shape values `phi` | only if `values` |
| Reference gradients `grad_ref` | only if `gradients` |
| Physical gradients `grad_phys` | only if `gradients` |
| JxW | only if `jxw` |
| Physical quadrature point `x_q` | not flag-gated; derivable from `map(p_q)`; may be cached |

Invariants that tests must check (for P1 triangle):

1. Partition of unity: `sum_i phi[i][q] == 1` for every `q` (within `1e-12`).
2. Nodal property: `phi[i](vertex_j) == delta_ij` for reference vertices.
3. Sum of reference gradients is zero: `sum_i grad_ref[i] == 0`.
4. `JxW[q] > 0` on any valid (counter-clockwise, non-degenerate) cell.
5. Physical gradients equal `J^{-T}` applied to reference gradients.
6. `JxW[q] == |det J| * w_q`.
7. Skipping: if a flag is false, the corresponding field is `None`.

Physical quadrature points are not part of `UpdateFlags` in this slice; callers
use the mapping. A future `update_quadrature_points`-style flag is additive.

### 4. Typed IDs and the `CellId` decision

Typed IDs remain index-based, exactly as in `src/grid.rs`:

```rust
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VertexId(pub usize);
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CellId(pub usize);
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BoundaryEdgeId(pub usize);
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DofId(pub usize);
```

**Decision: keep index-based IDs now.** They are valid only for the mesh/DoF
handler instance they came from. They are invalidated by refinement, renumbering,
or any structural change; no automatic remapping is provided. Persistent or
path-based IDs (e.g. refinement-tree paths) are **deferred to the adaptivity /
distributed seam** and are out of P0-P2. No parallel worker may assume
persistent IDs exist.

Operations that take IDs must validate bounds and return `Result` (see existing
`cell_points`, `boundary_edge_geometry`). IDs are `Copy` and carry no lifetime.

### 5. Documentation language

New repository-facing text is written in English. Existing French documentation
(`README.md`, `PLAN_IMPLEMENTATION.md`, `docs/*.md` in French) is left intact;
there is no unprompted mass translation. ADRs and RFCs follow the workspace
`docs_language` (English unless set otherwise). Code comments and repository
documentation are in English.

## Consequences

- Workers implement against the exact signatures above; deviations require an
  ADR/RFC entry.
- The base seam can be implemented and tested independently of grid/fe.
- Typed-ID code must not assume persistence; this is explicit in doc comments.

## Alternatives considered

- Generic `ReferenceCell<DIM>`: rejected (no second use yet; YAGNI).
- External `bitflags` crate for `UpdateFlags`: rejected (zero-dep policy).
- Persistent/path-based IDs now: rejected (adaptivity seam not stable; would
  over-constrain grid).
- Including normals in `UpdateFlags` now: rejected (no face-evaluation seam yet).