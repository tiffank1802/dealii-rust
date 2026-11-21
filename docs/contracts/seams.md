# Seam contracts

This document defines the seven acyclic seams of `dealii-rust`, their ownership,
dependencies, existing files, frozen contracts, and the rules that keep parallel
work from colliding. It complements ADR-0001 (parity scope), ADR-0002
(licensing/provenance), and ADR-0003 (concrete contracts).

## Acyclic dependency direction

```
base  <--  grid  <--  dofs  <--  assembly
base  <--  fe    <--  dofs  <--  assembly
base  <--  lac   <-----------  assembly
grid, fe, dofs, lac  <--  io
```

`base` has no dependencies. `grid` and `fe` depend only on `base`. `dofs`
depends on `grid` and `fe`. `lac` depends on `base`. `assembly` depends on
`grid`, `fe`, `dofs`, `lac`. `io` depends on `grid`, `fe`, `dofs`, `lac` (and
may read assembly outputs via data, not via seam dependency). No seam may depend
on a seam to its right in the merge order below.

## Seam table

| Seam | Owns | Depends on | Existing `src/*.rs` files | Frozen contract (summary) |
|---|---|---|---|---|
| `base` | `Point2`, `Point3`, `ReferenceCell`, `UpdateFlags`, `Error`, numeric helpers (`signed_double_area`, `triangle_determinant`) | (none) | `base.rs` | Rust 1.70, zero deps, no `unsafe`; concrete types only; `ReferenceCell` and `UpdateFlags` per ADR-0003 |
| `grid` | Vertices, cells, connectivity, boundary edges/vertices, meshes (`Mesh2`, tetrahedral mesh), refinement (uniform, red-green, non-conforming), IDs (`VertexId`, `CellId`, `BoundaryEdgeId`) | `base` | `grid.rs`, `tetra.rs`, `local_refinement.rs`, `nonconforming_refinement.rs`, `quad_refinement.rs`, `tetra_refinement.rs`, `refinement_history.rs` | Index-based IDs with documented invalidation (ADR-0003); validation returns `Result`; no persistent IDs in P0-P2 |
| `fe` | Reference shape functions, quadrature rules, mappings (affine triangle, bilinear quad, tetrahedron), `FEValues` evaluation, field interpolation/projection primitives | `base` | `fe.rs`, `quad.rs`, `quadrature.rs` | FEValues contract per ADR-0003; flags gate computation; JxW positive on valid cells; partition of unity and nodal property |
| `dofs` | Global DoF numbering, local-to-global maps, boundary DoFs, affine constraints (`AffineConstraints`) | `grid`, `fe` | `dofs.rs`, `constraints.rs` | DoFs are index-based (`DofId`); constraints expressed as `x = C y + g`; cycles detected |
| `lac` | COO/CSR storage, matrix-vector product, CG (with report), preconditioners (Jacobi), vector operations | `base` | `lac.rs` | Defined duplicate-summing on compression; SPD precondition for CG documented; errors explicit |
| `assembly` | Local/global assembly of diffusion-reaction, Neumann handling, error estimation, field evaluation/projection, integrators | `grid`, `fe`, `dofs`, `lac` | `elliptic.rs`, `elliptic_quad.rs`, `elliptic_tetra.rs`, `integrators.rs`, `numerics.rs`, `field.rs`, `quad_field.rs`, `tetra_field.rs`, `estimator.rs`, `neumann.rs`, `neumann_quad.rs`, `neumann3d.rs`, `neumann_core.rs` | Split by cell family to reduce conflicts; sequential COO then CSR; Dirichlet eliminates rows and columns |
| `io` | VTK export, Gmsh import (4.1/2.2), comparison harness I/O | `grid`, `fe`, `dofs`, `lac` | `vtk.rs`, `gmsh.rs`, `gmsh22.rs` | Limited, validated subsets; provenance per ADR-0002; no mandatory native deps |

## Ownership rules

1. **One owner per seam.** A change to a seam's public (`pub`) API requires an ADR or RFC entry before merge.
2. **Internals are `pub(crate)`.** Cross-seam access goes through the seam's `pub` API only.
3. **Additive-only within a milestone.** New items may be added; existing signatures are not changed inside a milestone. Breaking changes land in a dedicated PR with an ADR/RFC and a version bump note.
4. **Merge order:** `base` -> `grid`/`fe` -> `dofs`/`lac` -> `assembly` -> `io`. Respect this order to keep diffs reviewable and conflicts low.
5. **Split `assembly` by cell family.** P1 triangle, Q1 quadrilateral, and P1 tetrahedron assembly are separate sub-tracks; they must not edit shared assembly helpers in the same PR without coordination.
6. **No premature traits.** Concrete types first; extract a trait only after two implementations or a real extension need it (per plan section 3).
7. **No new external crates** in P0-P2 without a decision recorded in `docs/` (per plan section 5).

## Mapping to role-tracks

| Track | Owner seam(s) | Focus |
|---|---|---|
| A | decisions + `base` | ADRs, `ReferenceCell`, `UpdateFlags`, numeric conventions |
| B | `fe` | FEValues, quadrature, mappings, shape functions |
| C | `grid` + `dofs` | Meshes, refinement, DoF numbering, constraints |
| D | `lac` | Sparse algebra, solvers, preconditioners |
| E | validation | Tests, analytic checks, comparison harness, benchmarks |
| F | `assembly` + `io` | Integrators, assembly, VTK/Gmsh, examples |

Track E (validation) is cross-cutting and may touch any seam's tests but not its
implementation. Track A (decisions) owns the ADRs and the seam contract document.