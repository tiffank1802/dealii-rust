# ADR-0001: Parity scope for the dealii-rust reimplementation

- **Status:** Accepted
- **Date:** 2026-11-21
- **Deciders:** project maintainers

## Context

`dealii-rust` is an experimental, serial, safe Rust finite-element library inspired
by deal.II. deal.II is a large, mature, distributed, multi-element, multi-physics
library. Attempting full parity would dominate every other goal and would force
premature abstractions. The project has already shipped useful P0-P2 workflows
(P1 triangles, Q1 quadrilaterals, P1 tetrahedra, conforming refinement, basic
constraints) without claiming deal.II equivalence.

The maintainer has confirmed the parity target below. This ADR records it so
parallel workers share one definition of "done".

## Decision

The parity target is **scenario equivalence for the P0-P2 core workflows, plus
concept/API parity where it is cheap**, with **full parity explicitly out of
scope**.

Concretely "in scope":

- Scalar stationary Poisson (and symmetric diffusion-reaction) on conforming
  meshes of triangles, quadrilaterals and tetrahedra, with affine or
  isoparametric bilinear mappings.
- Dirichlet and Neumann boundary conditions, including non-homogeneous Dirichlet
  and pure-Neumann nullspace handling.
- Conforming (red) refinement with transfer of P1/Q1 fields, plus one-level
  non-conforming (hanging-node) refinement for triangles.
- DoF numbering, COO to CSR assembly, CG solve, VTK export, limited Gmsh import.
- Concept and API parity with deal.II where the cost is low (names, update-flag
  style, reference-cell topology, JxW convention), so that users familiar with
  deal.II can transfer knowledge.

Explicitly **out of scope for P0-P2** (and therefore not required for any
acceptance criterion in this milestone):

- Distributed/MPI execution, ghost data, partitioning.
- `hp` adaptivity and mixed-order constraints.
- Non-matching grids and mortar methods.
- Particles, particle-mesh coupling.
- Matrix-free operators and accelerators.
- All element families beyond P1/Q1 (higher order, H(div)/H(curl), etc.).
- Generic `const DIM` abstractions until two concrete uses justify them.
- Binary compatibility or class-for-class translation of the deal.II API.

These capabilities remain **optional** and may be added later under their own
seams and ADRs/RFCs; nothing in P0-P2 may depend on them.

## Consequences

**Positive**

- Abstractions stay concrete and small; the zero-dependency, `forbid(unsafe_code)`
  policy remains achievable.
- Acceptance criteria are testable with analytic checks and a single targeted
  deal.II comparison (Q1 unit square), instead of a full reference suite.
- Parallel workers do not block on capabilities that are not needed yet.

**Negative / accepted risks**

- Some deal.II idioms will not be reproduced; divergence must be documented in
  module docs rather than hidden.
- Later addition of deferred capabilities may require breaking API changes; those
  must follow the ownership rules in `docs/contracts/seams.md`.

**Re-evaluation trigger**

Re-open this decision (via RFC) when any of the following occurs:

1. A P0-P2 workflow cannot be completed without a deferred capability.
2. The deal.II comparison (T-004) reveals a correctness gap attributable to a
   deferred feature.
3. A second independent implementation needs an abstraction that the current
   concrete types cannot express.
4. The project moves beyond P2.
