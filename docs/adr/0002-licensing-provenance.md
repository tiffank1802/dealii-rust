# ADR-0002: Licensing and provenance rules for dealii-rust

- **Status:** Accepted
- **Date:** 2026-11-21
- **Deciders:** project maintainers

## Context

The upstream deal.II checkout used for reference (`.cache/dealii`, pinned commit
`1d230ad0be5ab2a3913cf4de8187301430721ed5`) is dual-licensed under
Apache-2.0 with LLVM Exception and LGPL-2.1-or-later. `dealii-rust` is licensed
`MIT OR Apache-2.0` (see `Cargo.toml` and `LICENSE-MIT` / `LICENSE-APACHE`).

deal.II is used as a *reference* only: comparison programs, conceptual guidance,
and validation of numerical results. It is not a dependency and no deal.II code
is linked or compiled into the crate. The repository already contains a
comparison harness (`comparison/dealii_q1/`) whose `reference.cc` is derived
from deal.II's API but is authored independently.

## Decision

1. **Upstream license recorded.** deal.II: `Apache-2.0 WITH LLVM-exception` OR
   `LGPL-2.1-or-later` (dual). `dealii-rust`: `MIT OR Apache-2.0`.
2. **No verbatim reuse.** No verbatim code, quadrature tables or constants, or
   test data may be copied from deal.II (or any third party) without an SPDX
   provenance record and a license check. Algorithms and mathematics are not
   copyrightable and may be re-implemented from first principles.
3. **SPDX headers for derived files.** Any file derived from or transcribed
   from an upstream source (including `comparison/dealii_q1/reference.cc` and
   any future fixture or script) must carry an SPDX license identifier in its
   header naming the source and the applicable license. Files with no such
   provenance are treated as original work under `MIT OR Apache-2.0`.
4. **Fixture and script audit.** All test fixtures, reference data, mesh files,
   and helper scripts must be audited for provenance and license compatibility.
   Each must record its origin, license, and (where applicable) the exact
   upstream commit. Unprovenanced data must be removed or replaced before it is
   used in CI.
5. **Wording discipline.** Until the deal.II Q1 comparison (T-004) actually runs
   on a machine with deal.II installed, results are described as "analytically
   checked", never "validated against deal.II". The comparison harness
   (`comparison/dealii_q1/README.md`) already states this; it is binding.

## Consequences

**Positive**

- The crate remains `MIT OR Apache-2.0` with no copyleft contamination.
- Provenance is auditable; reviewers can distinguish original from derived
  files at a glance.
- Comparison results cannot be overstated.

**Negative / accepted risks**

- Some convenient numerical tables (quadrature rules, reference values) must be
  derived or regenerated rather than transcribed, which costs effort.
- The audit is ongoing work; until it completes, fixtures are assumed
  unprovenanced and must not be promoted to reference data.

**Alternatives considered**

- Copying deal.II tables/constants verbatim: rejected (license risk, and it
  would make the crate's provenance unauditable).
- Re-licensing `dealii-rust` to match deal.II: rejected; the project has
  already committed to `MIT OR Apache-2.0`.
- Treating the comparison harness as deal.II-validated before it runs: rejected;
  wording stays "analytically checked".