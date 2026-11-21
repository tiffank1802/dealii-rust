# Targeted comparison with deal.II: P1 on a triangular mesh

**Status: harness added; deal.II run pending on a machine with deal.II.** The
comparator can be exercised without deal.II via `analytic_smoke.txt`.

This case compares on `[0,1]^2` with a conforming P1 triangle mesh: the unit
square split into two triangles (vertices `(0,0), (1,0), (0,1), (1,1)`; cells
`[0,1,3]` and `[0,3,2]`), matching Rust's `Mesh2::unit_square(1)`.

Quantities compared:

- the global stiffness matrix for the Laplacian;
- the load vector of the L2 projection for `f = x^2 + y^2`;
- the unconstrained L2 projection coefficients;
- the L2 error of that projection.

Both programs use a degree-5 simplex quadrature (`QGaussSimplex(2, 3)` in
deal.II, `triangle_degree_5` in Rust). The output is ordered by Rust vertex
index `(0,0), (1,0), (0,1), (1,1)`. The expected file has four lines of five
numbers (stiffness row plus load entry), then one line with four projection
coefficients and the L2 error: 25 finite numbers in total. The comparator
checks a maximum absolute difference of `1e-10` for each quantity.

This case does not validate deformed meshes, boundary conditions, or global
convergence.

## Running on a machine with deal.II >= 9.4 and CMake

From the repository root:

```sh
cmake -S comparison/dealii_p1_tri -B comparison/dealii_p1_tri/build -DDEAL_II_DIR=/path/to/deal.II
cmake --build comparison/dealii_p1_tri/build --config Release
mkdir -p .cache
comparison/dealii_p1_tri/build/p1_tri_reference > .cache/dealii_p1_tri_reference.txt
cargo run --release --example compare_p1_tri_reference -- .cache/dealii_p1_tri_reference.txt
```

Record the exact deal.II, CMake, C++ compiler, Rust and system versions, the
command and its output, to keep the comparison traceable. If differences
appear, check first the vertex/cell ordering, the quadrature and the compile
options; do not alter results to hide a discrepancy.

## Testing the comparator without deal.II

```sh
cargo run --example compare_p1_tri_reference -- comparison/dealii_p1_tri/analytic_smoke.txt
```