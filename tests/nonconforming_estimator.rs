//! RFC-0002 §7 validation suite for the nonconforming P1 residual indicator.
//!
//! The hanging-specific tests use crafted, constraint-satisfying vectors
//! (prolongation of a field sampled at original vertices) rather than solver
//! output, so every assertion is deterministic.
use dealii_rust::{
    base::Point2,
    estimator,
    grid::{
        CellId, Mesh2, NonconformingMesh, NonconformingRefinement, UniformRefinement, VertexId,
    },
};
use std::f64::consts::PI;

fn close(x: f64, y: f64, tolerance: f64) {
    assert!(
        (x - y).abs() <= tolerance,
        "{} != {} (tol {})",
        x,
        y,
        tolerance
    );
}

fn quadratic(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y
}

fn affine(p: Point2) -> f64 {
    1.0 + 2.0 * p.x - 3.0 * p.y
}

fn source(_p: Point2) -> f64 {
    -4.0
}

/// Nodal P1 values on a conforming mesh.
fn mesh_values<F: Fn(Point2) -> f64>(mesh: &Mesh2, field: F) -> Vec<f64> {
    mesh.vertices().iter().copied().map(field).collect()
}

/// Reconstructed values on a nonconforming level: prolong a field sampled on
/// original vertices so every midpoint (including sub-midpoints of an
/// unbalanced chain) equals the mean of its masters.
fn reconstructed_values<F: Fn(Point2) -> f64>(
    level: &NonconformingRefinement,
    field: F,
) -> Vec<f64> {
    let original: Vec<f64> = level.vertices()[..level.mesh().original_vertices]
        .iter()
        .copied()
        .map(field)
        .collect();
    level.prolongate_p1(&original).unwrap()
}

fn assert_hanging_rules(level: &NonconformingRefinement, values: &[f64]) {
    for h in level.hanging_nodes() {
        let expected = 0.5 * (values[h.masters[0].0] + values[h.masters[1].0]);
        close(values[h.vertex.0], expected, 1e-12);
    }
}

/// Locate the active cell containing both vertices, independent of numbering.
fn cell_with(level: &NonconformingRefinement, a: VertexId, b: VertexId) -> CellId {
    level
        .cells()
        .iter()
        .position(|cell| cell.contains(&a) && cell.contains(&b))
        .map(CellId)
        .unwrap_or_else(|| panic!("no active cell contains {:?} and {:?}", a, b))
}

fn assert_indicators_match(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (got, want) in actual.iter().zip(expected) {
        close(*got, *want, 1e-10);
    }
}

/// Marking every cell with `NonconformingRefinement::new` produces no
/// interfaces, so the new estimator must reproduce the frozen conforming
/// arithmetic exactly (per cell and in the global root-sum-square).
#[test]
fn nonconforming_estimator_matches_conforming_on_uniform_refinement() {
    let parent = Mesh2::unit_square(2).unwrap();
    let marked: Vec<CellId> = (0..parent.cells().len()).map(CellId).collect();
    let transition = NonconformingRefinement::new(&parent, &marked).unwrap();
    let mesh = transition.mesh();
    let uniform = UniformRefinement::new(&parent).unwrap();
    let reference = uniform.mesh();

    // The two refinements number their midpoints differently, so build each
    // solution vector from the coordinates of its own mesh.
    let values: Vec<f64> = mesh.vertices.iter().copied().map(quadratic).collect();
    let reference_values = mesh_values(reference, quadratic);

    let expected = estimator::poisson_dirichlet(reference, &reference_values, source).unwrap();
    let actual = estimator::poisson_dirichlet_nonconforming(mesh, &values, source).unwrap();

    assert_indicators_match(&actual.cell_indicators, &expected.cell_indicators);
    close(actual.global_indicator, expected.global_indicator, 1e-10);
}

/// Three uniformly refined levels: the residual indicator of `-Δu = -4` with
/// `u = x² + y²` scales like `h`, so each doubling of `n` halves `η` (the
/// documented band mirrors `tests/refinement.rs:170-174`).
#[test]
fn uniform_refinement_halves_global_indicator() {
    let mut mesh = Mesh2::unit_square(2).unwrap();
    let mut previous = f64::INFINITY;
    for _ in 0..3 {
        let level = NonconformingRefinement::new(&mesh, &[]).unwrap();
        let values = mesh_values(&mesh, quadratic);
        let eta =
            estimator::poisson_dirichlet_nonconforming(level.mesh(), &values, source).unwrap();

        assert!(eta.global_indicator.is_finite());
        assert!(eta.global_indicator > 0.0);
        if previous.is_finite() {
            let ratio = eta.global_indicator / previous;
            assert!(
                (0.25..0.65).contains(&ratio),
                "eta ratio = {} (previous = {}, current = {})",
                ratio,
                previous,
                eta.global_indicator
            );
        }
        previous = eta.global_indicator;
        mesh = UniformRefinement::new(&mesh).unwrap().mesh().clone();
    }
}

/// The primary RFC-0002 failure mode: a coarse/fine interface contributes a
/// real normal-derivative jump. Clearing the public split metadata makes
/// the same cells and values interface-blind; that estimate must be smaller.
#[test]
fn hanging_interface_jump_is_counted() {
    let parent = Mesh2::unit_square(2).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let mesh = level.mesh();
    assert!(!mesh.hanging.is_empty());

    // A curved field with alternating signs on the original vertices gives a
    // strong inter-cell gradient contrast; prolongation keeps every slave
    // continuous, so only the jump policy can separate the two estimates.
    let field = |p: Point2| (2.0 * PI * p.x).cos() * (2.0 * PI * p.y).cos();
    let values = reconstructed_values(&level, field);
    assert_hanging_rules(&level, &values);

    let aware = estimator::poisson_dirichlet_nonconforming(mesh, &values, source).unwrap();

    let mut blind_mesh = mesh.clone();
    blind_mesh.edge_splits.clear();
    blind_mesh.hanging.clear();
    blind_mesh.midpoint_parents.clear();
    let blind = estimator::poisson_dirichlet_nonconforming(&blind_mesh, &values, source).unwrap();

    assert!(
        aware.global_indicator > 1.01 * blind.global_indicator,
        "aware = {}, blind = {}",
        aware.global_indicator,
        blind.global_indicator
    );
}

/// Once both neighbours of the split edge are refined, the active cells form
/// a conforming mesh and the interface policy contributes nothing.
#[test]
fn hanging_interface_vanishes_after_both_sides_refined() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    assert_eq!(first.hanging_nodes().len(), 1);
    let h = first.hanging_nodes()[0];
    let neighbour = cell_with(&first, h.masters[0], h.masters[1]);
    let second = first.refine(&[neighbour]).unwrap();
    assert!(second.hanging_nodes().is_empty());

    let mesh = second.mesh();
    let conforming = Mesh2::new(mesh.vertices.clone(), mesh.cells.clone()).unwrap();
    let values: Vec<f64> = mesh.vertices.iter().copied().map(quadratic).collect();

    let expected = estimator::poisson_dirichlet(&conforming, &values, source).unwrap();
    let actual = estimator::poisson_dirichlet_nonconforming(mesh, &values, source).unwrap();
    assert_indicators_match(&actual.cell_indicators, &expected.cell_indicators);
    close(actual.global_indicator, expected.global_indicator, 1e-10);
}

/// An unbalanced two-level chain: the master of one slave is itself a slave.
/// The leaf-segment walk must handle the sub-midpoint and stay finite.
#[test]
fn unbalanced_chain_estimates_finite() {
    let parent = Mesh2::unit_square(2).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let h = first.hanging_nodes()[0];
    let fine_child = cell_with(&first, h.vertex, h.masters[0]);
    let second = first.refine(&[fine_child]).unwrap();

    let slaves: Vec<VertexId> = second
        .hanging_nodes()
        .iter()
        .map(|node| node.vertex)
        .collect();
    assert!(
        second
            .hanging_nodes()
            .iter()
            .any(|node| slaves.contains(&node.masters[0]) || slaves.contains(&node.masters[1])),
        "expected a cross-level hanging chain"
    );

    let values = reconstructed_values(&second, quadratic);
    assert_hanging_rules(&second, &values);

    let eta = estimator::poisson_dirichlet_nonconforming(second.mesh(), &values, source).unwrap();
    assert!(eta.global_indicator.is_finite());
    assert!(eta.global_indicator > 0.0);
    assert!(eta
        .cell_indicators
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0));
}

/// The estimator consumes the reconstructed vector; it rejects a value that
/// violates a hanging relation instead of silently repairing it.
#[test]
fn estimator_rejects_unreconstructed_hanging_values() {
    let parent = Mesh2::unit_square(2).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let values = reconstructed_values(&level, quadratic);
    assert!(estimator::poisson_dirichlet_nonconforming(level.mesh(), &values, source).is_ok());

    let mut unreconstructed = values.clone();
    unreconstructed[level.hanging_nodes()[0].vertex.0] += 0.1;
    assert!(
        estimator::poisson_dirichlet_nonconforming(level.mesh(), &unreconstructed, source).is_err()
    );
}

/// Invalid vectors, nonfinite sources and structurally invalid meshes are
/// rejected at the estimator seam (mirrors `tests/refinement.rs:181-191`).
#[test]
fn estimator_rejects_invalid_numerical_data() {
    let parent = Mesh2::unit_square(1).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let mesh: &NonconformingMesh = level.mesh();
    let zeros = vec![0.0; mesh.vertices.len()];

    assert!(estimator::poisson_dirichlet_nonconforming(mesh, &[0.0], source).is_err());
    assert!(estimator::poisson_dirichlet_nonconforming(
        mesh,
        &vec![f64::NAN; mesh.vertices.len()],
        source
    )
    .is_err());
    assert!(estimator::poisson_dirichlet_nonconforming(mesh, &zeros, |_| f64::NAN).is_err());
    assert!(estimator::poisson_dirichlet_nonconforming(mesh, &zeros, |_| f64::MAX).is_err());

    let mut invalid = mesh.clone();
    invalid.parent_cells.clear();
    assert!(estimator::poisson_dirichlet_nonconforming(&invalid, &zeros, source).is_err());
}

/// Bulk marking returns indices in the current active cell space, which
/// `refine` must accept directly; a zero fraction selects nothing.
#[test]
fn mark_bulk_on_nonconforming_indicators_returns_active_cell_ids() {
    let parent = Mesh2::unit_square(2).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0), CellId(6)]).unwrap();
    let values = reconstructed_values(&level, quadratic);
    let eta = estimator::poisson_dirichlet_nonconforming(level.mesh(), &values, source).unwrap();

    let marked = estimator::mark_bulk(&eta.cell_indicators, 0.5).unwrap();
    assert!(!marked.is_empty());
    for id in &marked {
        assert!(id.0 < level.cells().len(), "{:?} is not an active cell", id);
    }
    assert!(level.refine(&marked).is_ok());
    assert!(estimator::mark_bulk(&eta.cell_indicators, 0.0)
        .unwrap()
        .is_empty());
}

/// A globally affine P1 field has zero gradients everywhere; with `f = 0`
/// both the volume term and every interface jump vanish.
#[test]
fn affine_solution_has_zero_indicator() {
    let parent = Mesh2::unit_square(2).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let values = reconstructed_values(&level, affine);
    let eta = estimator::poisson_dirichlet_nonconforming(level.mesh(), &values, |_| 0.0).unwrap();

    for indicator in &eta.cell_indicators {
        assert!(*indicator <= 1e-12, "indicator = {}", indicator);
    }
    assert!(eta.global_indicator <= 1e-12);
}
