//! Two-level adaptive nonconforming P1 Poisson on a triangle mesh.
//!
//! The marker is deterministic: a patch of cells is refined from the conforming
//! parent, then one isolated cell is refined a second time. The solve uses the
//! hanging-node constraints composed across both levels, and a final selective
//! coarsening restores that cell from its four children. This is an additive
//! example; the one-level `poisson_hanging` example is untouched.
//!
//! A separate, independent section then drives the RFC-0002 estimator loop:
//! solve → estimate → mark the bulk fraction → refine, with a fixed cycle
//! budget. Only refinement is indicator-driven; see `run_estimator_loop` for
//! why coarsening cannot be driven by an indicator alone.
use dealii_rust::{
    base::Point2,
    estimator::{mark_bulk, poisson_dirichlet_nonconforming},
    grid::{CellId, Mesh2, NonconformingRefinement, VertexId},
    lac::{conjugate_gradient, CgOptions},
};

fn exact(p: Point2) -> f64 {
    p.x * p.x + p.y * p.y
}

fn centroid(vertices: &[Point2], cell: &[VertexId; 3]) -> (f64, f64) {
    (
        (vertices[cell[0].0].x + vertices[cell[1].0].x + vertices[cell[2].0].x) / 3.0,
        (vertices[cell[0].0].y + vertices[cell[1].0].y + vertices[cell[2].0].y) / 3.0,
    )
}

fn edge_is_split(mesh: &NonconformingRefinement, cell: &[VertexId; 3]) -> bool {
    let view = mesh.mesh();
    (0..3).any(|j| {
        let a = cell[j];
        let b = cell[(j + 1) % 3];
        let key = if a < b { (a, b) } else { (b, a) };
        view.edge_splits.contains_key(&key)
    })
}

/// Pick the cell closest to `(target_x, target_y)` whose edges have not been
/// split at an earlier level. Refining such a cell introduces new midpoints
/// that no neighbour uses, so its sibling group can be coarsened again.
fn isolated_cell(mesh: &NonconformingRefinement, target_x: f64, target_y: f64) -> Vec<CellId> {
    let mut best: Option<(usize, f64)> = None;
    for (i, cell) in mesh.cells().iter().enumerate() {
        if edge_is_split(mesh, cell) {
            continue;
        }
        let (x, y) = centroid(mesh.vertices(), cell);
        let distance = (x - target_x).powi(2) + (y - target_y).powi(2);
        if best.map_or(true, |(_, previous)| distance < previous) {
            best = Some((i, distance));
        }
    }
    best.map(|(i, _)| vec![CellId(i)]).unwrap_or_default()
}

struct SolveResult {
    values: Vec<f64>,
    free: usize,
    iterations: usize,
    residual: f64,
}

fn solve(mesh: &NonconformingRefinement) -> Result<SolveResult, Box<dyn std::error::Error>> {
    let (a, b) = mesh.assemble_poisson(|_| -4.0)?;
    let mut constraints = mesh.hanging_constraints()?;
    for (id, value) in mesh.interpolate_boundary(exact)? {
        constraints.add(id, value, &[])?;
    }
    let closed = constraints.close()?;
    let free = closed.free_dofs().len();
    let (a, b) = closed.condense(&a, &b)?;
    if free == 0 {
        return Ok(SolveResult {
            values: closed.reconstruct(&[])?,
            free,
            iterations: 0,
            residual: 0.0,
        });
    }
    let (reduced, report) = conjugate_gradient(&a, &b, CgOptions::default())?;
    Ok(SolveResult {
        values: closed.reconstruct(&reduced)?,
        free,
        iterations: report.iterations,
        residual: report.residual_norm,
    })
}

/// Bulk-marking fraction: select the largest indicators until their squared
/// sum reaches `ESTIMATOR_THETA` of the total (RFC-0002 §6). Example value,
/// not a convergence guarantee.
const ESTIMATOR_THETA: f64 = 0.3;
/// Absolute global-indicator target (`tol` in RFC-0002 §6) for this example.
/// The loop also stops on an empty mark or after `ESTIMATOR_MAX_CYCLES`.
const ESTIMATOR_TOLERANCE: f64 = 0.05;
/// Hard cycle budget: local nonconforming refinement need not reduce η
/// monotonically, so the demonstration is bounded even when `tol` is not met.
const ESTIMATOR_MAX_CYCLES: usize = 5;

/// Estimator-driven adaptive loop (RFC-0002 §6) on the current mesh.
///
/// Every cycle solves from scratch: assemble, add hanging constraints and
/// Dirichlet values, close, condense, solve the reduced system, then
/// reconstruct. The estimator requires that reconstructed vector — raw CG
/// output has unconstrained slave values and is rejected by the hanging-rule
/// check. The source `f = -4` matches the manufactured solution `u = x² + y²`.
///
/// Only refinement responds to η. Naive η-driven coarsening is unsafe:
/// `NonconformingRefinement::coarsen` restores a parent only when all four
/// children are marked, none has descendants, and no surviving midpoint
/// references a removed midpoint as slave or master
/// (`src/nonconforming_refinement.rs:598-616`). A scalar cell indicator cannot
/// establish any of those structural conditions, and the history offers only
/// level-pop `coarsen_last` (`src/nonconforming_refinement.rs:1170-1175`).
/// The deterministic coarsening demo at the end of `main` stays separate.
///
/// The estimator indexes `mesh.vertices` directly and needs no DofHandler;
/// `DofHandler::from_vertex_count` remains count-only.
fn run_estimator_loop(parent: &Mesh2) -> Result<(), Box<dyn std::error::Error>> {
    let mut transition = NonconformingRefinement::new(parent, &[])?;
    println!(
        "estimator loop: theta={} tolerance={} max_cycles={}",
        ESTIMATOR_THETA, ESTIMATOR_TOLERANCE, ESTIMATOR_MAX_CYCLES
    );
    let mut stop_reason = "cycle budget exhausted";
    for cycle in 0..ESTIMATOR_MAX_CYCLES {
        let solution = solve(&transition)?;
        let estimate =
            poisson_dirichlet_nonconforming(transition.mesh(), &solution.values, |_| -4.0)?;
        let marked = mark_bulk(&estimate.cell_indicators, ESTIMATOR_THETA)?;
        let cells = transition.mesh().cells.len();
        let vertices = transition.mesh().vertices.len();
        println!(
            "cycle {}: eta={:.6} marked={} cells={} vertices={}",
            cycle,
            estimate.global_indicator,
            marked.len(),
            cells,
            vertices
        );
        if marked.is_empty() {
            stop_reason = "empty mark";
            break;
        }
        if estimate.global_indicator <= ESTIMATOR_TOLERANCE {
            stop_reason = "tolerance reached";
            break;
        }
        transition = transition.refine(&marked)?;
        // All IDs from this level are stale after the structural change
        // (ADR-0003 §4). Discard the marks and re-solve from scratch on the new
        // mesh; no solution transfer is defined for this seam.
        drop(marked);
    }
    println!("estimator loop stopped: {}", stop_reason);
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = std::env::args().nth(1).map_or(Ok(4), |s| s.parse())?;
    let parent = Mesh2::unit_square(n)?;

    // Estimator-driven section (RFC-0002 §6), independent of the deterministic
    // two-level solve and the selective coarsening demonstration below.
    run_estimator_loop(&parent)?;

    // Level 1: a deterministic patch in the lower-left quadrant.
    let patch: Vec<CellId> = parent
        .cells()
        .iter()
        .enumerate()
        .filter(|(_, cell)| {
            let (x, y) = centroid(parent.vertices(), cell);
            x < 0.5 && y < 0.5
        })
        .map(|(i, _)| CellId(i))
        .collect();
    let level1 = NonconformingRefinement::new(&parent, &patch)?;

    // Level 2: refine one isolated cell away from the level-1 patch.
    let isolated = isolated_cell(&level1, 0.75, 0.75);
    let target = isolated
        .first()
        .ok_or("no unsplit cell available for a second level")?;
    let level2 = level1.refine(&isolated)?;

    let solution = solve(&level2)?;
    let mut max_hanging_error = 0.0;
    for node in level2.hanging_nodes() {
        let error = (solution.values[node.vertex.0]
            - 0.5 * (solution.values[node.masters[0].0] + solution.values[node.masters[1].0]))
            .abs();
        if error > max_hanging_error {
            max_hanging_error = error;
        }
    }
    println!(
        "level0: cells={} vertices={}",
        parent.cells().len(),
        parent.vertices().len()
    );
    println!(
        "level1: cells={} vertices={} hanging={}",
        level1.cells().len(),
        level1.vertices().len(),
        level1.hanging_nodes().len()
    );
    println!(
        "level2: cells={} vertices={} hanging={} free={} iterations={} residual={:.3e} max_hanging_error={:.3e}",
        level2.cells().len(),
        level2.vertices().len(),
        level2.hanging_nodes().len(),
        solution.free,
        solution.iterations,
        solution.residual,
        max_hanging_error,
    );

    // Selective coarsening: restore the parent of the four level-2 children.
    // All IDs come from `level2`; the returned mesh has its own numbering.
    let children: Vec<CellId> = level2.children(*target)?.to_vec();
    let coarsened = level2.coarsen(&children)?;
    println!(
        "coarsened: cells={} vertices={} hanging={}",
        coarsened.cells().len(),
        coarsened.vertices().len(),
        coarsened.hanging_nodes().len()
    );
    Ok(())
}
