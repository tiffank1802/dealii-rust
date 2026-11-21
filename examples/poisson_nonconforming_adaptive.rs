//! Two-level adaptive nonconforming P1 Poisson on a triangle mesh.
//!
//! The marker is deterministic: a patch of cells is refined from the conforming
//! parent, then one isolated cell is refined a second time. The solve uses the
//! hanging-node constraints composed across both levels, and a final selective
//! coarsening restores that cell from its four children. This is an additive
//! example; the one-level `poisson_hanging` example is untouched.
use dealii_rust::{
    base::Point2,
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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: usize = std::env::args().nth(1).map_or(Ok(4), |s| s.parse())?;
    let parent = Mesh2::unit_square(n)?;

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
