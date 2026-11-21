//! RFC-0001 §13 continuity and vertex-reuse checks for multi-level refine.
use dealii_rust::base::Point2;
use dealii_rust::constraints::{AffineConstraints, ClosedConstraints};
use dealii_rust::grid::{
    BoundaryEdgeId, CellId, Mesh2, NonconformingRefinement, NonconformingRefinementHistory,
    VertexId,
};
use dealii_rust::lac::CooBuilder;
use std::collections::{BTreeMap, BTreeSet};

fn near(a: f64, b: f64) {
    near_tol(a, b, 1e-13);
}

fn near_tol(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{} != {}", a, b);
}

/// Assert `p` lies on the closed segment `a`-`b`. Used for boundary ancestry:
/// a refined exterior edge must stay on the parent exterior edge it came from.
fn point_on_segment(p: Point2, a: Point2, b: Point2) {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let length_squared = dx * dx + dy * dy;
    assert!(length_squared > 0.0, "degenerate parent boundary edge");
    near_tol((p.x - a.x) * dy - (p.y - a.y) * dx, 0.0, 1e-12);
    let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / length_squared;
    assert!(
        (-1e-12..=1.0 + 1e-12).contains(&t),
        "point outside parent boundary segment: t={}",
        t
    );
}

fn key(a: VertexId, b: VertexId) -> (VertexId, VertexId) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

fn total_area(vertices: &[Point2], cells: &[[VertexId; 3]]) -> f64 {
    cells
        .iter()
        .map(|cell| {
            let p = cell.map(|v| vertices[v.0]);
            let cross =
                (p[1].x - p[0].x) * (p[2].y - p[0].y) - (p[1].y - p[0].y) * (p[2].x - p[0].x);
            0.5 * cross.abs()
        })
        .sum()
}

/// Locate an active cell containing both vertices (indices are deterministic
/// for `Mesh2::unit_square`, but locating by incidence keeps the test honest).
fn cell_with(r: &NonconformingRefinement, a: VertexId, b: VertexId) -> CellId {
    r.cells()
        .iter()
        .position(|cell| cell.contains(&a) && cell.contains(&b))
        .map(CellId)
        .unwrap_or_else(|| panic!("no active cell contains {:?} and {:?}", a, b))
}

/// A level-2 midpoint whose master is a level-1 midpoint must be constrained
/// to the mean of its masters, and both rules must close into one chain.
#[test]
fn multi_level_hanging_node_is_average_of_masters() {
    let parent = Mesh2::unit_square(2).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let h1 = first.hanging_nodes()[0];
    let fine_child = cell_with(&first, h1.vertex, h1.masters[0]);

    let second = first.refine(&[fine_child]).unwrap();

    // The level-1 slave is still hanging (its coarse neighbour is unmarked),
    // and a new slave sits on a sub-edge whose master is that level-1 vertex.
    let slaves: BTreeSet<VertexId> = second.hanging_nodes().iter().map(|h| h.vertex).collect();
    assert!(slaves.contains(&h1.vertex));
    let chain = second
        .hanging_nodes()
        .iter()
        .find(|h| slaves.contains(&h.masters[0]) || slaves.contains(&h.masters[1]))
        .expect("no cross-level hanging chain");

    for h in second.hanging_nodes() {
        assert_eq!(
            second.mesh().midpoint_parents[&h.vertex],
            (h.masters[0], h.masters[1])
        );
        let p = second.mesh().vertices[h.vertex.0];
        let a = second.mesh().vertices[h.masters[0].0];
        let b = second.mesh().vertices[h.masters[1].0];
        near(p.x, 0.5 * a.x + 0.5 * b.x);
        near(p.y, 0.5 * a.y + 0.5 * b.y);
    }

    // A prolonged level-0 P1 field satisfies the composed constraints exactly.
    let values: Vec<f64> = (0..second.mesh().original_vertices)
        .map(|i| (i as f64 * 0.37).sin())
        .collect();
    let fine = second.prolongate_p1(&values).unwrap();
    let closed = second.hanging_constraints().unwrap().close().unwrap();
    let reduced: Vec<f64> = closed.free_dofs().iter().map(|&i| fine[i]).collect();
    let reconstructed = closed.reconstruct(&reduced).unwrap();
    for (actual, expected) in reconstructed.iter().zip(&fine) {
        near(*actual, *expected);
    }
    near(
        fine[chain.vertex.0],
        0.5 * (fine[chain.masters[0].0] + fine[chain.masters[1].0]),
    );
}

/// Several intersecting interfaces in two levels must not emit the same
/// hanging slave (or the same master pair) twice.
#[test]
fn multi_level_no_duplicate_hanging_slaves() {
    let parent = Mesh2::unit_square(2).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0), CellId(6)]).unwrap();
    let second = first.refine(&[CellId(1), CellId(4)]).unwrap();

    assert!(!second.hanging_nodes().is_empty());
    let mut slaves = BTreeSet::new();
    let mut master_edges = BTreeSet::new();
    for h in second.hanging_nodes() {
        assert!(
            slaves.insert(h.vertex),
            "duplicate hanging slave {:?}",
            h.vertex
        );
        assert!(
            master_edges.insert(key(h.masters[0], h.masters[1])),
            "duplicate hanging master edge"
        );
        assert_eq!(
            second.mesh().edge_splits[&key(h.masters[0], h.masters[1])],
            h.vertex
        );
    }
    // Duplicate slaves would make `AffineConstraints::add` fail.
    assert!(second.hanging_constraints().unwrap().close().is_ok());
}

/// Reusing an already split edge keeps the level-1 midpoint id, adds no
/// duplicate vertex, and drops the interface once both neighbours are fine.
#[test]
fn already_split_edge_reuses_vertex() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    assert!(first.refine(&[]).is_err());
    assert!(first.refine(&[CellId(0), CellId(0)]).is_err());
    assert!(first.refine(&[CellId(first.cells().len())]).is_err());

    let h1 = first.hanging_nodes()[0];
    let edge = key(h1.masters[0], h1.masters[1]);
    let midpoint = first.mesh().edge_splits[&edge];
    let old_point: Point2 = first.mesh().vertices[midpoint.0];

    // Cell 4 is the coarse neighbour carried unchanged by the first level.
    let second = first.refine(&[CellId(4)]).unwrap();

    assert_eq!(second.mesh().edge_splits[&edge], midpoint);
    assert_eq!(second.mesh().vertices[midpoint.0], old_point);
    // Only the two previously unsplit edges of cell 1 get fresh midpoints.
    assert_eq!(
        second.mesh().vertices.len(),
        first.mesh().vertices.len() + 2
    );
    assert_eq!(second.mesh().cells.len(), 8);
    // Both sides of the diagonal are now refined: no T-junction remains.
    assert!(second.hanging_nodes().is_empty());
}

/// RFC-0001 §8: marking every child of every parent restores the original
/// parent mesh exactly — vertices, connectivity, area and no hanging nodes.
#[test]
fn coarsen_full_sibling_group_restores_parent() {
    let parent = Mesh2::unit_square(1).unwrap();
    let fine = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();
    let all: Vec<CellId> = (0..fine.cells().len()).map(CellId).collect();

    let coarse = fine.coarsen(&all).unwrap();

    assert_eq!(coarse.vertices(), parent.vertices());
    assert_eq!(coarse.cells(), parent.cells());
    assert!(coarse.hanging_nodes().is_empty());
    assert!(coarse.mesh().edge_splits.is_empty());
    assert!(coarse.boundary_edges() == parent.boundary_edges());
    near(
        total_area(coarse.vertices(), coarse.cells()),
        total_area(parent.vertices(), parent.cells()),
    );
    // The input transition is untouched, and the returned one is a usable
    // base for another level.
    assert_eq!(fine.cells().len(), 8);
    let again = coarse.refine(&[CellId(0)]).unwrap();
    assert_eq!(again.mesh().original_vertices, parent.vertices().len());
    assert!(again.hanging_constraints().unwrap().close().is_ok());
}

/// RFC-0001 §8 condition 1 and §11: a partial sibling group and an empty
/// selection are rejected.
#[test]
fn coarsen_rejects_partial_sibling_group() {
    let parent = Mesh2::unit_square(1).unwrap();
    let fine = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();

    assert!(fine.coarsen(&[]).is_err());
    assert!(fine.coarsen(&[CellId(0), CellId(1), CellId(2)]).is_err());
    assert!(fine
        .coarsen(&[CellId(0), CellId(0), CellId(1), CellId(2), CellId(3)])
        .is_err());
    assert!(fine.coarsen(&[CellId(8)]).is_err());
}

/// RFC-0001 §8 condition 2: a child that still has descendants cannot be
/// coarsened. Only leaves are active; restoring a parent whose child was
/// refined again would need two levels in one call, so mixing the level-2
/// children with the surviving level-1 children is rejected.
#[test]
fn coarsen_rejects_cell_with_descendants() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();
    let second = first.refine(&[CellId(0)]).unwrap();

    let mut marked: Vec<CellId> = second.children(CellId(0)).unwrap().to_vec();
    assert_eq!(marked.len(), 4);
    marked.extend_from_slice(&[CellId(4), CellId(5), CellId(6)]);

    assert!(second.coarsen(&marked).is_err());
}

/// RFC-0001 §8 conditions 3-4: a neighbour still refined across a shared
/// edge depends on its midpoint, so sibling groups cannot be restored
/// independently of each other.
#[test]
fn coarsen_rejects_neighbor_dependency() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();
    let second = first
        .refine(&[CellId(0), CellId(1), CellId(2), CellId(3)])
        .unwrap();

    let mut marked = Vec::new();
    marked.extend_from_slice(second.children(CellId(0)).unwrap());
    marked.extend_from_slice(second.children(CellId(1)).unwrap());
    assert!(second.coarsen(&marked).is_err());
}

/// RFC-0001 §8: the interface midpoint reused while refining the coarse
/// neighbour was introduced by the earlier transition, so it is not removed.
/// Coarsening the neighbour restores the level-1 mesh exactly, hanging node
/// included, instead of being rejected as a surviving vertex reference.
#[test]
fn coarsen_reused_interface_midpoint_restores_parent() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let shared = first.hanging_nodes()[0].vertex;

    let second = first.refine(&[CellId(4)]).unwrap();
    assert!(second.hanging_nodes().is_empty());
    assert!(second.cells().iter().any(|cell| cell.contains(&shared)));

    let kids: Vec<CellId> = second.children(CellId(4)).unwrap().to_vec();
    let restored = second.coarsen(&kids).unwrap();

    assert_eq!(restored.vertices(), first.vertices());
    assert_eq!(restored.cells(), first.cells());
    assert_eq!(restored.mesh().edge_splits, first.mesh().edge_splits);
    assert_eq!(restored.hanging_nodes(), first.hanging_nodes());
    assert!(restored.hanging_constraints().unwrap().close().is_ok());
}

/// Fix 3: after coarsening, `children` and `parent_cells` share one index
/// space, so a cell can be found through `children(parent_cell(cell))` even
/// when the parent-level group index differs from the new cell index.
#[test]
fn coarsen_children_and_parent_cells_share_index_space() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let second = first.refine(&[CellId(0), CellId(4)]).unwrap();
    let kids: Vec<CellId> = second.children(CellId(4)).unwrap().to_vec();

    let coarse = second.coarsen(&kids).unwrap();
    assert!(coarse.cells().len() > kids.len());
    for i in 0..coarse.cells().len() {
        let cell = CellId(i);
        let parent = coarse.parent_cell(cell).unwrap();
        assert!(
            coarse.children(parent).unwrap().contains(&cell),
            "cell {:?} missing from children of parent {:?}",
            cell,
            parent
        );
    }
}

/// RFC-0001 §8: one level per call must also work for independent refined
/// regions. Coarsening one of two disconnected sibling groups must leave a
/// transition the second group can be coarsened from in a separate call.
#[test]
fn coarsen_twice_on_disconnected_refined_regions() {
    let vertices = vec![
        Point2::new(0.0, 0.0).unwrap(),
        Point2::new(1.0, 0.0).unwrap(),
        Point2::new(0.0, 1.0).unwrap(),
        Point2::new(2.0, 0.0).unwrap(),
        Point2::new(3.0, 0.0).unwrap(),
        Point2::new(2.0, 1.0).unwrap(),
    ];
    let cells = vec![[0, 1, 2].map(VertexId), [3, 4, 5].map(VertexId)];
    let parent = Mesh2::new(vertices, cells).unwrap();
    let fine = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();

    let once = fine.coarsen(fine.children(CellId(1)).unwrap()).unwrap();
    let twice = once.coarsen(once.children(CellId(0)).unwrap()).unwrap();

    assert_eq!(twice.cells(), parent.cells());
    for i in 0..twice.cells().len() {
        let cell = CellId(i);
        let parent_id = twice.parent_cell(cell).unwrap();
        assert!(
            twice.children(parent_id).unwrap().contains(&cell),
            "cell {:?} missing from children of parent {:?}",
            cell,
            parent_id
        );
    }
}

/// RFC-0001 §1 validation seam: malformed caller-supplied meshes are rejected
/// at the boundary by `validate` and `try_new` instead of panicking later.
#[test]
fn mesh_validate_rejects_out_of_range_ids_and_nonfinite_coordinates() {
    let parent = Mesh2::unit_square(1).unwrap();
    let mut mesh = NonconformingRefinement::new(&parent, &[])
        .unwrap()
        .mesh()
        .clone();
    assert!(mesh.validate().is_ok());
    assert!(NonconformingRefinementHistory::try_new(mesh.clone()).is_ok());

    mesh.cells[0][0] = VertexId(usize::MAX);
    assert!(mesh.validate().is_err());
    assert!(NonconformingRefinementHistory::try_new(mesh.clone()).is_err());

    let mut mesh = NonconformingRefinement::new(&parent, &[])
        .unwrap()
        .mesh()
        .clone();
    mesh.vertices[0].y = f64::NAN;
    assert!(mesh.validate().is_err());
    assert!(NonconformingRefinementHistory::try_new(mesh).is_err());
}

/// ADR-0003 §4: coarsening invalidates every typed ID. Only IDs from the
/// returned mesh may be used afterwards; old midpoint/boundary-edge indices
/// are out of bounds in the compacted mesh.
#[test]
fn coarsen_invalidates_ids_documented() {
    let parent = Mesh2::unit_square(1).unwrap();
    let fine = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();
    let old_midpoint = *fine.mesh().edge_splits.values().next().unwrap();
    let old_boundary_edges = fine.boundary_edges().len();
    let all: Vec<CellId> = (0..fine.cells().len()).map(CellId).collect();

    let coarse = fine.coarsen(&all).unwrap();

    assert!(old_midpoint.0 >= coarse.vertices().len());
    assert!(BoundaryEdgeId(old_boundary_edges - 1).0 >= coarse.boundary_edges().len());
    assert!(coarse
        .parent_boundary_edge(BoundaryEdgeId(old_boundary_edges - 1))
        .is_err());
    for cell in coarse.cells() {
        for &v in cell {
            assert!(v.0 < coarse.vertices().len());
        }
    }
    for i in 0..coarse.cells().len() {
        assert!(coarse.parent_cell(CellId(i)).is_ok());
    }
    assert!(coarse.parent_boundary_edge(BoundaryEdgeId(0)).is_ok());
}

/// RFC-0001 §9: a history tracks the level across several refine/pop cycles
/// and returns to the exact retained mesh after each pop.
#[test]
fn history_multi_cycle_refine_coarsen() {
    let parent = Mesh2::unit_square(1).unwrap();
    let root = NonconformingRefinement::new(&parent, &[]).unwrap();
    let mut history = NonconformingRefinementHistory::new(root.mesh().clone());

    assert_eq!(history.level(), 0);

    // Cycle 1: refine cell 0, then its coarse neighbour (cell 4).
    history.refine(&[CellId(0)]).unwrap();
    assert_eq!(history.level(), 1);
    let level1 = history.mesh().clone();
    // Refining from the bare root reproduces the direct one-level result.
    let direct = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    assert_eq!(level1.cells, direct.mesh().cells);
    assert_eq!(level1.edge_splits, direct.mesh().edge_splits);
    assert_eq!(level1.hanging, direct.mesh().hanging);

    history.refine(&[CellId(4)]).unwrap();
    assert_eq!(history.level(), 2);
    assert_eq!(history.mesh().cells.len(), 8);
    assert!(history.mesh().hanging.is_empty());

    history.coarsen_last().unwrap();
    assert_eq!(history.level(), 1);
    assert_eq!(history.mesh().cells, level1.cells);
    assert_eq!(history.mesh().edge_splits, level1.edge_splits);
    assert_eq!(history.mesh().hanging, level1.hanging);

    // Cycle 2: a fresh refinement from the restored level, then pop twice.
    history.refine(&[CellId(4)]).unwrap();
    assert_eq!(history.level(), 2);
    history.coarsen_last().unwrap();
    assert_eq!(history.level(), 1);
    assert_eq!(history.mesh().cells, level1.cells);
    history.coarsen_last().unwrap();
    assert_eq!(history.level(), 0);
    assert_eq!(history.mesh().cells, root.mesh().cells);
    assert_eq!(history.mesh().edge_splits, root.mesh().edge_splits);
}

/// RFC-0001 §9: `mesh_at` exposes the retained view at each level and rejects
/// levels that were never reached or have been popped.
#[test]
fn history_mesh_at_level() {
    let parent = Mesh2::unit_square(2).unwrap();
    let root = NonconformingRefinement::new(&parent, &[]).unwrap();
    let mut history = NonconformingRefinementHistory::new(root.mesh().clone());

    let level0 = history.mesh().clone();
    history.refine(&[CellId(0)]).unwrap();
    let level1 = history.mesh().clone();
    history.refine(&[CellId(1), CellId(4)]).unwrap();
    let level2 = history.mesh().clone();

    assert_eq!(history.level(), 2);
    assert_eq!(history.mesh_at(0).unwrap().cells, level0.cells);
    assert_eq!(history.mesh_at(1).unwrap().cells, level1.cells);
    assert_eq!(history.mesh_at(2).unwrap().cells, level2.cells);
    assert_eq!(history.mesh_at(2).unwrap().vertices, level2.vertices);
    assert_eq!(history.mesh_at(2).unwrap().edge_splits, level2.edge_splits);
    assert_eq!(history.mesh_at(2).unwrap().hanging, level2.hanging);
    assert!(history.mesh_at(3).is_err());
    assert!(history.mesh_at(usize::MAX).is_err());
    // The earlier level is retained, not overwritten by the current view.
    assert_ne!(history.mesh_at(1).unwrap().cells.len(), level2.cells.len());

    history.coarsen_last().unwrap();
    assert!(history.mesh_at(2).is_err());
    assert_eq!(history.mesh_at(1).unwrap().cells, level1.cells);
}

/// RFC-0001 §6: a count-only handler takes its count from the nonconforming
/// mesh's vertex list; connectivity accessors are unavailable through it.
#[test]
fn count_only_dof_handler_uses_nonconforming_mesh_vertex_count() {
    use dealii_rust::dofs::DofHandler;

    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let second = first.refine(&[CellId(4)]).unwrap();
    let mesh = second.mesh();

    let dofs = DofHandler::from_vertex_count(mesh.vertices.len());
    assert_eq!(dofs.count(), mesh.vertices.len());
    assert!(dofs.cell_dofs(CellId(0)).is_err());
    assert!(dofs.try_boundary_dofs().is_err());
}

/// RFC-0001 §9/§11: rejected operations never append or remove a level.
#[test]
fn history_errors_are_non_mutating() {
    let parent = Mesh2::unit_square(1).unwrap();
    let root = NonconformingRefinement::new(&parent, &[]).unwrap();
    let mut history = NonconformingRefinementHistory::new(root.mesh().clone());
    let base = history.mesh().clone();

    assert!(history.refine(&[]).is_err());
    assert!(history.refine(&[CellId(usize::MAX)]).is_err());
    assert!(history.refine(&[CellId(0), CellId(0)]).is_err());
    assert_eq!(history.level(), 0);
    assert_eq!(history.mesh().cells, base.cells);
    assert_eq!(history.mesh().edge_splits, base.edge_splits);

    history.refine(&[CellId(0)]).unwrap();
    let level1 = history.mesh().clone();
    assert!(history.refine(&[CellId(level1.cells.len())]).is_err());
    assert_eq!(history.level(), 1);
    assert_eq!(history.mesh().cells, level1.cells);
    assert_eq!(history.mesh().edge_splits, level1.edge_splits);

    history.coarsen_last().unwrap();
    assert_eq!(history.level(), 0);
    assert!(history.coarsen_last().is_err());
    assert_eq!(history.level(), 0);
    assert_eq!(history.mesh().cells, base.cells);
    assert_eq!(history.mesh().edge_splits, base.edge_splits);
    assert!(history.mesh_at(1).is_err());
}

/// RFC-0001 A5 dropped the mesh-level `boundary_parents` field. The
/// representable invariant lives on the transition: every refined exterior
/// edge maps through `parent_boundary_edge` to an exterior edge of the
/// immediately preceding level, the children cover every parent exterior
/// edge, each child stays geometrically on its parent, and no interface edge
/// is exposed as a boundary edge.
#[test]
fn multi_level_boundary_edges_are_parent_boundary() {
    let parent = Mesh2::unit_square(1).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();

    // Level 1 -> root `Mesh2`: the origin of the retained boundary ids.
    let mut root_coverage = BTreeSet::new();
    for i in 0..first.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        let pid = first.parent_boundary_edge(id).unwrap();
        assert!(pid.0 < parent.boundary_edges().len());
        root_coverage.insert(pid);
        let ancestor = parent.boundary_edge_geometry(pid).unwrap();
        let child = first.boundary_edge_geometry(id).unwrap();
        point_on_segment(child.points[0], ancestor.points[0], ancestor.points[1]);
        point_on_segment(child.points[1], ancestor.points[0], ancestor.points[1]);
    }
    assert_eq!(root_coverage.len(), parent.boundary_edges().len());

    let second = first.refine(&[CellId(4)]).unwrap();
    assert_boundary_ancestry(&second, &first);

    // Independent path on a larger root, including a split interface.
    let base = Mesh2::unit_square(2).unwrap();
    let a = NonconformingRefinement::new(&base, &[CellId(0)]).unwrap();
    let b = a.refine(&[CellId(1), CellId(4)]).unwrap();
    assert_boundary_ancestry(&b, &a);
}

fn assert_boundary_ancestry(child: &NonconformingRefinement, parent: &NonconformingRefinement) {
    let mut covered = BTreeSet::new();
    let mut length_by_parent = BTreeMap::new();
    for i in 0..child.boundary_edges().len() {
        let id = BoundaryEdgeId(i);
        let pid = child.parent_boundary_edge(id).unwrap();
        assert!(pid.0 < parent.boundary_edges().len());
        covered.insert(pid);
        let ancestor = parent.boundary_edge_geometry(pid).unwrap();
        let edge = child.boundary_edge_geometry(id).unwrap();
        point_on_segment(edge.points[0], ancestor.points[0], ancestor.points[1]);
        point_on_segment(edge.points[1], ancestor.points[0], ancestor.points[1]);
        *length_by_parent.entry(pid.0).or_insert(0.0) += edge.length;
        let owners = child
            .cells()
            .iter()
            .filter(|cell| cell.contains(&edge.vertices[0]) && cell.contains(&edge.vertices[1]))
            .count();
        assert_eq!(owners, 1, "interface edge exposed as boundary");
    }
    assert_eq!(covered.len(), parent.boundary_edges().len());
    // Child exterior edges partition their parent edge: lengths must add up.
    for i in 0..parent.boundary_edges().len() {
        let ancestor = parent.boundary_edge_geometry(BoundaryEdgeId(i)).unwrap();
        let covered_length = length_by_parent.get(&i).copied().unwrap_or(0.0);
        near_tol(covered_length, ancestor.length, 1e-12);
    }
}

/// A two-level mesh with one hanging master that is itself hanging, used by
/// the composition tests below.
fn cross_level_chain() -> NonconformingRefinement {
    let parent = Mesh2::unit_square(2).unwrap();
    let first = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let h = first.hanging_nodes()[0];
    let fine_child = cell_with(&first, h.vertex, h.masters[0]);
    first.refine(&[fine_child]).unwrap()
}

/// RFC-0001 §7: `close` removes every constrained slave from the reduced set
/// and reports the remaining free original DoFs in increasing order.
#[test]
fn constraints_close_resolves_to_free_dofs() {
    let level = cross_level_chain();
    let closed = level.hanging_constraints().unwrap().close().unwrap();
    let constrained: BTreeSet<usize> = level.hanging_nodes().iter().map(|h| h.vertex.0).collect();
    let expected: Vec<usize> = (0..level.vertices().len())
        .filter(|i| !constrained.contains(i))
        .collect();
    assert_eq!(closed.size(), level.vertices().len());
    assert_eq!(closed.free_dofs(), expected.as_slice());
    assert!(closed.free_dofs().len() < closed.size());
    assert!(level.hanging_nodes().iter().any(|h| {
        constrained.contains(&h.masters[0].0) || constrained.contains(&h.masters[1].0)
    }));
}

/// RFC-0001 §7: reconstruction satisfies every original midpoint rule, both
/// without offsets and with nonzero offsets propagating through a chain.
#[test]
fn constraints_reconstruct_satisfies_original_rules() {
    let level = cross_level_chain();
    let plain = level.hanging_constraints().unwrap().close().unwrap();
    let reduced: Vec<f64> = (0..plain.free_dofs().len())
        .map(|j| ((j + 1) as f64 * 0.31).sin() + 0.05 * j as f64)
        .collect();
    let values = plain.reconstruct(&reduced).unwrap();
    for h in level.hanging_nodes() {
        near_tol(
            values[h.vertex.0],
            0.5 * (values[h.masters[0].0] + values[h.masters[1].0]),
            1e-12,
        );
    }

    let offsets: Vec<f64> = (0..level.hanging_nodes().len())
        .map(|i| 0.01 * (i + 1) as f64)
        .collect();
    let mut builder = AffineConstraints::new(level.vertices().len());
    for (&offset, h) in offsets.iter().zip(level.hanging_nodes()) {
        builder
            .add(
                h.vertex.0,
                offset,
                &[(h.masters[0].0, 0.5), (h.masters[1].0, 0.5)],
            )
            .unwrap();
    }
    let closed = builder.close().unwrap();
    let shifted = closed.reconstruct(&reduced).unwrap();
    for (&offset, h) in offsets.iter().zip(level.hanging_nodes()) {
        near_tol(
            shifted[h.vertex.0],
            0.5 * (shifted[h.masters[0].0] + shifted[h.masters[1].0]) + offset,
            1e-11,
        );
    }
}

/// RFC-0001 §7/§12: the same iterative `close` that resolves mesh chains
/// terminates on a 4000-link chain without recursion.
#[test]
fn constraints_deep_chain_still_terminates() {
    let n = 4000;
    let mut builder = AffineConstraints::new(n + 1);
    builder.add(n, 5.0, &[]).unwrap();
    for i in 1..n {
        builder.add(i, 0.0, &[(i - 1, 0.5), (n, 0.5)]).unwrap();
    }
    let closed = builder.close().unwrap();
    assert_eq!(closed.free_dofs(), &[0]);
    let values = closed.reconstruct(&[2.0]).unwrap();
    near_tol(values[1], 3.5, 1e-12);
    near_tol(values[n - 1], 5.0, 1e-12);
    near_tol(values[n], 5.0, 1e-12);
}

/// RFC-0001 §7: a fully constrained submesh has a 0x0 reduced system and
/// reconstructs from an empty reduced vector.
#[test]
fn fully_constrained_submesh_reconstructs_empty() {
    let parent = Mesh2::unit_square(1).unwrap();
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let n = level.vertices().len();
    let slaves: BTreeSet<usize> = level.hanging_nodes().iter().map(|h| h.vertex.0).collect();
    assert!(!slaves.is_empty());

    let mut builder = level.hanging_constraints().unwrap();
    for v in 0..n {
        if !slaves.contains(&v) {
            builder.add(v, 0.25 * v as f64, &[]).unwrap();
        }
    }
    let closed = builder.close().unwrap();
    assert!(closed.free_dofs().is_empty());
    let values = closed.reconstruct(&[]).unwrap();
    assert_eq!(values.len(), n);
    for h in level.hanging_nodes() {
        near_tol(
            values[h.vertex.0],
            0.5 * (values[h.masters[0].0] + values[h.masters[1].0]),
            1e-12,
        );
    }

    let mut identity = CooBuilder::new(n, n);
    for i in 0..n {
        identity.add(i, i, 1.0).unwrap();
    }
    let matrix = identity.into_csr().unwrap();
    let (reduced, rhs) = closed.condense(&matrix, &values).unwrap();
    assert_eq!((reduced.rows(), reduced.cols()), (0, 0));
    assert!(rhs.is_empty());
    assert!(closed.reconstruct(&[0.0]).is_err());
}

/// RFC-0001 §8: there is no explicit `restrict` operator. For nested P1
/// spaces the meaningful round trip is to prolongate from the retained parent
/// vertices, then restrict by nodal injection (fine values at the retained
/// parent indices) and inverse midpoint interpolation. The structural coarsen
/// step is checked separately by restoring the exact parent mesh.
#[test]
fn prolongate_then_restrict_recovers_parent_field() {
    let parent = Mesh2::unit_square(1).unwrap();
    let n = parent.vertices().len();
    let field: Vec<f64> = (0..n)
        .map(|i| (i as f64 * 0.7).sin() + 0.25 * i as f64)
        .collect();

    let fine_level = NonconformingRefinement::new(&parent, &[CellId(0), CellId(1)]).unwrap();
    let fine = fine_level.prolongate_p1(&field).unwrap();
    for i in 0..n {
        near_tol(fine[i], field[i], 1e-12);
    }
    for (&(a, b), &mid) in &fine_level.mesh().edge_splits {
        near_tol(2.0 * fine[mid.0] - fine[b.0], fine[a.0], 1e-12);
        near_tol(2.0 * fine[mid.0] - fine[a.0], fine[b.0], 1e-12);
    }

    let deeper = fine_level.refine(&[CellId(0)]).unwrap();
    let finer = deeper.prolongate_p1(&field).unwrap();
    for i in 0..n {
        near_tol(finer[i], field[i], 1e-12);
    }
    for h in deeper.hanging_nodes() {
        near_tol(
            finer[h.vertex.0],
            0.5 * (finer[h.masters[0].0] + finer[h.masters[1].0]),
            1e-12,
        );
    }

    let all: Vec<CellId> = (0..fine_level.cells().len()).map(CellId).collect();
    let restored = fine_level.coarsen(&all).unwrap();
    assert_eq!(restored.vertices(), parent.vertices());
    assert_eq!(restored.cells(), parent.cells());
}

/// RFC-0001 §11: `new` accepts an empty marking, but `refine`/`coarsen` do not,
/// and a rejected operation leaves the transition unchanged.
#[test]
fn empty_markings_are_rejected_by_refine_and_coarsen() {
    let parent = Mesh2::unit_square(1).unwrap();
    let root = NonconformingRefinement::new(&parent, &[]).unwrap();
    assert_eq!(root.cells().len(), parent.cells().len());
    assert!(root.refine(&[]).is_err());
    assert!(root.coarsen(&[]).is_err());

    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    let before_cells = level.cells().to_vec();
    let before_splits = level.mesh().edge_splits.clone();
    assert!(level.refine(&[]).is_err());
    assert!(level.coarsen(&[]).is_err());
    assert_eq!(level.cells(), before_cells.as_slice());
    assert_eq!(level.mesh().edge_splits, before_splits);
}

/// RFC-0001 §12: non-finite input cannot reach a midpoint. Geometry rejects
/// non-finite points at the mesh boundary, and interpolation rejects
/// non-finite fields and constraint coefficients.
#[test]
fn non_finite_midpoint_rejection() {
    assert!(Point2::new(f64::NAN, 0.0).is_err());
    assert!(Point2::new(0.0, f64::INFINITY).is_err());
    let bad = vec![
        Point2::new(0.0, 0.0).unwrap(),
        Point2::new(1.0, 0.0).unwrap(),
        Point2 {
            x: 0.0,
            y: f64::NAN,
        },
    ];
    assert!(Mesh2::new(bad, vec![[VertexId(0), VertexId(1), VertexId(2)]]).is_err());

    let level =
        NonconformingRefinement::new(&Mesh2::unit_square(1).unwrap(), &[CellId(0)]).unwrap();
    assert!(level.prolongate_p1(&[1.0, f64::NAN, 3.0, 4.0]).is_err());
    assert!(level
        .prolongate_p1(&[f64::INFINITY, 1.0, 2.0, 3.0])
        .is_err());
    let mut builder = AffineConstraints::new(3);
    assert!(builder.add(2, 0.0, &[(0, 0.5), (1, f64::NAN)]).is_err());
}

/// RFC-0001 §7/§12: malformed affine rules are rejected at `add` or `close`,
/// and a failed `add` leaves the builder unchanged.
#[test]
fn malformed_constraint_rules_are_rejected() {
    let mut builder = AffineConstraints::new(3);
    assert!(builder.add(3, 0.0, &[]).is_err());
    assert!(builder.add(0, f64::NAN, &[]).is_err());
    assert!(builder.add(0, 0.0, &[(3, 1.0)]).is_err());
    assert!(builder.add(0, 0.0, &[(1, f64::INFINITY)]).is_err());
    assert!(builder.add(0, 0.0, &[(1, 0.5), (1, 0.5)]).is_err());
    builder.add(0, 0.0, &[(1, 0.5), (2, 0.5)]).unwrap();
    assert!(builder.add(0, 0.0, &[(1, 1.0)]).is_err());
    assert_eq!(builder.close().unwrap().free_dofs(), &[1, 2]);

    let mut cycle = AffineConstraints::new(3);
    cycle.add(0, 0.0, &[(1, 1.0)]).unwrap();
    cycle.add(1, 0.0, &[(2, 1.0)]).unwrap();
    cycle.add(2, 0.0, &[(0, 1.0)]).unwrap();
    assert!(cycle.close().is_err());
    let mut self_dependency = AffineConstraints::new(1);
    self_dependency.add(0, 0.0, &[(0, 1.0)]).unwrap();
    assert!(self_dependency.close().is_err());

    let level =
        NonconformingRefinement::new(&Mesh2::unit_square(1).unwrap(), &[CellId(0)]).unwrap();
    let mut from_mesh = level.hanging_constraints().unwrap();
    for h in level.hanging_nodes() {
        assert!(from_mesh
            .add(
                h.vertex.0,
                0.0,
                &[(h.masters[0].0, 0.5), (h.masters[1].0, 0.5)]
            )
            .is_err());
    }
}

/// RFC-0001 §12: the `usize` accumulation guards (`checked_mul` in `new`/
/// `refine`/`unit_square`, `checked_add`/`checked_sub` in `coarsen`) cannot be
/// driven to overflow with a real allocation, so their input boundary is
/// tested with `usize::MAX` IDs. The floating-point guards in `close` are
/// exercised with offsets/coefficients that would overflow when expanded.
#[test]
fn cumulative_overflow_guards_reject_invalid_inputs() {
    let parent = Mesh2::unit_square(1).unwrap();
    assert!(Mesh2::unit_square(usize::MAX).is_err());
    assert!(NonconformingRefinement::new(&parent, &[CellId(usize::MAX)]).is_err());
    let level = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    assert!(level.refine(&[CellId(usize::MAX)]).is_err());
    assert!(level.coarsen(&[CellId(usize::MAX)]).is_err());
    assert!(level.refine(&[CellId(0), CellId(usize::MAX)]).is_err());
    assert!(level.coarsen(&[CellId(0), CellId(usize::MAX)]).is_err());

    let mut offsets = AffineConstraints::new(3);
    offsets.add(1, f64::MAX, &[(0, 2.0)]).unwrap();
    offsets.add(2, 0.0, &[(1, 2.0)]).unwrap();
    assert!(offsets.close().is_err());
    let mut weights = AffineConstraints::new(3);
    weights.add(1, 0.0, &[(0, f64::MAX)]).unwrap();
    weights.add(2, 0.0, &[(1, 2.0)]).unwrap();
    assert!(weights.close().is_err());
}

/// Keep refining the first fine cell: its first vertex is an original vertex,
/// so every level adds a midpoint on a sub-edge of the same interface and the
/// hanging relation remains a linear chain (no 2^level fan-out).
fn deep_single_cell_chain(levels: usize) -> NonconformingRefinement {
    let parent = Mesh2::unit_square(1).unwrap();
    let mut current = NonconformingRefinement::new(&parent, &[CellId(0)]).unwrap();
    for _ in 1..levels {
        current = current.refine(&[CellId(0)]).unwrap();
    }
    current
}

/// Column support of the closed constraint matrix, measured by reconstructing
/// each reduced basis vector: `support[s]` is the number of free DoFs that
/// original DoF `s` depends on, i.e. the row width of its expanded rule.
fn constraint_row_support(closed: &ClosedConstraints) -> Vec<usize> {
    let m = closed.free_dofs().len();
    let mut support = vec![0usize; closed.size()];
    for j in 0..m {
        let mut reduced = vec![0.0; m];
        reduced[j] = 1.0;
        let values = closed.reconstruct(&reduced).unwrap();
        for (slot, &value) in support.iter_mut().zip(&values) {
            if value.abs() > 1e-12 {
                *slot += 1;
            }
        }
    }
    support
}

/// RFC-0001 §13: worst-case constraint row width on a deep one-sided chain
/// stays bounded instead of doubling with each level.
#[test]
fn deep_refinement_row_width_bounded() {
    let shallow = deep_single_cell_chain(3);
    let deep = deep_single_cell_chain(6);
    assert!(deep.vertices().len() > shallow.vertices().len());
    assert!(deep.hanging_nodes().len() >= shallow.hanging_nodes().len());

    let shallow_support =
        constraint_row_support(&shallow.hanging_constraints().unwrap().close().unwrap());
    let deep_support =
        constraint_row_support(&deep.hanging_constraints().unwrap().close().unwrap());
    let shallow_max = shallow_support.iter().copied().max().unwrap_or(0);
    let deep_max = deep_support.iter().copied().max().unwrap_or(0);
    assert!(shallow_max >= 2 && deep_max >= 2);
    assert!(
        deep_max <= 8,
        "unexpected row width {} at depth 6",
        deep_max
    );
    assert!(
        deep_max <= shallow_max + 2,
        "row width grew from {} to {}",
        shallow_max,
        deep_max
    );
    let slaves: BTreeSet<usize> = deep.hanging_nodes().iter().map(|h| h.vertex.0).collect();
    assert!(slaves.iter().all(|&s| deep_support[s] >= 1));
}

/// RFC-0001 §13/§12: six refinement levels plus constraint close and
/// reconstruction must finish in a bounded time (no recursion or stack blowup).
#[test]
fn deep_refinement_completes_in_budget() {
    use std::time::{Duration, Instant};

    let start = Instant::now();
    let level = deep_single_cell_chain(6);
    let closed = level.hanging_constraints().unwrap().close().unwrap();
    let reduced: Vec<f64> = (0..closed.free_dofs().len())
        .map(|j| ((j + 1) as f64 * 0.19).sin())
        .collect();
    let values = closed.reconstruct(&reduced).unwrap();
    for h in level.hanging_nodes() {
        near_tol(
            values[h.vertex.0],
            0.5 * (values[h.masters[0].0] + values[h.masters[1].0]),
            1e-12,
        );
    }
    let elapsed = start.elapsed();
    assert!(
        elapsed < Duration::from_secs(5),
        "deep refinement took {:?}",
        elapsed
    );
}
