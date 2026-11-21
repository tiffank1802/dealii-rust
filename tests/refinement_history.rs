use dealii_rust::{
    base::Point2,
    grid::{CellId, Mesh2, RefinementHistory, VertexId},
};
fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y).unwrap()
}
#[test]
fn retains_levels_and_rolls_back_whole_levels_exactly() {
    let base = Mesh2::unit_square(3).unwrap();
    let mut h = RefinementHistory::new(base.clone());
    assert_eq!(h.level(), 0);
    assert_eq!(h.mesh_at(0).unwrap().cells(), base.cells());
    h.refine(&[CellId(0), CellId(5)]).unwrap();
    let level1 = h.mesh().clone();
    h.refine(&[CellId(0)]).unwrap();
    let level2 = h.mesh().clone();
    assert_eq!(h.level(), 2);
    assert_eq!(h.mesh_at(0).unwrap().vertices(), base.vertices());
    assert_eq!(h.mesh_at(1).unwrap().cells(), level1.cells());
    assert_eq!(h.mesh_at(2).unwrap().cells(), level2.cells());
    for i in 0..level2.cells().len() {
        let parent = h.transition(1).unwrap().parent_cell(CellId(i)).unwrap();
        assert!(parent.0 < level1.cells().len());
        assert!(h
            .transition(1)
            .unwrap()
            .children(parent)
            .unwrap()
            .contains(&CellId(i)));
    }
    h.coarsen_last().unwrap();
    assert_eq!(h.level(), 1);
    assert_eq!(h.mesh().vertices(), level1.vertices());
    assert_eq!(h.mesh().cells(), level1.cells());
    assert!(h.mesh_at(2).is_err());
    assert!(h.transition(1).is_err());
    h.coarsen_last().unwrap();
    assert_eq!(h.level(), 0);
    assert_eq!(h.mesh().cells(), base.cells());
    assert!(h.coarsen_last().is_err());
}
#[test]
fn transfer_then_injection_recovers_parent_at_multiple_levels() {
    let base = Mesh2::unit_square(2).unwrap();
    let values: Vec<_> = (0..base.vertices().len())
        .map(|i| ((i * i + 7) % 13) as f64)
        .collect();
    let mut h = RefinementHistory::new(base);
    h.refine(&[CellId(0)]).unwrap();
    let mid = h.transition(0).unwrap().prolongate_p1(&values).unwrap();
    h.refine(&[CellId(3)]).unwrap();
    let fine = h.transition(1).unwrap().prolongate_p1(&mid).unwrap();
    assert_eq!(h.restrict_injection(1, &fine).unwrap(), mid);
    assert_eq!(h.restrict_injection(0, &mid).unwrap(), values);
    h.coarsen_last().unwrap();
    assert_eq!(h.mesh().vertices().len(), mid.len());
    h.coarsen_last().unwrap();
    assert_eq!(h.mesh().vertices().len(), values.len());
}
#[test]
fn injection_discards_child_detail_and_is_not_projection() {
    let base = Mesh2::unit_square(1).unwrap();
    let mut h = RefinementHistory::new(base);
    h.refine(&[CellId(0)]).unwrap();
    let coarse = vec![0.; h.mesh_at(0).unwrap().vertices().len()];
    let mut child = h.transition(0).unwrap().prolongate_p1(&coarse).unwrap();
    let last = child.len() - 1;
    child[last] = 7.;
    let restricted = h.restrict_injection(0, &child).unwrap();
    assert_eq!(restricted, coarse);
    let reconstructed = h.transition(0).unwrap().prolongate_p1(&restricted).unwrap();
    assert_ne!(reconstructed, child);
}
#[test]
fn failed_operations_leave_levels_unchanged() {
    let base = Mesh2::unit_square(1).unwrap();
    let mut h = RefinementHistory::new(base.clone());
    assert!(h.refine(&[]).is_err());
    assert!(h.refine(&[CellId(999)]).is_err());
    assert!(h.refine(&[CellId(0), CellId(0)]).is_err());
    assert!(h.coarsen_last().is_err());
    assert_eq!(h.level(), 0);
    assert_eq!(h.mesh().cells(), base.cells());
    assert!(h.mesh_at(1).is_err());
    assert!(h.mesh_at(usize::MAX).is_err());
    assert!(h.transition(usize::MAX).is_err());
    h.refine(&[CellId(0)]).unwrap();
    let fine = h.mesh().clone();
    assert!(h.restrict_injection(0, &[]).is_err());
    assert!(h
        .restrict_injection(0, &vec![f64::NAN; fine.vertices().len()])
        .is_err());
    assert!(h
        .restrict_injection(1, &vec![0.; fine.vertices().len()])
        .is_err());
    assert!(h.refine(&[CellId(usize::MAX)]).is_err());
    assert_eq!(h.level(), 1);
    assert_eq!(h.mesh().cells(), fine.cells());
}
#[test]
fn unrepresentable_midpoint_does_not_modify_history() {
    let x = 1e16;
    let mesh = Mesh2::new(
        vec![p(x, x), p(x + 2., x), p(x, x + 2.)],
        vec![[VertexId(0), VertexId(1), VertexId(2)]],
    )
    .unwrap();
    let mut h = RefinementHistory::new(mesh.clone());
    assert!(h.refine(&[CellId(0)]).is_err());
    assert_eq!(h.level(), 0);
    assert_eq!(h.mesh().vertices(), mesh.vertices());
}
#[test]
fn can_refine_differently_after_backtracking() {
    let mut h = RefinementHistory::new(Mesh2::unit_square(3).unwrap());
    h.refine(&[CellId(0)]).unwrap();
    let a = h.mesh().clone();
    h.coarsen_last().unwrap();
    h.refine(&[CellId(17)]).unwrap();
    let b = h.mesh();
    assert_ne!(a.vertices(), b.vertices());
    assert_eq!(h.level(), 1);
    assert_eq!(
        h.transition(0).unwrap().children(CellId(17)).unwrap().len(),
        4
    );
}
