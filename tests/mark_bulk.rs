use dealii_rust::{estimator::mark_bulk, grid::CellId};
#[test]
fn marks_minimal_prefix_of_squared_indicators_with_stable_ties() {
    let values = [1., 2., 2., 0.];
    assert_eq!(mark_bulk(&values, 0.0).unwrap(), vec![]);
    assert_eq!(mark_bulk(&values, 4.0 / 9.0).unwrap(), vec![CellId(1)]);
    assert_eq!(mark_bulk(&values, 0.5).unwrap(), vec![CellId(1), CellId(2)]);
    assert_eq!(
        mark_bulk(&values, 1.).unwrap(),
        vec![CellId(1), CellId(2), CellId(0)]
    );
}
#[test]
fn handles_zero_extreme_scales_and_invalid_inputs() {
    assert!(mark_bulk(&[], 0.5).unwrap().is_empty());
    assert!(mark_bulk(&[0., 0.], 1.).unwrap().is_empty());
    assert_eq!(
        mark_bulk(&[f64::MAX, f64::MAX / 2.], 0.7).unwrap(),
        vec![CellId(0)]
    );
    assert_eq!(
        mark_bulk(&[f64::MIN_POSITIVE, f64::MIN_POSITIVE], 1.).unwrap(),
        vec![CellId(0), CellId(1)]
    );
    for theta in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(mark_bulk(&[1.], theta).is_err());
    }
    for bad in [f64::NAN, f64::INFINITY, -1.] {
        assert!(mark_bulk(&[bad], 0.5).is_err());
    }
}
