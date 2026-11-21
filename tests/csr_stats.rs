use dealii_rust::lac::CooBuilder;
#[test]
fn logical_csr_storage_counts_pointers_values_and_stored_zeros() {
    let mut b = CooBuilder::new(3, 3);
    b.add(0, 0, 3.).unwrap();
    b.add(1, 0, 2.).unwrap();
    b.add(1, 0, -2.).unwrap();
    b.add(2, 2, 1.).unwrap();
    let matrix = b.into_csr().unwrap();
    assert_eq!(matrix.stored_entries(), 3);
    assert_eq!(matrix.get(1, 0).unwrap(), 0.);
    let idx = std::mem::size_of::<usize>() as u128;
    let value = std::mem::size_of::<f64>() as u128;
    assert_eq!(matrix.logical_storage_bytes(), 7 * idx + 3 * value);
}
#[test]
fn empty_matrix_still_has_row_pointer() {
    let matrix = CooBuilder::new(0, 0).into_csr().unwrap();
    assert_eq!(matrix.stored_entries(), 0);
    assert_eq!(
        matrix.logical_storage_bytes(),
        std::mem::size_of::<usize>() as u128
    );
}
