use thistance::distance::{distance, DistOpts};
#[test]
fn test_euclidean() {
    let p = [1.0, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let val: f64 = 8.0;
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert_eq!(outcome.unwrap(), val)
}

#[test]
fn test_zeros_euclidean() {
    let p = [0.0];
    let q = [0.0];
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert!(outcome.is_ok())
}

#[test]
fn test_min_euclidean() {
    let p = [f64::MIN, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert!(!outcome.unwrap().is_normal())
}

#[test]
fn test_max_euclidean() {
    let p = [f64::MAX, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert!(!outcome.unwrap().is_normal())
}

#[test]
fn test_nan_euclidean() {
    let p = [f64::NAN, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert!(outcome.is_err())
}

#[test]
fn test_nan_err_euclidean() {
    let p = [f64::NAN, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new().nan(false));
    assert!(outcome.unwrap().is_nan())
}
