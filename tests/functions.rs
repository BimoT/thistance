use thistance::distance::{distance, DistOpts};

#[test]
fn test_euclidean() {
    let p = [1.0, 2.0, 3.0, 4.0];
    let q = [5.0, 6.0, 7.0, 8.0];
    let d = distance(&p, &q, "euclidean", &DistOpts::new());

    assert_eq!(d.unwrap(), 8.0)
}

#[test]
fn min_euclidean() {
    let p = vec![f64::MIN, 2.0, 3.0, 4.0];
    let q = vec![5.0, 6.0, 7.0, 8.0];

    // let outcome = euclidean(&p, &q, test_na).unwrap();
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());

    assert!(!outcome.unwrap().is_normal())
}

#[test]
fn max_euclidean() {
    let p = vec![f64::MAX, 2.0, 3.0, 4.0];
    let q = vec![5.0, 6.0, 7.0, 8.0];

    // let outcome = euclidean(&p, &q, test_na).unwrap();
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());

    assert!(!outcome.unwrap().is_normal())
}
#[test]
fn na_euclidean() {
    let p = vec![f64::NAN, 2.0, 3.0, 4.0];
    let q = vec![5.0, 6.0, 7.0, 8.0];
    // let test_na = true;

    // let outcome = euclidean(&p, &q, test_na);
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new());
    assert!(outcome.is_err())
}

#[test]
fn nan_euclidean_err() {
    let p = vec![f64::NAN, 2.0, 3.0, 4.0];
    let q = vec![5.0, 6.0, 7.0, 8.0];
    // let test_na = true;

    // let outcome = euclidean(&p, &q, test_na);
    let outcome = distance(&p, &q, "euclidean", &DistOpts::new().nan(false));
    assert!(outcome.unwrap().is_nan())
}
