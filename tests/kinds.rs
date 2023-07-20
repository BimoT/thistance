use std::str::FromStr;
use thistance::errors::DistanceError;

use thistance::distance_kinds::DistanceMeasures;

#[test]
fn bla() -> Result<(), DistanceError> {
    let euc = DistanceMeasures::Euclidean;
    let e = DistanceMeasures::from_str("euclidean")?;
    assert_eq!(e, euc);
    Ok(())
}
