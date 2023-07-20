use thiserror::Error;

#[derive(Error, Debug)]
pub enum DistanceError {
    #[error("The two input slices are required to be of equal length")]
    UnequalLengths,
    #[error("The distance cannot divide a value by zero")] // TODO: improve errormsg
    DivideByZero,
    #[error("One of the inputs has a length of zero and cannot be processed")]
    ZeroLength,
    #[error("Cannot parse the string `{0}` into a LogKind")]
    ParseLogKindError(String),
    #[error("Cannot parse the string `{0}` into a DistanceMeasure")]
    ParseDistKindError(String),
    #[error("There are NaN values in the input")]
    NaNError,
}
