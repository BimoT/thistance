use crate::utils::{logmatch, LogKind};

pub fn entropy(p: &[f64], unit: &LogKind) -> f64 {
    let log = |x| logmatch(unit, x);
    let res = p.iter().fold(0.0, |acc, &a| {
        if a > 0.0 {
            let accum = a * log(a);
            acc + accum
        } else {
            acc + 0.0
        }
    });
    -res
}

pub fn joint_entropy(joint_probabilities: &[f64], unit: &LogKind) -> f64 {
    entropy(joint_probabilities, unit)
}

pub fn conditional_entropy(
    joint_probabilities: &[f64],
    probabilities: &[f64],
    unit: &LogKind,
) -> f64 {
    // Using the chain rule: H(X | Y) = H(X,Y) - H(Y)
    // Note: it is important that the probabilities vector corresponds to Y
    joint_entropy(joint_probabilities, unit) - entropy(probabilities, unit)
}

pub fn mutual_information(x: &[f64], y: &[f64], xy: &[f64], unit: &LogKind) -> f64 {
    // Using the identity: I(X,Y) = H(X) + H(Y) - H(X, Y)
    (entropy(x, unit) + entropy(y, unit)) - joint_entropy(xy, unit)
}
