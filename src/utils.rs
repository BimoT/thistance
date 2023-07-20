use crate::errors::DistanceError;

use std::str::FromStr;

#[derive(Default, Clone, Copy)]
pub enum LogKind {
    #[default]
    Ln,
    Log2,
    Log10,
}

impl FromStr for LogKind {
    type Err = DistanceError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "log2" => Ok(LogKind::Log2),
            "log10" => Ok(LogKind::Log10),
            "log" | "ln" => Ok(LogKind::Ln),
            _ => Err(DistanceError::ParseLogKindError(s.to_string())),
        }
    }
}

pub(crate) fn logmatch(kind: &LogKind, x: f64) -> f64 {
    match kind {
        LogKind::Ln => x.ln(),
        LogKind::Log2 => x.log2(),
        LogKind::Log10 => x.log10(),
    }
}

/// Rudimentary function to calculate the sign of a float.
/// Returns 0 if it's 0, -1 if negative, 1 if positive
pub(crate) fn sign(x: f64) -> i8 {
    if x == 0.0 {
        0
    } else if x < 0.0 {
        -1
    } else {
        1
    }
}

/// Simple function to calculate the mean of a vector of floats
pub(crate) fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

pub fn print_dist_methods() {
    let v = vec![
        "euclidean",
        "manhattan",
        "minkowski",
        "chebyshev",
        "sorensen",
        "gower",
        "soergel",
        "kulczynski_d",
        "canberra",
        "lorentzian",
        "intersection",
        "non_intersection",
        "wave_hedges",
        "czekanowski",
        "motyka",
        "kulczynski_s",
        "tanimoto",
        "ruzicka",
        "inner_product",
        "harmonic_mean",
        "cosine",
        "kumar_hassebrook",
        "jaccard",
        "dice",
        "fidelity",
        "bhattacharyya",
        "hellinger",
        "matusita",
        "squared_chord",
        "squared_euclidean",
        "pearson_chi_sq",
        "neyman_chi_sq",
        "squared_chi_sq",
        "prob_sym_chi_sq",
        "divergence_sq",
        "clark_sq",
        "additive_symm_chi_sq",
        "kullback_leibler",
        "jeffreys",
        "k_divergence",
        "topsoe",
        "jensen_shannon",
        "jensen_difference",
        "taneja",
        "kumar_johnson",
        "avg",
    ];
    println!("{:?}", v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sign_positive() {
        let a: f64 = 5.0;
        let b: f64 = 0.0;
        let c: f64 = -5.0;

        assert_eq!(sign(a), 1);
        assert_eq!(sign(b), 0);
        assert_eq!(sign(c), -1);
    }
}
