use std::str::FromStr;

use crate::distance_functions;
use crate::distance_kinds::DistanceMeasures;
use crate::errors::DistanceError;
use crate::utils::LogKind;

#[derive(Default, Clone, Copy)]
/// A DistOpts object is used to control the behaviour of the `distance` function.
/// A quick way to create a DistOpts object with sensible defaults is by calling `DistOpts::new()`.
/// Each field of DistOpts can be modified, with the functions `DistOpts.nan()`,
/// `DistOpts.log()`, `DistOpts.eps()` and `DistOpts.mink()`
pub struct DistOpts {
    test_nan: bool,
    log_kind: LogKind,
    epsilon: f64,
    minkowski_coefficient: f64,
}

impl DistOpts {
    /// Creates a default DistOpts instance, with the following defaults:
    ///     test_nan: true,
    ///     log_kind: LogKind::Ln,
    ///     epsilon: 0.00001,
    ///     minkowski_coefficient: 2.0
    pub fn new() -> DistOpts {
        DistOpts {
            test_nan: true,
            log_kind: LogKind::Ln,
            epsilon: 0.00001,
            minkowski_coefficient: 2.0,
        }
    }

    /// Tries to make a DistOpts instance
    /// Acceptable values for log_unit_str are:
    ///     "ln", "log2", "log10"
    pub fn try_build(
        test_nan: bool,
        log_unit_str: &str,
        epsilon: f64,
        minkowski_coefficient: f64,
    ) -> Result<Self, DistanceError> {
        let unit = LogKind::from_str(log_unit_str)?;
        Ok(DistOpts {
            test_nan,
            log_kind: unit,
            epsilon,
            minkowski_coefficient,
        })
    }

    /// Sets the test_nan value of the DistOpts
    pub fn nan(&mut self, new_test_nan: bool) -> Self {
        self.test_nan = new_test_nan;
        self.clone()
    }

    /// Sets the epsilon value of the DistOpts
    pub fn eps(&mut self, new_epsilon: f64) -> Self {
        self.epsilon = new_epsilon;
        self.clone()
    }

    /// Sets the minkowski coefficient of the DistOpts
    pub fn mink(&mut self, new_minkowski_coefficient: f64) -> Self {
        self.minkowski_coefficient = new_minkowski_coefficient;
        self.clone()
    }

    /// Tries to set the LogKind of the DistOpts. Acceptable values are:
    ///     "ln", "log2", "log10"
    pub fn log(&mut self, new_log_kind: &str) -> Result<Self, DistanceError> {
        let l: LogKind = new_log_kind.parse()?;
        self.log_kind = l;
        Ok(self.clone())
    }
}

//
pub fn distance(
    p: &[f64],
    q: &[f64],
    method: &str,
    options: &DistOpts,
) -> Result<f64, DistanceError> {
    let distance_method = method.parse::<DistanceMeasures>()?;
    let test_nan = options.test_nan;
    let unit = options.log_kind;
    let epsilon = options.epsilon;
    let mink = options.minkowski_coefficient;

    #[rustfmt::skip]
    let res = match distance_method {
        DistanceMeasures::Euclidean         => distance_functions::euclidean(p, q, test_nan),
        DistanceMeasures::Manhattan         => distance_functions::manhattan(p, q, test_nan),
        DistanceMeasures::Minkowski         => distance_functions::minkowski(p, q, mink, test_nan),
        DistanceMeasures::Chebyshev         => distance_functions::chebyshev(p, q, test_nan),
        DistanceMeasures::Sorensen          => distance_functions::sorensen(p, q, test_nan),
        DistanceMeasures::Gower             => distance_functions::gower(p, q, test_nan),
        DistanceMeasures::Soergel           => distance_functions::soergel(p, q, test_nan),
        DistanceMeasures::KulczynskiD       => distance_functions::kulczynski_d(p, q, epsilon, test_nan),
        DistanceMeasures::KulczynskiS       => distance_functions::kulczynski_s(p, q, epsilon, test_nan),
        DistanceMeasures::Canberra          => distance_functions::canberra(p, q, test_nan),
        DistanceMeasures::Lorentzian        => distance_functions::lorentzian(p, q, &unit, test_nan),
        DistanceMeasures::Intersection      => distance_functions::intersection_dist(p, q, test_nan),
        DistanceMeasures::NonIntersection   => distance_functions::non_intersection_dist(p, q, test_nan),
        DistanceMeasures::WaveHedges        => distance_functions::wave_hedges(p, q, test_nan),
        DistanceMeasures::Czekanowski       => distance_functions::czekanowski(p, q, test_nan),
        DistanceMeasures::Motyka            => distance_functions::motyka(p, q, test_nan),
        DistanceMeasures::Tanimoto          => distance_functions::tanimoto(p, q, test_nan),
        DistanceMeasures::Ruzicka           => distance_functions::ruzicka(p, q, test_nan),
        DistanceMeasures::InnerProduct      => distance_functions::inner_product(p, q, test_nan),
        DistanceMeasures::HarmonicMean      => distance_functions::harmonic_mean_dist(p, q, test_nan),
        DistanceMeasures::Cosine            => distance_functions::cosine_dist(p, q, test_nan),
        DistanceMeasures::KumarHassebrook   => distance_functions::kumar_hassebrook(p, q, test_nan),
        DistanceMeasures::Jaccard           => distance_functions::jaccard(p, q, test_nan),
        DistanceMeasures::Dice              => distance_functions::dice_dist(p, q, test_nan),
        DistanceMeasures::Fidelity          => distance_functions::fidelity(p, q, test_nan),
        DistanceMeasures::Bhattacharyya     => distance_functions::bhattacharyya(p, q, &unit, epsilon, test_nan),
        DistanceMeasures::Hellinger         => distance_functions::hellinger(p, q, test_nan),
        DistanceMeasures::Matusita          => distance_functions::matusita(p, q, test_nan),
        DistanceMeasures::SquaredChord      => distance_functions::squared_chord(p, q, test_nan),
        DistanceMeasures::SquaredEuclidean  => distance_functions::squared_euclidean(p, q, test_nan),
        DistanceMeasures::PearsonChiSq      => distance_functions::pearson_chi_sq(p, q, epsilon, test_nan),
        DistanceMeasures::NeymanChiSq       => distance_functions::neyman_chi_sq(p, q, epsilon, test_nan),
        DistanceMeasures::SquaredChiSq      => distance_functions::squared_chi_sq(p, q, test_nan),
        DistanceMeasures::ProbSymmChiSq     => distance_functions::prob_symm_chi_sq(p, q, test_nan),
        DistanceMeasures::DivergenceSq      => distance_functions::divergence_sq(p, q, test_nan),
        DistanceMeasures::ClarkSq           => distance_functions::clark_sq(p, q, test_nan),
        DistanceMeasures::AdditiveSymmChiSq => distance_functions::additive_symm_chi_sq(p, q, test_nan),
        DistanceMeasures::KullbackLeibler   => distance_functions::kullback_leibler_distance(p, q, &unit, epsilon, test_nan),
        DistanceMeasures::Jeffreys          => distance_functions::jeffreys(p, q, &unit, epsilon, test_nan),
        DistanceMeasures::KDivergence       => distance_functions::k_divergence(p, q, &unit, test_nan),
        DistanceMeasures::Topsoe            => distance_functions::topsoe(p, q, &unit, test_nan),
        DistanceMeasures::JensenShannon     => distance_functions::jensen_shannon(p, q, &unit, test_nan),
        DistanceMeasures::JensenDifference  => distance_functions::jensen_difference(p, q, &unit, test_nan),
        DistanceMeasures::Taneja            => distance_functions::taneja(p, q, &unit, epsilon, test_nan),
        DistanceMeasures::KumarJohnson      => distance_functions::kumar_johnson(p, q, epsilon, test_nan),
        DistanceMeasures::Avg               => distance_functions::avg(p, q, test_nan),
    };
    res
}
