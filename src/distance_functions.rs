use crate::errors::DistanceError;
use crate::utils::{logmatch, sign, LogKind};

/// Calculates the Euclidean distance
pub fn euclidean(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q)
        .fold(0.0f64, |acc, (&a, &b)| {
            let diff = (a - b).abs().powi(2);
            acc + diff
        })
        // acc + (a - b).abs().powi(2))
        // .map(|(&a, &b)| (a - b).abs().powi(2))
        // .sum::<f64>()
        .sqrt();
    Ok(dist)
}

/// Calculates the Manhattan distance
pub fn manhattan(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q)
        .map(|(&a, &b)| (a - b).abs())
        .sum();
    Ok(dist)
}

/// Calculates the Minkowski distance
pub fn minkowski(p: &[f64], q: &[f64], n: f64, test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q)
        .fold(0.0f64, |acc, (&a, &b)| {
            let diff = (a - b).abs();
            acc + diff
        })
        .powf(1.0 / n);
    Ok(dist)
}

// TODO: remove unwrap asap!
// TODO: map followed by reduce? should be done more elegantly
/// Calculates the Chebyshev distance
pub fn chebyshev(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .map(|(&a, &b)| (a - b).abs())
        .reduce(f64::max)
        .unwrap();
    Ok(dist)
}

/// Calculates the Sorensen distance
pub fn sorensen(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let distances: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc1, acc2), (&a, &b)| {
                let diff = (a - b).abs();
                let sum = a + b;
                (acc1 + diff, acc2 + sum)
            });
    if distances.1 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = distances.0 / distances.1;
        Ok(dist)
    }
}

pub fn gower(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    let plen = p.len();
    if plen != q.len() {
        return Err(DistanceError::UnequalLengths);
    }
    if plen == 0 {
        return Err(DistanceError::ZeroLength);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let distbase: f64 = p
        .iter()
        .zip(q.iter())
        .map(|(&a, &b)| (a - b).abs())
        .sum();
    let dist = (1 / plen) as f64 * distbase;
    Ok(dist)
}

pub fn soergel(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let distbase: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc1, acc2), (&a, &b)| {
                let diff = (a - b).abs();
                let maxpoint = a.max(b);
                (acc1 + diff, acc2 + maxpoint)
            });
    if distbase.1 == 0.0 {
        Ok(0.0)
    } else {
        let dist = distbase.0 / distbase.1;
        Ok(dist)
    }
}

pub fn kulczynski_d(
    p: &[f64],
    q: &[f64],
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let distbase: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc1, acc2), (&a, &b)| {
                let diff = (a - b).abs();
                let minpoint = a.min(b);
                if minpoint == 0.0 {
                    (acc1 + diff, acc2 + epsilon)
                } else {
                    (acc1 + diff, acc2 + minpoint)
                }
            });

    if distbase.1 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = distbase.0 / distbase.1;
        Ok(dist)
    }
}

pub fn kulczynski_s(
    p: &[f64],
    q: &[f64],
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    let k_d = kulczynski_d(p, q, epsilon, test_nan)?;

    if k_d == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = 1.0 / k_d;
        Ok(dist)
    }
}

pub fn canberra(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |accum, (&a, &b)| {
            let diff: f64 = (a - b).abs();
            let sum: f64 = a + b;
            if (diff == 0.0) || (sum == 0.0) {
                accum + 0.0f64
            } else {
                accum + diff / sum
            }
        });

    Ok(dist)
}

pub fn lorentzian(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |accum, (&a, &b)| {
            let diff: f64 = (a - b).abs();
            let val = log(1.0 + diff);
            accum + val
        });
    Ok(dist)
}

pub fn intersection_dist(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |accum, (&a, &b)| {
            let minpoint = a.min(b);
            accum + minpoint
        });
    Ok(dist)
}

pub fn non_intersection_dist(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    let i_d = intersection_dist(p, q, test_nan)?;
    if i_d == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = 1.0 / i_d;
        Ok(dist)
    }
}

pub fn wave_hedges(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |accum, (&a, &b)| {
            let diff: f64 = (a - b).abs();
            // let max_point = if a >= b { a } else { b };
            let max_point = a.max(b);
            if (diff == 0.0) || max_point == 0.0 {
                accum + 0.0
            } else {
                accum + diff / max_point
            }
        });
    Ok(dist)
}

pub fn czekanowski(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dtup: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc1, acc2), (&a, &b)| {
                let diff: f64 = (a - b).abs();
                let sum = a + b; // TODO: check if this needs to be absolute or not
                (acc1 + diff, acc2 + sum)
            });

    if dtup.1 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = dtup.0 / dtup.1;
        Ok(dist)
    }
}

pub fn motyka(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dtup: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc1, acc2), (&a, &b)| {
                let sum = a + b;
                let minimum = a.min(b);

                (acc1 + minimum, acc2 + sum)
            });

    if dtup.1 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = 1.0 - (dtup.0 / dtup.1);
        Ok(dist)
    }
}

pub fn tanimoto(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    soergel(p, q, test_nan)
}

pub fn ruzicka(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    let soerg = soergel(p, q, test_nan)?;
    Ok(1.0 - soerg)
}

pub fn inner_product(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let prod = a * b;
            acc + prod
        });

    Ok(dist)
}

pub fn harmonic_mean_dist(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist_half: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let prod = a * b;
            let sum = a + b;
            if (prod == 0.0) || (sum == 0.0) {
                acc + 0.0
            } else {
                acc + prod / sum
            }
        });
    let dist: f64 = dist_half * 2.0;

    Ok(dist)
}

pub fn cosine_dist(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist_tup: (f64, f64, f64) = p.iter().zip(q.iter()).fold(
        (0.0f64, 0.0f64, 0.0f64),
        |(acc_prod, acc_asq, acc_bsq), (&a, &b)| {
            let prod = a * b;
            let a_square = a.powi(2);
            let b_square = b.powi(2);

            (acc_prod + prod, acc_asq + a_square, acc_bsq + b_square)
        },
    );

    let dist2 = (dist_tup.1).sqrt() * (dist_tup.2).sqrt();
    if dist2 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = dist_tup.0 / dist2;
        Ok(dist)
    }
}

pub fn kumar_hassebrook(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist_tup: (f64, f64, f64) = p.iter().zip(q.iter()).fold(
        (0.0f64, 0.0f64, 0.0f64),
        |(acc_prod, acc_asq, acc_bsq), (&a, &b)| {
            let prod = a * b;
            let a_square = a.powi(2);
            let b_square = b.powi(2);

            (acc_prod + prod, acc_asq + a_square, acc_bsq + b_square)
        },
    );

    let dist2 = dist_tup.1 + dist_tup.2 - dist_tup.0;
    if dist2 == 0.0 {
        Ok(0.0)
    } else {
        let dist = dist_tup.0 / dist2;
        Ok(dist)
    }
}

pub fn jaccard(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    let kh = kumar_hassebrook(p, q, test_nan)?;
    Ok(1.0 - kh)
}

pub fn dice_dist(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist_tup: (f64, f64, f64) = p.iter().zip(q.iter()).fold(
        (0.0f64, 0.0f64, 0.0f64),
        |(acc_diff, acc_asq, acc_bsq), (&a, &b)| {
            let diff_square = (a - b).powi(2);
            let a_square = a.powi(2);
            let b_square = b.powi(2);

            (
                acc_diff + diff_square,
                acc_asq + a_square,
                acc_bsq + b_square,
            )
        },
    );

    let dist2 = dist_tup.1 + dist_tup.2;
    if dist2 == 0.0 {
        Err(DistanceError::DivideByZero)
    } else {
        let dist = dist_tup.0 / dist2;
        Ok(dist)
    }
}

pub fn fidelity(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let prod = (a * b).sqrt();
            acc + prod
        });

    Ok(dist)
}

pub fn bhattacharyya(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let fid_value = fidelity(p, q, test_nan)?;
    let dist = {
        if fid_value == 0.0 {
            -log(fid_value + epsilon)
        } else {
            -log(fid_value)
        }
    };

    Ok(dist)
}

pub fn hellinger(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let fid = fidelity(p, q, test_nan)?;

    let dist = 2.0 * (1.0 - fid).sqrt();
    Ok(dist)
}

pub fn matusita(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let fid = fidelity(p, q, test_nan)?;

    let dist = 2.0 - (2.0 * fid).sqrt();
    Ok(dist)
}

pub fn squared_chord(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let diff = (a.sqrt() - b.sqrt()).powi(2);
            acc + diff
        });

    Ok(dist)
}

pub fn squared_euclidean(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let diff = (a - b).powi(2);
            acc + diff
        });

    Ok(dist)
}

pub fn pearson_chi_sq(
    p: &[f64],
    q: &[f64],
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let d = if b == 0.0 { epsilon } else { b };
            acc + ((a - b).powi(2) / d)
        });

    Ok(dist)
}

pub fn neyman_chi_sq(
    p: &[f64],
    q: &[f64],
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let d = if a == 0.0 { epsilon } else { a };
            acc + ((a - b).powi(2) / d)
        });

    Ok(dist)
}

pub fn squared_chi_sq(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let abdiff = (a - b).powi(2);
            let absum = a + b;

            if abdiff == 0.0 || absum == 0.0 {
                acc + 0.0
            } else {
                acc + (abdiff / absum)
            }
        });

    Ok(dist)
}

pub fn prob_symm_chi_sq(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = 2.0 * squared_chi_sq(p, q, test_nan)?;

    Ok(dist)
}

pub fn divergence_sq(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let abdiff = (a - b).powi(2);
            let absum = (a + b).powi(2);

            if abdiff == 0.0 || absum == 0.0 {
                acc + 0.0
            } else {
                acc + (abdiff / absum)
            }
        });

    Ok(dist)
}

pub fn clark_sq(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let abdiff = (a - b).abs();
            let absum = a + b;

            if abdiff == 0.0 || absum == 0.0 {
                acc + 0.0
            } else {
                acc + (abdiff / absum).powi(2)
            }
        })
        .sqrt();

    Ok(dist)
}

pub fn additive_symm_chi_sq(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let absum = a + b;
            let abprod = a * b;

            if abprod == 0.0 || absum == 0.0 {
                acc + 0.0
            } else {
                let difference = (a - b).powi(2) * (absum / abprod);
                acc + difference
            }
        });

    Ok(dist)
}

pub fn kullback_leibler_distance(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            if a == 0.0 && b == 0.0 {
                acc + 0.0
            } else {
                let abratio = if b == 0.0 { a / epsilon } else { a / b };
                if abratio == 0.0 {
                    acc + 0.0
                } else {
                    acc + (a * log(abratio))
                }
            }
        });

    Ok(dist)
}

pub fn jeffreys(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let abratio = if b == 0.0 { a / epsilon } else { a / b };

            let z = if abratio == 0.0 { epsilon } else { abratio };

            acc + ((a - b) * log(z))
        });

    Ok(dist)
}

pub fn k_divergence(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            if a == 0.0 && b == 0.0 {
                acc + 0.0
            } else {
                acc + log((2.0 * a) / (a + b))
            }
        });

    Ok(dist)
}

pub fn topsoe(p: &[f64], q: &[f64], unit: &LogKind, test_nan: bool) -> Result<f64, DistanceError> {
    // FIX: divide by zero error here, check jensen_shannon
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            if a == 0.0 && b == 0.0 {
                acc + 0.0
            } else {
                let absum = a + b;

                let left_part = a * log((2.0 * a) / absum);
                let right_part = b * log((2.0 * b) / absum);
                let sum = left_part + right_part;

                acc + sum
            }
        });

    Ok(dist)
}

pub fn jensen_shannon(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let absum = a + b;

            let left_part = if a == 0.0 || absum == 0.0 {
                0.0
            } else {
                a * log((2.0 * a) / absum)
                // match unit {
                //     LogKind::Ln => a * ((2.0 * a) / absum).ln(),
                //     LogKind::Log2 => a * ((2.0 * a) / absum).log2(),
                //     LogKind::Log10 => a * ((2.0 * a) / absum).log10(),
                // }
            };

            let right_part = if b == 0.0 || absum == 0.0 {
                0.0
            } else {
                b * log((2.0 * b) / absum)
                // match unit {
                //     LogKind::Ln => b * ((2.0 * b) / absum).ln(),
                //     LogKind::Log2 => b * ((2.0 * b) / absum).log2(),
                //     LogKind::Log10 => b * ((2.0 * b) / absum).log10(),
                // }
            };
            let sum = left_part + right_part;

            acc + sum
        })
        * 0.5;

    Ok(dist)
}

pub fn jensen_difference(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);
    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let absum = a + b;

            match (sign(a), sign(b), sign(absum)) {
                (0, 0, 0) => acc + 0.0,
                (0, 1, 1) => acc + ((b * log(b)) / 2.0) - ((absum / 2.0) * log(absum / 2.0)),
                (1, 0, 1) => acc + ((a * log(a)) / 2.0) - ((absum / 2.0) * log(absum / 2.0)),
                (0, 0, 1) => acc + -((absum / 2.0) * log(absum / 2.0)),
                (1, 1, 0) => acc + ((a * log(a)) + (b * log(b))) / 2.0,
                (1, 0, 0) => acc + (a * log(a)) / 2.0,
                (0, 1, 0) => acc + (b * log(b)) / 2.0,
                (_, _, _) => {
                    acc + ((a * log(a)) + (b * log(b))) / 2.0 - (absum / 2.0) * log(absum / 2.0)
                },
            }
        });

    Ok(dist)
}

pub fn taneja(
    p: &[f64],
    q: &[f64],
    unit: &LogKind,
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    // this is either ln, log2 or log10
    let log = |a| logmatch(unit, a);
    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            if a == 0.0 && b == 0.0 {
                acc + 0.0
            } else {
                let absum = a + b;
                let denominator = 2.0 * (a * b).sqrt();

                if denominator == 0.0 {
                    acc + (absum / 2.0) * log(absum / epsilon)
                } else {
                    acc + (absum / 2.0) * log(absum / denominator)
                }
            }
        });

    Ok(dist)
}

pub fn kumar_johnson(
    p: &[f64],
    q: &[f64],
    epsilon: f64,
    test_nan: bool,
) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let dist: f64 = p
        .iter()
        .zip(q.iter())
        .fold(0.0f64, |acc, (&a, &b)| {
            let divisor = 2.0 * (a * b).powf(1.5);
            let div = if divisor == 0.0 { epsilon } else { divisor };
            acc + (a.powi(2) - b.powi(2)).powi(2) / div
        });

    Ok(dist)
}

pub fn avg(p: &[f64], q: &[f64], test_nan: bool) -> Result<f64, DistanceError> {
    if p.len() != q.len() {
        return Err(DistanceError::UnequalLengths);
    }

    if test_nan {
        for (a, b) in p.iter().zip(q.iter()) {
            if a.is_nan() || b.is_nan() {
                return Err(DistanceError::NaNError);
            }
        }
    }

    let distbase: (f64, f64) =
        p.iter()
            .zip(q.iter())
            .fold((0.0f64, 0.0f64), |(acc, maxi), (&a, &b)| {
                let absdiff = (a - b).abs();
                if absdiff > maxi {
                    //HACK: i can't write `maxi = absdiff` so i have to take the difference first
                    let new_max_difference = absdiff - maxi;
                    (acc + absdiff, maxi + new_max_difference)
                } else {
                    (acc + absdiff, maxi + 0.0)
                }
            });

    let dist = (distbase.0 + distbase.1) / 2.0;

    Ok(dist)
}
