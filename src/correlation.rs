use crate::errors::DistanceError;
use crate::utils::mean;

pub enum CorrelationMethods {
    PearsonCentered,
    PearsonUncentered,
    SquaredPearson,
}

pub fn pearson_corr_centered(x: &[f64], y: &[f64]) -> Result<f64, DistanceError> {
    if x.len() != y.len() {
        return Err(DistanceError::UnequalLengths);
    }
    let mean_x = mean(x);
    let mean_y = mean(y);

    let corr_tup =
        x.iter()
            .zip(y.iter())
            .fold((0.0, 0.0, 0.0), |(prod, x_diff_sq, y_diff_sq), (&a, &b)| {
                let xd = a - mean_x;
                let xds = xd.powi(2);
                let yd = b - mean_y;
                let yds = yd.powi(2);
                let prod_of_diffs = xd * yd;

                (prod + prod_of_diffs, x_diff_sq + xds, y_diff_sq + yds)
            });

    let corr = corr_tup.0 / (corr_tup.1.sqrt() * corr_tup.2.sqrt());
    Ok(corr)
}

pub fn pearson_corr_uncentered(x: &[f64], y: &[f64]) -> Result<f64, DistanceError> {
    if x.len() != y.len() {
        return Err(DistanceError::UnequalLengths);
    }
    let mean_x = mean(x);
    let mean_y = mean(y);

    let corr_tup =
        x.iter()
            .zip(y.iter())
            .fold((0.0, 0.0, 0.0), |(prod, x_diff_sq, y_diff_sq), (&a, &b)| {
                let xd = a - mean_x;
                let xds = xd.powi(2);
                let yd = b - mean_y;
                let yds = yd.powi(2);
                let prod_of_diffs = a * b;

                (prod + prod_of_diffs, x_diff_sq + xds, y_diff_sq + yds)
            });

    let corr = corr_tup.0 / (corr_tup.1.sqrt() * corr_tup.2.sqrt());
    Ok(corr)
}

pub fn squared_pearson_corr(x: &[f64], y: &[f64]) -> Result<f64, DistanceError> {
    if x.len() != y.len() {
        return Err(DistanceError::UnequalLengths);
    }
    let mean_x = mean(x);
    let mean_y = mean(y);

    let corr_tup =
        x.iter()
            .zip(y.iter())
            .fold((0.0, 0.0, 0.0), |(prod, x_squared, y_squared), (&a, &b)| {
                let xd = a - mean_x;
                let yd = b - mean_y;
                let xs = a * a;
                let ys = b * b;
                let prod_of_diffs = xd * yd;

                (prod + prod_of_diffs, x_squared + xs, y_squared + ys)
            });

    let corr = corr_tup.0 / (corr_tup.1.sqrt() * corr_tup.2.sqrt());
    Ok(corr)
}
