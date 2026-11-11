///|
/// Threshold event convention: lower is <=; upper is >.
pub(all) enum ThresholdEvent {
  LowerBoth
  UpperBoth
  UpperEither
  FirstUpperSecondLower
  FirstLowerSecondUpper
}

///|
/// Wilson interval for a Bernoulli rate, conditional on iid observations.
pub struct RateInterval {
  successes : Int32
  trials : Int32
  observed : Double
  lower : Double
  upper : Double
}

///|
/// Caller chooses z (1.959963984540054 approximates a two-sided 95% interval).
/// No time-series dependence or fitted-parameter adjustment is implied.
pub fn wilson_interval(
  successes : Int32,
  trials : Int32,
  z~ : Double = 1.959963984540054,
) -> RateInterval raise CopulaError {
  if trials < 1 || successes < 0 || successes > trials {
    raise InvalidInput(
      "Bernoulli counts require 0<=successes<=trials and trials>0",
    )
  }
  finite(z, "Wilson z")
  if z <= 0.0 || z > 10.0 {
    raise InvalidInput("Wilson z must be in (0,10]")
  }
  let n = trials.to_double()
  let p = successes.to_double() / n
  let z2 = z * z
  let denominator = 1.0 + z2 / n
  let center = (p + z2 / (2.0 * n)) / denominator
  let half_width = z *
    (p * (1.0 - p) / n + z2 / (4.0 * n * n)).sqrt() /
    denominator
  {
    successes,
    trials,
    observed: p,
    lower: (center - half_width).max(0.0),
    upper: (center + half_width).min(1.0),
  }
}

///|
/// Observed threshold frequency and model forecast at a single threshold pair.
pub struct CalibrationPoint {
  threshold : Pair
  forecast : Double
  observed : RateInterval
  residual : Double
}

///|
fn matches_event(p : Pair, threshold : Pair, event : ThresholdEvent) -> Bool {
  match event {
    LowerBoth => p.x <= threshold.x && p.y <= threshold.y
    UpperBoth => p.x > threshold.x && p.y > threshold.y
    UpperEither => p.x > threshold.x || p.y > threshold.y
    FirstUpperSecondLower => p.x > threshold.x && p.y <= threshold.y
    FirstLowerSecondUpper => p.x <= threshold.x && p.y > threshold.y
  }
}

///|
/// Compare analytic forecasts with exact counts, including observed ties.
/// Wilson ranges describe observed Bernoulli frequencies only. In-sample model
/// fitting and rank transformations are not accounted for in their coverage.
pub fn Copula::calibration_curve(
  self : Copula,
  data : Array[Pair],
  thresholds : Array[Pair],
  event : ThresholdEvent,
  z~ : Double = 1.959963984540054,
) -> Array[CalibrationPoint] raise CopulaError {
  validate_uniforms(data, 2)
  validate_pairs(thresholds, 1)
  let result : Array[CalibrationPoint] = []
  for threshold in thresholds {
    let risk = self.risk(threshold.x, threshold.y)
    let forecast = match event {
      LowerBoth => risk.lower_both
      UpperBoth => risk.upper_both
      UpperEither => risk.upper_either
      FirstUpperSecondLower => risk.x_above_y_below
      FirstLowerSecondUpper => risk.x_below_y_above
    }
    let mut successes = 0
    for point in data {
      if matches_event(point, threshold, event) {
        successes = successes + 1
      }
    }
    let observed = wilson_interval(successes, data.length(), z~)
    result.push({
      threshold,
      forecast,
      observed,
      residual: observed.observed - forecast,
    })
  }
  result
}
