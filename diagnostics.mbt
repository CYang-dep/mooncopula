///|
/// Grid discrepancy is descriptive; no uncalibrated goodness-of-fit p-value.
pub struct GridDiagnostic {
  points : Int32
  maximum_absolute : Double
  root_mean_square : Double
  sample_scaled_mean_square : Double
}

///|
/// Compare an analytic CDF with a step empirical CDF at caller-chosen points.
/// The finite sample statistic depends on the chosen grid and rank policy.
pub fn Copula::grid_diagnostic(
  self : Copula,
  data : Array[Pair],
  grid : Array[Pair],
) -> GridDiagnostic raise CopulaError {
  let empirical = EmpiricalCopula::new(data)
  validate_pairs(grid, 1)
  let mut maximum = 0.0
  let mut squares = 0.0
  for point in grid {
    let residual = empirical.cdf(point.x, point.y) - self.cdf(point.x, point.y)
    maximum = maximum.max(residual.abs())
    squares = squares + residual * residual
  }
  let mean_square = squares / grid.length().to_double()
  {
    points: grid.length(),
    maximum_absolute: maximum,
    root_mean_square: mean_square.sqrt(),
    sample_scaled_mean_square: data.length().to_double() * mean_square,
  }
}

///|
/// Interior rectangular grid; limits bound both memory and diagnostic work.
pub fn probability_grid(divisions : Int32) -> Array[Pair] raise CopulaError {
  if divisions < 2 || divisions > 64 {
    raise BudgetExceeded("grid divisions must be in [2,64]")
  }
  let result : Array[Pair] = []
  let denominator = divisions.to_double()
  for i in 1..<divisions {
    for j in 1..<divisions {
      result.push({
        x: i.to_double() / denominator,
        y: j.to_double() / denominator,
      })
    }
  }
  result
}

///|
/// Univariate descriptive uniformity statistics, without a significance claim.
pub struct UniformDiagnostic {
  count : Int32
  mean : Double
  variance : Double
  ks_distance : Double
  cramer_von_mises : Double
}

///|
/// Finite-sample empirical-vs-uniform KS and Cramer-von Mises statistics.
/// Variance uses n, not n-1. Inputs may include exact endpoints.
pub fn uniform_diagnostic(
  values : Array[Double],
) -> UniformDiagnostic raise CopulaError {
  if values.length() < 2 || values.length() > 4096 {
    raise BudgetExceeded("uniform diagnostic needs 2 to 4096 values")
  }
  for value in values {
    probability(value, "diagnostic value")
  }
  let sorted = values.copy()
  sorted.sort()
  let n = values.length().to_double()
  let mut mean = 0.0
  let mut m2 = 0.0
  let mut ks = 0.0
  let mut cvm = 1.0 / (12.0 * n)
  for i in 0..<sorted.length() {
    let x = sorted[i]
    let delta = x - mean
    mean = mean + delta / (i + 1).to_double()
    m2 = m2 + delta * (x - mean)
    ks = ks.max((i + 1).to_double() / n - x).max(x - i.to_double() / n)
    let residual = x - (2.0 * i.to_double() + 1.0) / (2.0 * n)
    cvm = cvm + residual * residual
  }
  {
    count: values.length(),
    mean,
    variance: m2 / n,
    ks_distance: ks,
    cramer_von_mises: cvm,
  }
}

///|
/// Residual marginal diagnostics and remaining rank correlation.
pub struct RosenblattDiagnostic {
  first : UniformDiagnostic
  second : UniformDiagnostic
  residual_spearman : Double
}

///|
/// Transform observations then diagnose both marginal uniformity and dependence.
/// Parameters fitted on the same data invalidate ordinary null p-values;
/// this API intentionally returns descriptive statistics only.
pub fn Copula::residual_diagnostic(
  self : Copula,
  data : Array[Pair],
) -> RosenblattDiagnostic raise CopulaError {
  let transformed = self.rosenblatt(data)
  {
    first: uniform_diagnostic(transformed.map(fn(p) { p.x })),
    second: uniform_diagnostic(transformed.map(fn(p) { p.y })),
    residual_spearman: spearman(transformed),
  }
}
