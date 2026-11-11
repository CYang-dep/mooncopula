///|
/// Percentile summary uses the type-7 sample quantile and sample variance.
pub struct SampleSummary {
  count : Int32
  mean : Double
  standard_deviation : Double
  lower : Double
  median : Double
  upper : Double
}

///|
/// Summarize finite replicates. Coverage is a percentile range, not a guarantee.
pub fn sample_summary(
  values : Array[Double],
  coverage~ : Double = 0.95,
) -> SampleSummary raise CopulaError {
  finite(coverage, "coverage")
  if coverage <= 0.0 || coverage >= 1.0 {
    raise InvalidInput("coverage must be strictly between zero and one")
  }
  let margin = EmpiricalMargin::new(values)
  let mut mean = 0.0
  let mut m2 = 0.0
  for i in 0..<values.length() {
    let delta = values[i] - mean
    mean = mean + delta / (i + 1).to_double()
    m2 = m2 + delta * (values[i] - mean)
  }
  finite(mean, "replicate mean")
  finite(m2, "replicate sum of squares")
  let alpha = (1.0 - coverage) / 2.0
  {
    count: values.length(),
    mean,
    standard_deviation: (m2 / (values.length() - 1).to_double()).max(0.0).sqrt(),
    lower: margin.linear_quantile(alpha),
    median: margin.linear_quantile(0.5),
    upper: margin.linear_quantile(1.0 - alpha),
  }
}

///|
/// Copy paired rows by zero-based indices. The original array is never aliased.
pub fn resample_pairs(
  data : Array[Pair],
  indices : Array[Int32],
) -> Array[Pair] raise CopulaError {
  validate_pairs(data, 2)
  if indices.length() != data.length() {
    raise InvalidInput("a resampling plan must contain exactly n indices")
  }
  let result : Array[Pair] = []
  for index in indices {
    if index < 0 || index >= data.length() {
      raise InvalidInput("resampling index outside original sample")
    }
    result.push(data[index])
  }
  result
}

///|
/// Each raw-data replicate is reranked; no fixed original margins are reused.
pub struct BootstrapReplicate {
  model : Copula
  log_likelihood : Double
  at_boundary : Bool
  upper_joint_probability : Double
}

///|
/// Deterministic replicates, summaries and visible boundary-fit frequency.
pub struct BootstrapResult {
  replicates : Array[BootstrapReplicate]
  theta : SampleSummary
  upper_joint_probability : SampleSummary
  boundary_fits : Int32
}

///|
/// Resample raw paired observations using caller-supplied plans.
/// This library never invents an RNG seed or silently drops failed replicates.
/// Plan generation determines the resampling design; iid percentile inference
/// needs iid observations and sufficient replicates, which this API cannot prove.
pub fn bootstrap_raw(
  family : Family,
  raw : Array[Pair],
  plans : Array[Array[Int32]],
  threshold : Pair,
  rotation~ : Rotation = R0,
  options~ : FitOptions = default_fit_options(family),
  coverage~ : Double = 0.95,
) -> BootstrapResult raise CopulaError {
  validate_pairs(raw, 3)
  probability(threshold.x, "threshold u")
  probability(threshold.y, "threshold v")
  if plans.length() < 2 || plans.length() > 128 {
    raise BudgetExceeded("bootstrap needs 2 to 128 plans")
  }
  if plans.length() * raw.length() > 262144 {
    raise BudgetExceeded("bootstrap total draw budget exceeded")
  }
  // Validate the entire design before any expensive fitting.
  for plan in plans {
    if plan.length() != raw.length() {
      raise InvalidInput("each bootstrap plan must contain n indices")
    }
    for index in plan {
      if index < 0 || index >= raw.length() {
        raise InvalidInput("bootstrap index outside original sample")
      }
    }
  }
  finite(coverage, "coverage")
  if coverage <= 0.0 || coverage >= 1.0 {
    raise InvalidInput("coverage must be strictly between zero and one")
  }
  let replicates : Array[BootstrapReplicate] = []
  let parameters : Array[Double] = []
  let probabilities : Array[Double] = []
  let mut boundary_fits = 0
  for plan in plans {
    let sampled = resample_pairs(raw, plan)
    let uniforms = pseudo_observations(sampled)
    let fitted = fit(family, uniforms, rotation~, options~)
    if !fitted.converged {
      raise NumericalFailure("bootstrap fit exhausted its iteration budget")
    }
    let risk = fitted.model.risk(threshold.x, threshold.y).upper_both
    if fitted.at_boundary {
      boundary_fits = boundary_fits + 1
    }
    replicates.push({
      model: fitted.model,
      log_likelihood: fitted.log_likelihood,
      at_boundary: fitted.at_boundary,
      upper_joint_probability: risk,
    })
    parameters.push(fitted.model.theta)
    probabilities.push(risk)
  }
  {
    replicates,
    theta: sample_summary(parameters, coverage~),
    upper_joint_probability: sample_summary(probabilities, coverage~),
    boundary_fits,
  }
}
