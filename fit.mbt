///|
/// Fitting controls are validated on every call, including caller-built values.
pub(all) struct FitOptions {
  lower : Double
  upper : Double
  tolerance : Double
  iterations : Int32
}

///|
/// Defaults cover the positive branch; Frank's negative branch is explicit.
pub fn default_fit_options(family : Family) -> FitOptions {
  match family {
    Independence =>
      { lower: 0.0, upper: 0.0, tolerance: 1.0e-6, iterations: 80, }
    Clayton => { lower: 0.01, upper: 30.0, tolerance: 1.0e-6, iterations: 80, }
    Gumbel => { lower: 1.0, upper: 20.0, tolerance: 1.0e-6, iterations: 80, }
    Frank => { lower: 0.01, upper: 30.0, tolerance: 1.0e-6, iterations: 80, }
  }
}

///|
/// Result of coarse-grid bracketing plus golden-section local refinement.
pub struct FitResult {
  model : Copula
  log_likelihood : Double
  aic : Double
  bic : Double
  observations : Int32
  evaluations : Int32
  converged : Bool
  at_boundary : Bool
}

///|
/// Sum of log densities using compensated accumulation.
pub fn Copula::log_likelihood(
  self : Copula,
  data : Array[Pair],
) -> Double raise CopulaError {
  validate_uniforms(data, 2)
  let mut sum = 0.0
  let mut correction = 0.0
  for point in data {
    let y = self.log_density(point.x, point.y) - correction
    let next = sum + y
    correction = next - sum - y
    sum = next
  }
  finite(sum, "log likelihood")
  sum
}

///|
fn fit_score(
  family : Family,
  theta : Double,
  rotation : Rotation,
  data : Array[Pair],
) -> Double raise CopulaError {
  Copula::new(family, theta, rotation~).log_likelihood(data)
}

///|
/// Fit open-square data. Raw observations must first be transformed to ranks.
/// No global optimality or confidence-interval guarantee is implied.
pub fn fit(
  family : Family,
  data : Array[Pair],
  rotation~ : Rotation = R0,
  options~ : FitOptions = default_fit_options(family),
) -> FitResult raise CopulaError {
  validate_uniforms(data, 3)
  finite(options.lower, "lower bound")
  finite(options.upper, "upper bound")
  finite(options.tolerance, "fit tolerance")
  if options.tolerance < 1.0e-9 || options.tolerance > 0.01 {
    raise InvalidInput("fit tolerance must be in [1e-9,0.01]")
  }
  if options.iterations < 1 || options.iterations > 160 {
    raise BudgetExceeded("fit iterations must be in [1,160]")
  }
  ignore(Copula::new(family, options.lower, rotation~))
  ignore(Copula::new(family, options.upper, rotation~))
  if options.lower > options.upper {
    raise InvalidInput("reversed fit bounds")
  }
  match family {
    Frank =>
      if options.lower < 0.0 && options.upper > 0.0 {
        raise InvalidInput("Frank fit bounds must stay on one side of zero")
      }
    _ => ()
  }
  let mut best_theta = options.lower
  let mut best_score = fit_score(family, best_theta, rotation, data)
  let mut evaluations = 1
  let mut best_index = 0
  let width = options.upper - options.lower
  // Sampling all endpoints avoids dropping boundary optima during refinement.
  for i in 1..<=24 {
    let theta = options.lower + width * i.to_double() / 24.0
    let score = fit_score(family, theta, rotation, data)
    evaluations = evaluations + 1
    if score > best_score {
      best_theta = theta
      best_score = score
      best_index = i
    }
  }
  let mut lo = options.lower +
    width * (best_index - 1).max(0).to_double() / 24.0
  let mut hi = options.lower +
    width * (best_index + 1).min(24).to_double() / 24.0
  let ratio = 0.6180339887498949
  let mut a = hi - ratio * (hi - lo)
  let mut b = lo + ratio * (hi - lo)
  let mut fa = fit_score(family, a, rotation, data)
  let mut fb = fit_score(family, b, rotation, data)
  evaluations = evaluations + 2
  for _ in 0..<options.iterations {
    if hi - lo <= options.tolerance {
      break
    }
    if fa > fb {
      hi = b
      b = a
      fb = fa
      a = hi - ratio * (hi - lo)
      fa = fit_score(family, a, rotation, data)
    } else {
      lo = a
      a = b
      fa = fb
      b = lo + ratio * (hi - lo)
      fb = fit_score(family, b, rotation, data)
    }
    evaluations = evaluations + 1
  }
  if fa > best_score {
    best_theta = a
    best_score = fa
  }
  if fb > best_score {
    best_theta = b
    best_score = fb
  }
  let parameters = match family {
    Independence => 0.0
    _ => 1.0
  }
  {
    model: Copula::new(family, best_theta, rotation~),
    log_likelihood: best_score,
    aic: 2.0 * parameters - 2.0 * best_score,
    bic: @math.ln(data.length().to_double()) * parameters - 2.0 * best_score,
    observations: data.length(),
    evaluations,
    converged: hi - lo <= options.tolerance,
    at_boundary: (best_theta - options.lower).abs() <= options.tolerance ||
    (best_theta - options.upper).abs() <= options.tolerance,
  }
}
