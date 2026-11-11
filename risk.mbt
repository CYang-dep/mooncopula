///|
/// Joint event decomposition for uniform thresholds, evaluated in one model.
pub struct JointRisk {
  lower_both : Double
  upper_both : Double
  upper_either : Double
  x_above_y_below : Double
  x_below_y_above : Double
} derive(Debug)

///|
pub extend JointRisk with @debug.Debug::{to_repr}

///|
pub fn Copula::risk(
  self : Copula,
  u : Double,
  v : Double,
) -> JointRisk raise CopulaError {
  let c = self.cdf(u, v)
  {
    lower_both: c,
    upper_both: checked_probability(1.0 - u - v + c),
    upper_either: checked_probability(1.0 - c),
    x_above_y_below: checked_probability(v - c),
    x_below_y_above: checked_probability(u - c),
  }
}

///|
/// Continuous rectangle mass using inclusion-exclusion, with explicit order.
pub fn Copula::rectangle(
  self : Copula,
  lo_u : Double,
  hi_u : Double,
  lo_v : Double,
  hi_v : Double,
) -> Double raise CopulaError {
  probability(lo_u, "lo_u")
  probability(hi_u, "hi_u")
  probability(lo_v, "lo_v")
  probability(hi_v, "hi_v")
  if lo_u > hi_u || lo_v > hi_v {
    raise InvalidInput("reversed rectangle")
  }
  checked_probability(
    self.cdf(hi_u, hi_v) -
    self.cdf(lo_u, hi_v) -
    self.cdf(hi_u, lo_v) +
    self.cdf(lo_u, lo_v),
  )
}

///|
/// Expected waiting time under independent, identically distributed event periods.
/// Zero probability gives None, never an artificial finite return period.
pub fn return_period(
  event_probability : Double,
  period_length : Double,
) -> Double? raise CopulaError {
  probability(event_probability, "event probability")
  finite(period_length, "period length")
  if period_length <= 0.0 {
    raise InvalidInput("period length must be positive")
  }
  if event_probability == 0.0 {
    return None
  }
  let result = period_length / event_probability
  if result.is_inf() {
    raise NumericalFailure("return period overflow")
  }
  Some(result)
}

///|
/// P(V>v | U>u), a threshold event; differs from conditional(v|u) at a point.
pub fn Copula::conditional_exceedance(
  self : Copula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  probability(u, "u")
  probability(v, "v")
  if u == 1.0 {
    raise InvalidInput("conditioning event has zero probability")
  }
  checked_probability(self.risk(u, v).upper_both / (1.0 - u))
}

///|
/// Theoretical asymptotic same-tail coefficients; no inference from data.
pub fn Copula::tail_coefficients(self : Copula) -> (Double, Double) {
  let (lower, upper) = match self.family {
    Independence => (0.0, 0.0)
    Clayton => (@math.exp(-@math.ln(2.0) / self.theta), 0.0)
    Gumbel => (0.0, 2.0 - @math.exp(@math.ln(2.0) / self.theta))
    Frank => (0.0, 0.0)
  }
  match self.rotation {
    R0 => (lower, upper)
    R180 => (upper, lower)
    R90 | R270 => (0.0, 0.0)
  }
}
