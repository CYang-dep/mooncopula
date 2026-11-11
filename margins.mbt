///|
/// Owned sorted finite empirical margin. No parametric extrapolation.
pub struct EmpiricalMargin {
  priv sorted : Array[Double]
}

///|
pub fn EmpiricalMargin::new(
  values : Array[Double],
) -> EmpiricalMargin raise CopulaError {
  if values.length() < 2 {
    raise InvalidInput("margin requires at least two observations")
  }
  if values.length() > 4096 {
    raise BudgetExceeded("margin exceeds 4096 observations")
  }
  for x in values {
    finite(x, "margin observation")
  }
  let sorted = values.copy()
  sorted.sort()
  { sorted, }
}

///|
pub fn EmpiricalMargin::size(self : EmpiricalMargin) -> Int32 {
  self.sorted.length()
}

///|
/// Empirical CDF with upper-bound binary search; ties use <=.
pub fn EmpiricalMargin::cdf(
  self : EmpiricalMargin,
  x : Double,
) -> Double raise CopulaError {
  finite(x, "margin query")
  let mut lo = 0
  let mut hi = self.sorted.length()
  while lo < hi {
    let mid = lo + (hi - lo) / 2
    if self.sorted[mid] <= x {
      lo = mid + 1
    } else {
      hi = mid
    }
  }
  lo.to_double() / self.sorted.length().to_double()
}

///|
/// Left-continuous generalized inverse: x[ceil(n*p)-1], p=0 -> minimum.
/// Returning an observed value avoids imposing interpolation on discrete margins.
pub fn EmpiricalMargin::quantile(
  self : EmpiricalMargin,
  p : Double,
) -> Double raise CopulaError {
  probability(p, "margin probability")
  let index = ((p * self.sorted.length().to_double()).ceil().to_int() - 1).max(
    0,
  )
  self.sorted[index]
}

///|
/// Linear order-statistic quantile (type 7), separate from the step inverse.
pub fn EmpiricalMargin::linear_quantile(
  self : EmpiricalMargin,
  p : Double,
) -> Double raise CopulaError {
  probability(p, "margin probability")
  let position = p * (self.sorted.length() - 1).to_double()
  let lower = position.floor().to_int()
  let upper = (lower + 1).min(self.sorted.length() - 1)
  let fraction = position - lower.to_double()
  // Weighted sum avoids overflow in subtraction of opposite signed extremes.
  let value = (1.0 - fraction) * self.sorted[lower] +
    fraction * self.sorted[upper]
  if value.is_nan() || value.is_inf() {
    raise NumericalFailure("quantile overflow")
  }
  value
}

///|
/// Map simulated dependence uniforms back to observed marginal units.
pub fn restore_margins(
  data : Array[Pair],
  x_margin : EmpiricalMargin,
  y_margin : EmpiricalMargin,
) -> Array[Pair] raise CopulaError {
  validate_uniforms(data, 1)
  let result : Array[Pair] = []
  for point in data {
    result.push({
      x: x_margin.quantile(point.x),
      y: y_margin.quantile(point.y),
    })
  }
  result
}

///|
/// Full raw-units scenario workflow using input observations as empirical margins.
pub fn Copula::scenarios(
  self : Copula,
  raw : Array[Pair],
  uniforms : Array[Pair],
) -> Array[Pair] raise CopulaError {
  validate_pairs(raw, 2)
  let x_margin = EmpiricalMargin::new(raw.map(fn(p) { p.x }))
  let y_margin = EmpiricalMargin::new(raw.map(fn(p) { p.y }))
  restore_margins(self.simulate(uniforms), x_margin, y_margin)
}
