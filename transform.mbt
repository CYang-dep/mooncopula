///|
/// Bounded inversion result exposes numerical residual instead of claiming exactness.
pub struct Inversion {
  value : Double
  residual : Double
  iterations : Int32
} derive(Debug)

///|
pub extend Inversion with @debug.Debug::{to_repr}

///|
/// Invert conditional CDF by monotone bisection (up to 80 iterations).
/// Endpoint q is supported, conditioning u remains strictly interior.
pub fn Copula::inverse_conditional(
  self : Copula,
  u : Double,
  q : Double,
  tolerance~ : Double = 1.0e-10,
  max_iterations~ : Int32 = 80,
) -> Inversion raise CopulaError {
  interior(u, "conditioning u")
  probability(q, "conditional quantile")
  finite(tolerance, "tolerance")
  if tolerance < 1.0e-13 || tolerance > 1.0e-3 {
    raise InvalidInput("tolerance must be in [1e-13,1e-3]")
  }
  if max_iterations < 1 || max_iterations > 80 {
    raise BudgetExceeded("inverse iterations must be in [1,80]")
  }
  if q == 0.0 || q == 1.0 {
    return { value: q, residual: 0.0, iterations: 0, }
  }
  let mut lo = 0.0
  let mut hi = 1.0
  let mut value = 0.5
  let mut residual = 1.0
  let mut used = 0
  for i in 0..<max_iterations {
    value = (lo + hi) / 2.0
    let actual = self.conditional(u, value)
    residual = (actual - q).abs()
    used = i + 1
    if residual <= tolerance {
      break
    }
    if actual < q {
      lo = value
    } else {
      hi = value
    }
  }
  { value, residual, iterations: used, }
}

///|
/// Rosenblatt transform (u,v) -> (u,H(v|u)). Does not mutate input.
pub fn Copula::rosenblatt(
  self : Copula,
  data : Array[Pair],
) -> Array[Pair] raise CopulaError {
  validate_uniforms(data, 1)
  let result : Array[Pair] = []
  for point in data {
    result.push({ x: point.x, y: self.conditional(point.x, point.y), })
  }
  result
}

///|
/// Inverse Rosenblatt map. Input consists of caller-provided open-square uniforms.
/// Enforces residual <=1e-8; concentrated/extreme cases may fail explicitly.
pub fn Copula::simulate(
  self : Copula,
  uniforms : Array[Pair],
) -> Array[Pair] raise CopulaError {
  validate_uniforms(uniforms, 1)
  let result : Array[Pair] = []
  for point in uniforms {
    let inverse = self.inverse_conditional(point.x, point.y)
    if inverse.residual > 1.0e-8 {
      raise NumericalFailure(
        "conditional inversion did not achieve residual tolerance",
      )
    }
    interior(inverse.value, "simulated v")
    result.push({ x: point.x, y: inverse.value, })
  }
  result
}
