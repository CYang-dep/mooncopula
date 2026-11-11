///|
/// Explicit failures; numerical singularities are not silently clipped.
pub(all) suberror CopulaError {
  InvalidInput(String)
  BudgetExceeded(String)
  NumericalFailure(String)
} derive(Debug)

///|
/// Supported bivariate analytic families.
pub(all) enum Family {
  Independence
  Clayton
  Gumbel
  Frank
} derive(Debug)

///|
/// Rotation of a base copula on the unit square.
pub(all) enum Rotation {
  R0
  R90
  R180
  R270
} derive(Debug)

///|
/// Construction validates parameters. Fields cannot be changed by consumers.
pub struct Copula {
  family : Family
  theta : Double
  rotation : Rotation
} derive(Debug)

///|
/// Owned pair of finite observations. Validation is performed by operations.
pub(all) struct Pair {
  x : Double
  y : Double
} derive(Debug)

///|
fn finite(x : Double, label : String) -> Unit raise CopulaError {
  if x.is_nan() || x.is_inf() {
    raise InvalidInput(label + " must be finite")
  }
}

///|
fn probability(x : Double, label : String) -> Unit raise CopulaError {
  finite(x, label)
  if x < 0.0 || x > 1.0 {
    raise InvalidInput(label + " must be in [0,1]")
  }
}

///|
fn interior(x : Double, label : String) -> Unit raise CopulaError {
  probability(x, label)
  if x <= 0.0 || x >= 1.0 {
    raise InvalidInput(label + " must be strictly between 0 and 1")
  }
}

///|
fn validate_pairs(
  data : Array[Pair],
  minimum : Int32,
) -> Unit raise CopulaError {
  if data.length() < minimum {
    raise InvalidInput("insufficient paired observations")
  }
  if data.length() > 4096 {
    raise BudgetExceeded("at most 4096 paired observations")
  }
  for item in data {
    finite(item.x, "x")
    finite(item.y, "y")
  }
}

///|
pub fn validate_uniforms(
  data : Array[Pair],
  minimum : Int32,
) -> Unit raise CopulaError {
  validate_pairs(data, minimum)
  for item in data {
    interior(item.x, "u")
    interior(item.y, "v")
  }
}

///|
pub fn checked_probability(x : Double) -> Double raise CopulaError {
  if x.is_nan() || x.is_inf() || x < -1.0e-10 || x > 1.0000000001 {
    raise NumericalFailure("probability outside numerical tolerance")
  }
  x.max(0.0).min(1.0)
}

///|
/// Create an analytic model in its documented numerically bounded domain.
/// Clayton: [0.01,30], Gumbel: [1,20], Frank: [-30,-0.01] or [0.01,30].
/// Independence requires theta=0. Use it for the exact independence limit.
pub fn Copula::new(
  family : Family,
  theta : Double,
  rotation~ : Rotation = R0,
) -> Copula raise CopulaError {
  finite(theta, "theta")
  let valid = match family {
    Independence => theta == 0.0
    Clayton => theta >= 0.01 && theta <= 30.0
    Gumbel => theta >= 1.0 && theta <= 20.0
    Frank => theta.abs() >= 0.01 && theta.abs() <= 30.0
  }
  if !valid {
    raise InvalidInput("theta outside supported family range")
  }
  { family, theta, rotation, }
}

///|
/// Exact independent model.
pub fn independence() -> Copula {
  { family: Independence, theta: 0.0, rotation: R0, }
}

///|
pub extend CopulaError with @debug.Debug::{to_repr}

///|
pub extend Family with @debug.Debug::{to_repr}

///|
pub extend Rotation with @debug.Debug::{to_repr}

///|
pub extend Copula with @debug.Debug::{to_repr}

///|
pub extend Pair with @debug.Debug::{to_repr}
