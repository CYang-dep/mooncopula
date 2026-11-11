///|
fn log_add(a : Double, b : Double) -> Double {
  let m = a.max(b)
  m + @math.ln(@math.exp(a - m) + @math.exp(b - m))
}

///|
/// log(u^-theta+v^-theta-1), evaluated without overflowing powers.
fn clayton_sum(t : Double, u : Double, v : Double) -> Double {
  let a = -t * @math.ln(u)
  let b = -t * @math.ln(v)
  let m = a.max(b)
  m + @math.ln(@math.exp(a - m) + @math.exp(b - m) - @math.exp(-m))
}

///|
fn gumbel_sum(t : Double, u : Double, v : Double) -> Double {
  log_add(t * @math.ln(-@math.ln(u)), t * @math.ln(-@math.ln(v)))
}

///|
/// Stable Frank denominator; avoids subtracting nearly equal unit values.
fn frank_den(t : Double, u : Double, v : Double) -> Double {
  @math.exp(-t) - @math.exp(-t * u) - @math.exp(-t * v) * -@math.expm1(-t * u)
}

///|
fn base_cdf(c : Copula, u : Double, v : Double) -> Double {
  match c.family {
    Independence => u * v
    Clayton => @math.exp(-clayton_sum(c.theta, u, v) / c.theta)
    Gumbel => @math.exp(-@math.exp(gumbel_sum(c.theta, u, v) / c.theta))
    Frank =>
      // D+A*B and D have the same sign. The ratio stays positive.
      -@math.ln(frank_den(c.theta, u, v) / @math.expm1(-c.theta)) / c.theta
  }
}

///|
/// Closed-square CDF with exact margins and explicit finite validation.
pub fn Copula::cdf(
  self : Copula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  probability(u, "u")
  probability(v, "v")
  if u == 0.0 || v == 0.0 {
    return 0.0
  }
  if u == 1.0 {
    return v
  }
  if v == 1.0 {
    return u
  }
  let value = match self.rotation {
    R0 => base_cdf(self, u, v)
    R90 => v - base_cdf(self, 1.0 - u, v)
    R180 => u + v - 1.0 + base_cdf(self, 1.0 - u, 1.0 - v)
    R270 => u - base_cdf(self, u, 1.0 - v)
  }
  checked_probability(value)
}

///|
fn base_log_density(c : Copula, u : Double, v : Double) -> Double {
  let t = c.theta
  let lu = @math.ln(u)
  let lv = @math.ln(v)
  match c.family {
    Independence => 0.0
    Clayton =>
      @math.ln_1p(t) -
      (1.0 + t) * (lu + lv) -
      (2.0 + 1.0 / t) * clayton_sum(t, u, v)
    Gumbel => {
      let ls = gumbel_sum(t, u, v)
      let z = @math.exp(ls / t)
      -z -
      lu -
      lv +
      (t - 1.0) * (@math.ln(-lu) + @math.ln(-lv)) +
      (2.0 / t - 2.0) * ls +
      @math.ln_1p((t - 1.0) / z)
    }
    Frank =>
      @math.ln((-t * @math.expm1(-t)).abs()) -
      t * (u + v) -
      2.0 * @math.ln(frank_den(t, u, v).abs())
  }
}

///|
/// Log density is preferred for fitting, including near singular boundaries.
pub fn Copula::log_density(
  self : Copula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  interior(u, "u")
  interior(v, "v")
  let (a, b) = match self.rotation {
    R0 => (u, v)
    R90 => (1.0 - u, v)
    R180 => (1.0 - u, 1.0 - v)
    R270 => (u, 1.0 - v)
  }
  let result = base_log_density(self, a, b)
  if result.is_nan() || result.is_inf() {
    raise NumericalFailure("non-finite log density")
  }
  result
}

///|
/// Density on the open square; boundary singularities are rejected.
pub fn Copula::density(
  self : Copula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  let value = @math.exp(self.log_density(u, v))
  if value.is_inf() || value.is_nan() {
    raise NumericalFailure("density overflow; use log_density")
  }
  value
}

///|
fn base_conditional(c : Copula, u : Double, v : Double) -> Double {
  let t = c.theta
  match c.family {
    Independence => v
    Clayton =>
      @math.exp(
        -(t + 1.0) * @math.ln(u) - (1.0 + 1.0 / t) * clayton_sum(t, u, v),
      )
    Gumbel => {
      let ls = gumbel_sum(t, u, v)
      @math.exp(
        -@math.exp(ls / t) +
        (1.0 / t - 1.0) * ls +
        (t - 1.0) * @math.ln(-@math.ln(u)) -
        @math.ln(u),
      )
    }
    Frank => @math.exp(-t * u) * @math.expm1(-t * v) / frank_den(t, u, v)
  }
}

///|
/// P(V<=v | U=u) = derivative C(u,v)/du, for 0<u<1 and 0<=v<=1.
pub fn Copula::conditional(
  self : Copula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  interior(u, "conditioning u")
  probability(v, "v")
  if v == 0.0 || v == 1.0 {
    return v
  }
  let value = match self.rotation {
    R0 => base_conditional(self, u, v)
    R90 => base_conditional(self, 1.0 - u, v)
    R180 => 1.0 - base_conditional(self, 1.0 - u, 1.0 - v)
    R270 => 1.0 - base_conditional(self, u, 1.0 - v)
  }
  checked_probability(value)
}
