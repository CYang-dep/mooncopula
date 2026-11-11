///|
/// Owned empirical distribution on the unit square. Construction copies input.
pub struct EmpiricalCopula {
  priv points : Array[Pair]
}

///|
/// Already ranked uniforms; use from_raw for finite raw measurements.
pub fn EmpiricalCopula::new(
  data : Array[Pair],
) -> EmpiricalCopula raise CopulaError {
  validate_uniforms(data, 2)
  { points: data.copy(), }
}

///|
pub fn EmpiricalCopula::from_raw(
  data : Array[Pair],
  policy~ : TiePolicy = Midrank,
) -> EmpiricalCopula raise CopulaError {
  EmpiricalCopula::new(pseudo_observations(data, policy~))
}

///|
pub fn EmpiricalCopula::size(self : EmpiricalCopula) -> Int32 {
  self.points.length()
}

///|
/// Return a defensive copy of pseudo-observations.
pub fn EmpiricalCopula::observations(self : EmpiricalCopula) -> Array[Pair] {
  self.points.copy()
}

///|
/// Empirical joint CDF, both coordinates use <=. Margins are discrete.
pub fn EmpiricalCopula::cdf(
  self : EmpiricalCopula,
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  probability(u, "u")
  probability(v, "v")
  let mut count = 0
  for p in self.points {
    if p.x <= u && p.y <= v {
      count = count + 1
    }
  }
  count.to_double() / self.points.length().to_double()
}

///|
/// Probability of (lo_u,hi_u] x (lo_v,hi_v], including exact discrete ties.
pub fn EmpiricalCopula::rectangle(
  self : EmpiricalCopula,
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
  let mut count = 0
  for p in self.points {
    if p.x > lo_u && p.x <= hi_u && p.y > lo_v && p.y <= hi_v {
      count = count + 1
    }
  }
  count.to_double() / self.points.length().to_double()
}

///|
/// Smoothed checkerboard CDF uses ordinal (untied) ranks only.
/// Each rank cell spreads mass uniformly over its n-by-n grid cell.
pub fn checkerboard_cdf(
  raw : Array[Pair],
  u : Double,
  v : Double,
) -> Double raise CopulaError {
  probability(u, "u")
  probability(v, "v")
  let data = pseudo_observations(raw, policy=RejectTies)
  let n = data.length().to_double()
  let mut total = 0.0
  for point in data {
    let rx = point.x * (n + 1.0)
    let ry = point.y * (n + 1.0)
    let fx = (n * u - rx + 1.0).max(0.0).min(1.0)
    let fy = (n * v - ry + 1.0).max(0.0).min(1.0)
    total = total + fx * fy
  }
  checked_probability(total / n)
}
