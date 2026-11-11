///|
/// Exact finite-sample concordance counts, including coordinate-specific ties.
pub struct Concordance {
  concordant : Int32
  discordant : Int32
  x_only_ties : Int32
  y_only_ties : Int32
  both_ties : Int32
  tau_b : Double
} derive(Debug)

///|
pub extend Concordance with @debug.Debug::{to_repr}

///|
/// Kendall tau-b; fails when either coordinate is constant.
pub fn concordance(data : Array[Pair]) -> Concordance raise CopulaError {
  validate_pairs(data, 2)
  let mut concordant = 0
  let mut discordant = 0
  let mut x_only_ties = 0
  let mut y_only_ties = 0
  let mut both_ties = 0
  for i in 0..<data.length() {
    for j in 0..<i {
      let dx = data[i].x.compare(data[j].x)
      let dy = data[i].y.compare(data[j].y)
      if dx == 0 && dy == 0 {
        both_ties = both_ties + 1
      } else if dx == 0 {
        x_only_ties = x_only_ties + 1
      } else if dy == 0 {
        y_only_ties = y_only_ties + 1
      } else if dx == dy {
        concordant = concordant + 1
      } else {
        discordant = discordant + 1
      }
    }
  }
  let untied = concordant + discordant
  let denominator = ((untied + x_only_ties).to_double() *
  (untied + y_only_ties).to_double()).sqrt()
  if denominator == 0.0 {
    raise InvalidInput("tau-b undefined for constant coordinate")
  }
  {
    concordant,
    discordant,
    x_only_ties,
    y_only_ties,
    both_ties,
    tau_b: (concordant - discordant).to_double() / denominator,
  }
}

///|
/// Spearman rho is Pearson correlation of midranks, including tied data.
pub fn spearman(data : Array[Pair]) -> Double raise CopulaError {
  let ranked = pseudo_observations(data)
  let mut sx = 0.0
  let mut sy = 0.0
  for point in ranked {
    sx = sx + point.x
    sy = sy + point.y
  }
  let mx = sx / ranked.length().to_double()
  let my = sy / ranked.length().to_double()
  let mut xx = 0.0
  let mut yy = 0.0
  let mut xy = 0.0
  for point in ranked {
    let dx = point.x - mx
    let dy = point.y - my
    xx = xx + dx * dx
    yy = yy + dy * dy
    xy = xy + dx * dy
  }
  if xx == 0.0 || yy == 0.0 {
    raise InvalidInput("rho undefined for constant coordinate")
  }
  xy / (xx * yy).sqrt()
}

///|
/// Finite-threshold tail probabilities, not asymptotic tail estimates.
pub struct TailDiagnostic {
  threshold : Double
  lower_joint : Double
  upper_joint : Double
  lower_coefficient : Double
  upper_coefficient : Double
  observations : Int32
} derive(Debug)

///|
pub extend TailDiagnostic with @debug.Debug::{to_repr}

///|
/// threshold in (0,0.5]; lower uses <=q, upper uses >1-q.
pub fn empirical_tails(
  data : Array[Pair],
  threshold : Double,
) -> TailDiagnostic raise CopulaError {
  validate_uniforms(data, 2)
  interior(threshold, "tail threshold")
  if threshold > 0.5 {
    raise InvalidInput("tail threshold exceeds 0.5")
  }
  let mut lower = 0
  let mut upper = 0
  for p in data {
    if p.x <= threshold && p.y <= threshold {
      lower = lower + 1
    }
    if p.x > 1.0 - threshold && p.y > 1.0 - threshold {
      upper = upper + 1
    }
  }
  let lower_joint = lower.to_double() / data.length().to_double()
  let upper_joint = upper.to_double() / data.length().to_double()
  {
    threshold,
    lower_joint,
    upper_joint,
    lower_coefficient: lower_joint / threshold,
    upper_coefficient: upper_joint / threshold,
    observations: data.length(),
  }
}
