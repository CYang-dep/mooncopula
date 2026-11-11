///|
/// An explicit model candidate, including Frank branch and rotation.
pub(all) struct Candidate {
  label : String
  family : Family
  rotation : Rotation
  options : FitOptions
}

///|
/// Ranked model fit plus normalized AIC weight.
pub struct RankedModel {
  label : String
  fit : FitResult
  delta_aic : Double
  weight : Double
}

///|
/// A compact standard model set; negative dependence is represented explicitly.
pub fn standard_candidates() -> Array[Candidate] {
  let result : Array[Candidate] = [
    {
      label: "independence",
      family: Independence,
      rotation: R0,
      options: default_fit_options(Independence),
    },
  ]
  for (label, family) in [("clayton", Clayton), ("gumbel", Gumbel)] {
    for
      (suffix, rotation) in [
        ("", R0),
        ("-90", R90),
        ("-180", R180),
        ("-270", R270),
      ] {
      result.push({
        label: label + suffix,
        family,
        rotation,
        options: default_fit_options(family),
      })
    }
  }
  result.push({
    label: "frank-positive",
    family: Frank,
    rotation: R0,
    options: default_fit_options(Frank),
  })
  result.push({
    label: "frank-negative",
    family: Frank,
    rotation: R0,
    options: { lower: -30.0, upper: -0.01, tolerance: 1.0e-6, iterations: 80, },
  })
  result
}

///|
/// AIC comparison on identical data. Failure/nonconvergence is explicit;
/// unsuccessful candidates are never silently discarded.
pub fn compare_models(
  candidates : Array[Candidate],
  data : Array[Pair],
) -> Array[RankedModel] raise CopulaError {
  validate_uniforms(data, 3)
  if candidates.length() < 1 || candidates.length() > 32 {
    raise BudgetExceeded("candidate count must be in [1,32]")
  }
  let scored : Array[(String, FitResult)] = []
  for i in 0..<candidates.length() {
    let candidate = candidates[i]
    if candidate.label.length() == 0 || candidate.label.length() > 80 {
      raise InvalidInput("candidate label must have 1 to 80 characters")
    }
    for j in 0..<i {
      if candidates[j].label == candidate.label {
        raise InvalidInput("duplicate candidate label")
      }
    }
    let fitted = fit(
      candidate.family,
      data,
      rotation=candidate.rotation,
      options=candidate.options,
    )
    if !fitted.converged {
      raise NumericalFailure(
        "candidate optimizer did not converge: " + candidate.label,
      )
    }
    scored.push((candidate.label, fitted))
  }
  scored.sort_by(fn(a, b) {
    let cmp = a.1.aic.compare(b.1.aic)
    if cmp == 0 {
      a.0.compare(b.0)
    } else {
      cmp
    }
  })
  let best = scored[0].1.aic
  let mut total = 0.0
  for item in scored {
    total = total + @math.exp(-0.5 * (item.1.aic - best))
  }
  scored.map(fn(item) {
    let delta = item.1.aic - best
    {
      label: item.0,
      fit: item.1,
      delta_aic: delta,
      weight: @math.exp(-0.5 * delta) / total,
    }
  })
}

///|
/// Fit a specific family from raw paired observations using midrank uniforms.
pub fn fit_raw(
  family : Family,
  raw : Array[Pair],
  rotation~ : Rotation = R0,
) -> FitResult raise CopulaError {
  fit(family, pseudo_observations(raw), rotation~)
}
