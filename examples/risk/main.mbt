///|
/// Compare dependence families, check event calibration and resample raw rows.
fn main raise {
  let raw : Array[@mooncopula.Pair] = [
    { x: 1.0, y: 3.0, },
    { x: 2.0, y: 1.0, },
    { x: 3.0, y: 4.0, },
    { x: 4.0, y: 2.0, },
    { x: 5.0, y: 6.0, },
    { x: 6.0, y: 5.0, },
    { x: 7.0, y: 8.0, },
    { x: 8.0, y: 7.0, },
  ]
  let uniforms = @mooncopula.pseudo_observations(raw)
  let ranking = @mooncopula.compare_models(
    @mooncopula.standard_candidates(),
    uniforms,
  )
  println("risk: comparable candidate likelihoods")
  for candidate in ranking {
    println(
      "\{candidate.label}: AIC=\{candidate.fit.aic}; weight=\{candidate.weight}",
    )
  }
  let best = ranking[0].fit.model
  let calibration = best.calibration_curve(
    uniforms,
    [{ x: 0.75, y: 0.75, }],
    @mooncopula.UpperBoth,
  )
  println(
    "upper_joint_forecast=\{calibration[0].forecast}; observed=\{calibration[0].observed.observed}",
  )
  let residual = best.residual_diagnostic(uniforms)
  println("residual_rho=\{residual.residual_spearman}")
  // Four fixed plans demonstrate mechanics only, not adequate inferential coverage.
  let plans = [
    [0, 1, 2, 3, 4, 5, 6, 7],
    [0, 0, 2, 3, 4, 5, 7, 7],
    [1, 2, 2, 3, 5, 6, 6, 7],
    [0, 1, 3, 3, 4, 4, 5, 7],
  ]
  let bootstrap = @mooncopula.bootstrap_raw(@mooncopula.Clayton, raw, plans, {
    x: 0.75,
    y: 0.75,
  })
  println(
    "bootstrap_replicates=\{bootstrap.replicates.length()}; theta_median=\{bootstrap.theta.median}",
  )
}
