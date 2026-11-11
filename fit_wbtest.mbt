///|
fn design_uniforms(n : Int32) -> Array[Pair] {
  Array::makei(n, fn(i) {
    {
      x: (i.to_double() + 0.5) / n.to_double(),
      y: ((i * 37 % n).to_double() + 0.5) / n.to_double(),
    }
  })
}

///|
test "bounded likelihood fit recovers simulated Clayton dependence" {
  let truth = Copula::new(Clayton, 2.5)
  let data = truth.simulate(design_uniforms(127))
  let result = fit(Clayton, data)
  assert_true(result.converged)
  assert_true((result.model.theta - 2.5).abs() < 0.4)
  assert_true(
    result.log_likelihood >= Copula::new(Clayton, 1.0).log_likelihood(data),
  )
  close(result.aic, 2.0 - 2.0 * result.log_likelihood, 1.0e-12)
}

///|
test "fit exposes boundaries iteration exhaustion and invalid branch" {
  let data = Copula::new(Clayton, 4.0).simulate(design_uniforms(61))
  let options = { lower: 0.1, upper: 1.0, tolerance: 1.0e-9, iterations: 1, }
  let limited = fit(Clayton, data, options~)
  assert_true(!limited.converged)
  assert_true(limited.at_boundary)
  let bad = try {
    ignore(
      fit(Frank, data, options={
        lower: -2.0,
        upper: 2.0,
        tolerance: 1.0e-6,
        iterations: 60,
      }),
    )
    false
  } catch {
    _ => true
  }
  assert_true(bad)
}
