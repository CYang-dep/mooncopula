///|
test "model comparison has normalized deterministic weights" {
  let data = Copula::new(Clayton, 3.0).simulate(design_uniforms(127))
  let selected = compare_models(standard_candidates(), data)
  assert_eq(selected.length(), 11)
  let mut sum = 0.0
  for i in 0..<selected.length() {
    sum = sum + selected[i].weight
    if i > 0 {
      assert_true(selected[i].fit.aic >= selected[i - 1].fit.aic)
    }
  }
  close(sum, 1.0, 1.0e-12)
  assert_true(selected[0].fit.log_likelihood > 0.0)
}
