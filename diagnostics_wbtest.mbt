///|
test "grid discrepancy has independently counted residuals" {
  let data = [{ x: 0.25, y: 0.25, }, { x: 0.75, y: 0.75, }]
  let d = independence().grid_diagnostic(data, [{ x: 0.5, y: 0.5, }])
  close(d.maximum_absolute, 0.25, 1.0e-14)
  close(d.root_mean_square, 0.25, 1.0e-14)
  close(d.sample_scaled_mean_square, 0.125, 1.0e-14)
  assert_eq(probability_grid(3).length(), 4)
}

///|
test "uniform statistics match closed form order statistics" {
  let d = uniform_diagnostic([0.25, 0.75])
  close(d.mean, 0.5, 1.0e-14)
  close(d.variance, 0.0625, 1.0e-14)
  close(d.ks_distance, 0.25, 1.0e-14)
  close(d.cramer_von_mises, 1.0 / 24.0, 1.0e-14)
  let data = [{ x: 0.2, y: 0.7, }, { x: 0.5, y: 0.3, }, { x: 0.8, y: 0.5, }]
  let r = independence().residual_diagnostic(data)
  close(r.residual_spearman, -0.5, 1.0e-14)
}
