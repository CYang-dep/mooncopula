///|
test "tau-b and rho distinguish ties and reversals" {
  let data = [
    { x: 1.0, y: 4.0, },
    { x: 2.0, y: 3.0, },
    { x: 3.0, y: 2.0, },
    { x: 4.0, y: 1.0, },
  ]
  close(concordance(data).tau_b, -1.0, 1.0e-15)
  close(spearman(data), -1.0, 1.0e-15)
  let tied = [{ x: 1.0, y: 1.0, }, { x: 1.0, y: 2.0, }, { x: 2.0, y: 2.0, }]
  let c = concordance(tied)
  assert_eq(c.concordant, 1)
  assert_eq(c.x_only_ties, 1)
  assert_eq(c.y_only_ties, 1)
  close(c.tau_b, 0.5, 1.0e-15)
}

///|
test "tail threshold empirical counts use explicit inequalities" {
  let data = [
    { x: 0.1, y: 0.2, },
    { x: 0.9, y: 0.8, },
    { x: 0.3, y: 0.6, },
    { x: 0.2, y: 0.2, },
  ]
  let t = empirical_tails(data, 0.2)
  close(t.lower_joint, 0.5, 1.0e-15)
  close(t.upper_joint, 0.0, 1.0e-15)
}
