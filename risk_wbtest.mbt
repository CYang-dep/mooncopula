///|
test "independent event algebra and return-period assumptions" {
  let c = independence()
  let r = c.risk(0.8, 0.7)
  close(r.upper_both, 0.06, 1.0e-14)
  close(r.upper_either, 0.44, 1.0e-14)
  close(
    r.lower_both + r.upper_both + r.x_above_y_below + r.x_below_y_above,
    1.0,
    1.0e-14,
  )
  close(c.conditional_exceedance(0.8, 0.7), 0.3, 1.0e-14)
  close(c.rectangle(0.2, 0.5, 0.3, 0.9), 0.18, 1.0e-14)
  match return_period(0.1, 2.0) {
    Some(x) => close(x, 20.0, 1.0e-14)
    None => fail("missing period")
  }
  assert_true(return_period(0.0, 1.0) is None)
}

///|
test "rotated asymptotic tails" {
  let (l, u) = Copula::new(Clayton, 2.0).tail_coefficients()
  close(l, 0.5.sqrt(), 1.0e-14)
  close(u, 0.0, 1.0e-14)
  let (a, b) = Copula::new(Clayton, 2.0, rotation=R180).tail_coefficients()
  close(a, 0.0, 1.0e-14)
  close(b, l, 1.0e-14)
}
