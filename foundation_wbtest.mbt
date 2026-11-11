///|
test "constructors enforce finite bounded parameters" {
  let c = Copula::new(Clayton, 2.0)
  assert_eq(c.theta, 2.0)
  assert_eq(independence().theta, 0.0)
  let bad = try {
    ignore(Copula::new(Gumbel, 0.5))
    false
  } catch {
    _ => true
  }
  assert_true(bad)
  let nan = try {
    ignore(Copula::new(Frank, 0.0 / 0.0))
    false
  } catch {
    _ => true
  }
  assert_true(nan)
  validate_uniforms([{ x: 0.3, y: 0.6, }], 1)
  assert_eq(checked_probability(1.0 + 1.0e-12), 1.0)
}
