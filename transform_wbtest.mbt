///|
test "conditional inverse and Rosenblatt round trips across families rotations" {
  for f in [Clayton, Gumbel, Frank] {
    for r in [R0, R90, R180, R270] {
      let c = Copula::new(f, 2.0, rotation=r)
      for i in 1..<16 {
        let u = i.to_double() / 16.0
        let q = 1.0 - u
        let inv = c.inverse_conditional(u, q)
        assert_true(inv.residual <= 1.0e-10)
        close(c.conditional(u, inv.value), q, 1.0e-10)
      }
      let data = [{ x: 0.2, y: 0.4, }, { x: 0.7, y: 0.8, }]
      let mapped = c.simulate(data)
      let restored = c.rosenblatt(mapped)
      close(restored[0].y, data[0].y, 1.0e-10)
      close(restored[1].y, data[1].y, 1.0e-10)
    }
  }
}

///|
test "inverse iteration limit remains observable" {
  let c = Copula::new(Gumbel, 2.0)
  let inv = c.inverse_conditional(0.3, 0.123, max_iterations=1)
  assert_eq(inv.iterations, 1)
  assert_true(inv.residual > 1.0e-5)
  assert_eq(c.inverse_conditional(0.3, 0.0).value, 0.0)
  assert_eq(c.inverse_conditional(0.3, 1.0).value, 1.0)
}
