///|
fn close(a : Double, b : Double, tolerance : Double) -> Unit raise {
  assert_true((a - b).abs() <= tolerance)
}

///|
test "analytic reference values and exact margins" {
  let c = Copula::new(Clayton, 2.0)
  close(c.cdf(0.5, 0.5), 1.0 / 7.0.sqrt(), 1.0e-12)
  close(c.density(0.5, 0.5), 1.481003649342278, 1.0e-12)
  for family in [Independence, Clayton, Gumbel, Frank] {
    let t = match family {
      Independence => 0.0
      _ => 2.0
    }
    let model = Copula::new(family, t)
    assert_eq(model.cdf(1.0, 0.3), 0.3)
    assert_eq(model.cdf(0.0, 0.3), 0.0)
    for i in 1..<20 {
      let u = i.to_double() / 20.0
      let v = 0.43
      let h = 1.0e-5
      close(
        (model.cdf(u + h, v) - model.cdf(u - h, v)) / (2.0 * h),
        model.conditional(u, v),
        1.0e-7,
      )
      close(
        (model.conditional(u, v + h) - model.conditional(u, v - h)) / (2.0 * h),
        model.density(u, v),
        1.0e-7,
      )
    }
  }
}

///|
test "rotations preserve margins and positive rectangle masses" {
  for r in [R0, R90, R180, R270] {
    let c = Copula::new(Clayton, 3.0, rotation=r)
    for i in 0..<10 {
      for j in 0..<10 {
        let u = i.to_double() / 10.0
        let v = j.to_double() / 10.0
        let mass = c.cdf(u + 0.1, v + 0.1) -
          c.cdf(u, v + 0.1) -
          c.cdf(u + 0.1, v) +
          c.cdf(u, v)
        assert_true(mass >= -1.0e-12)
      }
    }
  }
}
