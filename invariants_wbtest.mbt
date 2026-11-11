///|
test "all families rotations respect Frechet bounds monotonicity and nonnegative cell mass" {
  for
    (family, theta) in [
      (Independence, 0.0),
      (Clayton, 0.01),
      (Clayton, 30.0),
      (Gumbel, 1.0),
      (Gumbel, 20.0),
      (Frank, -30.0),
      (Frank, -0.01),
      (Frank, 0.01),
      (Frank, 30.0),
    ] {
    for rotation in [R0, R90, R180, R270] {
      let c = Copula::new(family, theta, rotation~)
      let mut mass = 0.0
      for i in 0..<10 {
        for j in 0..<10 {
          let u = i.to_double() / 10.0
          let v = j.to_double() / 10.0
          let p = c.cdf(u, v)
          assert_true(p >= (u + v - 1.0).max(0.0) - 1.0e-9)
          assert_true(p <= u.min(v) + 1.0e-9)
          let rectangle = c.rectangle(u, u + 0.1, v, v + 0.1)
          assert_true(rectangle >= 0.0)
          mass = mass + rectangle
        }
      }
      close(mass, 1.0, 1.0e-8)
    }
  }
}

///|
test "all invalid input classes propagate typed errors" {
  let bad_uniform = [{ x: 0.0, y: 0.5, }, { x: 0.5, y: 0.5, }]
  let c = independence()
  assert_true(
    try {
      ignore(c.simulate(bad_uniform))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(Copula::new(Frank, 0.0))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(ranks([0.0 / 0.0]))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(c.cdf(1.0 / 0.0, 0.5))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(probability_grid(65))
      false
    } catch {
      CopulaError::BudgetExceeded(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(wilson_interval(2, 1))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(sample_summary([1.0, 2.0], coverage=1.0))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
}

///|
test "resampling rejects invalid plans and fitting exhaustion without dropping observations" {
  let raw = [{ x: 1.0, y: 1.2, }, { x: 2.0, y: 1.8, }, { x: 3.0, y: 2.2, }]
  assert_true(
    try {
      ignore(resample_pairs(raw, [-1, 0, 1]))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(
        bootstrap_raw(Clayton, raw, [[0, 1, 2], [1, 2, 3]], { x: 0.8, y: 0.8, }),
      )
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(
        bootstrap_raw(
          Clayton,
          raw,
          [[0, 1, 2], [1, 2, 0]],
          { x: 0.8, y: 0.8, },
          options={
            lower: 0.01,
            upper: 30.0,
            tolerance: 1.0e-9,
            iterations: 1,
          },
        ),
      )
      false
    } catch {
      CopulaError::NumericalFailure(_) => true
      _ => false
    },
  )
}

///|
test "diagnostics reject constant correlation and unsupported grids" {
  let raw = [{ x: 0.3, y: 0.2, }, { x: 0.3, y: 0.8, }]
  assert_true(
    try {
      ignore(independence().residual_diagnostic(raw))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(independence().grid_diagnostic(raw, [{ x: 1.01, y: 0.5, }]))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
  assert_true(
    try {
      ignore(checkerboard_cdf(raw, 0.5, 0.5))
      false
    } catch {
      CopulaError::InvalidInput(_) => true
      _ => false
    },
  )
}

///|
test "empirical inverse satisfies discrete quantile inequalities with ties" {
  let margin = EmpiricalMargin::new([1.0, 1.0, 3.0, 7.0, 7.0])
  for i in 1..<=100 {
    let p = i.to_double() / 100.0
    let q = margin.quantile(p)
    assert_true(margin.cdf(q) >= p - 1.0e-14)
    assert_true(margin.cdf(q - 0.001) < p + 1.0e-14)
  }
}
