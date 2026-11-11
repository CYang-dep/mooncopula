///|
test "Wilson rate interval matches published score interval formula" {
  let w = wilson_interval(5, 10)
  close(w.lower, 0.236593090512564, 1.0e-12)
  close(w.upper, 0.763406909487436, 1.0e-12)
  close(wilson_interval(0, 10).lower, 0.0, 1.0e-14)
  close(wilson_interval(10, 10).upper, 1.0, 1.0e-14)
}

///|
test "calibration counts honor exact threshold ties for all event types" {
  let data = [
    { x: 0.5, y: 0.5, },
    { x: 0.8, y: 0.2, },
    { x: 0.7, y: 0.9, },
    { x: 0.1, y: 0.8, },
  ]
  let thresholds = [{ x: 0.5, y: 0.5, }]
  let c = independence()
  for
    (event, count, forecast) in [
      (LowerBoth, 1, 0.25),
      (UpperBoth, 1, 0.25),
      (UpperEither, 3, 0.75),
      (FirstUpperSecondLower, 1, 0.25),
      (FirstLowerSecondUpper, 1, 0.25),
    ] {
    let result = c.calibration_curve(data, thresholds, event)
    assert_eq(result[0].observed.successes, count)
    close(result[0].forecast, forecast, 1.0e-14)
  }
}
