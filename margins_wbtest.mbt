///|
test "empirical margins preserve discrete quantile definitions and ownership" {
  let values = [30.0, 10.0, 20.0, 20.0]
  let margin = EmpiricalMargin::new(values)
  values[0] = 999.0
  assert_eq(margin.cdf(20.0), 0.75)
  assert_eq(margin.quantile(0.75), 20.0)
  assert_eq(margin.quantile(1.0), 30.0)
  assert_eq(margin.quantile(0.0), 10.0)
  assert_eq(margin.linear_quantile(0.5), 20.0)
  let points = [
    { x: 10.0, y: 100.0, },
    { x: 20.0, y: 200.0, },
    { x: 30.0, y: 300.0, },
  ]
  let scenarios = independence().scenarios(points, [
    { x: 0.2, y: 0.8, },
    { x: 0.8, y: 0.2, },
  ])
  assert_eq(scenarios[0].x, 10.0)
  assert_eq(scenarios[0].y, 300.0)
  assert_eq(scenarios[1].x, 30.0)
}
