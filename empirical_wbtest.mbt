///|
test "empirical distribution owns input and preserves interval ties" {
  let points = [{ x: 0.2, y: 0.3, }, { x: 0.6, y: 0.7, }]
  let empirical = EmpiricalCopula::new(points)
  points[0] = { x: 0.9, y: 0.9, }
  assert_eq(empirical.cdf(0.4, 0.4), 0.5)
  assert_eq(empirical.rectangle(0.2, 1.0, 0.0, 1.0), 0.5)
  let exported = empirical.observations()
  exported[0] = { x: 0.9, y: 0.9, }
  assert_eq(empirical.cdf(0.4, 0.4), 0.5)
}

///|
test "checkerboard has exact uniform margins unlike step empirical CDF" {
  let raw = [
    { x: 1.0, y: 4.0, },
    { x: 2.0, y: 1.0, },
    { x: 3.0, y: 3.0, },
    { x: 4.0, y: 2.0, },
  ]
  for i in 0..<=20 {
    let u = i.to_double() / 20.0
    close(checkerboard_cdf(raw, u, 1.0), u, 1.0e-12)
    close(checkerboard_cdf(raw, 1.0, u), u, 1.0e-12)
  }
  close(checkerboard_cdf(raw, 0.5, 0.5), 0.25, 1.0e-12)
}
