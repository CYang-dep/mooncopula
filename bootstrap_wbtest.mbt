///|
test "percentile summaries match hand computed type7 and sample deviation" {
  let s = sample_summary([1.0, 2.0, 3.0, 4.0], coverage=0.5)
  close(s.mean, 2.5, 1.0e-14)
  close(s.standard_deviation, (5.0 / 3.0).sqrt(), 1.0e-14)
  close(s.lower, 1.75, 1.0e-14)
  close(s.upper, 3.25, 1.0e-14)
}

///|
test "paired plans preserve ownership and rerank reproducibly" {
  let data = [
    { x: 10.0, y: 12.0, },
    { x: 20.0, y: 18.0, },
    { x: 30.0, y: 25.0, },
    { x: 40.0, y: 35.0, },
  ]
  let copied = resample_pairs(data, [3, 0, 3, 1])
  close(copied[0].x, 40.0, 0.0)
  copied[0] = { x: -1.0, y: -1.0, }
  close(data[3].x, 40.0, 0.0)
  let result = bootstrap_raw(Independence, data, [[0, 1, 2, 3], [3, 1, 1, 0]], {
    x: 0.8,
    y: 0.9,
  })
  assert_eq(result.replicates.length(), 2)
  close(result.theta.mean, 0.0, 0.0)
  close(result.upper_joint_probability.mean, 0.02, 1.0e-14)
}
