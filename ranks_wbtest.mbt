///|
test "midranks ordinal and rejected ties" {
  let values = [4.0, 1.0, 4.0, 2.0]
  let mid = ranks(values)
  close(mid[0], 0.7, 1.0e-15)
  close(mid[2], 0.7, 1.0e-15)
  assert_eq(values[0], 4.0)
  let ord = ranks(values, policy=Ordinal)
  close(ord[0], 0.6, 1.0e-15)
  close(ord[2], 0.8, 1.0e-15)
  let failed = try {
    ignore(ranks(values, policy=RejectTies))
    false
  } catch {
    _ => true
  }
  assert_true(failed)
  let data = [{ x: 3.0, y: 9.0, }, { x: 1.0, y: 6.0, }, { x: 3.0, y: 6.0, }]
  let (tx, ty) = tie_counts(data)
  assert_eq(tx, 1)
  assert_eq(ty, 1)
  let pseudo = pseudo_observations(data)
  validate_uniforms(pseudo, 2)
}
