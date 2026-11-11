///|
/// Midrank policy keeps ties equal; Ordinal is stable input-order tie breaking.
pub(all) enum TiePolicy {
  Midrank
  Ordinal
  RejectTies
} derive(Debug)

///|
pub extend TiePolicy with @debug.Debug::{to_repr}

///|
/// One-based ranks divided by n+1, ensuring every output is interior.
pub fn ranks(
  values : Array[Double],
  policy~ : TiePolicy = Midrank,
) -> Array[Double] raise CopulaError {
  if values.length() < 1 {
    raise InvalidInput("empty rank input")
  }
  if values.length() > 4096 {
    raise BudgetExceeded("rank input exceeds 4096")
  }
  for value in values {
    finite(value, "rank input")
  }
  let order = Array::makei(values.length(), fn(i) { i })
  order.sort_by(fn(a, b) {
    let cmp = values[a].compare(values[b])
    if cmp == 0 {
      a.compare(b)
    } else {
      cmp
    }
  })
  let result = Array::make(values.length(), 0.0)
  let denominator = (values.length() + 1).to_double()
  let mut start = 0
  while start < order.length() {
    let mut end = start + 1
    while end < order.length() && values[order[end]] == values[order[start]] {
      end = end + 1
    }
    match policy {
      RejectTies => {
        if end - start > 1 {
          raise InvalidInput("ties forbidden by policy")
        }
        result[order[start]] = (start + 1).to_double() / denominator
      }
      Midrank => {
        let rank = (start + 1 + end).to_double() / 2.0 / denominator
        for k in start..<end {
          result[order[k]] = rank
        }
      }
      Ordinal =>
        for k in start..<end {
          result[order[k]] = (k + 1).to_double() / denominator
        }
    }
    start = end
  }
  result
}

///|
/// Convert finite paired raw observations into paired open-square ranks.
pub fn pseudo_observations(
  data : Array[Pair],
  policy~ : TiePolicy = Midrank,
) -> Array[Pair] raise CopulaError {
  validate_pairs(data, 2)
  let xs = data.map(fn(p) { p.x })
  let ys = data.map(fn(p) { p.y })
  let u = ranks(xs, policy~)
  let v = ranks(ys, policy~)
  Array::makei(data.length(), fn(i) { { x: u[i], y: v[i], } })
}

///|
/// Count excess observations tied with an earlier occurrence, per coordinate.
pub fn tie_counts(data : Array[Pair]) -> (Int32, Int32) raise CopulaError {
  validate_pairs(data, 1)
  let mut tx = 0
  let mut ty = 0
  for i in 0..<data.length() {
    let mut x_seen = false
    let mut y_seen = false
    for j in 0..<i {
      if data[i].x == data[j].x {
        x_seen = true
      }
      if data[i].y == data[j].y {
        y_seen = true
      }
    }
    if x_seen {
      tx = tx + 1
    }
    if y_seen {
      ty = ty + 1
    }
  }
  (tx, ty)
}
