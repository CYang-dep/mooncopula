///|
/// Synthetic component lifetimes; caller-provided uniforms ensure reproducibility.
fn main raise {
  let raw : Array[@mooncopula.Pair] = [
    { x: 100.0, y: 130.0, },
    { x: 180.0, y: 200.0, },
    { x: 240.0, y: 210.0, },
    { x: 300.0, y: 360.0, },
    { x: 400.0, y: 410.0, },
    { x: 500.0, y: 470.0, },
  ]
  let model = @mooncopula.Copula::new(@mooncopula.Clayton, 2.0)
  let uniforms : Array[@mooncopula.Pair] = [
    { x: 0.1, y: 0.2, },
    { x: 0.4, y: 0.8, },
    { x: 0.7, y: 0.3, },
    { x: 0.9, y: 0.6, },
  ]
  let dependent = model.simulate(uniforms)
  let scenarios = model.scenarios(raw, uniforms)
  let residual = model.rosenblatt(dependent)
  assert_true((residual[0].y - uniforms[0].y).abs() < 1.0e-8)
  println("reliability: dependent lifetime scenarios")
  for i in 0..<scenarios.length() {
    println("scenario[\{i}]=\{scenarios[i].x},\{scenarios[i].y}")
  }
  let conditional = model.inverse_conditional(0.25, 0.9)
  println(
    "conditional_q90=\{conditional.value}; residual=\{conditional.residual}",
  )
}
