///|
/// Synthetic paired rainfall and river-flow observations, not a real forecast.
fn main raise {
  let raw : Array[@mooncopula.Pair] = [
    { x: 18.0, y: 30.0, },
    { x: 26.0, y: 35.0, },
    { x: 12.0, y: 20.0, },
    { x: 40.0, y: 60.0, },
    { x: 32.0, y: 38.0, },
    { x: 21.0, y: 32.0, },
    { x: 45.0, y: 55.0, },
    { x: 37.0, y: 51.0, },
    { x: 16.0, y: 25.0, },
    { x: 28.0, y: 45.0, },
    { x: 48.0, y: 66.0, },
    { x: 35.0, y: 49.0, },
  ]
  let uniforms = @mooncopula.pseudo_observations(raw)
  let result = @mooncopula.fit(@mooncopula.Gumbel, uniforms)
  assert_true(result.converged)
  let risk = result.model.risk(0.8, 0.8)
  let empirical = @mooncopula.EmpiricalCopula::new(uniforms)
  let diagnostic = result.model.grid_diagnostic(
    uniforms,
    @mooncopula.probability_grid(10),
  )
  println("environment: Gumbel dependence fitted")
  println("theta=\{result.model.theta}")
  println("joint_exceedance=\{risk.upper_both}")
  println("empirical_exceedance=\{empirical.rectangle(0.8, 1.0, 0.8, 1.0)}")
  println("grid_rmse=\{diagnostic.root_mean_square}")
  match @mooncopula.return_period(risk.upper_both, 1.0) {
    Some(p) => println("return_period=\{p}")
    None => println("return_period=infinite")
  }
}
