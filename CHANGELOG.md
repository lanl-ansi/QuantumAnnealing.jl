QuantumAnnealing.jl Change Log
==============================

### Staged
- Remove DifferentialEquations.jl as a dependency in favor of OrdiaryDiffEq.jl
- Specify algorithm for `simulate_de` to avoid crashing issues
- Update tests to remove lower bounds for `simulate_de` functions

### v0.2.1
- Reduce memory usage in `simulate` and `simulate_de`
- Add Gibbs distribution (#34)
- Fix kwargs bug in `simulate_de` (#39)
- Fix analytic two spin model for non-integer `s` values (#40)

### v0.2.0
- Add a generic Magnus expansion solver for any order
- Update hard-coded Magnus expansion solver to support orders 1 through 4
- Improve default initial number of steps in `simulate`
- Change default solver order from 2 to 4 (breaking)
- Update d-wave simulation tools to use adaptive solvers (breaking)
- Reversed coefficient ordering of `get_quadratic_coefficients` (breaking)
- Remove export of a variety of internal helper functions (breaking)

### v0.1.0
- Add variant of `solve_de` with adaptive solve tolerance
- Add tools for working with classical Ising models (#8)
- Add data processing tools (#6)
- Fix units in dwave annealing schedules (#17)

### v0.0.1
- Initial release
