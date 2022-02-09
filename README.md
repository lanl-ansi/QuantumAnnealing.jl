# QuantumAnnealing.jl

Tools for the simulation and execution of quantum annealing algorithms.


## Quick Start

Install the package,
```
] add QuantumAnnealing
```

Load the package and build a two spin ferromagnetic Ising model for simulation,
```
using QuantumAnnealing

ising_model = Dict((1,) => 0.1, (1,2) => -1.0)
```

Perform a basic simulation with an annealing time of `2.0` and the trigonometric annealing schedule,
```
ρ = simulate(ising_model, 2.0, AS_CIRCULAR)
print_z_state_probabilities(ρ)
```

Increase the annealing time to approach the adiabatic limit,
```
ρ = simulate(ising_model, 10.0, AS_CIRCULAR)
print_z_state_probabilities(ρ)
```

Change the annealing schedule and observer different state probabilities,
```
ρ = simulate(ising_model, 10.0, AS_QUADRATIC)
print_z_state_probabilities(ρ)
```

# License
This software is provided under a BSD-ish license with a "modifications must be indicated" clause.  See the `LICENSE.md` file for the full text. This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
