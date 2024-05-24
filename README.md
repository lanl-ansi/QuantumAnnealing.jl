# QuantumAnnealing.jl

Status:
[![CI](https://github.com/lanl-ansi/QuantumAnnealing.jl/workflows/CI/badge.svg)](https://github.com/lanl-ansi/QuantumAnnealing.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/lanl-ansi/QuantumAnnealing.jl/branch/main/graph/badge.svg?token=0MYSS2hWWH)](https://codecov.io/gh/lanl-ansi/QuantumAnnealing.jl)
[![Documentation](https://github.com/lanl-ansi/QuantumAnnealing.jl/actions/workflows/documentation.yml/badge.svg)](https://lanl-ansi.github.io/QuantumAnnealing.jl/stable/)
</p>


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

Perform a basic simulation with an annealing time of `1.0` and a linear annealing schedule,
```
ρ = simulate(ising_model, 1.0, AS_LINEAR)
print_z_state_probabilities(ρ)
```

Increase the annealing time to approach the adiabatic limit,
```
ρ = simulate(ising_model, 5.0, AS_LINEAR)
print_z_state_probabilities(ρ)
```

Change the annealing schedule and observe different state probabilities,
```
ρ = simulate(ising_model, 5.0, AS_QUADRATIC)
print_z_state_probabilities(ρ)
```

For a more comprehensive introduction see the [package documentation](https://lanl-ansi.github.io/QuantumAnnealing.jl/stable/).

# Citation
If QuantumAnnealing.jl proves useful in your work, please consider citing it using the following BibTeX citation for the pre-print, available [here](https://arxiv.org/abs/2404.14501):

```
@misc{morrell2024quantumannealing,
      title={QuantumAnnealing: A Julia Package for Simulating Dynamics of Transverse Field Ising Models}, 
      author={Zachary Morrell and Marc Vuffray and Sidhant Misra and Carleton Coffrin},
      year={2024},
      eprint={2404.14501},
      archivePrefix={arXiv},
      primaryClass={quant-ph}
}
```

# License
This software is provided under a BSD-ish license with a "modifications must be indicated" clause.  See the `LICENSE.md` file for the full text. This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
