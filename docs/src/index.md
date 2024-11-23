# QuantumAnnealing Documentation

```@meta
CurrentModule = QuantumAnnealing
```

## Overview

QuantumAnnealing is a Julia package for simulation of quantum annealing protocols.
Due to the complexity of modeling quantum systems QuantumAnnealing is not expected to scale to systems with more than 20 qubits.
QuantumAnnealing also provides tools for emulating the quantum annealing protocols that are implemented in hardware by D-Wave Systems.


## Installation

The latest stable release of QuantumAnnealing can be installed using the Julia package manager with

```julia
] add QuantumAnnealing
```

For the current development version, "checkout" this package with

```julia
] add QuantumAnnealing#master
```

Test that the package works by running

```julia
] test QuantumAnnealing
```


## What is Quantum Annealing?

The objective of QuantumAnnealing is to solve ODEs arising in dynamic quantum systems.
Specifically it solves the Schrödinger equation with a time dependent Hamiltonian $H(t)$, acting over a set of $n$ qubits in natural units, as follows,
```math
i \frac{d}{dt}\left|\Psi(t)\right\rangle = H(t)\left|\Psi(t)\right\rangle
```
For $t \in [0, T]$ and for the initial condition $\left|\Psi_{0}\right\rangle$.  Note that we are taking $\hbar = 1$.

Currently the implementation focuses on solving Hamiltonians of the so-called Transverse-Field Ising model, which has the following form,
```math
H(t) = A(t) \left( \sum_i \hat{\sigma}^x_i \right) + B(t) \left( \sum_i  h_i \hat{\sigma}^z_i + \sum_{i,j} J_{ij} \hat{\sigma}^z_i \hat{\sigma}^z_j \right)
```
where the $h,J$ parameters are used to encode an Ising model of interest. It is generally assumed that $A(0) \gg B(0)$ and $A(T) \ll B(T)$, so that there is a smooth transition between the $\hat{\sigma}^x$ and $\hat{\sigma}^z$ terms, the so-called annealing process. The default initial state of this system is  $\left|\Psi_{0}\right\rangle =\bigotimes^{n}\frac{1}{\sqrt{2}}\left(\left|\uparrow\right\rangle -\left|\downarrow\right\rangle \right)$, which corresponds to the state of minimum energy when $B(0)=0$ and $A(0) > 0$.

The motivation of this model is if the evolution of the $H(t)$ is slow enough (i.e. adiabatic) then the final state of this dynamical system will be the lowest energy states of the input Ising model (i.e., $h,J$), which can encode a variety of hard computational tasks of practical interest.

The structure of the annealing functions $A(t),B(t)$ can have a dramatic impact on the final states of the dynamical system and hence exploring different versions of these functions is of general interest.

### Simulation of Quantum Annealing

In the [`simulate`](@ref) function is used to perform the ODE simulation in QuantumAnnealing.  Its primary arguments are:
- The Ising model of interest (i.e. $h,J$)
- The annealing time (i.e. $T$)
- The annealing schedule (i.e. $A(t),B(t)$)
The output of the `simulate` function is a density matrix, $\rho$, representing the complete quantum state in the z-basis spanned by $\left|\uparrow\right\rangle =\left(\begin{array}{c}
1\\
0
\end{array}\right)$ and $\left|\downarrow\right\rangle =\left(\begin{array}{c}
0\\
1
\end{array}\right)$.
The function [`z_measure_probabilities`](@ref) can be used to project the density matrix into probabilities of state vectors on the z-basis.

Users are encouraged to explore their own annealing functions but some canonical ones are provided for convenience: [`AS_LINEAR`](@ref), [`AS_QUADRATIC`](@ref), [`AS_CIRCULAR`](@ref), [`AS_DW_QUADRATIC`](@ref).


## Emulation of D-Wave Systems Hardware

QuantumAnnealing includes support for emulating the computations conducted by
the quantum annealing hardware produced by [D-Wave Systems](https://www.dwavesys.com/).

!!! info
    QuantumAnnealing's [`simulate`](@ref) function implements an idealized
    closed-quantum system, while real-world hardware is exposed to the environment.
    Discrepancies between QuantumAnnealing and hardware experiments are expected
    due to the effects of open-quantum systems.

### Ising Model Specification

QuantumAnnealing supports reading Ising models in the [`bqpjson`](https://github.com/lanl-ansi/bqpjson) format.
These data files can be generated from a specific D-Wave hardware graph using the [`DWIG`](https://github.com/lanl-ansi/dwig) tool and can be executed on hardware using [`DWISC`](https://github.com/lanl-ansi/dwisc) tool.
Tools for performing classical optimization of `bqpjson` encoded Ising models are available in [`ising-solvers`](https://github.com/lanl-ansi/ising-solvers).

QuantumAnnealing provides the [`read_bqpjson`](@ref) function for reading `bqpjson` files into the Ising models that can be executed with [`simulate`](@ref).


### Annealing Schedule Specification

As part of [D-Wave Systems' documentation](https://docs.dwavesys.com/docs/latest/doc_physical_properties.html) annealing schedules are provided as `xlsx` files
with four columns `s`, `A(s)`, `B(s)`, `C (normalized)`.
The [`read_dwave_annealing_schedule`](@ref) parses this data (in `csv` format) and converts it into the conventions used by QuantumAnnealing.

Users are strongly encouraged to download annealing schedules for specific
D-Wave Systems devices of interest, however a canonical D-Wave annealing
schedule, [`AS_DW_QUADRATIC`](@ref), is included with QuantumAnnealing for
preliminary testing and debugging.

!!! info
    D-Wave System's uses a different Hamiltonian convention than
    QuantumAnnealing. Specifically D-Wave System assumes $A(t) \leq 0.0$.
    In this case the default initial state of this system is
    $\left|\Psi_{0}\right\rangle =\bigotimes^{n}\frac{1}{\sqrt{2}}\left(\left|\uparrow\right\rangle + \left|\downarrow\right\rangle \right)$, which corresponds to the state of
    minimum energy when $B(0)=0$ and $A(0) < 0$.  QuantumAnnealing manages this
    change for the user when working with D-Wave System annealing schedules.

Although the _global_ annealing schedule is fixed in D-Wave hardware an
[_annealing schedule parameter_](https://docs.dwavesys.com/docs/latest/c_qpu_annealing.html#varying-the-global-anneal-schedule) can be used to modify how the schedule is 
executed in time by specifying a list of control points. QuantumAnnealing 
provides the function [`annealing_protocol_dwave`](@ref) to apply these changes
to an annealing schedule inside of QuantumAnnealing.

### D-Wave Simulation

Given an Ising model specification in the `bqpjson` format and a suitable annealing schedule, [`simulate_bqpjson`](@ref) will perform a complete simulation of a D-Wave hardware execution.
This function emulates the workflow of running [`DWISC`](https://github.com/lanl-ansi/dwisc) on real hardware.


## Implementation Conventions


### Qubits

QuantumAnnealing assumes that the qubits under consideration are numbered from 1-to-$n$, which we denote as $q_1, q_2, \dots, q_n$.
The state of quantum systems is given by a vector of $2^n$ complex numbers indexed from 0-to-$(2^n-1)$.
Each of these complex numbers are mapped to the qubit states using the following mathematical formula,
```math
state\_index = \sum^n_{i = 1} 2^{i-1} q_i
```
In this encoding, qubit 1 (i.e. $q_1$) is the least significant digit and qubit $n$ (i.e. $q_n$) is the most significant digit.
For example, in a 3 qibit system the a state value of 4 is equivalent to the binary array in Julia `[0,0,1]`.
The translation functions [`int_to_binary`](@ref) and [`binary_to_int`](@ref) performs these conversions following the code's conventions.

When printing qubit states as strings using the bra-ket notation (i.e. `|xyz⟩`) the least significant qubit is presented as the right most value.
For example, in a 3 qubit system the state index 4 is presented as `|100⟩` in the bra-ket notation.
The translation function [`binary_to_braket`](@ref) performs this encoding.

When translating binary 0/1 qubit values into up/down spin values the following conversion is used $\sigma_i = (-1)^{q_i}$, specifically $0 \rightarrow 1$ (up) and $1 \rightarrow -1$ (down).
This conversion is done so that the following property holds, $mod(q_1 + q_2, 2) = \sigma_1 \sigma_2$.
The function [`binary_to_spin`](@ref) performs this translation.
When presenting up/down qubit states as strings using the bra-ket notation the values are presented as $\uparrow$/$\downarrow$ using the helper function [`spin_to_braket`](@ref).
For example, in a 3 qubit system the state index 4 is presented as `|↓↑↑⟩` in the up/down bra-ket notation.


### Annealing Schedules

QuantumAnnealing assumes that annealing functions are univariate and specified over the domain from 0.0-to-1.0.  The default initial state of the system assumes that $A(0) > 0$ and $B(0) = 0$ and that the user would like to converge a ground state of the input Ising model.


### Ising Models

The sign conventions of this implementation ensure that:
- Hamiltonians with positive fields, i.e. $h>0$, are minimized by $\left|\downarrow\right\rangle$ states in the adiabatic limit.
- Hamiltonians with negative couplings, i.e. $J<0$, are ferromagnetic and are minimized by $\left|\uparrow \uparrow \right\rangle$ and $\left|\downarrow \downarrow \right\rangle$ states in the adiabatic limit.


### Units

The underlying mathematics of QuantumAnnealing assumes a natural unit system. However, in practice it is convent to define the annealing schedule in terms of gigahertz (GHz) and in this case after applying a suitable Plank constant the annealing time has units of nanoseconds.
The Ising model parameters are not assumed to be given any particular units but it is common to use values of $h,J$ in the range of -1.0-to-1.0.

