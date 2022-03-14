var documenterSearchIndex = {"docs":
[{"location":"api/#Base-Functions","page":"Library","title":"Base Functions","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [QuantumAnnealing]\nPages   = [\"base.jl\"]\nOrder   = [:function]\nPrivate  = true","category":"page"},{"location":"api/#QuantumAnnealing._check_ising_model_ids-Tuple{Dict}","page":"Library","title":"QuantumAnnealing._check_ising_model_ids","text":"checks that an Ising model qubit ids are in the range 1-to-n and returns the value of n.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.binary_to_braket-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.binary_to_braket","text":"converts a binary state vector into a bra-ket notation string note: reverses qubit order for presentation\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.binary_to_int-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.binary_to_int","text":"converts a binary state vector into an integer id following the package conventions valid ints are from 0-to-2^n-1\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.binary_to_spin-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.binary_to_spin","text":"converts a binary state vector (0/1) into an spin state vector (-1/1)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.initial_state_default-Tuple{Any}","page":"Library","title":"QuantumAnnealing.initial_state_default","text":"ground state of sumi A(0) Xi where A(0) > 0 and B(0) = 0\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.int_to_binary-Tuple{Int64}","page":"Library","title":"QuantumAnnealing.int_to_binary","text":"converts a integer id into a binary state vector following the package conventions valid ints are from 0-to-2^n-1 pad should be the total number qubits in the system\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.int_to_spin-Tuple{Int64}","page":"Library","title":"QuantumAnnealing.int_to_spin","text":"converts a integer id into a spin state vector following the package conventions valid ints are from 0-to-2^n-1 pad should be the total number qubits in the system\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.print_z_state_probabilities-Tuple{Matrix}","page":"Library","title":"QuantumAnnealing.print_z_state_probabilities","text":"given a 2^n vector of probably values, prints each value and its associated state vector. limit is used to limit the total number of states that are printed. sort is used to re-order the states by most likely instead of the default which is numerical order from 0-to-(2^n-1)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.spin_to_binary-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.spin_to_binary","text":"converts a spin state vector (-1/1) into an binary state vector (0/1)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.spin_to_braket-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.spin_to_braket","text":"converts a spin state vector (-1/1) into bra-ket notation (↓/↑) note: reverses qubit order for presentation\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.spin_to_int-Tuple{Vector}","page":"Library","title":"QuantumAnnealing.spin_to_int","text":"converts a spin state vector into an integer id following the package conventions valid ints are from 0-to-2^n-1\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.z_measure_probabilities-Tuple{Matrix}","page":"Library","title":"QuantumAnnealing.z_measure_probabilities","text":"given a 2^n-by-2^n density matrix returns the probably of seeing each of the 2^n states on z-basis.\n\n\n\n\n\n","category":"method"},{"location":"api/#Simulate-Functions","page":"Library","title":"Simulate Functions","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [QuantumAnnealing]\nPages   = [\"simulate.jl\"]\nOrder   = [:function]\nPrivate  = true","category":"page"},{"location":"api/#QuantumAnnealing._bernoulli_factorial-Tuple{Any}","page":"Library","title":"QuantumAnnealing._bernoulli_factorial","text":"computes the n-th Bernoulli number divided by n-th factorial number\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._get_quadratic_coefficients-Tuple{Any, Any, Any}","page":"Library","title":"QuantumAnnealing._get_quadratic_coefficients","text":"converts an arbitrary function f(s) into a quadratic form based on interpolation between two extreme points s0 and s1\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._hamiltonian_commutator-Tuple{Vector, Vector}","page":"Library","title":"QuantumAnnealing._hamiltonian_commutator","text":"given two hamiltonians of orders 0-to-(n-1) and 0-to-(m-1), applies the matrixcommutator operation\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._hamiltonian_eval-Tuple{Real, Vector}","page":"Library","title":"QuantumAnnealing._hamiltonian_eval","text":"given a hamiltonian of orders 0-to-(n-1), multiplies all orders by a scalar\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._hamiltonian_integrate-Tuple{Vector}","page":"Library","title":"QuantumAnnealing._hamiltonian_integrate","text":"given a hamiltonian of orders 0-to-(n-1), returns an integrated hamiltonian of 0-to-n\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._hamiltonian_scalar-Tuple{Real, Vector}","page":"Library","title":"QuantumAnnealing._hamiltonian_scalar","text":"given a hamiltonian of orders 0-to-(n-1), multiplies all orders by a scalar\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._hamiltonian_sum-Tuple{Vector}","page":"Library","title":"QuantumAnnealing._hamiltonian_sum","text":"given a list of order-based hamiltonians (which are also vectors), returns the sum of these\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._matrix_commutator-Tuple{Any, Any}","page":"Library","title":"QuantumAnnealing._matrix_commutator","text":"given to matricies, applies the commutator operation\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing._shift_quadratic_coefficients-NTuple{4, Any}","page":"Library","title":"QuantumAnnealing._shift_quadratic_coefficients","text":"shifts a quadatric function by a value of x\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.hamiltonian_transverse_ising-Tuple{Dict, AnnealingSchedule, Real}","page":"Library","title":"QuantumAnnealing.hamiltonian_transverse_ising","text":"Function to build the transverse field Ising model hamiltonian at a given unitless timestep s.\n\nArguments: isingmodel - ising model represented as a dictionary.  The qubits               and couplings are represented as tuples, and the weights               are numbers.               For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2) annealingschedule - The annealing schedule, of the form given by the struct s - the imaginary timestep. This should usually be in the range from 0.0-to-1.0\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.simulate-Tuple{Dict, Real, AnnealingSchedule}","page":"Library","title":"QuantumAnnealing.simulate","text":"Main function for performing quantum annealing simulation via a Magnus Expansion (second order). Noise can be simulated by running multiple times with randomized constant fields.\n\nArguments: isingmodel - ising model represented as a dictionary.  The qubits               and couplings are represented as tuples, and the weights               are numbers.               For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2) annealingschedule - The annealing schedule, of the form given by the struct\n\nParameters: initialstate - Initial state vector. Defaults to uniform superposition state on n qubits constantfieldx - vector of constant biases in the X basis on each qubit. Default is zeros(n) constantfieldz - vector of constant biases in the Z basis on each qubit. Default is zeros(n) The parameters `meantolandmax_tolspecify the desired simulation accuracy. Thesilence` parameter can be used to suppress the progress log.\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.simulate_flexible_order-Tuple{Dict, Real, AnnealingSchedule, Int64, Int64}","page":"Library","title":"QuantumAnnealing.simulate_flexible_order","text":"an any-order magnus expansion solver with a fixed number of time steps\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [QuantumAnnealing]\nPages   = [\"simulate_de.jl\"]\nOrder   = [:function]\nPrivate  = true","category":"page"},{"location":"api/#QuantumAnnealing.simulate_de-NTuple{4, Any}","page":"Library","title":"QuantumAnnealing.simulate_de","text":"simulator that uses the DifferentialEquations.jl differential equation solver to solve the Schrodinger Equation, rather than the Magnus Expansion.  This simulator can run into issues at higher anneal times.  The API is the same as for the simulate function.\n\nArguments: isingmodel - ising model represented as a dictionary.  The qubits               and couplings are represented as tuples, and the weights               are numbers.               For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2) annealingschedule - The annealing schedule, of the form given by the struct reltol - the relative tolerance that will be used in the ODE simulation\n\nParameters: initialstate - Initial state vector. Defaults to uniform superposition state on n qubits constantfieldx - vector of constant biases in the X basis on each qubit. Default is zeros(n) constantfield_z - vector of constant biases in the Z basis on each qubit. Default is zeros(n)\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.simulate_de-Tuple{Dict, Real, AnnealingSchedule}","page":"Library","title":"QuantumAnnealing.simulate_de","text":"A simplified interface to the simulate_de routine that determines suitable tolerance values to ensure a high accuracy simulation. The parameters mean_tol and max_tol specify the desired simulation accuracy. The silence parameter can be used to suppress the progress log.\n\n\n\n\n\n","category":"method"},{"location":"api/#Annealing-Schedules","page":"Library","title":"Annealing Schedules","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"AnnealingSchedule\nAS_LINEAR\nAS_QUADRATIC\nAS_CIRCULAR\nAS_DW_QUADRATIC","category":"page"},{"location":"api/#QuantumAnnealing.AnnealingSchedule","page":"Library","title":"QuantumAnnealing.AnnealingSchedule","text":"A data structure containing two uni-variate functions defined in the domain of 0.0-to-1.0. The A function drives the evolution of the X basis and the B function drives the evolution of the Z basis. init_default is a uni-variate function that given an integer n builds a initial state vector for an n qubit system.\n\n\n\n\n\n","category":"type"},{"location":"api/#QuantumAnnealing.AS_LINEAR","page":"Library","title":"QuantumAnnealing.AS_LINEAR","text":"An AnnealingSchedule implementing a simple linear form\n\n\n\n\n\n","category":"constant"},{"location":"api/#QuantumAnnealing.AS_QUADRATIC","page":"Library","title":"QuantumAnnealing.AS_QUADRATIC","text":"An AnnealingSchedule implementing a simple quadratic form\n\n\n\n\n\n","category":"constant"},{"location":"api/#QuantumAnnealing.AS_CIRCULAR","page":"Library","title":"QuantumAnnealing.AS_CIRCULAR","text":"An AnnealingSchedule implementing uniform circular motion with an analytical solution on a single qubit\n\n\n\n\n\n","category":"constant"},{"location":"api/#QuantumAnnealing.AS_DW_QUADRATIC","page":"Library","title":"QuantumAnnealing.AS_DW_QUADRATIC","text":"An AnnealingSchedule approximating those used in hardware by D-Wave Systems\n\nNOTE: users are strongly encouraged to download annealing schedules for specific D-Wave Systems devices and load them using parse_dwave_annealing_schedule.\n\n\n\n\n\n","category":"constant"},{"location":"api/#Ising-Functions","page":"Library","title":"Ising Functions","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [QuantumAnnealing]\nPages   = [\"ising.jl\"]\nOrder   = [:function]\nPrivate  = true","category":"page"},{"location":"api/#QuantumAnnealing.compute_ising_energy_levels-Tuple{Dict}","page":"Library","title":"QuantumAnnealing.compute_ising_energy_levels","text":"given an Ising model computes a mapping from energy values to collections of state integers\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.compute_ising_state_energies-Tuple{Dict}","page":"Library","title":"QuantumAnnealing.compute_ising_state_energies","text":"given an Ising model computes a mapping from state integers to energy values\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.eval_ising_state_energy-Tuple{Vector, Dict}","page":"Library","title":"QuantumAnnealing.eval_ising_state_energy","text":"given a state vector of spin values and an Ising model computes the energy of that spin configuration\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.print_ising_energy_levels-Tuple{Dict}","page":"Library","title":"QuantumAnnealing.print_ising_energy_levels","text":"given an Ising model prints state strings by ascending energy levels. limit is used the stop the printing after a number of energy levels have been presented (limit <= 0 will print all states).\n\n\n\n\n\n","category":"method"},{"location":"api/#D-Wave-Functions","page":"Library","title":"D-Wave Functions","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [QuantumAnnealing]\nPages   = [\"dwave.jl\"]\nOrder   = [:function]\nPrivate  = true","category":"page"},{"location":"api/#QuantumAnnealing.annealing_protocol_dwave-Tuple{AnnealingSchedule}","page":"Library","title":"QuantumAnnealing.annealing_protocol_dwave","text":"Function to modify an existing annealing schedule to use a customized annealing schedule (asch).  These parameters are the same as those used in a dwisc call or a dwave schedule. Inputs: annealingschedule - annealingschedule\n\nParameters: asch - This is the annealing-schedule parameter.  This is a list of tuples of the form        [(s₀,seffective₀), (s₀,seffective₁), ..., (sₙ,s_effectiveₙ)].\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.initial_state_default_dwave-Tuple{Any}","page":"Library","title":"QuantumAnnealing.initial_state_default_dwave","text":"ground state of sumi A(0)Xi where A(0) < 0 and B(0) = 0\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.read_bqpjson-Tuple{String}","page":"Library","title":"QuantumAnnealing.read_bqpjson","text":"This function reads in a bqpjson file and generates the ising model dictionary  inputs: bqpjson::String - a bqpjson file (v1.0.0) that can be run on D-Wave hardware outputs: n - number of qubits ising_model::Dict{Tuple => Float64} - Dictionary of qubits and couplings to weights mapping:Dict{Int => Int} - mapping of the qubits in the bqpjson file to simulation qubits\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.read_dwave_annealing_schedule-Tuple{Any}","page":"Library","title":"QuantumAnnealing.read_dwave_annealing_schedule","text":"function to take a CSV of DWave annealing schedule values and convert it into an annealing schedule usable by the simulator. valid values for interpolation are :none, :linear, :quadratic\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.simulate_bqpjson-NTuple{4, Any}","page":"Library","title":"QuantumAnnealing.simulate_bqpjson","text":"function that allows for simulation from a bqpjson data file\n\n\n\n\n\n","category":"method"},{"location":"api/#QuantumAnnealing.simulate_bqpjson_noisy-NTuple{4, Any}","page":"Library","title":"QuantumAnnealing.simulate_bqpjson_noisy","text":"function that allows for simulation with x and z noise from a bqpjson data file. The x_bias and z_bias parameters provide vectors of noise realizations.\n\n\n\n\n\n","category":"method"},{"location":"#QuantumAnnealing-Documentation","page":"Home","title":"QuantumAnnealing Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumAnnealing","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing is a Julia package for simulation of quantum annealing protocols. Due to the complexity of modeling quantum systems QuantumAnnealing is not expected to scale to systems with more than 20 qubits. QuantumAnnealing also provides tools for emulating the quantum annealing protocols that are implemented in hardware by D-Wave Systems.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The latest stable release of QuantumAnnealing can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add QuantumAnnealing","category":"page"},{"location":"","page":"Home","title":"Home","text":"For the current development version, \"checkout\" this package with","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add QuantumAnnealing#master","category":"page"},{"location":"","page":"Home","title":"Home","text":"Test that the package works by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"] test QuantumAnnealing","category":"page"},{"location":"#What-is-Quantum-Annealing?","page":"Home","title":"What is Quantum Annealing?","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The objective of QuantumAnnealing is to solve ODEs arising in dynamic quantum systems. Specifically it solves the Schrödinger equation with a time dependent Hamiltonian H(t), acting over a set of n qubits in natural units, as follows,","category":"page"},{"location":"","page":"Home","title":"Home","text":"i fracddtleftPsi(t)rightrangle = H(t)leftPsi(t)rightrangle","category":"page"},{"location":"","page":"Home","title":"Home","text":"For t in 0 T and for the initial condition leftPsi_0rightrangle.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently the implementation focuses on solving Hamiltonians of the so-called Transverse-Field Ising model, which has the following form,","category":"page"},{"location":"","page":"Home","title":"Home","text":"H(t) = A(t) left( sum_i hatsigma^x_i right) + B(t) left( sum_i  h_i hatsigma^z_i + sum_ij J_ij hatsigma^z_i hatsigma^z_j right)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where the hJ parameters are used to encode an Ising model of interest. It is generally assumed that A(0) gg B(0) and A(T) ll B(T), so that there is a smooth transition between the hatsigma^x and hatsigma^z terms, the so-called annealing process. The default initial state of this system is  leftPsi_0rightrangle =bigotimes^nfrac1sqrt2left(leftuparrowrightrangle -leftdownarrowrightrangle right), which corresponds to the state of minimum energy when B(0)=0 and A(0)  0.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The motivation of this model is if the evolution of the H(t) is slow enough (i.e. adiabatic) then the final state of this dynamical system will be the lowest energy states of the input Ising model (i.e., hJ), which can encode a variety of hard computational tasks of practical interest.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The structure of the annealing functions A(t)B(t) can have a dramatic impact on the final states of the dynamical system and hence exploring different versions of these functions is of general interest.","category":"page"},{"location":"#Simulation-of-Quantum-Annealing","page":"Home","title":"Simulation of Quantum Annealing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In the simulate function is used to perform the ODE simulation in QuantumAnnealing.  Its primary arguments are:","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Ising model of interest (i.e. hJ)\nThe annealing time (i.e. T)\nThe annealing schedule (i.e. A(t)B(t))","category":"page"},{"location":"","page":"Home","title":"Home","text":"The output of the simulate function is a density matrix, rho, representing the complete quantum state in the z-basis spanned by leftuparrowrightrangle =left(beginarrayc 1\n0 endarrayright) and leftdownarrowrightrangle =left(beginarrayc 0\n1 endarrayright). The function z_measure_probabilities can be used to project the density matrix into probabilities of state vectors on the z-basis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Users are encouraged to explore their own annealing functions but some canonical ones are provided for convenience: AS_LINEAR, AS_QUADRATIC, AS_CIRCULAR, AS_DW_QUADRATIC.","category":"page"},{"location":"#Emulation-of-D-Wave-Systems-Hardware","page":"Home","title":"Emulation of D-Wave Systems Hardware","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing includes support for emulating the computations conducted by the quantum annealing hardware produced by D-Wave Systems.","category":"page"},{"location":"","page":"Home","title":"Home","text":"info: Info\nQuantumAnnealing's simulate function implements an idealized closed-quantum system, while real-world hardware is exposed to the environment. Discrepancies between QuantumAnnealing and hardware experiments are expected due to the effects of open-quantum systems.","category":"page"},{"location":"#Ising-Model-Specification","page":"Home","title":"Ising Model Specification","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing supports reading Ising models in the bqpjson format. These data files can be generated from a specific D-Wave hardware graph using the DWIG tool and can be executed on hardware using DWISC tool. Tools for performing classical optimization of bqpjson encoded Ising models are available in ising-solvers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing provides the read_bqpjson function for reading bqpjson files into the Ising models that can be executed with simulate.","category":"page"},{"location":"#Annealing-Schedule-Specification","page":"Home","title":"Annealing Schedule Specification","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As part of D-Wave Systems' documentation annealing schedules are provided as xlsx files with four columns s, A(s), B(s), C (normalized). The parse_dwave_annealing_schedule parses this data (in csv format) and converts it into the conventions used by QuantumAnnealing.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Users are strongly encouraged to download annealing schedules for specific D-Wave Systems devices of interest, however a canonical D-Wave annealing schedule, AS_DW_QUADRATIC, is included with QuantumAnnealing for preliminary testing and debugging.","category":"page"},{"location":"","page":"Home","title":"Home","text":"info: Info\nD-Wave System's uses a different Hamiltonian convention than QuantumAnnealing. Specifically D-Wave System assumes A(t) leq 00. In this case the default initial state of this system is leftPsi_0rightrangle =bigotimes^nfrac1sqrt2left(leftuparrowrightrangle + leftdownarrowrightrangle right), which corresponds to the state of minimum energy when B(0)=0 and A(0)  0.  QuantumAnnealing manages this change for the user when working with D-Wave System annealing schedules.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Although the global annealing schedule is fixed in D-Wave hardware an annealing schedule parameter can be used to modify how the schedule is  executed in time by specifying a list of control points. QuantumAnnealing  provides the function dwave_annealing_protocol to apply these changes to an annealing schedule inside of QuantumAnnealing.","category":"page"},{"location":"#D-Wave-Simulation","page":"Home","title":"D-Wave Simulation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Given an Ising model specification in the bqpjson format and a suitable annealing schedule, simulate_bqpjson will perform a complete simulation of a D-Wave hardware execution. This function emulates the workflow of running DWISC on real hardware.","category":"page"},{"location":"#Implementation-Conventions","page":"Home","title":"Implementation Conventions","text":"","category":"section"},{"location":"#Qubits","page":"Home","title":"Qubits","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing assumes that the qubits under consideration are numbered from 1-to-n, which we denote as q_1 q_2 dots q_n. The state of quantum systems is given by a vector of 2^n complex numbers indexed from 0-to-(2^n-1). Each of these complex numbers are mapped to the qubit states using the following mathematical formula,","category":"page"},{"location":"","page":"Home","title":"Home","text":"state_index = sum^n_i = 1 2^i-1 q_i","category":"page"},{"location":"","page":"Home","title":"Home","text":"In this encoding, qubit 1 (i.e. q_1) is the least significant digit and qubit n (i.e. q_n) is the most significant digit. For example, in a 3 qibit system the a state value of 4 is equivalent to the binary array in Julia [0,0,1]. The translation functions int2binary and binary2int performs these conversions following the code's conventions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When printing qubit states as strings using the bra-ket notation (i.e. |xyz⟩) the least significant qubit is presented as the right most value. For example, in a 3 qubit system the state index 4 is presented as |100⟩ in the bra-ket notation. The translation function binary2braket performs this encoding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When translating binary 0/1 qubit values into up/down spin values the following conversion is used sigma_i = (-1)^q_i, specifically 0 rightarrow 1 (up) and 1 rightarrow -1 (down). This conversion is done so that the following property holds, mod(q_1 + q_2 2) = sigma_1 sigma_2. The function binary2spin performs this translation. When presenting up/down qubit states as strings using the bra-ket notation the values are presented as uparrow/downarrow using the helper function spin2braket. For example, in a 3 qubit system the state index 4 is presented as |↓↑↑⟩ in the up/down bra-ket notation.","category":"page"},{"location":"#Annealing-Schedules","page":"Home","title":"Annealing Schedules","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"QuantumAnnealing assumes that annealing functions are univariate and specified over the domain from 0.0-to-1.0.  The default initial state of the system assumes that A(0)  0 and B(0) = 0 and that the user would like to converge a ground state of the input Ising model.","category":"page"},{"location":"#Ising-Models","page":"Home","title":"Ising Models","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The sign conventions of this implementation ensure that:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Hamiltonians with positive fields, i.e. h0, are minimized by leftdownarrowrightrangle states in the adiabatic limit.\nHamiltonians with negative couplings, i.e. J0, are ferromagnetic and are minimized by leftuparrow uparrow rightrangle and leftdownarrow downarrow rightrangle states in the adiabatic limit.","category":"page"},{"location":"#Units","page":"Home","title":"Units","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The underlying mathematics of QuantumAnnealing assumes a natural unit system. However, in practice it is convent to define the annealing schedule in terms of gigahertz (GHz) and in this case after applying a suitable Plank constant the annealing time has units of nanoseconds. The Ising model parameters are not assumed to be given any particular units but it is common to use values of hJ in the range of -1.0-to-1.0.","category":"page"}]
}
