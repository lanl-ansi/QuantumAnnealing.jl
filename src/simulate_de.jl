
"""
simulator that uses the DifferentialEquations.jl differential equation solver
to solve the Schrodinger Equation, rather than the Magnus Expansion.  This
simulator can run into issues at higher anneal times.  The API is the same as
for the simulate function.

Arguments:
ising_model - ising model represented as a dictionary.  The qubits
              and couplings are represented as tuples, and the weights
              are numbers.
              For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2)
annealing_schedule - The annealing schedule, of the form given by the struct

Parameters:
initial_state - Initial state vector. Defaults to uniform superposition state on n qubits
constant_field_x - vector of constant biases in the X basis on each qubit. Default is zeros(n)
constant_field_z - vector of constant biases in the Z basis on each qubit. Default is zeros(n)
"""
function simulate_de(ising_model, annealing_time, annealing_schedule; initial_state=nothing, constant_field_x=nothing, constant_field_z=nothing)
    n = _check_ising_model_ids(ising_model)

    if initial_state == nothing
        initial_state = annealing_schedule.init_default(n)
    end

    if constant_field_x == nothing
        constant_field_x = zeros(n)
    end

    if constant_field_z == nothing
        constant_field_z = zeros(n)
    end

    const_x_component = sum_x(n, constant_field_x)
    const_z_component = sum_z(n, constant_field_z)

    H(s) = transverse_ising_hamiltonian(ising_model, annealing_schedule, s) + const_x_component + const_z_component
    schrod_eq(state, time, s) = -im * time * H(s) * state

    s_range = (0.0, 1.0)
    prob = DifferentialEquations.ODEProblem(schrod_eq, initial_state, s_range, annealing_time)
    sol = DifferentialEquations.solve(prob)

    state = sol(1)
    density = state * state'

    return density
end

