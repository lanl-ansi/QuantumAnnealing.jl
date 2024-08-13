
"""
simulator that uses the OrdinaryDiffEq.jl differential equation solver
to solve the Schrodinger Equation, rather than the Magnus Expansion.  This
simulator can run into issues at higher anneal times.  The API is the same as
for the simulate function.

Arguments:
ising_model - ising model represented as a dictionary.  The qubits
              and couplings are represented as tuples, and the weights
              are numbers.
              For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2)
annealing_schedule - The annealing schedule, of the form given by the struct
reltol - the relative tolerance that will be used in the ODE simulation

Parameters:
initial_state - Initial state vector. Defaults to uniform superposition state on n qubits
constant_field_x - vector of constant biases in the X basis on each qubit. Default is zeros(n)
constant_field_z - vector of constant biases in the Z basis on each qubit. Default is zeros(n)
"""
function simulate_de(ising_model, annealing_time, annealing_schedule, reltol; abstol=reltol/1e3, initial_state=nothing, constant_field_x=nothing, constant_field_z=nothing, kwargs...)
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

    const_x_component = _sum_X(n, constant_field_x)
    const_z_component = _sum_Z(n, constant_field_z)

    H(s) = hamiltonian_transverse_ising(ising_model, annealing_schedule, s) + const_x_component + const_z_component
    schrod_eq(state, time, s) = -im * time * H(s) * state

    s_range = (0.0, 1.0)
    prob = OrdinaryDiffEq.ODEProblem(schrod_eq, initial_state, s_range, annealing_time)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(); abstol=abstol, reltol=reltol, saveat=[1.0], kwargs...)


    state = sol(1)
    density = state * state'

    return density
end



"""
A simplified interface to the `simulate_de` routine that determines
suitable tolerance values to ensure a high accuracy simulation.
The parameters `mean_tol` and `max_tol` specify the desired simulation accuracy.
The `silence` parameter can be used to suppress the progress log.
"""
function simulate_de(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule; reltol=1e-4, mean_tol=1e-6, max_tol=1e-4, iteration_limit=100, silence=false, kwargs...)
    start_time = time()
    mean_delta = mean_tol + 1.0
    max_delta = max_tol + 1.0

    if !silence
        println()
        println("iter |    reltol    |    max(Δ)    |    mean(Δ)   |")
    end

    ρ_prev = simulate_de(ising_model, annealing_time, annealing_schedule, reltol; kwargs...)

    iteration = 1
    while mean_delta >= mean_tol || max_delta >= max_tol
        reltol /= 10.0

        ρ = simulate_de(ising_model, annealing_time, annealing_schedule, reltol;  kwargs...)

        ρ_delta = abs.(ρ .- ρ_prev)
        mean_delta = sum(ρ_delta)/length(ρ_delta)
        max_delta = maximum(ρ_delta)

        !silence && Printf.@printf("%4d | %e | %e | %e |\n", iteration, reltol, max_delta, mean_delta)

        ρ_prev = ρ
        iteration += 1
        if iteration > iteration_limit
            error("iteration limit reached in simulate function without reaching convergence criteria")
        end
    end

    if !silence
        println("")
        println("\033[1mconverged\033[0m")
        Printf.@printf("   iterations........: %d\n", iteration-1)
        Printf.@printf("   simulation reltol.: %e\n", reltol)
        Printf.@printf("   maximum difference: %e <= %e\n", max_delta, max_tol)
        Printf.@printf("   mean difference...: %e <= %e\n", mean_delta, mean_tol)
        Printf.@printf("   runtime (seconds).: %f\n", time()-start_time)
        println("")
    end

    return ρ_prev
end


