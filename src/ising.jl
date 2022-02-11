### Helper Functions for Classical Ising Models ###

"""
given a state vector of spin values and an Ising model
computes the energy of that spin configuration
"""
function eval_ising_state_energy(spin_state::Vector, ising_model::Dict)
    energy = 0.0
    for (ids,v) in ising_model
        val = v
        for qid in ids
            val *= spin_state[qid]
        end
        energy += val
    end
    return energy
end

"""
given an Ising model computes a mapping from state integers to energy values
"""
function compute_ising_state_energies(ising_model::Dict)
    n = _check_ising_model_ids(ising_model)

    state_energy = Dict{Int,Float64}()
    for state_int in 0:(2^n-1)
        spin_vector = binary2spin(int2binary(state_int, pad=n))
        energy = eval_ising_state_energy(spin_vector, ising_model)
        state_energy[state_int] = energy
    end

    return state_energy
end

"""
given an Ising model computes a mapping from energy values to collections of state integers
"""
function compute_ising_energy_levels(ising_model::Dict)
    state_energies = compute_ising_state_energies(ising_model)

    energy_levels = Dict{Float64,Set{Int}}()
    for (state_id, energy) in state_energies
        if !haskey(energy_levels, energy)
            energy_levels[energy] = Set{Int}()
        end
        push!(energy_levels[energy], state_id)
    end

    return energy_levels
end

"""
given an Ising model prints state strings by ascending energy levels
"""
function print_ising_energy_levels(ising_model::Dict)
    n = 1

    for (k,v) in ising_model
        for qid in k
            if qid > n
                n = qid
            end
        end
    end

    energy_levels = compute_ising_energy_levels(ising_model)

    energies = sort(collect(keys(energy_levels)))

    for energy in energies
        state_ids = energy_levels[energy]
        println("\033[1menergy: $(energy)\033[0m")
        for state_id in sort(collect(state_ids))
            state = binary2spin(int2binary(state_id, pad=n))
            state_string = spin2braket(state)
            println("   $(state_string)")
        end
    end
end
