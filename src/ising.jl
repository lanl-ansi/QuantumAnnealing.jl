### Helper Functions for Classical Ising Models ###


function eval_ising_state_energy(state::Vector, ising_model::Dict)
    energy = 0.0
    for (ids,v) in ising_model
        val = v
        for qid in ids
            val *= state[qid]
        end
        energy += val
    end
    return energy
end

function compute_state_energies(ising_model::Dict)
    n = _check_ising_model_ids(ising_model)

    state_energy = Dict{Int,Float64}()
    for state_int in 0:(2^n-1)
        spin_vector = binary2spin(int2binary(state_int, pad=n))
        energy = eval_ising_state_energy(spin_vector, ising_model)
        state_energy[state_int] = energy
    end

    return state_energy
end

function print_state_enrgies(ising_model::Dict)
    state_energies = compute_state_energies(ising_model)

    energy_levels = Dict{Float64,Set{Int}}()
    for (state_id, energy) in state_energies
        if !haskey(energy_levels, energy)
            energy_levels[energy] = Set{Int}()
        end
        push!(energy_levels[energy], state_id)
    end

    energies = sort(collect(keys(energy_levels)))

    n = round(Int, log2(length(state_energies)))


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