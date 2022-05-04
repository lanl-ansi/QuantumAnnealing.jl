### Foundations of the QuantumAnnealing Package ###

"""
A data structure containing two uni-variate functions defined in the domain of 0.0-to-1.0.
The `A` function drives the evolution of the X basis and
the `B` function drives the evolution of the Z basis.
`init_default` is a uni-variate function that given an integer n builds a
initial state vector for an n qubit system.
"""
struct AnnealingSchedule
    A::Function
    B::Function
    init_default::Function
end

"""
A short hand AnnealingSchedule constructor that uses the initial_state_default,
which is the most common case for the conventions of this implementation.
"""
AnnealingSchedule(A,B) = AnnealingSchedule(A, B, initial_state_default)

#predefining Pauli Matrices
const _IMAT = SparseArrays.sparse([1,2], [1,2], [1.0+0im;1.0+0im])
const _XMAT = SparseArrays.sparse([1,2], [2,1], [1.0+0im;1.0+0im])
const _YMAT = SparseArrays.sparse([1,2], [2,1], [-im;im])
const _ZMAT = SparseArrays.sparse([1,2], [1,2], [1.0+0im;-1.0+0im])


"""
ground state of sum_i A(0) X_i where A(0) > 0 and B(0) = 0
"""
function initial_state_default(n)
    return complex(foldl(kron,[[1;-1] for i in 1:n]) ./ (2^(n/2)))
end


"""
checks that an Ising model qubit ids are in the range 1-to-n and returns the
value of n.
"""
function _check_ising_model_ids(ising_model::Dict)
    qubit_ids = Set{Int}()
    for (k,v) in ising_model
        for qid in k
            push!(qubit_ids, qid)
        end
    end

    qid_min = minimum(qubit_ids)
    qid_max = maximum(qubit_ids)

    if qid_min < 1
        error("ising model qubit id of $(qid_min) found but only values greater than or equal to 1 are supported.")
    end

    for i in 1:qid_max
        if !(i in qubit_ids)
            @warn "qubit id of $(i) is not used, renumber qubit ids from 1-to-n to increase performance"
        end
    end

    return qid_max
end

"""
converts a integer id into a binary state vector following the package conventions
valid ints are from 0-to-2^n-1
pad should be the total number qubits in the system
"""
function int_to_binary(x::Int; pad=0)
    return digits(x, base=2, pad=pad)
end

"""
converts a binary state vector into an integer id following the package conventions
valid ints are from 0-to-2^n-1
"""
function binary_to_int(states::Vector)
    return sum(v * (1<<(i-1)) for (i,v) in enumerate(states))
end

"""
converts a binary state vector (0/1) into an spin state vector (-1/1)
"""
function binary_to_spin(states::Vector)
    return [v == 0 ? 1 : -1 for v in states]
end

"""
converts a spin state vector into an integer id following the package conventions
valid ints are from 0-to-2^n-1
"""
function spin_to_int(spin::Vector)
    return binary_to_int(spin_to_binary(spin))
end

"""
converts a integer id into a spin state vector following the package conventions
valid ints are from 0-to-2^n-1
pad should be the total number qubits in the system
"""
function int_to_spin(x::Int; pad=0)
    return binary_to_spin(int_to_binary(x, pad=pad))
end

"""
converts a spin state vector (-1/1) into an binary state vector (0/1)
"""
function spin_to_binary(spin::Vector)
    return [i == 1 ? 0 : 1 for i in spin]
end


"""
converts a binary state vector into a bra-ket notation string
note: reverses qubit order for presentation
"""
function binary_to_braket(states::Vector)
    return "|$(join(reverse(states)))⟩"
end

"""
converts a spin state vector (-1/1) into bra-ket notation (↓/↑)
note: reverses qubit order for presentation
"""
function spin_to_braket(states::Vector)
    return "|$(join([s < 0 ? "↓" : "↑" for s in reverse(states)]))⟩"
end


"""
given a 2^n-by-2^n density matrix returns the probably of seeing each of the
2^n states on z-basis.
"""
function z_measure_probabilities(density::Matrix)
    return [real(density[i,i]) for i in 1:(size(density)[1])]
end


"""
given a 2^n vector of probably values, prints each value and its associated
state vector. `limit` is used to limit the total number of states that are
printed. `sort` is used to re-order the states by most likely instead of the
default which is numerical order from 0-to-(2^n-1)
"""
function print_z_state_probabilities(density::Matrix; limit=50, sort=false)
    probs = z_measure_probabilities(density)
    n = ceil(Int, log2(length(probs)))

    prob_order = enumerate(probs)
    if sort
        prob_order = Base.sort(collect(prob_order), by=(x) -> x[2], rev=true)
    end

    i = 0
    for (state_id,pr) in prob_order
        state = int_to_spin(state_id-1, pad=n)
        state_string = spin_to_braket(state)
        prob_string = rpad(round(pr, digits=6),8, " ")
        println("$(prob_string) $(state_string)")
        i += 1
        if limit > 0 && i >= limit
            println("first $(limit) of $(length(probs)) states shown")
            break
        end
    end
end

function distribution_gibbs(H::Function, s::Real, β::Real)
    tr(A) = sum([A[i,i] for i=1:(size(A)[1])])
    exp_ham = exp(-β * Matrix(H(s)))
    exp_ham = exp_ham ./ tr(exp_ham)
    return exp_ham
end
