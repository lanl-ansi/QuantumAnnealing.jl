
"""
An AnnealingSchedule implementing uniform circular motion with an analytical
solution on a single qubit
"""
const AS_CIRCULAR = AnnealingSchedule(s -> cos(π/2*s), s -> sin(π/2*s))

"An AnnealingSchedule implementing a simple linear form"
const AS_LINEAR = AnnealingSchedule(s -> 1.0-s, s -> s)

"An AnnealingSchedule implementing a simple quadratic form"
const AS_QUADRATIC = AnnealingSchedule(s -> (1.0-s)^2, s -> s^2)

"""
converts an arbitrary function `f(s)` into a quadratic form based on
interpolation between two extreme points `s0` and `s1`
"""
function _get_quadratic_coefficients(f, s0, s1)
    sm = (s0+s1)/2.0

    b = f.([s0, sm, s1])
    A = [
        1 s0 s0^2;
        1 sm sm^2;
        1 s1 s1^2
    ]

    x = A\b

    return x
end

"shifts a quadatric function by a value of x"
function _shift_quadratic_coefficients(x, c0, c1, c2)
    c0_shift = c0 + c1*x + c2*x^2
    c1_shift = c1 + 2*c2*x
    c2_shift = c2

    return [c0_shift, c1_shift, c2_shift]
end


function _lie_bracket(A, B)
    return -im * (A*B - B*A)
end

function _tensor_sum_single_qubit(mat, n::Int)
    return sum([foldl(kron,[j == i ? mat : _IMAT for i in n:-1:1]) for j in 1:n])
end

function _tensor_sum_single_qubit(mat, n::Int, weights::Vector)
    return sum([foldl(kron,[j == i ? weights[j] * mat : _IMAT for i in n:-1:1]) for j in 1:n])
end

function _sum_X(n::Int)
    return _tensor_sum_single_qubit(_XMAT, n)
end

function _sum_X(n::Int, w::Vector)
    return _tensor_sum_single_qubit(_XMAT, n, w)
end

function _sum_Y(n::Int)
    return _tensor_sum_single_qubit(_YMAT, n)
end

function _sum_Y(n::Int, w::Vector)
    return _tensor_sum_single_qubit(_YMAT, n, w)
end

function _sum_Z(n::Int)
    return _tensor_sum_single_qubit(_ZMAT, n)
end

function _sum_Z(n::Int, w::Vector)
    return _tensor_sum_single_qubit(_ZMAT, n, w)
end

function _kron_Z(n::Int, t, w::Real)
    Ivec = [1;1]
    Zvec = [1;-1]
    matvec = [k in t ? Zvec : Ivec for k in n:-1:1]
    return SparseArrays.spdiagm(complex(w * foldl(kron, matvec)))
end

"""
Function to build the transverse field Ising model hamiltonian at a given
unitless timestep `s`.

Arguments:
ising_model - ising model represented as a dictionary.  The qubits
              and couplings are represented as tuples, and the weights
              are numbers.
              For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2)
annealing_schedule - The annealing schedule, of the form given by the struct
s - the imaginary timestep. This should usually be in the range from 0.0-to-1.0
"""
function hamiltonian_transverse_ising(ising_model::Dict, annealing_schedule::AnnealingSchedule, s::Real)
    n = _check_ising_model_ids(ising_model)

    x_component = _sum_X(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component += _kron_Z(n, tup, w)
    end

    return annealing_schedule.A(s) * x_component + annealing_schedule.B(s) * z_component
end


"a Magnus expansion implementation simulate specialized to two quadratic functions of orders 1 to 4"
function simulate_magnus_optimized(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule, steps::Int, order::Int; initial_state=nothing, state_steps=nothing)
    if !(1 <= order && order <= 4)
        error("simulate_magnus_optimized only supports orders from 1-to-4, given $(order)")
    end
    if steps < 2
        error("at least two simulation steps are required by simulate_magnus_optimized, given $(steps)")
    end
    if steps < annealing_time/1000
        @warn("the number of simulation steps ($(steps)) is small relative to the annealing time ($(round(Int, annealing_time))), this may cause numerical issues.")
    end

    n = _check_ising_model_ids(ising_model)

    if initial_state == nothing
        initial_state = annealing_schedule.init_default(n)
    end

    track_states = !(state_steps == nothing)

    t0 = 0
    s0 = 0

    R0 = initial_state * initial_state'

    x_component = _sum_X(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component = z_component + _kron_Z(n, tup, w)
    end

    H_parts = _H_parts(x_component, z_component, order)

    s_steps = range(0, 1, length=steps)
    R_current = R0
    U = foldl(kron, [_IMAT for i = 1:n])
    U_next = similar(U)
    if track_states
        push!(state_steps, R_current)
    end

    for i in 1:(steps-1)
        s0 = s_steps[i]
        s1 = s_steps[i+1]

        Ω_list = _Ω_list_optimized(annealing_time, s0, s1, annealing_schedule, H_parts, order)

        #for (i,Ωi) in enumerate(Ω_list)
        #   println("Ω_$(i)")
        #   display(Matrix(Ωi))
        #end

        U_next = LinearAlgebra.exp!(Matrix(sum(Ω_list)))
        U = U_next * U

        if track_states
            R_current = U * R0 * U'
            push!(state_steps, R_current)
        end
    end

    return U * R0 * U'
end


function _H_parts(x_component, z_component, order::Int)
    @assert(1 <= order && order <= 4)

    parts = Dict{Any,Any}(
        (1,) => x_component,
        (2,) => z_component
    )

    if order >= 2
        parts[(2,1)] = _lie_bracket(x_component, z_component)
    end
    if order >= 3
        parts[(3,1)] = _lie_bracket(x_component, parts[(2,1)])
        parts[(3,2)] = _lie_bracket(parts[(2,1)], z_component)
    end
    if order >= 4
        parts[(4,1)] = _lie_bracket(x_component, parts[(3,1)])
        parts[(4,2)] = _lie_bracket(x_component, parts[(3,2)])
        parts[(4,3)] = _lie_bracket(z_component, parts[(3,2)])
    end

    return parts
end


function _Ω_list_optimized(annealing_time::Real, s0::Real, s1::Real, annealing_schedule::AnnealingSchedule, H_parts::Dict, order::Int)
    @assert(1 <= order && order <= 4)
    δs = s1 - s0
    δst = annealing_time*δs

    aqc = _get_quadratic_coefficients(annealing_schedule.A, s0, s1)
    bqc = _get_quadratic_coefficients(annealing_schedule.B, s0, s1)

    aqc = _shift_quadratic_coefficients(s0, aqc...)
    bqc = _shift_quadratic_coefficients(s0, bqc...)

    aqc = [aqc[1], δs*aqc[2], δs^2*aqc[3]]
    bqc = [bqc[1], δs*bqc[2], δs^2*bqc[3]]

    Ω1 = -im*δst*(_integral_1(aqc)*H_parts[(1,)] + _integral_1(bqc)*H_parts[(2,)])
    Ω_list = [Ω1]

    if order >= 2
        Ω2 = -im*δst^2/2*(_integral_21(aqc,bqc)*H_parts[(2,1)])
        push!(Ω_list, Ω2)
    end
    if order >= 3
        Ω3 = -im*δst^3/6*(_integral_31(aqc,bqc)*H_parts[(3,1)] + _integral_31(bqc,aqc)*H_parts[(3,2)])
        push!(Ω_list, Ω3)
    end
    if order >= 4
        Ω4 = -im*δst^4/6*(_integral_41(aqc,bqc)*H_parts[(4,1)] + _integral_42(aqc,bqc)*H_parts[(4,2)] + _integral_41(bqc,aqc)*H_parts[(4,3)])
        push!(Ω_list, Ω4)
    end

    return Ω_list
end


function _integral_1(u::Vector)
    return u[1] + u[2]/2 + u[3]/3
end

function _integral_21(u::Vector, v::Vector)
    return u[2]*v[1]/6 + u[3]*v[1]/6 -
        u[1]*v[2]/6 + u[3]*v[2]/30 -
        u[1]*v[3]/6 - u[2]*v[3]/30
end

function _integral_31(u::Vector, v::Vector)
    return u[2]^2*v[1]/40 - u[1]*u[3]*v[1]/60 +
        u[2]*u[3]*v[1]/24 + 5*u[3]^2*v[1]/252 -
        u[1]*u[2]*v[2]/40 - u[1]*u[3]*v[2]/30 +
        u[2]*u[3]*v[2]/840 + u[3]^2*v[2]/360 +
        u[1]^2*v[3]/60 - u[1]*u[2]*v[3]/120 -
        u[2]^2*v[3]/840 - 5*u[1]*u[3]*v[3]/252 -
        u[2]*u[3]*v[3]/360
end

function _integral_41(u::Vector, v::Vector)
    return -u[1]^2*u[2]*v[1]/120 - u[1]*u[2]^2*v[1]/120 -
        u[2]^3*v[1]/840 - u[1]^2*u[3]*v[1]/120 -
        19*u[1]*u[2]*u[3]*v[1]/1260 - u[2]^2*u[3]*v[1]/360 -
        17*u[1]*u[3]^2*v[1]/2520 - u[2]*u[3]^2*v[1]/504 -
        u[3]^3*v[1]/2520 + u[1]^3*v[2]/120 +
        u[1]^2*u[2]*v[2]/120 + u[1]*u[2]^2*v[2]/840 +
        11*u[1]^2*u[3]*v[2]/2520 - u[1]*u[2]*u[3]*v[2]/1260 -
        u[2]^2*u[3]*v[2]/2520 - u[1]*u[3]^2*v[2]/720 -
        u[2]*u[3]^2*v[2]/2016 - u[3]^3*v[2]/7920 +
        u[1]^3*v[3]/120 + 3*u[1]^2*u[2]*v[3]/280 +
        u[1]*u[2]^2*v[3]/280 + u[2]^3*v[3]/2520 +
        17*u[1]^2*u[3]*v[3]/2520 + 17*u[1]*u[2]*u[3]*v[3]/5040 +
        u[2]^2*u[3]*v[3]/2016 + u[1]*u[3]^2*v[3]/2520 +
        u[2]*u[3]^2*v[3]/7920
end

function _integral_42(u::Vector, v::Vector)
    return u[1]*u[2]*v[1]^2/60 + u[2]^2*v[1]^2/120 +
        u[1]*u[3]*v[1]^2/60 + 19*u[2]*u[3]*v[1]^2/1260 +
        17*u[3]^2*v[1]^2/2520 - u[1]^2*v[1]*v[2]/60 +
        u[2]^2*v[1]*v[2]/420 + 2*u[1]*u[3]*v[1]*v[2]/315 +
        2*u[2]*u[3]*v[1]*v[2]/315 + 17*u[3]^2*v[1]*v[2]/5040 -
        u[1]^2*v[2]^2/120 - u[1]*u[2]*v[2]^2/420 +
        u[1]*u[3]*v[2]^2/1260 + u[2]*u[3]*v[2]^2/1260 +
        u[3]^2*v[2]^2/2016 - u[1]^2*v[1]*v[3]/60 -
        2*u[1]*u[2]*v[1]*v[3]/315 - u[2]^2*v[1]*v[3]/1260 +
        u[2]*u[3]*v[1]*v[3]/1680 + u[3]^2*v[1]*v[3]/1260 -
        19*u[1]^2*v[2]*v[3]/1260 - 2*u[1]*u[2]*v[2]*v[3]/315 -
        u[2]^2*v[2]*v[3]/1260 - u[1]*u[3]*v[2]*v[3]/1680 +
        u[3]^2*v[2]*v[3]/3960 - 17*u[1]^2*v[3]^2/2520 -
        17*u[1]*u[2]*v[3]^2/5040 - u[2]^2*v[3]^2/2016 -
        u[1]*u[3]*v[3]^2/1260 - u[2]*u[3]*v[3]^2/3960
end


"computes the n-th Bernoulli number divided by n-th factorial number"
function _bernoulli_factorial(n)
    if n == 1
        return -0.5
    end

    v = Vector{Rational{BigInt}}(undef, n+1)
    for m in 0:n
        v[m + 1] = 1//(m+1)
        for j in m:-1:1
            v[j] = j*(v[j] - v[j+1])
        end
    end

    return Float64(v[1]/factorial(big(n)))
end

"given to matrices, applies the commutator operation"
function _matrix_commutator(a,b)
    return a*b - b*a
end

"given two hamiltonians of orders 0-to-(n-1) and 0-to-(m-1), applies the _matrix_commutator operation"
function _hamiltonian_commutator(h1::Vector, h2::Vector)
    max_index = (length(h1)-1)+(length(h2)-1)+1
    zeros = 0*(h1[1]*h2[1])
    h_prod = [deepcopy(zeros) for i in 1:max_index]

    for (i,m1) in enumerate(h1)
        for (j,m2) in enumerate(h2)
            pow = (i-1)+(j-1)
            h_prod[pow+1] += _matrix_commutator(m1,m2)
        end
    end

    return h_prod
end

"given a hamiltonian of orders 0-to-(n-1), multiplies all orders by a scalar"
function _hamiltonian_eval(x::Real, h_list::Vector)
    return sum(h*x^(i-1) for (i,h) in enumerate(h_list))
end

"given a hamiltonian of orders 0-to-(n-1), multiplies all orders by a scalar"
function _hamiltonian_scalar(x::Real, h_list::Vector)
    return [x*h for h in h_list]
end

"given a hamiltonian of orders 0-to-(n-1), returns an integrated hamiltonian of 0-to-n"
function _hamiltonian_integrate(h::Vector)
    ih = [i==1 ? (0*h[1]) : h[i-1] for i in 1:length(h)+1]

    for (i,h) in enumerate(ih)
        if i != 1
            ih[i] = h./(i-1)
        end
    end

    return ih
end


"given a list of order-based hamiltonians (which are also vectors), returns the sum of these"
function _hamiltonian_sum(h_list::Vector)
    max_index = maximum(length(h) for h in h_list)
    zeros = 0*(h_list[1][1])
    h_sum = [deepcopy(zeros) for i in 1:max_index]

    for h in h_list
        for (i,m) in enumerate(h)
            h_sum[i] += m
        end
    end

    return h_sum
end


# https://iopscience.iop.org/article/10.1088/2399-6528/aab291
function _Ω_list_generic(δs::Real, H::Vector, order::Int)
    @assert(order >= 1)

    Ω_list = [_hamiltonian_integrate(H)]
    S_list = Dict{Tuple{Int64,Int64},Any}()

    for n in 2:order
        S_list[(1,n)] = _hamiltonian_commutator(Ω_list[n-1], H)
        for j in 2:n-1
            S_list[(j,n)] = _hamiltonian_sum([
                _hamiltonian_commutator(Ω_list[m], S_list[(j-1,n-m)])
            for m in 1:n-j])
        end

        Ω_n = _hamiltonian_sum([
            _hamiltonian_scalar(_bernoulli_factorial(j),_hamiltonian_integrate(S_list[(j,n)]))
        for j in 1:n-1])

        push!(Ω_list, Ω_n)
    end

    Ω_list = [_hamiltonian_eval(δs, Ω_n) for Ω_n in Ω_list]

    return Ω_list
end


"""
a generic Magnus expansion solver with a fixed number of time steps
"""
function simulate_magnus_generic(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule, steps::Int, order::Int; initial_state=nothing, constant_field_x=nothing, constant_field_z=nothing, state_steps=nothing)
    @warn("this generic magnus expansion solver is not optimized, runtime overheads for high orders are significant", maxlog=1)
    if steps < 2
        error("at least two steps are required by simulate, given $(steps)")
    end
    if order > 10
        @warn("magnus expansion orders above 10 can produce numerical stability issues, given $(order)", maxlog=1)
    end
    if steps < annealing_time/1000
        @warn("the number of simulation steps ($(steps)) is small relative to the annealing time ($(round(Int, annealing_time))), this may cause numerical issues.")
    end

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

    track_states = !(state_steps == nothing)

    t0 = 0
    s0 = 0

    R0 = initial_state * initial_state'

    x_component = _sum_X(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component = z_component + _kron_Z(n, tup, w)
    end

    s_steps = range(0, 1, length=steps)
    R_current = R0
    U = foldl(kron, [_IMAT for i = 1:n])

    if track_states
        push!(state_steps, R_current)
    end

    # explore use of https://github.com/JuliaSymbolics/Symbolics.jl
    for i in 1:(steps-1)
        s0 = s_steps[i]
        s1 = s_steps[i+1]

        aqc = _get_quadratic_coefficients(annealing_schedule.A, s0, s1)
        bqc = _get_quadratic_coefficients(annealing_schedule.B, s0, s1)

        aqc = _shift_quadratic_coefficients(s0, aqc...)
        bqc = _shift_quadratic_coefficients(s0, bqc...)

        constant_x = _sum_X(n, constant_field_x)
        constant_z = _sum_Z(n, constant_field_z)

        H = -im*annealing_time*[
            aqc[1] * x_component + bqc[1] * z_component + constant_x + constant_z,
            aqc[2] * x_component + bqc[2] * z_component,
            aqc[3] * x_component + bqc[3] * z_component,
        ]

        #Ω_list = _Ω_list_generic(s1 - s0, [Matrix(h) for h in H], order)
        Ω_list = _Ω_list_generic(s1 - s0, H, order)
        #for (i,Ωi) in enumerate(Ω_list)
        #   println("Ω_$(i)")
        #   display(Matrix(Ωi))
        #end

        U_next = LinearAlgebra.exp!(Matrix(sum(Ω_list)))
        U = U_next * U

        if track_states
            R_current = U * R0 * U'
            push!(state_steps, R_current)
        end
    end

    return U * R0 * U'
end



"""
Main function for performing quantum annealing simulation via a Magnus Expansion (second order).
Noise can be simulated by running multiple times with randomized constant fields.

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
The parameters `mean_tol` and `max_tol` specify the desired simulation accuracy.
The `silence` parameter can be used to suppress the progress log.
"""
function simulate(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule; steps=0, order=4, mean_tol=1e-6, max_tol=1e-4, iteration_limit=100, silence=false, state_steps=nothing, kwargs...)
    start_time = time()
    mean_delta = mean_tol + 1.0
    max_delta = max_tol + 1.0

    simulate_method = simulate_magnus_optimized
    if order > 4 || haskey(kwargs, :constant_field_x) || haskey(kwargs, :constant_field_z)
        simulate_method = simulate_magnus_generic
    end

    if steps <= 1
        steps = max(2, round(Int, annealing_time/100))
    end

    if !silence
        println()
        println("\033[1mmethod: $(simulate_method)  order: $(order)  steps: $(steps)\033[0m")
        println("iter |  steps  |    max(Δ)    |    mean(Δ)   |")
    end

    ρ_prev = simulate_method(ising_model, annealing_time, annealing_schedule, steps, order; kwargs...)

    iteration = 1
    while mean_delta >= mean_tol || max_delta >= max_tol
        steps *= 2

        if state_steps != nothing
            empty!(state_steps)
        end

        ρ = simulate_method(ising_model, annealing_time, annealing_schedule, steps, order; state_steps=state_steps, kwargs...)

        ρ_delta = abs.(ρ .- ρ_prev)
        mean_delta = sum(ρ_delta)/length(ρ_delta)
        max_delta = maximum(ρ_delta)

        !silence && Printf.@printf("%4d | %7d | %e | %e |\n", iteration, steps, max_delta, mean_delta)

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
        Printf.@printf("   simulation steps..: %d\n", steps)
        Printf.@printf("   maximum difference: %e <= %e\n", max_delta, max_tol)
        Printf.@printf("   mean difference...: %e <= %e\n", mean_delta, mean_tol)
        Printf.@printf("   runtime (seconds).: %f\n", time()-start_time)
        println("")
    end

    return ρ_prev
end

