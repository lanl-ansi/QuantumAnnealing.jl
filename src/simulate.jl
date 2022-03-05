

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
function get_function_coefficients(f, s0, s1)
    smid = (s0+s1)/2.0

    b = f.([s0, smid, s1])
    A = [(s0)^2 (s0) 1
         (smid)^2 (smid) 1
         (s1)^2 (s1) 1]

    x = A\b

    return x
end

function lie_bracket(A, B)
    return -im * (A*B - B*A)
end

function tensor_sum_single_qubit(mat, n::Int)
    return sum([foldl(kron,[j == i ? mat : IMAT for i in n:-1:1]) for j in 1:n])
end

function tensor_sum_single_qubit(mat, n::Int, weights::Vector)
    return sum([foldl(kron,[j == i ? weights[j] * mat : IMAT for i in n:-1:1]) for j in 1:n])
end

function sum_x(n::Int)
    return tensor_sum_single_qubit(XMAT, n)
end

function sum_x(n::Int, w::Vector)
    return tensor_sum_single_qubit(XMAT, n, w)
end

function sum_y(n::Int)
    return tensor_sum_single_qubit(YMAT, n)
end

function sum_y(n::Int, w::Vector)
    return tensor_sum_single_qubit(YMAT, n, w)
end

function sum_z(n::Int)
    return tensor_sum_single_qubit(ZMAT, n)
end

function sum_z(n::Int, w::Vector)
    return tensor_sum_single_qubit(ZMAT, n, w)
end

function zizj_vectorized(n::Int, i::Int, j::Int, J_val)
    Ivec = [1.0+0im,1.0+0im]
    Zvec = [1.0+0im,-1.0+0im]
    matvec = [(k == i || k == j) ? Zvec : Ivec for k in n:-1:1]
    return J_val * foldl(kron, matvec)
end

function sum_zizj(n, J::Dict)
    sum_zs = zeros(2^n)
    for (i,j) in keys(J)
        J_val = J[(i,j)]
        sum_zs +=  zizj_vectorized(n, i, j, J_val)
    end
    return SparseArrays.spdiagm(sum_zs)
end

function sum_z_tup(n, tup, w)
    Ivec = [1;1]
    Zvec = [1;-1]
    matvec = [k in tup ? Zvec : Ivec for k in n:-1:1]
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
function transverse_ising_hamiltonian(ising_model::Dict, annealing_schedule::AnnealingSchedule, s::Real)
    n = _check_ising_model_ids(ising_model)

    x_component = sum_x(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component = z_component + sum_z_tup(n, tup, w)
    end

    return annealing_schedule.A(s) * x_component + annealing_schedule.B(s) * z_component
end


function integral_1_sched(a_2, a_1, a_0, s0, δ)
    #TODO: Possibly change to integrate from s0 to s1 to allow
    #      for a more fine-grained integration with As and Bs
    #      from spline
    #δ = s1 - s0
    left = [(2*(3*s0^2+3*s0*δ+δ^2)) (3*(2*s0+δ)) 6] * δ / 6
    right = [a_2, a_1, a_0]
    return (left * right)[1]
end

function integral_2_sched(a_2, a_1, a_0, b_2, b_1, b_0, s0, δ)
    #TODO: see integral_1_sched todo.
    #returns the magnus expansions double integral for 
    #(A(s1)B(s2) - A(s2)B(s1))
    #δ = s0 - s1
    left = [a_2 a_1 a_0]
    right = [b_2,b_1,b_0]
    integ_mat = [(0) (s0^2 + s0*δ + δ^2/5) (2*s0 + δ);
                 (-(s0^2 + s0*δ + δ^2/5)) (0) (1);
                 -(2*s0+δ) (-1) (0)] * δ^3 / 6
    return (left * integ_mat * right)[1]
end




"""
Main function for performing quantum annealing simulation via a Magnus Expansion.
Noise can be simulated by running multiple times with randomized constant fields.

Arguments:
ising_model - ising model represented as a dictionary.  The qubits
              and couplings are represented as tuples, and the weights
              are numbers.
              For Example: im = Dict((1,) => 1, (2,) => 0.5, (1,2) => 2)
annealing_schedule - The annealing schedule, of the form given by the struct
steps - number of iterations for the Magnus Expansion

Parameters:
initial_state - Initial state vector. Defaults to uniform superposition state on n qubits
constant_field_x - vector of constant biases in the X basis on each qubit. Default is zeros(n)
constant_field_z - vector of constant biases in the Z basis on each qubit. Default is zeros(n)
"""
function simulate(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule, steps::Int; initial_state=nothing, constant_field_x=nothing, constant_field_z=nothing, state_steps=nothing)
    if steps < 2
        error("at least two steps are required by simulate, given $(steps)")
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

    ηs = ones(n)
    hs = zeros(n)

    x_component = sum_x(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component = z_component + sum_z_tup(n, tup, w)
    end

    xz_bracket = lie_bracket(x_component, z_component)

    constant_component = sum_x(n, constant_field_x) + sum_z(n, constant_field_z)
    constant_bracket_x = lie_bracket(x_component, constant_component)
    constant_bracket_z = lie_bracket(z_component, constant_component)

    s_steps = range(0, 1, length=steps)
    R_current = R0
    U = foldl(kron, [IMAT for i = 1:n])

    if track_states
        push!(state_steps, R_current)
    end

    for i in 1:(steps-1)
        s0 = s_steps[i]
        s1 = s_steps[i+1]
        δs = s1 - s0

        a_2, a_1, a_0 = get_function_coefficients(annealing_schedule.A, s0, s1)
        b_2, b_1, b_0 = get_function_coefficients(annealing_schedule.B, s0, s1)

        integA = integral_1_sched(a_2, a_1, a_0, s0, δs)
        integB = integral_1_sched(b_2, b_1, b_0, s0, δs)
        #println(integA, " ", integB)
        integ2 = integral_2_sched(a_2, a_1, a_0, b_2, b_1, b_0, s0, δs)

        integ2A = a_1*δs^3/6 + a_2*s0*δs^3/3 + a_2*δs^4/6
        integ2B = b_1*δs^3/6 + b_2*s0*δs^3/3 + b_2*δs^4/6

        Ω1Sched = integA * x_component + integB * z_component
        Ω1Const = δs * constant_component
        Ω1 = Ω1Sched + Ω1Const

        #Ω2Sched = integ2 * xz_bracket 
        #Ω2Const = integ2A * constant_bracket_x + integ2B * constant_bracket_z
        #Ω2 = (Ω2Sched + Ω2Const)/2

        #println(Matrix(-im .* annealing_time.*Ω1))
        #U_next = exp(Matrix(-im * (annealing_time*Ω1 + (annealing_time^2)*Ω2)))
        U_next = exp(Matrix(-im * (annealing_time*Ω1)))
        U = U_next * U

        if track_states
            R_current = U * R0 * U'
            push!(state_steps, R_current)
        end
    end

    return U * R0 * U'
end



"""
A simplified interface to the core `simulate` routine that determines a
suitable number of steps to ensure a high accuracy simulation.
The parameters `mean_tol` and `max_tol` specify the desired simulation accuracy.
The `silence` parameter can be used to suppress the progress log.
"""
function simulate(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule; steps=2, mean_tol=1e-6, max_tol=1e-4, iteration_limit=100, silence=false, state_steps=nothing, kwargs...)
    start_time = time()
    mean_delta = mean_tol + 1.0
    max_delta = max_tol + 1.0

    if !silence
        println()
        println("iter |  steps  |    max(Δ)    |    mean(Δ)   |")
    end

    ρ_prev = simulate(ising_model, annealing_time, annealing_schedule, steps; kwargs...)

    iteration = 1
    while mean_delta >= mean_tol || max_delta >= max_tol
        steps *= 2

        if state_steps != nothing
            empty!(state_steps)
        end

        ρ = simulate(ising_model, annealing_time, annealing_schedule, steps; state_steps=state_steps, kwargs...)

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






# multiplication of polynomiuals
# integration of polynomiuals
# commuator of hamilatoan
# hamiltonian DS

function _poly_product(p1::Vector{<:Real}, p2::Vector{<:Real})
    p_prod = Dict{Int,Float64}()
    for (i,c1) in enumerate(p1)
        for (j,c2) in enumerate(p2)
            pow = (i-1)+(j-1)
            if !haskey(p_prod, pow+1)
                p_prod[pow+1] = 0.0
            end
            p_prod[pow+1] += c1 * c2
        end
    end
    max_index = maximum(keys(p_prod))
    p = [0.0 for i in 1:max_index]
    for (k,v) in p_prod
        p[k] = v
    end
    return p
end
#_poly_product([1,1], [1,1])


function _poly_integrate(p::Vector{<:Real})
    ip = [0.0 for i in 1:length(p)+1]

    for (i,c) in enumerate(p)
        ip[i+1] = c/(i)
    end

    return ip
end
#_poly_integrate([1,1])


function _poly_eval(x::Real, p::Vector{<:Real})
    return sum(v*x^(i-1) for (i,v) in enumerate(p))
end
# _poly_eval(2, [1,1,1])
# _poly_eval(2, _poly_integrate([1,1,1]))


function _matrix_commutator(a::Matrix, b::Matrix)
    # if cancels return zero (TBD)
    return a*b - b*a
end

function _hamiltonian_commutator(h1::Vector, h2::Vector)
    h = []
    for (p1,m1) in h1
        for (p2,m2) in h2
            pp = _poly_product(p1,p2)
            mc = _matrix_commutator(m1,m2)
            # push only if MC is non-zero (if any value above abs tol)
            push!(h, (poly=pp, matrix=mc))
        end
    end
    return h
end

function _hamiltonian_integrate(h::Vector)
    hi = [(poly=_poly_integrate(p), matrix=m) for (p,m) in h]
    return hi
end

function _hamiltonian_eval(x::Real, h::Vector)
    val = sum(_poly_eval(x,p) .* m for (p,m) in h)
    return val
end

# H that is a poly for X and Z
function _magnus_generator(A::Vector, order::Int)
    #h = []
    Ω1 = _hamiltonian_integrate(A)
    return [Ω1]
end


#=
two_spin_X = Complex{Float64}[0.0 1.0 1.0 0.0; 1.0  0.0 0.0 1.0; 1.0 0.0  0.0 1.0; 0.0 1.0 1.0 0.0]
two_spin_Z = Complex{Float64}[2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0; 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0]

a_2, a_1, a_0 = get_function_coefficients(AS_CIRCULAR.A, 0.1, 0.2)
b_2, b_1, b_0 = get_function_coefficients(AS_CIRCULAR.B, 0.1, 0.2)

H = [(poly=[a_0,a_1,a_2], matrix=-im.*two_spin_X), (poly=[b_0,b_1,b_2], matrix=-im.*two_spin_Z)]
_magnus_generator(H, 1)


two_spin_model = Dict((1,2) => 2)
ρ = QuantumAnnealing.simulate_tmp(two_spin_model, 10.0, AS_CIRCULAR, 2)
ρ = QuantumAnnealing.simulate_tmp(two_spin_model, 1.0, AS_CIRCULAR, mean_tol=1e-4, max_tol=1e-2)
=#

function simulate_tmp(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule, steps::Int; initial_state=nothing, constant_field_x=nothing, constant_field_z=nothing, state_steps=nothing)
    if steps < 2
        error("at least two steps are required by simulate, given $(steps)")
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

    ηs = ones(n)
    hs = zeros(n)

    x_component = sum_x(n)
    z_component = SparseArrays.spzeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component = z_component + sum_z_tup(n, tup, w)
    end

    xz_bracket = lie_bracket(x_component, z_component)

    constant_component = sum_x(n, constant_field_x) + sum_z(n, constant_field_z)
    constant_bracket_x = lie_bracket(x_component, constant_component)
    constant_bracket_z = lie_bracket(z_component, constant_component)

    s_steps = range(0, 1, length=steps)
    R_current = R0
    U = foldl(kron, [IMAT for i = 1:n])

    if track_states
        push!(state_steps, R_current)
    end

    for i in 1:(steps-1)
        s0 = s_steps[i]
        s1 = s_steps[i+1]
        δs = s1 - s0

        a_2, a_1, a_0 = get_function_coefficients(annealing_schedule.A, s0, s1)
        b_2, b_1, b_0 = get_function_coefficients(annealing_schedule.B, s0, s1)

        a_2_shift = a_2
        a_1_shift = a_1 + 2*a_2*s0
        a_0_shift = a_0 + a_1*s0 + a_2*s0^2

        b_2_shift = b_2
        b_1_shift = b_1 + 2*b_2*s0
        b_0_shift = b_0 + b_1*s0 + b_2*s0^2

        # println(a_2, " ", a_1, " ", a_0)
        # println(a_2_shift, " ", a_1_shift, " ", a_0_shift)
        # println()
        
        a_int_eval = _poly_eval(δs, _poly_integrate([a_0_shift,a_1_shift,a_2_shift]))
        b_int_eval = _poly_eval(δs, _poly_integrate([b_0_shift,b_1_shift,b_2_shift]))
        

        #println(a_int_eval, " ", b_int_eval)

        H = [
            (poly=[a_0_shift,a_1_shift,a_2_shift], matrix=Matrix(-im*annealing_time .* x_component)),
            (poly=[b_0_shift,b_1_shift,b_2_shift], matrix=Matrix(-im*annealing_time .* z_component))
        ]

        Ω_list = _magnus_generator(H, 1)
        #display(Ω_list[1])
        Ω1 = _hamiltonian_eval(δs, Ω_list[1])

        #println(Ω1)

        # integA = integral_1_sched(a_2, a_1, a_0, s0, δs)
        # integB = integral_1_sched(b_2, b_1, b_0, s0, δs)
        # integ2 = integral_2_sched(a_2, a_1, a_0, b_2, b_1, b_0, s0, δs)

        # integ2A = a_1*δs^3/6 + a_2*s0*δs^3/3 + a_2*δs^4/6
        # integ2B = b_1*δs^3/6 + b_2*s0*δs^3/3 + b_2*δs^4/6

        # Ω1Sched = integA * x_component + integB * z_component
        # Ω1Const = δs * constant_component
        # Ω1 = Ω1Sched + Ω1Const

        # Ω2Sched = integ2 * xz_bracket 
        # Ω2Const = integ2A * constant_bracket_x + integ2B * constant_bracket_z
        # Ω2 = (Ω2Sched + Ω2Const)/2

        #U_next = exp(Matrix(-im * (annealing_time*Ω1 + (annealing_time^2)*Ω2)))
        #println(Matrix(Ω1))
        U_next = exp(Matrix(Ω1))
        U = U_next * U

        if track_states
            R_current = U * R0 * U'
            push!(state_steps, R_current)
        end
    end

    return U * R0 * U'
end


function simulate_tmp(ising_model::Dict, annealing_time::Real, annealing_schedule::AnnealingSchedule; steps=2, mean_tol=1e-6, max_tol=1e-4, iteration_limit=100, silence=false, state_steps=nothing, kwargs...)
    start_time = time()
    mean_delta = mean_tol + 1.0
    max_delta = max_tol + 1.0

    if !silence
        println()
        println("iter |  steps  |    max(Δ)    |    mean(Δ)   |")
    end

    ρ_prev = simulate(ising_model, annealing_time, annealing_schedule, steps; kwargs...)

    iteration = 1
    while mean_delta >= mean_tol || max_delta >= max_tol
        steps *= 2

        if state_steps != nothing
            empty!(state_steps)
        end

        ρ = simulate_tmp(ising_model, annealing_time, annealing_schedule, steps; state_steps=state_steps, kwargs...)

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


