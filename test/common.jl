### Helper Functions for Tests ###

I = [1 0; 0 1]
X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]



### An analytical solution to the following 1 qubit model
### using the AS_CIRCULAR annealing schedule

one_spin_model = Dict((1,) => 1)

function one_spin_H(s)
    return cos(s*π/2)*X + sin(s*π/2)*Z
end

function one_spin_ρ(T; s=1.0)
    r = 2.0 * T
    ω0 = π/2
    ω1 = sqrt(4.0 * T^2 + π^2/4.0)

    a = -((r^2.0) + (ω0^2.0)*cos(ω1*s))*cos(ω0*s) - ω0*ω1*sin(ω0*s)*sin(ω1*s)
    b = -r*ω0*(1.0 - cos(ω1*s))
    c = -((r^2.0) + (ω0^2.0)*cos(ω1*s))*sin(ω0*s) + ω0*ω1*cos(ω0*s)*sin(ω1*s)

    ψs = [a; b; c;] ./ (ω1^2.0)

    return (I + ψs[1]*X + ψs[2]*Y + ψs[3]*Z)/2.0
end



### An analytical solution to the following 2 qubit model
### using the AS_CIRCULAR annealing schedule

two_spin_model = Dict((1,2) => 2)

function two_spin_H(s)
    X = Complex{Float64}[0.0 1.0 1.0 0.0; 1.0  0.0 0.0 1.0; 1.0 0.0  0.0 1.0; 0.0 1.0 1.0 0.0]
    Z = Complex{Float64}[2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0; 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0]
    return cos(s*π/2)*X + sin(s*π/2)*Z
end

function two_spin_ρ(t; s=1.0)
    s0 = 1/2*cos(π/4*s*sqrt(1+64*t^2/π^2))*sqrt(1-sin(π/2*s)) +
        (8im*t*sqrt(1-sin(π/2*s))/π + sqrt(1+sin(π/2*s)))*sin(π/4*s*sqrt(1+64*t^2/π^2)) / (2*sqrt(1+64*t^2/π^2))
    s1 = -(cos(π/4*s*sqrt(1+64*t^2/π^2))*(1+sin(π/2*s)) + 
        ((-cos(π/2*s) + 8im*t*(1+sin(π/2*s))/π)*(sin(π/4*s*sqrt(1+64*t^2/π^2))))/(sqrt(1+64*t^2/π^2))
        )/(2*sqrt(1+sin(π/2*s)))
    s2 = s1
    s3 = s0

    sv = [s0, s1, s2, s3]

    return sv*sv'
end



s_100 = range(0, 1, length=100)
s_10000 = range(0, 1, length=10000)

AS_CIRCULAR_pwc_csv_100 = read_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:none, initial_state=initial_state_default)
AS_CIRCULAR_pwc_csv_1000 = read_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:none, initial_state=initial_state_default)

AS_CIRCULAR_pwl_csv_100 = read_dwave_annealing_schedule("data/trig_sched_100.csv", initial_state=initial_state_default)
AS_CIRCULAR_pwl_csv_1000 = read_dwave_annealing_schedule("data/trig_sched_1000.csv", initial_state=initial_state_default)

AS_CIRCULAR_pwq_csv_100 = read_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:quadratic, initial_state=initial_state_default)
AS_CIRCULAR_pwq_csv_1000 = read_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:quadratic, initial_state=initial_state_default)


function sum_zizj(n, J::Dict)
    sum_zs = zeros(2^n)
    for (i,j) in keys(J)
        J_val = J[(i,j)]
        sum_zs += zizj_vectorized(n, i, j, J_val)
    end
    return diagm(sum_zs)
end

function zizj_vectorized(n::Int, i::Int, j::Int, J_val)
    Ivec = [1.0+0im,1.0+0im]
    Zvec = [1.0+0im,-1.0+0im]
    matvec = [(k == i || k == j) ? Zvec : Ivec for k in n:-1:1]
    return J_val * foldl(kron, matvec)
end

function ising_hamiltonian(ising_model)
    n = QuantumAnnealing._check_ising_model_ids(ising_model)

    z_component = zeros(2^n, 2^n)
    for (tup,w) in ising_model
        z_component += QuantumAnnealing._kron_Z(n, tup, w)
    end

    return z_component
end
