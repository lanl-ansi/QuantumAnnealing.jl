### Helper Functions for Tests ###

I = [1 0; 0 1]
X = [0 1; 1 0]
Y = [0 -im; im 0]
Z = [1 0; 0 -1]

function single_spin_h(s; T=1.0)
    r = 2.0 * T
    ω0 = π/2
    ω1 = sqrt(4.0 * T^2 + π^2/4.0)

    return r*cos(ω0*s)*[1.0;0.0;0.0] + r*sin(ω0*s)*[0.0;0.0;1.0]
end

function single_spin_ψ(s; T=1.0)
    r = 2.0 * T
    ω0 = π/2
    ω1 = sqrt(4.0 * T^2 + π^2/4.0)

    a = -((r^2.0) + (ω0^2.0)*cos(ω1*s))*cos(ω0*s) - ω0*ω1*sin(ω0*s)*sin(ω1*s)
    b = -r*ω0*(1.0 - cos(ω1*s))
    c = -((r^2.0) + (ω0^2.0)*cos(ω1*s))*sin(ω0*s) + ω0*ω1*cos(ω0*s)*sin(ω1*s)
    return [a; b; c;] ./ (ω1^2.0)
end

function single_spin_H(s)
    return cos(s*π/2)*X + sin(s*π/2)*Z
end

function single_spin_ρ(s; T=1.0)
    ψs = single_spin_ψ(s, T=T)
    return (I + ψs[1]*X + ψs[2]*Y + ψs[3]*Z)/2.0
end



two_spin_model = Dict((1,2) => 2)

function two_spin_H(s)
    X = Complex{Float64}[0.0 1.0 1.0 0.0; 1.0 0.0 0.0 1.0; 1.0 0.0 0.0 1.0; 0.0 1.0 1.0 0.0]
    Z = Complex{Float64}[2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0; 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0]
    return cos(s*π/2)*X + sin(s*π/2)*Z
end

function two_spin_porbs(t)
    t = 2*t/π
    p1 = sin(1/4*π*sqrt(1 + 16*t^2))^2/(2 + 32*t^2)
    p2 = (1 + 32*t^2 + cos(1/2*π*sqrt(1 + 16*t^2)))/(4 + 64*t^2)
    return [p1,p2,p2,p1]
end

function two_spin_ρ(s; t=1.0)
    t = 2*t/π

    s0 = (
        cosh(π/4*s*sqrt(-1-16*t^2 + 0im)) * sqrt((1+16*t^2)*(-1+sin((π*s)/2)) + 0im) +
        (4im * t * sqrt(1-sin((π*s)/2)) + sqrt(1+sin((π*s)/2)))*sinh(π/4*s*sqrt(-1-16*t^2 + 0im))
    )/(2*sqrt(-1-16*t^2 + 0im))
    s1 = exp(-π/4*s*sqrt(-1-16*t^2 + 0im))*(
        (-1+exp(π/2*s*sqrt(-1-16*t^2 + 0im)))*sqrt(1+sin((π*s)/2))*(sec((π*s)/2)-tan((π*s)/2)) -
        (
            4im * (-1+exp(π/2*s*sqrt(-1-16*t^2 + 0im))) * t +
            (1+exp(π/2*s*sqrt(-1-16*t^2 + 0im)))*sqrt(-1-16*t^2 + 0im)
        )*sqrt(1-sin((π*s)/2))*(sec((π*s)/2)+tan((π*s)/2))
    )/(4*sqrt(-1-16*t^2 + 0im))
    s2 = s1
    s3 = s0

    sv = [s0, s1, s2, s3]

    return sv*sv'
end

single_spin_model = Dict((1,) => 1)
single_spin_analytic_ρ = single_spin_ρ(1.0)

two_spin_analytic_ρ = two_spin_ρ(1.0-1e-6)


s_100 = range(0, 1, length=100)
s_10000 = range(0, 1, length=10000)

AS_CIRCULAR_pwc_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:none, initial_state=default_initial_state)
AS_CIRCULAR_pwc_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:none, initial_state=default_initial_state)

AS_CIRCULAR_pwl_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", initial_state=default_initial_state)
AS_CIRCULAR_pwl_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", initial_state=default_initial_state)

AS_CIRCULAR_pwq_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:quadratic, initial_state=default_initial_state)
AS_CIRCULAR_pwq_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:quadratic, initial_state=default_initial_state)

