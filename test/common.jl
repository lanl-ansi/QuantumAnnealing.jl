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


single_spin_model = Dict((1,) => 1)
single_spin_analytic_ρ = single_spin_ρ(1.0)
single_spin_analytic_prob = real(tr(single_spin_analytic_ρ * [0 0; 0 1]))

s_100 = range(0, 1, length=100)
s_10000 = range(0, 1, length=10000)

AS_CIRCULAR_pwc_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:none, initial_state=default_initial_state)
AS_CIRCULAR_pwc_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:none, initial_state=default_initial_state)

AS_CIRCULAR_pwl_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", initial_state=default_initial_state)
AS_CIRCULAR_pwl_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", initial_state=default_initial_state)

AS_CIRCULAR_pwq_csv_100 = parse_dwave_annealing_schedule("data/trig_sched_100.csv", interpolation=:quadratic, initial_state=default_initial_state)
AS_CIRCULAR_pwq_csv_1000 = parse_dwave_annealing_schedule("data/trig_sched_1000.csv", interpolation=:quadratic, initial_state=default_initial_state)

