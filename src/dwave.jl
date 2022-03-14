### Helper Functions for Running Simulations of D-Wave Hardware ###


"""
ground state of sum_i A(0)X_i where A(0) < 0 and B(0) = 0
"""
function initial_state_default_dwave(n)
    return complex(ones(2^n)./(2^(n/2)))
end


"""
An AnnealingSchedule approximating those used in hardware by D-Wave Systems

NOTE: users are strongly encouraged to download annealing schedules for specific
D-Wave Systems devices and load them using `parse_dwave_annealing_schedule`.
"""
const AS_DW_QUADRATIC = AnnealingSchedule(
    (s) -> begin
        if s >= 0.69
            return 0
        else
            return (6.366401*((1.449275)^2*s^2 + (-2.898551)*s + 1.0)*(2.0*π))/-2.0
        end
    end,
    (s) -> (14.55571*(0.85*s^2 + 0.15*s + 0.0)*(2.0*π))/2.0,
    initial_state_default_dwave
)


# construct a piecewise constant function
function _calc_constant_pwp(x_values, y_values)
    @assert(length(x_values) == length(y_values))

    pwp = []
    for (x,y) in zip(x_values, y_values)
        poly = (x_min=x, coefficients=[y])
        push!(pwp, poly)
    end

    return pwp
end


# construct a piecewise linear function
function _calc_linear_pwp(x_values, y_values)
    @assert(length(x_values) == length(y_values))

    points = collect(zip(x_values, y_values))

    pwp = []
    for i in 1:(length(points)-1)
        x1, y1 = points[i]
        x2, y2 = points[i+1]

        m = (y2 - y1)/(x2 - x1)
        b = y1 - m * x1

        poly = (x_min=x1, coefficients=[b,m])
        push!(pwp, poly)
    end

    return pwp
end


# construct a spline-based piecewise quadratic function
function _calc_quadratic_pwp(x_values, y_values)
    @assert(length(x_values) == length(y_values))

    n = length(x_values) - 1
    A = zeros(3n, 3n)
    b = zeros(3n)
    for j = 1:n
        A[j,j] = x_values[j]^2
        A[j,j+n] = x_values[j]
        A[j,j+2n] = 1
        b[j] = y_values[j]

        A[j+n,j] = x_values[j+1]^2
        A[j+n,j+n] = x_values[j+1]
        A[j+n,j+2n] = 1
        b[j+n] = y_values[j+1]
    end

    for j = 1:n-1
        A[j+2n,j] = 2*x_values[j+1]
        A[j+2n,j+1] = -2*x_values[j+1]
        A[j+2n,j+n] = 1
        A[j+2n,j+n+1] = -1
        b[j+2n] = 0
    end
    A[3n,n] = 1
    b[3n] = 0

    coefficient_matrix = reshape(A\b, (n,3))

    pwp = []
    for i in 1:(length(x_values)-1)
        a,b,c = coefficient_matrix[i,:]

        poly = (x_min=x_values[i], coefficients=[c,b,a])
        push!(pwp, poly)
    end

    return pwp
end

# evaluator for a generic piecewise polynomial function
function _eval_pwp_function(x, pwp_points)
    if x <= pwp_points[1].x_min
        pwp = pwp_points[1]
    elseif x >= pwp_points[end].x_min
        pwp = pwp_points[end]
    else
        index = 1
        while x >= pwp_points[index].x_min
            index += 1
        end
        pwp = pwp_points[index - 1]
    end

    value = sum(v*x^(i-1) for (i,v) in enumerate(pwp.coefficients))

    return value
end


"""
function to take a CSV of DWave annealing schedule values and convert it into
an annealing schedule usable by the simulator.
valid values for interpolation are :none, :linear, :quadratic
"""
function read_dwave_annealing_schedule(infile; header=1, delim=',', interpolation=:linear, initial_state=initial_state_default_dwave)
    s_values = Float64[]
    a_values = Float64[]
    b_values = Float64[]

    csv = CSV.File(infile, header=header, delim=delim)

    if length(csv) < 2
        error("The CSV annealing schedule requires at least 2 rows, given $(length(csv))")
    end

    for (i,row) in enumerate(csv)
        if length(row) < 3
            error("The CSV annealing schedule requires at least 3 values in each row, given $(length(row)) on row $(i)")
        end

        push!(s_values, row[1])
        push!(a_values, row[2])
        push!(b_values, row[3])
    end

    # change from GHz to natural units (annealing time in nanoseconds)
    a_values = a_values .* (2.0*π)
    b_values = b_values .* (2.0*π)

    # rescale and swap sign based on D-Wave hamiltonian convention
    # https://docs.dwavesys.com/docs/latest/c_qpu_annealing.html
    a_values = a_values ./ -2.0
    b_values = b_values ./ 2.0

    if interpolation == :none
        a_pwp = _calc_constant_pwp(s_values, a_values)
        b_pwp = _calc_constant_pwp(s_values, b_values)

    elseif interpolation == :linear
        a_pwp = _calc_linear_pwp(s_values, a_values)
        b_pwp = _calc_linear_pwp(s_values, b_values)

    elseif interpolation == :quadratic
        a_pwp = _calc_quadratic_pwp(s_values, a_values)
        b_pwp = _calc_quadratic_pwp(s_values, b_values)

    else
        error("interpolation of type \":$(interpolation)\" is not supported by parse_dwave_annealing_schedule")
    end

    return AnnealingSchedule(
        (s) -> _eval_pwp_function(s, a_pwp),
        (s) -> _eval_pwp_function(s, b_pwp),
        initial_state
    )
end


"""
Function to modify an existing annealing schedule to use a customized
annealing schedule (asch).  These parameters are the same as those
used in a dwisc call or a dwave schedule.
Inputs:
annealing_schedule - annealing_schedule

Parameters:
asch - This is the annealing-schedule parameter.  This is a list of tuples of the form
       [(s₀,s_effective₀), (s₀,s_effective₁), ..., (sₙ,s_effectiveₙ)].
"""
function annealing_protocol_dwave(annealing_schedule::AnnealingSchedule; asch=[(0,0) (1,1)])
    asch_slopes = zeros(length(asch)-1)
    for i in 1:(length(asch)-1)
        s0,s_eff_0 = asch[i]
        s1,s_eff_1 = asch[i+1]
        asch_slopes[i] = (s_eff_1 - s_eff_0)/(s1 - s0)
    end

    #branchless piecewise function using linear interpolation from y = m*(x-x0) + y0
    function asch_func(s)
        return sum([(asch_slopes[i]*(s-asch[i][1]) + asch[i][2]) * (asch[i][1] <= s < asch[i+1][1]) for i = 1:(length(asch)-1)]) + ((s == asch[end][1])*asch[end][2])
    end

    new_annealing_schedule = AnnealingSchedule(
        s -> annealing_schedule.A(asch_func(s)),
        s -> annealing_schedule.B(asch_func(s))
    )
    return new_annealing_schedule
end


"""
function that allows for simulation from a bqpjson data file
"""
function simulate_bqpjson(infile, outfile, annealing_time, annealing_schedule; simulated_num_reads=1e17, scale=1.0, kwargs...)
    ising_model, qubit_ids = read_bqpjson(infile)
    n = length(qubit_ids)
    for (k,v) in ising_model
        ising_model[k] = v*scale
    end

    ρ = simulate(ising_model, annealing_time, annealing_schedule; kwargs...)

    write_dwisc(outfile, ρ, ising_model, qubit_ids, simulated_num_reads=simulated_num_reads, annealing_time=annealing_time)
end


"""
function that allows for simulation with x and z noise from a bqpjson data file.
The `x_bias` and `z_bias` parameters provide vectors of noise realizations.
"""
function simulate_bqpjson_noisy(infile, outfile, annealing_time, annealing_schedule; simulated_num_reads=1e17, scale=1.0, x_bias::Vector{<:Any}=[], z_bias::Vector{<:Any}=[], kwargs...)
    if length(x_bias) > 0 && length(z_bias) > 0 && length(x_bias) != length(z_bias)
        error("x_bias and z_bias require the same number of parameters given, $(length(x_bias)) and $(length(z_bias)) respectively")
    end

    if length(x_bias) > 0 && length(z_bias) == 0
        z_bias = [0.0 for i in length(x_bias)]
    end

    if length(z_bias) > 0 && length(x_bias) == 0
        x_bias = [0.0 for i in length(z_bias)]
    end

    x_bias = x_bias .* scale
    z_bias = z_bias .* scale

    ising_model, qubit_ids = read_bqpjson(infile)
    n = length(qubit_ids)
    for (k,v) in ising_model
        ising_model[k] = v*scale
    end

    accumulator = zeros(2^n,2^n)
    for shot in 1:length(x_bias)
        x_field = x_bias[shot]
        z_field = z_bias[shot]

        ρ = simulate(ising_model, annealing_time, annealing_schedule, constant_field_x=[x_field], constant_field_z=[z_field]; kwargs...)

        accumulator = accumulator + ρ
    end
    mean_ρ = accumulator ./ length(x_bias)

    write_dwisc(outfile, mean_ρ, ising_model, qubit_ids, simulated_num_reads=simulated_num_reads, annealing_time=annealing_time)
end


"""
This function reads in a bqpjson file and generates the ising model dictionary 
inputs:
bqpjson::String - a bqpjson file (v1.0.0) that can be run on D-Wave hardware
outputs:
n - number of qubits
ising_model::Dict{Tuple => Float64} - Dictionary of qubits and couplings to weights
mapping:Dict{Int => Int} - mapping of the qubits in the bqpjson file to simulation qubits
"""
function read_bqpjson(bqpjson_file::String)
    bqp_data = JSON.parsefile(bqpjson_file)

    # maps qubit ids in 1-to-n to bqpjson ids
    qubit_ids = sort(bqp_data["variable_ids"])

    # maps bqpjson ids to qubit ids
    mapping = Dict(qid => i for (i, qid) in enumerate(qubit_ids))

    ising_model = Dict{Any,Float64}()

    for lt in bqp_data["linear_terms"]
        i = lt["id"]
        ising_model[(mapping[i],)] = lt["coeff"]
    end

    for qt in bqp_data["quadratic_terms"]
        i = qt["id_tail"]
        j = qt["id_head"]
        ising_model[(mapping[i],mapping[j])] = qt["coeff"]
    end

    return ising_model, qubit_ids
end


function write_dwisc(outfile::String, ρ, ising_model, qubit_ids; simulated_num_reads=1e17, annealing_time=1000)
    dwisc_data = Dict()
    dwisc_data["variable_ids"] = qubit_ids
    dwisc_data["metadata"] = Dict(
        "sim_annealing_time" => annealing_time/1000
    )

    solutions = dwisc_data["solutions"] = Any[]

    probs = z_measure_probabilities(ρ)

    n = length(qubit_ids)

    for state_int in 0:(2^n-1)
        prob = probs[state_int+1]
        spin_vector = int_to_spin(state_int, pad=n)
        energy = eval_ising_state_energy(spin_vector, ising_model)

        sol_data = Dict(
            "energy" => energy,
            "prob" => prob, 
            "num_occurrences" => round(Int, simulated_num_reads * prob),
            "solution" => spin_vector
        )
        push!(solutions, sol_data)
    end

    sort!(solutions, by=(x) -> x["prob"], rev=true)

    json_string = JSON.json(dwisc_data)

    open(outfile,"w") do io
        write(io, json_string)
    end
end
