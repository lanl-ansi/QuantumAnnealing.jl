using DelimitedFiles

# add a factor of 2x and sign inversion based on D-Wave Hamiltonian convention
# https://docs.dwavesys.com/docs/latest/c_qpu_annealing.html
A = (s) -> -2.0*cos(π/2*s)/(2.0*π)
B = (s) ->  2.0*sin(π/2*s)/(2.0*π)

header = ["s", "A(s)", "B(s)"]

open("trig_sched_100.csv", "w") do io
    table = Any[header]
    for s in range(0, 1, length=100)
        push!(table, [s, A(s), B(s)])
    end
    writedlm(io, table, ',')
end

open("trig_sched_1000.csv", "w") do io
    table = Any[header]
    for s in range(0, 1, length=1000)
        push!(table, [s, A(s), B(s)])
    end
    writedlm(io, table, ',')
end
