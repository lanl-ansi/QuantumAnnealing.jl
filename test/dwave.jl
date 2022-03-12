
@testset "simulate bqpjson" begin

    @testset "1 qubit, function schedule" begin
        annealing_time = 1000.0
        annealing_schedule = AS_CIRCULAR
        #steps = 2500

        ising_intended = Dict((1,) => 1)

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, silence=true)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate(ising_intended, annealing_time, annealing_schedule, silence=true)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[2], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "1 qubit, dwave schedule" begin
        annealing_time = 1000.0
        annealing_schedule = AS_DW_QUADRATIC

        ising_intended = Dict((1,) => 1)

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, silence=true)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate(ising_intended, annealing_time, annealing_schedule, silence=true)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[2], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "2 qubit, function schedule" begin
        annealing_time = 1000.0
        annealing_schedule = AS_CIRCULAR

        # the ising model that is encoded in bqpjson_2q.json
        ising_intended = Dict((1,) => -1, (2,) => -1, (1,2) => -1)

        bqpjson_file = "data/bqpjson_2q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, silence=true)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate(ising_intended, annealing_time, annealing_schedule, silence=true)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[1], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [304, 308]
    end

    @testset "1 qubit, function schedule, z noise" begin
        annealing_time = 100.0
        annealing_schedule = AS_CIRCULAR
        numshots = 10

        # random numbers between -0.75 and 0.75
        z_bias = (Random.rand(numshots) .- 0.5)*1.5

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson_noisy(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, z_bias=z_bias, silence=true)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        @test dwisc_data["solutions"][1]["prob"] > 0.70 # true w.h.p.
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "1 qubit, function schedule, x and z noise" begin
        annealing_time = 100.0
        annealing_schedule = AS_CIRCULAR
        numshots = 10

        # random numbers between -0.75 and 0.75
        x_bias = (Random.rand(numshots) .- 0.5)*1.5
        z_bias = (Random.rand(numshots) .- 0.5)*1.5

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson_noisy(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, x_bias=x_bias, z_bias=z_bias, silence=true)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        @test dwisc_data["solutions"][1]["prob"] > 0.70 # true w.h.p.
        @test dwisc_data["variable_ids"] == [100]
    end

end
