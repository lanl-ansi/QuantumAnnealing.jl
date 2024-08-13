
@testset "simulate, 1 qubit" begin
    @testset "function schedule, default anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR, 100, 4)
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "analytical solution, adaptive steps" begin
        ρ = simulate(one_spin_model, 1.0, AS_CIRCULAR, mean_tol=1e-7, max_tol=1e-7, silence=true)
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "function schedule, fast anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 0.5, AS_CIRCULAR, 1000, 4)
        @test isapprox(one_spin_ρ(0.5), ρ)
    end

    @testset "function schedule, slow anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 2.0, AS_CIRCULAR, 1000, 4)
        @test isapprox(one_spin_ρ(2.0), ρ)
    end

    @testset "fractional field value" begin
        ρ_target = [0.420186+0.0im -0.409634+0.275372im; -0.409634-0.275372im 0.579814+2.77556e-17im]
        ρ = simulate_magnus_optimized(Dict((1,) => 0.5), 1.0, AS_CIRCULAR, 100, 4)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "field value above 1.0" begin
        ρ_target = [0.291065-2.77556e-17im 0.114524+0.43958im; 0.114524-0.43958im 0.708935+5.55112e-17im]
        ρ = simulate_magnus_optimized(Dict((1,) => 1.5), 1.0, AS_CIRCULAR, 100, 4)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule, constant terms" begin
        ρ_target = [0.0578906+1.38778e-17im -0.165069-0.165202im; -0.165069+0.165202im 0.942109+0.0im]
        ρ = simulate_magnus_generic(one_spin_model, 1.0, AS_CIRCULAR, 100, 4, constant_field_x = [1], constant_field_z = [1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule, different initial state" begin
        ρ_target = [0.326006+0.0im -0.413095-0.221536im; -0.413095+0.221536im 0.673994+0.0im]
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR, 100, 4, initial_state = [0,1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_LINEAR), default anneal time" begin
        ρ_target = [0.422382+2.77556e-17im -0.278818+0.40772im; -0.278818-0.40772im 0.577618+2.77556e-17im]
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_LINEAR, 100, 4)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_QUADRATIC), default anneal time" begin
        ρ_target = [0.489037+0.0im -0.393381+0.308433im; -0.393381-0.308433im 0.510963+5.55112e-17im]
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_QUADRATIC, 100, 4)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_DW_QUADRATIC), default anneal time" begin
        ρ_target = [0.0162536+0.0im 0.121897-0.0336245im; 0.121897+0.0336245im  0.983746+2.77556e-17im]
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_DW_QUADRATIC, 100, 4)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end


    @testset "function schedule, too few steps, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR, 10, 4)

        # NOTE, atol required due to too few iterations
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-5)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
    end

    @testset "csv schedule pwq, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR_pwq_csv_1000, 100, 4)

        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "csv schedule pwq, low resolution, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR_pwq_csv_100, 100, 4)

        # NOTE, atol required due to pwq approx_IMATion in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-7)
    end

    @testset "csv schedule pwl, low resolution, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR_pwl_csv_100, 100, 4)

        # NOTE, atol required due to pwl approx_IMATion in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-4)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-5)
    end

    @testset "csv schedule pwc, low resolution, analytical solution" begin
        ρ = simulate_magnus_optimized(one_spin_model, 1.0, AS_CIRCULAR_pwc_csv_100, 100, 4)

        # NOTE, atol required due to pwc approx_IMATion in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-2)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-3)
    end

    @testset "reverse annealing test" begin
        ρ_target = [0.88389+0.0im -0.197484-0.252247im; -0.197484+0.252247im 0.11611+0.0im]

        asch = [(0.0, 1.0), (0.5, 0.5), (1.0, 1.0)]
        mod_asch = annealing_protocol_dwave(AS_LINEAR, asch=asch)

        ρ = simulate_magnus_optimized(Dict((1,) => 0.1), 5.0, mod_asch, 100, 4, initial_state = [0,1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "collect probability trajectory, nonadaptive" begin
        ρ_list = []
        steps=100
        ρ = simulate_magnus_optimized(one_spin_model, 1, AS_CIRCULAR, steps, 4, state_steps=ρ_list)

        @test isapprox(one_spin_ρ(1.0), ρ)
        @test length(ρ_list) == steps
        @test isapprox(ρ_list[1], (initial_state_default(1) * initial_state_default(1)'))
        @test isapprox(ρ_list[steps], one_spin_ρ(1.0))
    end

    @testset "collect probability trajectory, adaptive" begin
        ρ_list = []
        ρ = simulate(one_spin_model, 1.0, AS_CIRCULAR, mean_tol=1e-7, max_tol=1e-7, silence=true, state_steps=ρ_list)

        @test isapprox(one_spin_ρ(1.0), ρ)
        @test length(ρ_list) == 64
        @test isapprox(ρ_list[1], (initial_state_default(1) * initial_state_default(1)'))
        @test isapprox(ρ_list[64], one_spin_ρ(1.0))
    end
end


@testset "simulate, 2 qubit" begin
    @testset "function schedule, default anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(two_spin_model, 1.0, AS_CIRCULAR, 100, 4)
        @test isapprox(two_spin_ρ(1.0), ρ)
    end

    @testset "function schedule, fast anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(two_spin_model, 0.5, AS_CIRCULAR, 100, 4)
        @test isapprox(two_spin_ρ(0.5), ρ)
    end

    @testset "function schedule, slow anneal time, analytical solution" begin
        ρ = simulate_magnus_optimized(two_spin_model, 2.0, AS_CIRCULAR, 100, 4)
        @test isapprox(two_spin_ρ(2.0), ρ)
    end
end

@testset "analytic models, unitary check" begin
    @testset "one qubit model" begin
        s_vals = 0.0:0.5:1.0
        for s in s_vals
            ρ = one_spin_ρ(100.0, s=s)
            @test tr(ρ) == 1
        end
    end

    @testset "two qubit model" begin
        s_vals = 0.0:0.5:1.0
        for s in s_vals
            ρ = two_spin_ρ(100.0, s=s)
            @test tr(ρ) == 1
        end
    end
end

@testset "simulate, generic magnus expansion" begin
    @testset "1 qubit, adaptive, orders 1 to 8" begin
        for i in 1:8
            ρ = simulate(one_spin_model, 1.0, AS_CIRCULAR, order=i, max_tol=1e-9, silence=true)
            @test isapprox(one_spin_ρ(1.0), ρ)
        end
    end

    @testset "2 qubit, adaptive, orders 1 to 8" begin
        for i in 1:8
            ρ = simulate(two_spin_model, 1.0, AS_CIRCULAR, order=i, max_tol=1e-9, silence=true)
            @test isapprox(two_spin_ρ(1.0), ρ)
        end
    end

    @testset "5 qubit, fixed_order solver, from 1 to 4" begin
        ising_model = Dict((1,) => -1, (1,2) => -1, (1,3) => -1, (1,4) => -1, (1,5) => -1, (2,3) => 1, (4,5) => 1)

        for i in 1:4
            ρ_target = simulate_magnus_optimized(ising_model, 2.0, AS_CIRCULAR, 2, i)
            ρ = simulate_magnus_generic(ising_model, 2.0, AS_CIRCULAR, 2, i)
            @test isapprox(ρ_target, ρ)
        end

    end

    @testset "5 qubit, hardcoded forth order solver, function schedules" begin
        ising_model = Dict((1,) => -1, (1,2) => -1, (1,3) => -1, (1,4) => -1, (1,5) => -1, (2,3) => 1, (4,5) => 1)

        ρ_target = simulate_magnus_optimized(ising_model, 2.0, AS_LINEAR, 2, 4)
        ρ = simulate_magnus_generic(ising_model, 2.0, AS_LINEAR, 2, 4)
        @test isapprox(ρ_target, ρ)

        ρ_target = simulate_magnus_optimized(ising_model, 2.0, AS_QUADRATIC, 2, 4)
        ρ = simulate_magnus_generic(ising_model, 2.0, AS_QUADRATIC, 2, 4)
        @test isapprox(ρ_target, ρ)

        ρ_target = simulate_magnus_optimized(ising_model, 2.0, AS_DW_QUADRATIC, 2, 4)
        ρ = simulate_magnus_generic(ising_model, 2.0, AS_DW_QUADRATIC, 2, 4)
        @test isapprox(ρ_target, ρ)
    end
end


@testset "simulate, multi-qubit" begin

    @testset "2 qubit, function schedules (AS_CIRCULAR, AS_LINEAR, AS_QUADRATIC), near adiabatic limit" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 100.0
        steps = 100

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_CIRCULAR, steps, 4)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_LINEAR, steps, 4)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_QUADRATIC, steps, 4)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit, function schedules (AS_DW_QUADRATIC), near adiabatic limit" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 10.0
        steps = 1000

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_DW_QUADRATIC, steps, 4)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit, csv schedule, near adiabatic limit" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        ρ = simulate_magnus_optimized(ising_model, 100.0, AS_CIRCULAR_pwl_csv_1000, 100, 4)

        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)
    end

    @testset "2 qubit, function schedule, fractional values, near adiabatic limit" begin
        ising_model = Dict((1,) => 0.5, (1,2) => -0.75)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        ρ = simulate_magnus_optimized(ising_model, 50.0, AS_CIRCULAR, 100, 4)

        # NOTE, Non-1.0 due to annealing time not being in the fill adiabatic limit
        @test real(tr(ρ*min_ρ)) >= 0.999

        @test isapprox(min_ρ, ρ, atol=1e-1)
        @test !isapprox(min_ρ, ρ, atol=1e-2)
    end

    @testset "2 qubit probability trajectory, nonadaptive" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        n = QuantumAnnealing._check_ising_model_ids(ising_model)

        annealing_time = 10.0
        steps = 1000

        ρ_list = []
        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_DW_QUADRATIC, steps, 4, state_steps=ρ_list)
        @test real(ρ[4,4]) >= 0.9999
        @test length(ρ_list) == steps
        @test isapprox(ρ_list[1], (initial_state_default_dwave(n) * initial_state_default_dwave(n)'))
        @test real(ρ_list[steps][4,4]) >= 0.9999
        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit probability trajectory, adaptive" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        n = QuantumAnnealing._check_ising_model_ids(ising_model)

        annealing_time = 10.0

        ρ_list = []
        ρ = simulate(ising_model, annealing_time, AS_DW_QUADRATIC, silence=true, state_steps=ρ_list)
        @test real(ρ[4,4]) >= 0.9999
        @test length(ρ_list) == 512
        @test isapprox(ρ_list[1], (initial_state_default_dwave(n) * initial_state_default_dwave(n)'))
        @test real(ρ_list[end][4,4]) >= 0.9999
        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "3 qubit, degenerate, function schedules (AS_CIRCULAR, AS_QUADRATIC), near adiabatic limit" begin
        # ring of disagrees => 6 ground states
        ising_model = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))

        annealing_time = 100.0
        steps = 500

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_CIRCULAR, steps, 4)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0, atol=1e-7)
        end

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_QUADRATIC, steps, 4)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0)
        end
    end

    @testset "3 qubit, degenerate, function schedules (AS_DW_QUADRATIC), near adiabatic limit" begin
        # ring of disagrees => 6 ground states
        ising_model = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)

        H = ising_hamiltonian(ising_model)
        evals,evecs = eigen(Matrix(H))

        annealing_time = 100.0
        steps = 4000

        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_DW_QUADRATIC, steps, 4)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0)
        end
    end

    @testset "3 qubit, csv schedule, near adiabatic limit" begin
        #ring of disagrees => 6 ground states
        ising_model = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)
        annealing_time = 100.0
        steps = 500

        ρ_target = simulate_magnus_optimized(ising_model, annealing_time, AS_CIRCULAR, steps, 4)
        ρ = simulate_magnus_optimized(ising_model, annealing_time, AS_CIRCULAR_pwq_csv_1000, steps, 4)

        @test isapprox(ρ_target, ρ, atol=1e-7)
        @test !isapprox(ρ_target, ρ, atol=1e-8)
    end
end


@testset "2 qubit, print z state probabilities" begin
    mktemp() do path,io
        out = stdout
        err = stderr
        redirect_stdout(io)
        redirect_stderr(io)

        ρ = simulate_magnus_optimized(two_spin_model, 1.0, AS_CIRCULAR, 100, 4)
        print_z_state_probabilities(ρ)
        print_z_state_probabilities(ρ, sort=true)
        print_z_state_probabilities(ρ, sort=true, limit=2)

        flush(io)
        redirect_stdout(out)
        redirect_stderr(err)
    end
end


@testset "simulate_de" begin

    @testset "1 qubit, function schedule, analytical solution" begin
        ρ = simulate_de(one_spin_model, 1.0, AS_CIRCULAR, 1e-6)
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "1 qubit, function schedule, analytical solution with kwargs" begin
        ρ = simulate_de(one_spin_model, 1.0, AS_CIRCULAR, 1e-6, saveat=[1])
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "2 qubit, function schedule, analytical solution" begin
        ρ = simulate_de(two_spin_model, 1.0, AS_CIRCULAR, 1e-6)
        @test isapprox(two_spin_ρ(1.0), ρ)
    end

    @testset "1 qubit, function schedule, analytical solution, with kwargs" begin
        ρ = simulate_de(one_spin_model, 1.0, AS_CIRCULAR, 1e-6, saveat=[1])
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "1 qubit, function schedule, constant terms" begin
        ρ = simulate_magnus_generic(one_spin_model, 1.0, AS_CIRCULAR, 100, 4, constant_field_x=[1], constant_field_z=[1])
        ρ_de = simulate_de(one_spin_model, 1.0, AS_CIRCULAR, 1e-6, constant_field_x=[1], constant_field_z=[1])
        @test isapprox(ρ, ρ_de, atol=1e-7)
    end

    @testset "1 qubit, csv schedule, analytical solution" begin
        ρ = simulate_de(one_spin_model, 1.0, AS_CIRCULAR_pwl_csv_1000, 1e-6)
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
    end

    @testset "2 qubit, function schedule" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_magnus_optimized(ising_model, 1.0, AS_CIRCULAR, 100, 4)
        ρ_de = simulate_de(ising_model, 1.0, AS_CIRCULAR, 1e-6)
        @test isapprox(ρ, ρ_de, atol=1e-7)
    end

    @testset "2 qubit, function schedule, long annealing time" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_magnus_optimized(ising_model, 1000.0, AS_CIRCULAR, 4000, 4)
        ρ_de = simulate_de(ising_model, 1000.0, AS_CIRCULAR, 1e-8)
        @test isapprox(ρ, ρ_de, atol=1e-6)
    end

    @testset "2 qubit, function schedule, adaptive tolerance" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_magnus_optimized(ising_model, 100.0, AS_CIRCULAR, 500, 4)
        ρ_de = simulate_de(ising_model, 100.0, AS_CIRCULAR, silence=true)
        @test isapprox(ρ, ρ_de, atol=1e-6)
    end
end
