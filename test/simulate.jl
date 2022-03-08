
@testset "simulate, 1 qubit" begin
    @testset "function schedule, default anneal time, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR, 100)
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "analytical solution, adaptive steps" begin
        ρ = simulate(one_spin_model, 1.0, AS_CIRCULAR, mean_tol=1e-7, max_tol=1e-7, silence=true)
        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "function schedule, fast anneal time, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 0.5, AS_CIRCULAR, 1000)
        @test isapprox(one_spin_ρ(0.5), ρ)
    end

    @testset "function schedule, slow anneal time, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 2.0, AS_CIRCULAR, 1000)
        @test isapprox(one_spin_ρ(2.0), ρ)
    end

    @testset "fractional field value" begin
        ρ_target = [0.420186+0.0im -0.409634+0.275372im; -0.409634-0.275372im 0.579814+2.77556e-17im]
        ρ = simulate_o2(Dict((1,) => 0.5), 1.0, AS_CIRCULAR, 100)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "field value above 1.0" begin
        ρ_target = [0.291065-2.77556e-17im 0.114524+0.43958im; 0.114524-0.43958im 0.708935+5.55112e-17im]
        ρ = simulate_o2(Dict((1,) => 1.5), 1.0, AS_CIRCULAR, 100)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule, constant terms" begin
        ρ_target = [0.0578906+1.38778e-17im -0.165069-0.165202im; -0.165069+0.165202im 0.942109+0.0im]
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR, 100, constant_field_x = [1], constant_field_z = [1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule, different initial state" begin
        ρ_target = [0.326006+0.0im -0.413095-0.221536im; -0.413095+0.221536im 0.673994+0.0im]
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR, 100, initial_state = [0,1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_LINEAR), default anneal time" begin
        ρ_target = [0.422382+2.77556e-17im -0.278818+0.40772im; -0.278818-0.40772im 0.577618+2.77556e-17im]
        ρ = simulate_o2(one_spin_model, 1.0, AS_LINEAR, 100)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_QUADRATIC), default anneal time" begin
        ρ_target = [0.489037+0.0im -0.393381+0.308433im; -0.393381-0.308433im 0.510963+5.55112e-17im]
        ρ = simulate_o2(one_spin_model, 1.0, AS_QUADRATIC, 100)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "function schedule (AS_DW_QUADRATIC), default anneal time" begin
        ρ_target = [0.0162536+0.0im 0.121897-0.0336245im; 0.121897+0.0336245im  0.983746+2.77556e-17im]
        ρ = simulate_o2(one_spin_model, 1.0, AS_DW_QUADRATIC, 100)

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end


    @testset "function schedule, too few steps, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR, 10)

        # NOTE, atol required due to too few iterations
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-5)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
    end

    @testset "csv schedule pwq, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR_pwq_csv_1000, 100)

        @test isapprox(one_spin_ρ(1.0), ρ)
    end

    @testset "csv schedule pwq, low resolution, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR_pwq_csv_100, 100)

        # NOTE, atol required due to pwq approximation in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-7)
    end

    @testset "csv schedule pwl, low resolution, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR_pwl_csv_100, 100)

        # NOTE, atol required due to pwl approximation in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-4)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-5)
    end

    @testset "csv schedule pwc, low resolution, analytical solution" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR_pwc_csv_100, 100)

        # NOTE, atol required due to pwc approximation in the schedule file
        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-2)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-3)
    end

    @testset "reverse annealing test" begin
        ρ_target = [0.88389+0.0im -0.197484-0.252247im; -0.197484+0.252247im 0.11611+0.0im]

        asch = [(0.0, 1.0), (0.5, 0.5), (1.0, 1.0)]
        mod_asch = dwave_annealing_protocol(AS_LINEAR, asch=asch)

        ρ = simulate_o2(Dict((1,) => 0.1), 5.0, mod_asch, 100, initial_state = [0,1])

        # NOTE, atol required due to too few digits in target
        @test isapprox(ρ_target, ρ, atol=1e-6)
    end

    @testset "collect probability trajectory, nonadaptive" begin
        ρ_list = []
        steps=100
        ρ = simulate_o2(one_spin_model, 1, AS_CIRCULAR, steps, state_steps=ρ_list)

        @test isapprox(one_spin_ρ(1.0), ρ)
        @test length(ρ_list) == steps
        @test isapprox(ρ_list[1], (default_initial_state(1) * default_initial_state(1)'))
        @test isapprox(ρ_list[steps], one_spin_ρ(1.0))
    end

    @testset "collect probability trajectory, adaptive" begin
        ρ_list = []
        ρ = simulate(one_spin_model, 1.0, AS_CIRCULAR, mean_tol=1e-7, max_tol=1e-7, silence=true, state_steps=ρ_list)

        @test isapprox(one_spin_ρ(1.0), ρ)
        @test length(ρ_list) == 64
        @test isapprox(ρ_list[1], (default_initial_state(1) * default_initial_state(1)'))
        @test isapprox(ρ_list[64], one_spin_ρ(1.0))
    end
end


@testset "simulate, 2 qubit" begin
    @testset "function schedule, default anneal time, analytical solution" begin
        ρ = simulate_o2(two_spin_model, 1.0, AS_CIRCULAR, 100)
        @test isapprox(two_spin_ρ(1.0), ρ)
    end

    @testset "function schedule, fast anneal time, analytical solution" begin
        ρ = simulate_o2(two_spin_model, 0.5, AS_CIRCULAR, 100)
        @test isapprox(two_spin_ρ(0.5), ρ)
    end

    @testset "function schedule, slow anneal time, analytical solution" begin
        ρ = simulate_o2(two_spin_model, 2.0, AS_CIRCULAR, 100)
        @test isapprox(two_spin_ρ(2.0), ρ)
    end
end


@testset "simulate, any-order magnus expansion" begin
    @testset "1 qubit, adaptive, orders 1 to 8" begin
        at = 1.0

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 1, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 2, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 3, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 4, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 5, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 6, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 7, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)

        ρ = simulate(one_spin_model, at, AS_CIRCULAR, 8, max_tol=1e-9, silence=true)
        @test isapprox(one_spin_ρ(at), ρ)
    end

    @testset "2 qubit, adaptive, orders 1 to 8" begin
        at = 1.0

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 1, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 2, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 3, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 4, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 5, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 6, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 7, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)

        ρ = simulate(two_spin_model, at, AS_CIRCULAR, 8, max_tol=1e-9, silence=true)
        @test isapprox(two_spin_ρ(at), ρ)
    end

    @testset "5 qubit, hardcoded first order solver" begin
        ising_model = Dict((1,) => -1, (1,2) => -1, (1,3) => -1, (1,4) => -1, (1,5) => -1, (2,3) => 1, (4,5) => 1)

        ρ_target = simulate_o1(ising_model, 2.0, AS_CIRCULAR, 2)
        ρ = simulate(ising_model, 2.0, AS_CIRCULAR, 2, 1)
        @test isapprox(ρ_target, ρ)
    end

    @testset "5 qubit, hardcoded second order solver" begin
        ising_model = Dict((1,) => -1, (1,2) => -1, (1,3) => -1, (1,4) => -1, (1,5) => -1, (2,3) => 1, (4,5) => 1)

        ρ_target = simulate_o2(ising_model, 2.0, AS_CIRCULAR, 2)
        ρ = simulate(ising_model, 2.0, AS_CIRCULAR, 2, 2)
        @test isapprox(ρ_target, ρ)
    end

    @testset "5 qubit, hardcoded second order solver, function schedules" begin
        ising_model = Dict((1,) => -1, (1,2) => -1, (1,3) => -1, (1,4) => -1, (1,5) => -1, (2,3) => 1, (4,5) => 1)

        ρ_target = simulate_o2(ising_model, 2.0, AS_LINEAR, 2)
        ρ = simulate(ising_model, 2.0, AS_LINEAR, 2, 2)
        @test isapprox(ρ_target, ρ)

        ρ_target = simulate_o2(ising_model, 2.0, AS_QUADRATIC, 2)
        ρ = simulate(ising_model, 2.0, AS_QUADRATIC, 2, 2)
        @test isapprox(ρ_target, ρ)

        ρ_target = simulate_o2(ising_model, 2.0, AS_DW_QUADRATIC, 2)
        ρ = simulate(ising_model, 2.0, AS_DW_QUADRATIC, 2, 2)
        @test isapprox(ρ_target, ρ)
    end
end


@testset "simulate, multi-qubit" begin

    @testset "2 qubit, function schedules (AS_CIRCULAR, AS_LINEAR, AS_QUADRATIC), near adiabatic limit" begin
        n = 2
        h = ones(n)
        J = Dict((1,2) => -1)
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 100.0
        steps = 100

        ρ = simulate_o2(ising_model, annealing_time, AS_CIRCULAR, steps)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)

        ρ = simulate_o2(ising_model, annealing_time, AS_LINEAR, steps)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)

        ρ = simulate_o2(ising_model, annealing_time, AS_QUADRATIC, steps)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit, function schedules (AS_DW_QUADRATIC), near adiabatic limit" begin
        n = 2
        h = ones(n)
        J = Dict((1,2) => -1)
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 10.0
        steps = 1000

        ρ = simulate_o2(ising_model, annealing_time, AS_DW_QUADRATIC, steps)
        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit, csv schedule, near adiabatic limit" begin
        n = 2
        h = ones(2)
        J = Dict((1,2) => -1)
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        ρ = simulate_o2(ising_model, 100.0, AS_CIRCULAR_pwl_csv_1000, 100)

        @test real(ρ[4,4]) >= 0.9999

        @test isapprox(min_ρ, ρ, atol=1e-2)
        @test !isapprox(min_ρ, ρ, atol=1e-3)
    end

    @testset "2 qubit, function schedule, fractional values, near adiabatic limit" begin
        n = 2
        h = ones(2)
        J = Dict((1,2) => -0.76)
        ising_model = Dict((1,) => 0.5, (1,2) => -0.75)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        ρ = simulate_o2(ising_model, 50.0, AS_CIRCULAR, 100)

        # NOTE, Non-1.0 due to annealing time not being in the fill adiabatic limit
        @test real(tr(ρ*min_ρ)) >= 0.999

        @test isapprox(min_ρ, ρ, atol=1e-1)
        @test !isapprox(min_ρ, ρ, atol=1e-2)
    end

    @testset "2 qubit probability trajectory, nonadaptive" begin
        n = 2
        h = ones(n)
        J = Dict((1,2) => -1)
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 10.0
        steps = 1000

        ρ_list = []
        ρ = simulate_o2(ising_model, annealing_time, AS_DW_QUADRATIC, steps, state_steps=ρ_list)
        @test real(ρ[4,4]) >= 0.9999
        @test length(ρ_list) == steps
        @test isapprox(ρ_list[1], (default_dwave_initial_state(n) * default_dwave_initial_state(n)'))
        @test real(ρ_list[steps][4,4]) >= 0.9999
        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "2 qubit probability trajectory, adaptive" begin
        n = 2
        h = ones(n)
        J = Dict((1,2) => -1)
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))
        min_vec=evecs[:,1]
        min_ρ = min_vec * min_vec'

        annealing_time = 10.0

        ρ_list = []
        ρ = simulate(ising_model, annealing_time, AS_DW_QUADRATIC, silence=true, state_steps=ρ_list)
        @test real(ρ[4,4]) >= 0.9999
        @test length(ρ_list) == 1024
        @test isapprox(ρ_list[1], (default_dwave_initial_state(n) * default_dwave_initial_state(n)'))
        @test real(ρ_list[1024][4,4]) >= 0.9999
        @test isapprox(min_ρ, ρ, atol=1e-3)
        @test !isapprox(min_ρ, ρ, atol=1e-4)
    end

    @testset "3 qubit, degenerate, function schedules (AS_CIRCULAR, AS_QUADRATIC), near adiabatic limit" begin
        # ring of disagrees => 6 ground states
        n = 3
        h = zeros(n)
        J = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)
        ising_model = J

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))

        annealing_time = 100.0
        steps = 100

        ρ = simulate_o2(ising_model, annealing_time, AS_CIRCULAR, steps)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0, atol=1e-8)
        end

        ρ = simulate_o2(ising_model, annealing_time, AS_QUADRATIC, steps)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0)
        end
    end

    @testset "3 qubit, degenerate, function schedules (AS_DW_QUADRATIC), near adiabatic limit" begin
        # ring of disagrees => 6 ground states
        n = 3
        h = zeros(n)
        J = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)
        ising_model = J

        H = sum_z(n,h) + sum_zizj(n,J)
        evals,evecs = eigen(Matrix(H))

        annealing_time = 100.0
        steps = 100

        ρ = simulate_o2(ising_model, annealing_time, AS_DW_QUADRATIC, steps)
        for i = 1:6
            min_vec = evecs[:,i]
            min_ρ = min_vec * min_vec'
            @test isapprox(real(tr(ρ*min_ρ)), 1.0/6.0, atol=1e-3)
            @test !isapprox(real(tr(ρ*min_ρ)), 1.0/6.0, atol=1e-4)
        end
    end

    @testset "3 qubit, csv schedule, near adiabatic limit" begin
        #ring of disagrees => 6 ground states
        n = 3
        h = zeros(n)
        J = Dict((1,2) => 1, (1,3) => 1, (2,3) => 1)
        ising_model = J

        annealing_time = 100.0
        steps = 100

        ρ_target = simulate_o2(ising_model, annealing_time, AS_CIRCULAR, steps)
        ρ = simulate_o2(ising_model, annealing_time, AS_CIRCULAR_pwq_csv_1000, steps)

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

        ρ = simulate_o2(two_spin_model, 1.0, AS_CIRCULAR, 100)
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

    @testset "2 qubit, function schedule, analytical solution" begin
        ρ = simulate_de(two_spin_model, 1.0, AS_CIRCULAR, 1e-6)
        @test isapprox(two_spin_ρ(1.0), ρ)
    end

    @testset "1 qubit, function schedule, constant terms" begin
        ρ = simulate_o2(one_spin_model, 1.0, AS_CIRCULAR, 100, constant_field_x=[1], constant_field_z=[1])
        ρ_de = simulate_de(one_spin_model, 1.0, AS_CIRCULAR, 1e-6, constant_field_x=[1], constant_field_z=[1])

        @test isapprox(ρ, ρ_de, atol=1e-8)
        @test !isapprox(ρ, ρ_de, atol=1e-9)
    end

    @testset "1 qubit, csv schedule, analytical solution" begin
        ρ = simulate_de(one_spin_model, 1.0, AS_CIRCULAR_pwl_csv_1000, 1e-6)

        @test isapprox(one_spin_ρ(1.0), ρ, atol=1e-6)
        @test !isapprox(one_spin_ρ(1.0), ρ, atol=1e-7)
    end

    @testset "2 qubit, function schedule" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_o2(ising_model, 1.0, AS_CIRCULAR, 100)
        ρ_de = simulate_de(ising_model, 1.0, AS_CIRCULAR, 1e-6)

        @test isapprox(ρ, ρ_de, atol=1e-8)
        @test !isapprox(ρ, ρ_de, atol=1e-9)
    end

    @testset "2 qubit, function schedule, long annealing time" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_o2(ising_model, 1000.0, AS_CIRCULAR, 4000)
        ρ_de = simulate_de(ising_model, 1000.0, AS_CIRCULAR, 1e-7)

        @test isapprox(ρ, ρ_de, atol=1e-6)
        @test !isapprox(ρ, ρ_de, atol=1e-7)
    end

    @testset "2 qubit, function schedule, adaptive tolerance" begin
        ising_model = Dict((1,) => -0.1, (2,) => -1, (1,2) => -1)
        ρ = simulate_o2(ising_model, 100.0, AS_CIRCULAR, 1000)
        ρ_de = simulate_de(ising_model, 100.0, AS_CIRCULAR, silence=true)

        @test isapprox(ρ, ρ_de, atol=1e-6)
        @test !isapprox(ρ, ρ_de, atol=1e-7)
    end
end


@testset "simulate bqpjson" begin

    @testset "1 qubit, function schedule" begin
        annealing_time = 10000.0
        annealing_schedule = AS_CIRCULAR
        steps = 1000

        ising_intended = Dict((1,) => 1)

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, steps)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate_o2(ising_intended, annealing_time, annealing_schedule, steps)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[2], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "1 qubit, dwave schedule" begin
        annealing_time = 10000.0
        annealing_schedule = AS_DW_QUADRATIC
        steps = 1000

        ising_intended = Dict((1,) => 1)

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, steps)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate_o2(ising_intended, annealing_time, annealing_schedule, steps)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[2], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "2 qubit, function schedule" begin
        annealing_time = 10000.0
        annealing_schedule = AS_CIRCULAR
        steps = 1000

        # the ising model that is encoded in bqpjson_2q.json
        ising_intended = Dict((1,) => -1, (2,) => -1, (1,2) => -1)

        bqpjson_file = "data/bqpjson_2q.json"
        dwisc_file = "tmp.json"

        simulate_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, steps)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        ρ = simulate_o2(ising_intended, annealing_time, annealing_schedule, steps)

        @test dwisc_data["solutions"][1]["prob"] >= 0.99
        @test isapprox(z_measure_probabilities(ρ)[1], dwisc_data["solutions"][1]["prob"])
        @test dwisc_data["variable_ids"] == [304, 308]
    end


    @testset "1 qubit, function schedule, z noise" begin
        annealing_time = 100.0
        annealing_schedule = AS_CIRCULAR
        steps = 100
        numshots = 10

        # random numbers between -0.75 and 0.75
        z_bias = (Random.rand(numshots) .- 0.5)*1.5

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_noisy_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, steps, z_bias=z_bias)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        @test dwisc_data["solutions"][1]["prob"] > 0.70 # true w.h.p.
        @test dwisc_data["variable_ids"] == [100]
    end

    @testset "1 qubit, function schedule, x and z noise" begin
        annealing_time = 100.0
        annealing_schedule = AS_CIRCULAR
        steps = 100
        numshots = 10

        # random numbers between -0.75 and 0.75
        x_bias = (Random.rand(numshots) .- 0.5)*1.5
        z_bias = (Random.rand(numshots) .- 0.5)*1.5

        bqpjson_file = "data/bqpjson_1q.json"
        dwisc_file = "tmp.json"

        simulate_noisy_bqpjson(bqpjson_file, dwisc_file, annealing_time, annealing_schedule, steps, x_bias=x_bias, z_bias=z_bias)

        dwisc_data = JSON.parsefile(dwisc_file)
        rm(dwisc_file)

        @test dwisc_data["solutions"][1]["prob"] > 0.70 # true w.h.p.
        @test dwisc_data["variable_ids"] == [100]
    end
end
