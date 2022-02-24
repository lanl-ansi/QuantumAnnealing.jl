
@testset "state encoding and transformation" begin

    @testset "int2binary" begin
        @test int2binary(0) == []
        @test int2binary(1) == [1]
        @test int2binary(2) == [0,1]

        @test int2binary(0, pad=3) == [0, 0, 0]
        @test int2binary(1, pad=3) == [1, 0, 0]
        @test int2binary(2, pad=3) == [0, 1, 0]
        @test int2binary(4, pad=3) == [0, 0, 1]

        @test int2binary(4, pad=4) == [0, 0, 1, 0]
    end

    @testset "int2spin" begin
        @test int2spin(0) == []
        @test int2spin(1) == [-1]
        @test int2spin(2) == [1,-1]

        @test int2spin(0, pad=3) == [ 1,  1,  1]
        @test int2spin(1, pad=3) == [-1,  1,  1]
        @test int2spin(2, pad=3) == [ 1, -1,  1]
        @test int2spin(4, pad=3) == [ 1,  1, -1]

        @test int2spin(4, pad=4) == [1, 1, -1, 1]
    end

    @testset "int2binary/binary2int" begin
        for i in 0:16
            @test binary2int(int2binary(i, pad=10)) == i
        end
    end

    @testset "int2spin/spin2int" begin
        for i in 0:16
            @test spin2int(int2spin(i, pad=10)) == i
        end
    end

    @testset "binary2spin" begin
        @test binary2spin([0, 0, 0]) == [ 1,  1,  1]
        @test binary2spin([1, 0, 0]) == [-1,  1,  1]
        @test binary2spin([0, 1, 0]) == [ 1, -1,  1]
        @test binary2spin([0, 0, 1]) == [ 1,  1, -1]
    end

    @testset "spin2binary" begin
        @test spin2binary([ 1,  1,  1]) == [0, 0, 0]
        @test spin2binary([-1,  1,  1]) == [1, 0, 0]
        @test spin2binary([ 1, -1,  1]) == [0, 1, 0]
        @test spin2binary([ 1,  1, -1]) == [0, 0, 1]
    end

    @testset "binary2braket" begin
        @test binary2braket(int2binary(0, pad=3)) == "|000⟩"
        @test binary2braket(int2binary(1, pad=3)) == "|001⟩"
        @test binary2braket(int2binary(2, pad=3)) == "|010⟩"
        @test binary2braket(int2binary(4, pad=3)) == "|100⟩"
    end

    @testset "binary2braket" begin
        @test spin2braket(int2spin(0, pad=3)) == "|↑↑↑⟩"
        @test spin2braket(int2spin(1, pad=3)) == "|↑↑↓⟩"
        @test spin2braket(int2spin(2, pad=3)) == "|↑↓↑⟩"
        @test spin2braket(int2spin(4, pad=3)) == "|↓↑↑⟩"
    end

end


@testset "ising energy computations" begin

    @testset "single qubit, single state" begin
        @test isapprox(eval_ising_state_energy([1],Dict((1,) => 1)), 1)
        @test isapprox(eval_ising_state_energy([-1],Dict((1,) => 1)), -1)
    end

    @testset "single qubit, all states" begin
        energies = compute_ising_state_energies(Dict((1,) => 1))
        @test isapprox(energies[0], 1)
        @test isapprox(energies[1], -1)
    end

    @testset "two qubit, single state" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)
        @test isapprox(eval_ising_state_energy([1, 1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([1, -1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([-1, 1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([-1, -1], ising_model), -3)
    end

    @testset "two qubit, all state energies" begin
        energies = compute_ising_state_energies(Dict((1,) => 1, (1,2) => -1))
        @test isapprox(energies[0], 0)
        @test isapprox(energies[1], 0)
        @test isapprox(energies[2], 2)
        @test isapprox(energies[3], -2)
    end

    @testset "two qubit, energy levels" begin
        energy_levels = compute_ising_energy_levels(Dict((1,2) => -1))
        @test energy_levels[1].energy == -1.0
        @test energy_levels[1].states == Set([0,3])
        @test energy_levels[2].energy == 1.0
        @test energy_levels[2].states == Set([1,2])
    end

    @testset "two qubit, print energy levels" begin
        mktemp() do path,io
            out = stdout
            err = stderr
            redirect_stdout(io)
            redirect_stderr(io)

            print_ising_energy_levels(Dict((1,2) => -1))
            print_ising_energy_levels(Dict((1,2) => -1), limit=0)
            print_ising_energy_levels(Dict((1,2) => -1), limit=1)
            print_ising_energy_levels(Dict((1,2) => -1), limit=10)

            flush(io)
            redirect_stdout(out)
            redirect_stderr(err)
        end
    end

end


@testset "transverse ising hamiltonian" begin

    @testset "single spin, analytical solution" begin
        @test all(isapprox(single_spin_H(s), transverse_ising_hamiltonian(single_spin_model, AS_CIRCULAR, s)) for s in s_100)
    end

    @testset "two spin, analytical solution" begin
        @test all(isapprox(two_spin_H(s), transverse_ising_hamiltonian(two_spin_model, AS_CIRCULAR, s)) for s in s_100)
    end

    @testset "single spin, linear schedule" begin
        annealing_schedule = AS_LINEAR

        H_00 = transverse_ising_hamiltonian(single_spin_model, annealing_schedule, 0.0)
        H_05 = transverse_ising_hamiltonian(single_spin_model, annealing_schedule, 0.5)
        H_10 = transverse_ising_hamiltonian(single_spin_model, annealing_schedule, 1.0)

        @test isapprox(H_00, [0 1; 1 0])
        @test isapprox(H_05, [0.5 0.5; 0.5 -0.5])
        @test isapprox(H_10, [1 0; 0 -1])
    end

    @testset "two spin, linear schedule" begin
        annealing_schedule = AS_LINEAR

        H_00 = transverse_ising_hamiltonian(two_spin_model, annealing_schedule, 0.0)
        H_05 = transverse_ising_hamiltonian(two_spin_model, annealing_schedule, 0.5)
        H_10 = transverse_ising_hamiltonian(two_spin_model, annealing_schedule, 1.0)

        @test isapprox(H_00, [0.0 1.0 1.0 0.0; 1.0  0.0 0.0 1.0; 1.0 0.0  0.0 1.0; 0.0 1.0 1.0 0.0])
        @test isapprox(H_05, [1.0 0.5 0.5 0.0; 0.5 -1.0 0.0 0.5; 0.5 0.0 -1.0 0.5; 0.0 0.5 0.5 1.0])
        @test isapprox(H_10, [2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0; 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0])
    end

end


@testset "csv annealing schedules" begin

    @testset "piecewise constant" begin
        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwc_csv_100.A(s) for s in s_100]
        @test isapprox(maximum(abs.(deltas)), 0.0, atol=1e-15)

        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwc_csv_100.A(s) for s in s_10000]
        @test isapprox(maximum(abs.(deltas)), 0.015708868493240165, atol=1e-5)
        @test isapprox(sum(abs.(deltas))/length(s_10000), 0.004986646839878622, atol=1e-5)
    end

    @testset "piecewise linear" begin
        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwl_csv_100.A(s) for s in s_100]
        @test isapprox(maximum(abs.(deltas)), 0.0, atol=1e-15)

        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwl_csv_100.A(s) for s in s_10000]
        @test isapprox(maximum(abs.(deltas)), 3.146450816005064e-5, atol=1e-10)
        @test isapprox(sum(abs.(deltas))/length(s_10000), 1.335316012137937e-5, atol=1e-10)
    end

    @testset "piecewise quadratic" begin
        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwq_csv_100.A(s) for s in s_100]
        @test isapprox(maximum(abs.(deltas)), 0.0, atol=1e-15)

        deltas = [AS_CIRCULAR.A(s) - AS_CIRCULAR_pwq_csv_100.A(s) for s in s_10000]
        @test isapprox(maximum(abs.(deltas)), 2.5623974451993714e-7, atol=1e-10)
        @test isapprox(sum(abs.(deltas))/length(s_10000), 1.6639955388651617e-7, atol=1e-10)
    end

    @testset "d-wave asch annealing schedule modification" begin
        annealing_schedule = AS_CIRCULAR

        ss = 0.0:0.001:1.0
        asch = [(0.0,0.0), (1.0,1.0)]
        mod_asch = dwave_annealing_protocol(annealing_schedule, asch=asch)
        @test isapprox(mod_asch.A.(ss), annealing_schedule.A.(ss))
        @test isapprox(mod_asch.B.(ss), annealing_schedule.B.(ss))

        ss_reg = [0.0, 0.25, 0.5, 0.5, 0.5, 0.75, 1.0]
        ss_mod = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
        asch = [(0.0,0.0), (0.4,0.5), (0.6,0.5), (1.0,1.0)]
        mod_asch = dwave_annealing_protocol(annealing_schedule, asch=asch)
        @test isapprox(mod_asch.A.(ss_mod), annealing_schedule.A.(ss_reg))
        @test isapprox(mod_asch.B.(ss_mod), annealing_schedule.B.(ss_reg))
    end

end



