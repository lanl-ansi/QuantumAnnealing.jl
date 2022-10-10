
@testset "state encoding and transformation" begin

    @testset "int2binary" begin
        @test int_to_binary(0) == []
        @test int_to_binary(1) == [1]
        @test int_to_binary(2) == [0,1]

        @test int_to_binary(0, pad=3) == [0, 0, 0]
        @test int_to_binary(1, pad=3) == [1, 0, 0]
        @test int_to_binary(2, pad=3) == [0, 1, 0]
        @test int_to_binary(4, pad=3) == [0, 0, 1]

        @test int_to_binary(4, pad=4) == [0, 0, 1, 0]
    end

    @testset "int2spin" begin
        @test int_to_spin(0) == []
        @test int_to_spin(1) == [-1]
        @test int_to_spin(2) == [1,-1]

        @test int_to_spin(0, pad=3) == [ 1,  1,  1]
        @test int_to_spin(1, pad=3) == [-1,  1,  1]
        @test int_to_spin(2, pad=3) == [ 1, -1,  1]
        @test int_to_spin(4, pad=3) == [ 1,  1, -1]

        @test int_to_spin(4, pad=4) == [1, 1, -1, 1]
    end

    @testset "int2binary/binary2int" begin
        for i in 0:16
            @test binary_to_int(int_to_binary(i, pad=10)) == i
        end
    end

    @testset "int2spin/spin2int" begin
        for i in 0:16
            @test spin_to_int(int_to_spin(i, pad=10)) == i
        end
    end

    @testset "binary2spin" begin
        @test binary_to_spin([0, 0, 0]) == [ 1,  1,  1]
        @test binary_to_spin([1, 0, 0]) == [-1,  1,  1]
        @test binary_to_spin([0, 1, 0]) == [ 1, -1,  1]
        @test binary_to_spin([0, 0, 1]) == [ 1,  1, -1]
    end

    @testset "spin2binary" begin
        @test spin_to_binary([ 1,  1,  1]) == [0, 0, 0]
        @test spin_to_binary([-1,  1,  1]) == [1, 0, 0]
        @test spin_to_binary([ 1, -1,  1]) == [0, 1, 0]
        @test spin_to_binary([ 1,  1, -1]) == [0, 0, 1]
    end

    @testset "binary2braket" begin
        @test binary_to_braket(int_to_binary(0, pad=3)) == "|000⟩"
        @test binary_to_braket(int_to_binary(1, pad=3)) == "|001⟩"
        @test binary_to_braket(int_to_binary(2, pad=3)) == "|010⟩"
        @test binary_to_braket(int_to_binary(4, pad=3)) == "|100⟩"
    end

    @testset "binary2braket" begin
        @test spin_to_braket(int_to_spin(0, pad=3)) == "|↑↑↑⟩"
        @test spin_to_braket(int_to_spin(1, pad=3)) == "|↑↑↓⟩"
        @test spin_to_braket(int_to_spin(2, pad=3)) == "|↑↓↑⟩"
        @test spin_to_braket(int_to_spin(4, pad=3)) == "|↓↑↑⟩"
    end

end


@testset "ising energy computations" begin

    @testset "1 qubit, single state" begin
        @test isapprox(eval_ising_state_energy([1],Dict((1,) => 1)), 1)
        @test isapprox(eval_ising_state_energy([-1],Dict((1,) => 1)), -1)
    end

    @testset "1 qubit, all states" begin
        energies = compute_ising_state_energies(Dict((1,) => 1))
        @test isapprox(energies[0], 1)
        @test isapprox(energies[1], -1)
    end

    @testset "2 qubit, single state" begin
        ising_model = Dict((1,) => 1, (2,) => 1, (1,2) => -1)
        @test isapprox(eval_ising_state_energy([1, 1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([1, -1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([-1, 1], ising_model), 1)
        @test isapprox(eval_ising_state_energy([-1, -1], ising_model), -3)
    end

    @testset "2 qubit, all state energies" begin
        energies = compute_ising_state_energies(Dict((1,) => 1, (1,2) => -1))
        @test isapprox(energies[0], 0)
        @test isapprox(energies[1], 0)
        @test isapprox(energies[2], 2)
        @test isapprox(energies[3], -2)
    end

    @testset "2 qubit, energy levels" begin
        energy_levels = compute_ising_energy_levels(Dict((1,2) => -1))
        @test energy_levels[1].energy == -1.0
        @test energy_levels[1].states == Set([0,3])
        @test energy_levels[2].energy == 1.0
        @test energy_levels[2].states == Set([1,2])
    end

    @testset "2 qubit, print energy levels" begin
        mktemp() do path,io
            out = stdout
            err = stderr
            redirect_stdout(io)
            redirect_stderr(io)

            print_ising_energy_levels(two_spin_model)
            print_ising_energy_levels(two_spin_model, limit=0)
            print_ising_energy_levels(two_spin_model, limit=1)
            print_ising_energy_levels(two_spin_model, limit=10)
            flush(io)
            redirect_stdout(out)
            redirect_stderr(err)
        end
    end

end


@testset "transverse ising hamiltonian" begin

    @testset "1 qubit, analytical solution" begin
        @test all(isapprox(one_spin_H(s), hamiltonian_transverse_ising(one_spin_model, AS_CIRCULAR, s)) for s in s_100)
    end

    @testset "2 qubit, analytical solution" begin
        @test all(isapprox(two_spin_H(s), hamiltonian_transverse_ising(two_spin_model, AS_CIRCULAR, s)) for s in s_100)
    end

    @testset "1 qubit, linear schedule" begin
        annealing_schedule = AS_LINEAR

        H_00 = hamiltonian_transverse_ising(one_spin_model, annealing_schedule, 0.0)
        H_05 = hamiltonian_transverse_ising(one_spin_model, annealing_schedule, 0.5)
        H_10 = hamiltonian_transverse_ising(one_spin_model, annealing_schedule, 1.0)

        @test isapprox(H_00, [0 1; 1 0])
        @test isapprox(H_05, [0.5 0.5; 0.5 -0.5])
        @test isapprox(H_10, [1 0; 0 -1])
    end

    @testset "2 qubit, linear schedule" begin
        annealing_schedule = AS_LINEAR

        H_00 = hamiltonian_transverse_ising(two_spin_model, annealing_schedule, 0.0)
        H_05 = hamiltonian_transverse_ising(two_spin_model, annealing_schedule, 0.5)
        H_10 = hamiltonian_transverse_ising(two_spin_model, annealing_schedule, 1.0)

        @test isapprox(H_00, [0.0 1.0 1.0 0.0; 1.0  0.0 0.0 1.0; 1.0 0.0  0.0 1.0; 0.0 1.0 1.0 0.0])
        @test isapprox(H_05, [1.0 0.5 0.5 0.0; 0.5 -1.0 0.0 0.5; 0.5 0.0 -1.0 0.5; 0.0 0.5 0.5 1.0])
        @test isapprox(H_10, [2.0 0.0 0.0 0.0; 0.0 -2.0 0.0 0.0; 0.0 0.0 -2.0 0.0; 0.0 0.0 0.0 2.0])
    end

    @testset "boltzmann sampling, 1 qubit" begin
        @test isapprox(z_measure_probabilities(distribution_gibbs(one_spin_H, 0, 1)), [0.5, 0.5])
        @test isapprox(z_measure_probabilities(distribution_gibbs(one_spin_H, 1, 100)), [0, 1])
    end
end

@testset "Analytic Density Matrices" begin
    @testset "One Spin" begin
        one_spin_s_0 = one_spin_ρ(2, s=0)
        one_spin_s_0123 = one_spin_ρ(2, s=0.123)
        one_spin_s_05 = one_spin_ρ(2, s=0.5)
        one_spin_s_1 = one_spin_ρ(2, s=1)

        @test isapprox(sum([one_spin_s_0[i,i] for i in 1:2]), 1)
        @test isapprox(sum([one_spin_s_0123[i,i] for i in 1:2]), 1)
        @test isapprox(sum([one_spin_s_05[i,i] for i in 1:2]), 1)
        @test isapprox(sum([one_spin_s_1[i,i] for i in 1:2]), 1)
    end

    @testset "Two Spin" begin
        two_spin_s_0 = two_spin_ρ(2, s=0)
        two_spin_s_0123 = two_spin_ρ(2, s=0.123)
        two_spin_s_05 = two_spin_ρ(2, s=0.5)
        two_spin_s_1 = two_spin_ρ(2, s=1)

        @test isapprox(sum([two_spin_s_0[i,i] for i in 1:4]), 1)
        @test isapprox(sum([two_spin_s_0123[i,i] for i in 1:4]), 1)
        @test isapprox(sum([two_spin_s_05[i,i] for i in 1:4]), 1)
        @test isapprox(sum([two_spin_s_1[i,i] for i in 1:4]), 1)
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
        mod_asch = annealing_protocol_dwave(annealing_schedule, asch=asch)
        @test isapprox(mod_asch.A.(ss), annealing_schedule.A.(ss))
        @test isapprox(mod_asch.B.(ss), annealing_schedule.B.(ss))

        ss_reg = [0.0, 0.25, 0.5, 0.5, 0.5, 0.75, 1.0]
        ss_mod = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
        asch = [(0.0,0.0), (0.4,0.5), (0.6,0.5), (1.0,1.0)]
        mod_asch = annealing_protocol_dwave(annealing_schedule, asch=asch)
        @test isapprox(mod_asch.A.(ss_mod), annealing_schedule.A.(ss_reg))
        @test isapprox(mod_asch.B.(ss_mod), annealing_schedule.B.(ss_reg))
    end

end


