@testset "magnus expansion method" begin

    @testset "magnus generator subroutines" begin
        mc = QuantumAnnealing._matrix_commutator([1 2; 3 4], [5 6; 7 8])
        @test isapprox([-4 -12; 12 4], mc)

        he = QuantumAnnealing._hamiltonian_eval(2, [[1],[2],[3],[4]])
        @test isapprox([49], he)

        hs = QuantumAnnealing._hamiltonian_scalar(2, [[1],[2],[3],[4]])
        @test isapprox([[2],[4],[6],[8]], hs)

        hi = QuantumAnnealing._hamiltonian_integrate([[1],[2],[3],[4]])
        @test isapprox([[0], [1], [1], [1], [1]], hi)

        hs = QuantumAnnealing._hamiltonian_sum([[1,2,3,4],[5,6]])
        @test isapprox([6,8,3,4], hs)

        hc = QuantumAnnealing._hamiltonian_commutator([[1 2; 3 4],[1 2; 3 4],[1 2; 3 4],[1 2; 3 4]],[[5 6; 7 8],[9 10; 11 12]])
        @test isapprox([[-4 -12; 12 4], [-12 -36; 36 12], [-12 -36; 36 12], [-12 -36; 36 12], [-8 -24; 24 8]], hc)

        bernoulli_fact_numbers = [-0.5, 0.08333333333333333, 0.0, -0.001388888888888889, 0.0, 3.306878306878307e-5]
        for (i,v) in enumerate(bernoulli_fact_numbers)
            @test isapprox(QuantumAnnealing._bernoulli_factorial(i), v)
        end
    end

    @testset "Ω computations, orders 1 to 4" begin
        ising_model = Dict((1,2) => -1, (1,) => 1/2, (2,) => -5/7)
        n = 2
        annealing_time = 2.0
        annealing_schedule = AS_LINEAR
        order = 4

        s0 = 0.0
        s1 = 1.0
        δs = s1 - s0

        x_component = QuantumAnnealing._sum_X(n)
        z_component = zeros(2^n, 2^n)
        for (tup,w) in ising_model
            z_component += QuantumAnnealing._kron_Z(n, tup, w)
        end

        H_parts = QuantumAnnealing._H_parts(x_component, z_component, order)
        Ω_list1 = QuantumAnnealing._Ω_list_optimized(annealing_time, s0, s1, annealing_schedule, H_parts, order)

        aqc = QuantumAnnealing._get_quadratic_coefficients(annealing_schedule.A, s0, s1)
        bqc = QuantumAnnealing._get_quadratic_coefficients(annealing_schedule.B, s0, s1)

        aqc = QuantumAnnealing._shift_quadratic_coefficients(s0, aqc...)
        bqc = QuantumAnnealing._shift_quadratic_coefficients(s0, bqc...)

        H = -im*annealing_time*[
            aqc[1] * x_component + bqc[1] * z_component,
            aqc[2] * x_component + bqc[2] * z_component,
            aqc[3] * x_component + bqc[3] * z_component,
        ]

        Ω_list2 = QuantumAnnealing._Ω_list_generic(δs, H, order)

        for i in 1:order
            @test isapprox(Ω_list1[i], Ω_list2[i])
        end
    end
end
