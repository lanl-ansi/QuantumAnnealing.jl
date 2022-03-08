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

    @testset "Ω computations, orders 1 to 3" begin
        ising_model = Dict((1,2) => -1, (1,) => 1/2, (2,) => -5/7)
        n = 2
        annealing_time = 2.0
        annealing_schedule = AS_LINEAR
        order = 3

        s0 = 0.0
        s1 = 1.0
        δs = s1 - s0

        x_component = sum_x(n)
        z_component = zeros(2^n, 2^n)
        for (tup,w) in ising_model
            z_component = z_component + sum_z_tup(n, tup, w)
        end


        a_2, a_1, a_0 = get_function_coefficients(annealing_schedule.A, s0, s1)
        b_2, b_1, b_0 = get_function_coefficients(annealing_schedule.B, s0, s1)

        Ω_list1 = QuantumAnnealing._Ω_list(annealing_time, s0, s1, [a_2, a_1, a_0], [b_2, b_1, b_0], x_component, z_component, order)

        a_2_shift = a_2
        a_1_shift = a_1 + 2*a_2*s0
        a_0_shift = a_0 + a_1*s0 + a_2*s0^2

        b_2_shift = b_2
        b_1_shift = b_1 + 2*b_2*s0
        b_0_shift = b_0 + b_1*s0 + b_2*s0^2

        H = -im*annealing_time*[
            a_0_shift * x_component + b_0_shift * z_component,
            a_1_shift * x_component + b_1_shift * z_component,
            a_2_shift * x_component + b_2_shift * z_component,
        ]

        Ω_list_tmp = QuantumAnnealing._magnus_generator(H, order)
        Ω_list2 = [QuantumAnnealing._hamiltonian_eval(δs, Ωi) for Ωi in Ω_list_tmp]

        for i in 1:order
            @test isapprox(Ω_list1[i], Ω_list2[i])
        end
    end
end
