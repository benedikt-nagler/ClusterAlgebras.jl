@testset "denominator_vector" begin
    @testset "initial cluster" begin
        s = Seed(Quiver(:A, 2))
        # Initial variables have d-vector −eᵢ
        @test denominator_vector(s, 1) == [-1, 0]
        @test denominator_vector(s, 2) == [0, -1]
    end

    @testset "A_2 orbit: five cluster variables" begin
        # A_2 cluster algebra has 5 cluster variables with d-vectors:
        # −e₁, −e₂  (initial)  and  e₁, e₂, e₁+e₂  (from mutation)
        q = Quiver(:A, 2)
        s0 = Seed(q)

        s1  = mutate(s0, 1)          # x₁ → (x₂+1)/x₁;  new x₁ has denom x₁ → e₁
        s12 = mutate(s1,  2)         # x₂ → (x₁+x₂+1)/(x₁x₂)           → e₁+e₂
        s121 = mutate(s12, 1)        # x₁ → (x₁+1)/x₂                   → e₂
        # After one more mutation we cycle back; collect all seen dvecs
        all_dvecs = Set([
            denominator_vector(s0,   1),
            denominator_vector(s0,   2),
            denominator_vector(s1,   1),
            denominator_vector(s12,  2),
            denominator_vector(s121, 1),
        ])
        expected = Set([[-1,0],[0,-1],[1,0],[0,1],[1,1]])
        @test all_dvecs == expected
    end

    @testset "first mutation denom" begin
        # Generic: after mutating vertex k, new variable has d-vector eₖ
        for n in 2:4
            q = Quiver(:A, n)
            s = Seed(q)
            s1 = mutate(s, 1)
            @test denominator_vector(s1, 1) == [1, zeros(Int, n-1)...]
        end
    end

    @testset "invalid vertex" begin
        s = Seed(Quiver(:A, 2))
        @test_throws ClusterAlgebraError denominator_vector(s, 0)
        @test_throws ClusterAlgebraError denominator_vector(s, 3)
    end
end
