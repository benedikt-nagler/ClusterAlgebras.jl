@testset "Periodicity" begin

    # ─── is_bipartite ────────────────────────────────────────────────────────

    @testset "is_bipartite" begin
        @test is_bipartite(Quiver(:A, 2))                    # 1 → 2
        @test !is_bipartite(Quiver(:A, 3))                   # linear 1 → 2 → 3
        @test is_bipartite(Quiver([(1, 2), (3, 2)]))         # alternating A₃
        markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
        @test !is_bipartite(markov)                          # oriented 3-cycle
    end

    # ─── mutation_period ─────────────────────────────────────────────────────

    @testset "mutation_period: A₂ pentagon" begin
        s = Seed(Quiver(:A, 2))
        # Pentagon recurrence: the alternating sequence returns after
        # h + 2 = 5 applications of [1, 2] (10 single mutations).
        h = RootSystem(:A, 2).coxeter_number
        @test mutation_period(s, [1, 2]) == h + 2 == 5
    end

    @testset "mutation_period: involutive sequence" begin
        s = Seed(Quiver(:A, 3))
        @test mutation_period(s, [1, 1]) == 1     # μ₁μ₁ = id
    end

    @testset "mutation_period: non-periodic (Markov)" begin
        markov = Seed(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
        @test mutation_period(markov, [1, 2]; max_period=3) === nothing
    end

    @testset "mutation_period: errors" begin
        s = Seed(Quiver(:A, 2))
        @test_throws InvalidVertex mutation_period(s, [1, 7])
        @test_throws InvalidArgument mutation_period(s, Int[])
        @test_throws InvalidArgument mutation_period(s, [1]; max_period=0)
    end

    # ─── y_system: Zamolodchikov periodicity, period h + 2 ───────────────────

    @testset "y_system: A₂ (pentagon, period 5)" begin
        q = Quiver(:A, 2)
        h = RootSystem(:A, 2).coxeter_number
        result = y_system(q)
        @test result.period == h + 2 == 5
        @test length(result.y_trajectory) == result.period + 1
        # Initial entry is the initial y-variables; intermediate Y-seeds are
        # genuinely distinct until the first return.
        @test result.y_trajectory[1] == y_variables(extend(Seed(q)))
        @test allunique(result.y_trajectory[1:result.period])
    end

    @testset "y_system: A₃ bipartite (period 6)" begin
        q = Quiver([(1, 2), (3, 2)])                          # sources 1,3 → sink 2
        h = RootSystem(:A, 3).coxeter_number
        @test y_system(q).period == h + 2 == 6
    end

    @testset "y_system: A₄ bipartite (period 7)" begin
        q = Quiver([(1, 2), (3, 2), (3, 4)])                  # sources 1,3; sinks 2,4
        h = RootSystem(:A, 4).coxeter_number
        @test y_system(q).period == h + 2 == 7
    end

    @testset "y_system: D₄ bipartite (period 8)" begin
        q = Quiver([(1, 4), (2, 4), (3, 4)])                  # leaves → center
        h = RootSystem(:D, 4).coxeter_number
        @test y_system(q).period == h + 2 == 8
    end

    @testset "y_system: principal seed input" begin
        s = extend(Seed(Quiver(:A, 2)))
        @test y_system(s).period == 5
    end

    @testset "y_system: errors" begin
        @test_throws InvalidArgument y_system(Quiver(:A, 3))  # linear, not bipartite
        @test_throws InvalidArgument y_system(Quiver(:A, 2); max_period=0)
    end

    @testset "y_system: cutoff returns nothing" begin
        result = y_system(Quiver([(1, 2), (3, 2)]); max_period=3)
        @test result.period === nothing
        @test length(result.y_trajectory) == 4                # initial + 3 half-steps
    end
end
