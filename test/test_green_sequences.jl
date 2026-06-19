@testset "Green sequences" begin

    # ─── Predicates ──────────────────────────────────────────────────────────

    @testset "is_green / is_red / is_all_red" begin
        s = extend(Seed(Quiver(:A, 2)))

        # Initial seed: C = I₂, both vertices green.
        @test is_green(s, 1)
        @test is_green(s, 2)
        @test !is_red(s, 1)
        @test !is_all_red(s)

        # After μ₁: c₁ = (−1,0) red, c₂ = (0,1) green.
        s1 = mutate(s, 1)
        @test is_red(s1, 1)
        @test is_green(s1, 2)
        @test !is_all_red(s1)

        # After μ₁,μ₂: both red.
        s12 = mutate(s1, 2)
        @test is_red(s12, 1)
        @test is_red(s12, 2)
        @test is_all_red(s12)
    end

    # ─── A1: single MGS [1] ──────────────────────────────────────────────────

    @testset "A1" begin
        s = Seed(Quiver(:A, 1))
        seqs = maximal_green_sequences(s)
        @test length(seqs) == 1
        @test seqs[1] == [1]
    end

    # ─── A2: exactly two MGS ─────────────────────────────────────────────────

    @testset "A2" begin
        s = Seed(Quiver(:A, 2))
        seqs = maximal_green_sequences(s)

        # Verify count.
        @test length(seqs) == 2

        # Both known MGS must be present.
        @test [1, 2]    ∈ seqs
        @test [2, 1, 2] ∈ seqs
    end

    # ─── Structural checks (A3) ───────────────────────────────────────────────

    @testset "A3 structural" begin
        s = Seed(Quiver(:A, 3))
        seqs = maximal_green_sequences(s)

        # There must be at least one MGS.
        @test !isempty(seqs)

        for seq in seqs
            # Replay the sequence and check every step mutates a green vertex.
            sp = extend(s)
            for k in seq
                @test is_green(sp, k)
                sp = mutate(sp, k)
            end
            # Terminal seed must be all-red.
            @test is_all_red(sp)
        end
    end

    # ─── Auto-extend of TrivialCoefficients seed ─────────────────────────────

    @testset "auto-extend" begin
        s_trivial  = Seed(Quiver(:A, 2))
        s_extended = extend(s_trivial)
        @test maximal_green_sequences(s_trivial)  == maximal_green_sequences(s_extended)
    end

    # ─── D4: finite type, all sequences valid ─────────────────────────────────

    @testset "D4 structural" begin
        s = Seed(Quiver(:D, 4))
        seqs = maximal_green_sequences(s)

        @test !isempty(seqs)
        for seq in seqs
            sp = extend(s)
            for k in seq
                @test is_green(sp, k)
                sp = mutate(sp, k)
            end
            @test is_all_red(sp)
        end
    end

end
