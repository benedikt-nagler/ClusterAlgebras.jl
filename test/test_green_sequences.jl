@testset "Green sequences" begin

    # ─── Predicates ──────────────────────────────────────────────────────────

    @testset "is_green / is_red / is_all_red" begin
        s = extend(Seed(Quiver(:A, 2)))

        # Initial seed: C = I₂, both vertices green.
        @test is_green(s, 1)
        @test is_green(s, 2)
        @test !is_red(s, 1)
        @test !is_all_red(s)

        # After μ₁: C = [-1 1; 0 1] — c₁ = (−1,0) red, c₂ = (1,1) green.
        s1 = mutate(s, 1)
        @test is_red(s1, 1)
        @test is_green(s1, 2)
        @test !is_all_red(s1)

        # After μ₁,μ₂: C = [0 -1; 1 -1] — c₁ = (0,1) green, c₂ = (−1,−1) red.
        s12 = mutate(s1, 2)
        @test is_green(s12, 1)
        @test is_red(s12, 2)
        @test !is_all_red(s12)

        # After μ₁,μ₂,μ₁: C = [0 -1; -1 0] — both red.
        s121 = mutate(s12, 1)
        @test is_red(s121, 1)
        @test is_red(s121, 2)
        @test is_all_red(s121)
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

        # Both known MGS must be present.  For B = [0 1; -1 0] (arrow 1 → 2)
        # the sink-first sequence [2,1] has length 2 and the source-first
        # sequence [1,2,1] has length 3.
        @test [2, 1]    ∈ seqs
        @test [1, 2, 1] ∈ seqs
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

    # ─── Sign tracking (green_sequence_signs) ──────────────────────────────────

    @testset "green_sequence_signs" begin
        s = Seed(Quiver(:A, 2))
        seq = first(maximal_green_sequences(s))    # an actual A₂ MGS

        signs = green_sequence_signs(s, seq)
        @test size(signs) == (length(seq) + 1, 2)  # (length+1) × n_mutable
        @test all(==(:green), signs[1, :])         # start: all green
        @test all(==(:red),   signs[end, :])       # end: all red (MGS theorem)
        # sign-coherent: no seed shows a :mixed c-vector
        @test !any(==(:mixed), signs)

        # matches the step-by-step per-vertex predicates
        sp = extend(s)
        @test signs[1, 1] == (is_green(sp, 1) ? :green : :red)

        # accepts an already-extended seed too, same answer
        @test green_sequence_signs(extend(s), seq) == signs

        # a non-maximal prefix need not end all-red — [1] on A₂ leaves v2 green
        partial = green_sequence_signs(s, [1])
        @test size(partial) == (2, 2)
        @test partial[2, 2] == :green

        # every MGS of A₃ starts all-green, ends all-red
        for seq3 in maximal_green_sequences(Seed(Quiver(:A, 3)))
            m = green_sequence_signs(Seed(Quiver(:A, 3)), seq3)
            @test all(==(:green), m[1, :])
            @test all(==(:red), m[end, :])
        end
    end

end
