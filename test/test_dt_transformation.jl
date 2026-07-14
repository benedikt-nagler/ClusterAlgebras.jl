using Test
using ClusterAlgebras
using AbstractAlgebra

@testset "DT transformation" begin

    # Compose the DT map with itself: substitute values into a rational y-image.
    _subst(f, vals) =
        evaluate(numerator(f), vals) // evaluate(denominator(f), vals)

    # ─── A₂ convention pin ───────────────────────────────────────────────────
    # Pins the exact sign/direction conventions (FZ-IV Prop. 3.9 y-mutation,
    # C-columns = c-vectors, σ from C_final = −P_σ) that the DDP layer relies on.

    @testset "A₂ explicit: seq = [2, 1]" begin
        s = extend(Seed(Quiver(:A, 2)))
        y = y_variables(s)
        r = dt_transformation(s; seq = [2, 1])

        @test r.sequence == [2, 1]
        @test r.sigma == [1, 2]
        @test cmatrix(r.seed) == [-1 0; 0 -1]
        # DT: y₁ ↦ 1/(y₁(1+y₂)),  y₂ ↦ (1+y₁+y₁y₂)/y₂
        @test r.y_images == [inv(y[1] * (1 + y[2])),
                             (1 + y[1] + y[1] * y[2]) // y[2]]
    end

    @testset "A₂ second MGS: seq = [1, 2, 1]" begin
        s = extend(Seed(Quiver(:A, 2)))
        r1 = dt_transformation(s; seq = [2, 1])
        r2 = dt_transformation(s; seq = [1, 2, 1])

        @test r2.sigma == [2, 1]
        @test cmatrix(r2.seed) == [0 -1; -1 0]
        @test r2.y_images == r1.y_images       # Keller invariance, explicitly
    end

    # ─── Keller invariance (primary oracle) ──────────────────────────────────
    # The DT map is independent of the chosen maximal green sequence; the final
    # C-matrix is −P_σ and the final quiver is the σ-relabeled initial quiver.

    @testset "Keller invariance: $type" for type in [(:A, 2), (:A, 3), (:B, 2)]
        s = extend(Seed(Quiver(type...)))
        n = s.quiver.n_mutable
        B0 = s.quiver.B[1:n, 1:n]
        sequences = maximal_green_sequences(s)
        @test length(sequences) >= 2

        results = [dt_transformation(s; seq) for seq in sequences]
        for r in results
            σ = r.sigma
            P = zeros(Int, n, n)
            for j in 1:n; P[σ[j], j] = 1; end
            @test cmatrix(r.seed) == -P
            @test all(r.seed.quiver.B[i, j] == B0[σ[i], σ[j]]
                      for i in 1:n, j in 1:n)
            @test is_all_red(r.seed)
            @test r.y_images == results[1].y_images
        end
    end

    # ─── Periodicity (secondary oracle) ──────────────────────────────────────
    # In type A₂ the DT transformation has order 5 (pentagon). Iterate by
    # substituting y_images into itself - NOT by replaying the raw sequence.

    @testset "A₂ DT has order 5" begin
        s = extend(Seed(Quiver(:A, 2)))
        y0 = y_variables(s)
        r = dt_transformation(s)

        current = y0
        returns = Int[]
        for t in 1:6
            current = [_subst(im, current) for im in r.y_images]
            current == y0 && push!(returns, t)
        end
        @test returns == [5]
    end

    # ─── Auto-extend and Quiver method ───────────────────────────────────────

    @testset "auto-extend / Quiver method" begin
        q = Quiver(:A, 2)
        r0 = dt_transformation(extend(Seed(q)); seq = [2, 1])
        r1 = dt_transformation(Seed(q); seq = [2, 1])
        r2 = dt_transformation(q; seq = [2, 1])

        @test r1.sigma == r2.sigma == r0.sigma
        @test string.(r1.y_images) == string.(r2.y_images) == string.(r0.y_images)
    end

    # ─── Errors ──────────────────────────────────────────────────────────────

    @testset "errors" begin
        s = extend(Seed(Quiver(:A, 2)))
        @test_throws InvalidArgument dt_transformation(s; seq = [1, 1])  # not green
        @test_throws InvalidArgument dt_transformation(s; seq = [1])    # not maximal
        @test_throws InvalidVertex dt_transformation(s; seq = [7])

        # The Markov quiver admits no MGS (theorem); bounded search finds none.
        markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
        @test_throws InvalidArgument dt_transformation(markov; max_length = 5)

        # DT is attached to the initial seed: mutated principal seeds rejected.
        @test_throws InvalidArgument dt_transformation(mutate(s, 1))
    end
end
