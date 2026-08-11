using Test
using ClusterAlgebras
using AbstractAlgebra: gen, base_ring

@testset "KS dilog / Ω(γ) export" begin

    # ─── A₂: ordered charges ─────────────────────────────────────────────────
    # A₂ has exactly two maximal green sequences, of lengths 2 and 3.  The short
    # one visits the simple charges; the long one inserts e₁+e₂ in the middle
    # (the BPS chamber with the bound state).
    @testset "A₂ ordered c-vectors" begin
        q   = Quiver(:A, 2)
        mgs = sort(maximal_green_sequences(extend(Seed(q))); by = length)
        @test length(mgs) == 2
        short, long = mgs

        cs = ordered_c_vectors(q, short)
        cl = ordered_c_vectors(q, long)
        @test sort(cs) == [[0, 1], [1, 0]]
        @test sort(cl) == [[0, 1], [1, 0], [1, 1]]
        @test cl[2] == [1, 1]                    # bound state in the middle
        @test cs == reverse([cl[1], cl[3]])      # simples in opposite order

        # Seed and Quiver entry points agree; trivial seeds are extended.
        @test ordered_c_vectors(Seed(q), short) == cs
    end

    # ─── Pentagon identity (pins order + skew-form conventions) ──────────────
    @testset "A₂ pentagon identity" begin
        q   = Quiver(:A, 2)
        mgs = sort(maximal_green_sequences(extend(Seed(q))); by = length)
        w1  = quantum_dilog_word(q, mgs[1])
        w2  = quantum_dilog_word(q, mgs[2])
        p1  = ks_dilog_product(w1; truncation_degree = 6)
        p2  = ks_dilog_product(w2; truncation_degree = 6)
        @test p1 == p2
        @test !isempty(p1)
    end

    # ─── A₃: full wall-crossing invariance ───────────────────────────────────
    @testset "A₃ wall-crossing invariance" begin
        q   = Quiver(:A, 3)
        # 9 = number of maximal chains of the Tamari lattice T₄ (the oriented
        # exchange graph of linearly oriented A₃).  The "8" once quoted on the
        # roadmap was wrong.
        mgs = maximal_green_sequences(extend(Seed(q)))
        @test length(mgs) == 9
        products = [ks_dilog_product(quantum_dilog_word(q, seq);
                                     truncation_degree = 4) for seq in mgs]
        @test all(p == products[1] for p in products)
    end

    # ─── Charges are positive roots; Ω(γ) = 1 ────────────────────────────────
    @testset "charges are positive roots, Ω = 1" begin
        for (type, n) in ((:A, 2), (:A, 3), (:B, 2), (:G, 2))
            q     = Quiver(type, n)
            roots = Set(RootSystem(type, n).positive_roots)
            for seq in maximal_green_sequences(extend(Seed(q)))
                w = quantum_dilog_word(q, seq)
                @test all(γ in roots for γ in w.charges)
                @test all(==(1), values(omega(w)))
                @test sum(values(omega(w))) == length(seq)
            end
        end
    end

    # ─── Word structure ──────────────────────────────────────────────────────
    @testset "QuantumDilogWord fields" begin
        q   = Quiver(:A, 2)
        seq = sort(maximal_green_sequences(extend(Seed(q))); by = length)[1]
        w   = quantum_dilog_word(q, seq)
        @test w.sequence == seq
        @test w.skew == -w.skew'                 # skew form from B
        @test length(w.charges) == length(seq)
        @test occursin("E(ŷ^", sprint(show, w))
    end

    # ─── Error paths ─────────────────────────────────────────────────────────
    @testset "errors" begin
        q   = Quiver(:A, 2)
        mgs = sort(maximal_green_sequences(extend(Seed(q))); by = length)
        # Mutating the same vertex twice: it is red after the first mutation.
        @test_throws InvalidArgument ordered_c_vectors(q, [mgs[1][1], mgs[1][1]])
        # A green but non-maximal prefix.
        @test_throws InvalidArgument quantum_dilog_word(q, mgs[1][1:1])
        # Vertex out of range.
        @test_throws InvalidVertex ordered_c_vectors(q, [5])
        # Truncation degree must be positive.
        w = quantum_dilog_word(q, mgs[1])
        @test_throws InvalidArgument ks_dilog_product(w; truncation_degree = 0)
    end

    # ─── Skew-symmetrizable quivers: Λ = −D·B and the q^{d_k} weighting ───────
    # Λ = Bᵀ would be a skew form only for skew-symmetric B (for B₂ it is
    # [0 -2; 1 0], for G₂ [0 -3; 1 0], neither skew).  −D·B is skew for every
    # skew-symmetrizable B and equals Bᵀ when d ≡ 1, so the skew-symmetric
    # conventions above are untouched.
    @testset "skew form is −D·B, back-compatible on d ≡ 1" begin
        for (letter, rank) in ((:A, 2), (:A, 3), (:B, 2), (:G, 2), (:C, 3))
            q   = Quiver(letter, rank)
            n   = q.n_mutable
            Bm  = q.B[1:n, 1:n]
            seq = first(maximal_green_sequences(extend(Seed(q))))
            w   = quantum_dilog_word(q, seq)
            @test w.skew == -permutedims(w.skew)                  # always skew
            @test w.skew == [-q.d[i] * Bm[i, j] for i in 1:n, j in 1:n]
            @test w.weights == [q.d[k] for k in seq]
            if all(==(1), q.d)
                @test w.skew == permutedims(Bm)                   # the old form
                @test all(==(1), w.weights)
            end
        end
        # The weighting is not decorative: it is non-trivial exactly here, and
        # dropping it breaks the invariance oracles below already for B₂.
        for (letter, rank) in ((:B, 2), (:G, 2))
            q = Quiver(letter, rank)
            @test !all(==(1),
                       quantum_dilog_word(q,
                           first(maximal_green_sequences(extend(Seed(q))))).weights)
        end
    end

    # ─── B₂ hexagon and G₂ octagon identities ────────────────────────────────
    # The non-simply-laced analogues of the A₂ pentagon: the two maximal green
    # sequences have lengths 2 and 4 (B₂) resp. 2 and 6 (G₂), and the weighted
    # products agree.
    @testset "B₂ hexagon / G₂ octagon identities" begin
        for (letter, rank, lengths, trunc) in ((:B, 2, [2, 4], 8),
                                               (:G, 2, [2, 6], 8))
            q   = Quiver(letter, rank)
            mgs = sort(maximal_green_sequences(extend(Seed(q))); by = length)
            @test length.(mgs) == lengths
            products = [ks_dilog_product(quantum_dilog_word(q, seq);
                                         truncation_degree = trunc) for seq in mgs]
            @test all(p == products[1] for p in products)
            @test !isempty(products[1])
        end
    end

    # ─── C₃: full wall-crossing invariance, skew-symmetrizable ───────────────
    @testset "C₃ wall-crossing invariance" begin
        q   = Quiver(:C, 3)
        mgs = maximal_green_sequences(extend(Seed(q)))
        @test length(mgs) == 14
        products = [ks_dilog_product(quantum_dilog_word(q, seq);
                                     truncation_degree = 4) for seq in mgs]
        @test all(p == products[1] for p in products)
    end
    # ─── Shifted / inverted factors: the refined-Ω generalization ────────────
    # A state with a nontrivial refined index Ω(γ, y) contributes shifted
    # dilogarithm factors with integer exponents.  The two conventions are pinned
    # here: a shift is the substitution ŷ^γ → v^s ŷ^γ, and an exponent is a
    # genuine inverse in the truncated quantum torus.
    @testset "shifted and inverted dilogarithm factors" begin
        Λ = [0 2; -2 0]
        γ = [[1, 0]]

        # 𝔼(ŷ^γ)·𝔼(ŷ^γ)^{-1} = 1
        w = QuantumDilogWord([γ[1], γ[1]], Λ, [1, 1], [1, 1], [0, 0], [1, -1])
        prod = ks_dilog_product(w; truncation_degree = 6)
        @test prod == Dict([0, 0] => one(first(values(prod))))

        # a shift is a substitution: the coefficient of ŷ^{nγ} picks up v^{ns}
        for s in (-2, 1, 3)
            plain   = ks_dilog_product(QuantumDilogWord([γ[1]], Λ, [1]);
                                       truncation_degree = 6)
            shifted = ks_dilog_product(
                QuantumDilogWord([γ[1]], Λ, [1], [1], [s], [1]);
                truncation_degree = 6)
            @test keys(shifted) == keys(plain)
            for (α, c) in plain
                n = α[1]                                     # α = n·γ, γ = (1,0)
                @test shifted[α] == c * parent(c)(gen(base_ring(parent(c))))^(n * s)
            end
        end

        # a positive exponent is repeated multiplication
        sq = ks_dilog_product(QuantumDilogWord([γ[1]], Λ, [1], [1], [0], [2]);
                              truncation_degree = 6)
        twice = ks_dilog_product(QuantumDilogWord([γ[1], γ[1]], Λ, [1, 1]);
                                 truncation_degree = 6)
        @test sq == twice
    end

    # ─── The adjoint action and its classical limit ──────────────────────────
    # The product itself has no v → 1 limit (its coefficients have poles there);
    # the conjugation does, and at v = 1 it is exactly the Poisson automorphism
    # X_μ ↦ X_μ (1 + X_γ)^{Λ(γ,μ)} that the DDP jump uses.
    @testset "adjoint action and classical limit" begin
        Λ = [0 2; -2 0]
        w = QuantumDilogWord([[1, 0]], Λ, [1])

        @test_throws ClusterAlgebras.InvalidArgument ks_classical_limit(
            ks_dilog_product(w; truncation_degree = 4))

        ad = ks_classical_limit(ks_dilog_adjoint(w, [0, 1]; truncation_degree = 4))
        # Λ(γ, μ) = (1,0)·Λ·(0,1) = 2, so X_μ ↦ X_μ (1 + X_γ)²
        expected = Dict([0, 1] => 1 // 1, [1, 1] => 2 // 1, [2, 1] => 1 // 1)
        @test ad == Dict{Vector{Int}, Rational{BigInt}}(expected)

        # a cycle is unchanged by its own wall: Λ(γ, γ) = 0
        self = ks_classical_limit(ks_dilog_adjoint(w, [1, 0]; truncation_degree = 4))
        @test self == Dict{Vector{Int}, Rational{BigInt}}([1, 0] => 1 // 1)
    end

end
