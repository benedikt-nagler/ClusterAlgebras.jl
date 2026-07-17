using Test
using ClusterAlgebras
using Random

@testset "permute_vertices" begin
    q = Quiver(:A, 3)

    @testset "identity" begin
        @test permute_vertices(q, [1, 2, 3]).B == q.B
    end

    @testset "relabeling matches B′[σi, σj] = B[i, j]" begin
        σ = [3, 1, 2]
        p = permute_vertices(q, σ)
        for i in 1:3, j in 1:3
            @test p.B[σ[i], σ[j]] == q.B[i, j]
        end
    end

    @testset "labels travel with their vertices" begin
        p = permute_vertices(q, [2, 3, 1])
        @test p.labels[2] == q.labels[1]
        @test p.labels[3] == q.labels[2]
        @test p.labels[1] == q.labels[3]
    end

    @testset "composition and inverse" begin
        σ = [2, 3, 1]
        τ = invperm(σ)
        @test permute_vertices(permute_vertices(q, σ), τ).B == q.B
    end

    @testset "symmetrizers travel with their vertices" begin
        # B₂: skew-symmetrizable with d = [1, 2].
        b = Quiver([0 2; -1 0], 2, [1, 2])
        p = permute_vertices(b, [2, 1])
        @test p.d == [2, 1]
    end

    @testset "invalid permutations throw" begin
        @test_throws InvalidArgument permute_vertices(q, [1, 2])
        @test_throws InvalidArgument permute_vertices(q, [1, 1, 2])
        # Vertex 3 is frozen here, so no permutation may move it into 1:2.
        qf = Quiver([0 1 1; -1 0 1; -1 -1 0], 2)
        @test_throws InvalidArgument permute_vertices(qf, [1, 3, 2])
        @test permute_vertices(qf, [2, 1, 3]).n_mutable == 2
    end
end

@testset "canonical_form" begin
    @testset "invariant under relabeling" begin
        rng = MersenneTwister(20260717)
        for _ in 1:40
            q = random_quiver(rng, 5)
            σ = randperm(rng, 5)
            @test canonical_form(permute_vertices(q, σ)).B == canonical_form(q).B
        end
    end

    @testset "idempotent" begin
        rng = MersenneTwister(11)
        for _ in 1:20
            q = random_quiver(rng, 4)
            c = canonical_form(q)
            @test canonical_form(c).B == c.B
        end
    end

    @testset "canonical form is isomorphic to its input" begin
        rng = MersenneTwister(12)
        for _ in 1:20
            q = random_quiver(rng, 4)
            @test is_isomorphic(canonical_form(q), q)
        end
    end

    @testset "canonical_permutation realises canonical_form" begin
        rng = MersenneTwister(13)
        for _ in 1:20
            q = random_quiver(rng, 4)
            @test permute_vertices(q, canonical_permutation(q)).B == canonical_form(q).B
        end
    end

    @testset "frozen vertices stay in the frozen block" begin
        qf = Quiver([0 1 2 1; -1 0 1 0; -2 -1 0 1; -1 0 -1 0], 2)
        c = canonical_form(qf)
        @test c.n_mutable == 2
        @test c.n_frozen == 2
        @test is_isomorphic(c, qf)
    end

    @testset "edge cases" begin
        @test canonical_form(Quiver(zeros(Int, 0, 0))).B == zeros(Int, 0, 0)
        @test canonical_form(Quiver(zeros(Int, 1, 1))).B == zeros(Int, 1, 1)
    end
end

@testset "is_isomorphic" begin
    q = Quiver(:A, 3)

    @testset "reflexive and symmetric" begin
        @test is_isomorphic(q, q)
        p = permute_vertices(q, [2, 3, 1])
        @test is_isomorphic(q, p)
        @test is_isomorphic(p, q)
    end

    @testset "ignores labels" begin
        relabeled = Quiver(q.B, q.n_mutable, q.d, ["x", "y", "z"])
        @test is_isomorphic(q, relabeled)
    end

    @testset "separates genuinely different quivers" begin
        # A₃ (path) vs the 3-cycle: same rank and arrow count, different structure.
        cycle = Quiver([0 1 -1; -1 0 1; 1 -1 0])
        @test !is_isomorphic(q, cycle)
        # Different rank, and different arrow multiplicities.
        @test !is_isomorphic(q, Quiver(:A, 2))
        @test !is_isomorphic(Quiver([0 1; -1 0]), Quiver([0 2; -2 0]))
    end

    @testset "respects symmetrizers" begin
        # Same exchange matrix shape, different symmetrizer data.
        b2 = Quiver([0 2; -1 0], 2, [1, 2])
        @test is_isomorphic(b2, b2)
        @test !is_isomorphic(b2, Quiver([0 2; -2 0]))
    end

    @testset "isomorphism is finer than mutation equivalence" begin
        # A mutation of A₃ lands in the same mutation class but need not be
        # isomorphic to it.
        m = mutate(q, 2)
        @test is_finite_type(m)
        @test !is_isomorphic(m, q)
    end

    @testset "distinguishes orientations of the same underlying graph" begin
        # Linear A₃ orientations: 1→2→3 vs 1→2←3 are not isomorphic as quivers.
        path  = Quiver([0 1 0; -1 0 1; 0 -1 0])
        sink  = Quiver([0 1 0; -1 0 -1; 0 1 0])
        @test !is_isomorphic(path, sink)
    end

    @testset "agrees with brute force over all relabelings" begin
        # The definition, tested directly: at rank 4 all 24 relabelings are cheap
        # to try, so the canonical form has to agree with exhaustive search.
        all_perms = [collect(p) for p in Iterators.product(fill(1:4, 4)...) if isperm(collect(p))]
        rng = MersenneTwister(99)
        for _ in 1:30
            q1 = random_quiver(rng, 4)
            q2 = random_quiver(rng, 4)
            brute = any(σ -> permute_vertices(q1, σ).B == q2.B, all_perms)
            @test is_isomorphic(q1, q2) == brute
        end
    end
end
