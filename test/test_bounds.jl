using Test
using ClusterAlgebras

# The Markov quiver: the oriented 3-cycle with double arrows, the standard
# example of a mutation class containing no acyclic seed.
const _MARKOV = [0 2 -2; -2 0 2; 2 -2 0]

# A₃ with principal coefficients realized as frozen vertices, so that the
# exchange polynomials carry their y's.  B̃ = [B; I] padded to a square quiver.
function _principal_quiver(B::Matrix{Int})
    n  = size(B, 1)
    Bt = vcat(B, [i == j ? 1 : 0 for i in 1:n, j in 1:n])
    return Quiver(hcat(Bt, zeros(Int, 2n, n)), n)
end

# Every cluster variable reachable from `s` within `depth` mutations.
function _reachable_variables(s::Seed, depth::Int)
    n        = s.quiver.n_mutable
    vars     = Set(s.cluster[1:n])
    frontier = [s]
    for _ in 1:depth
        next = typeof(s)[]
        for t in frontier, k in 1:n
            t2 = mutate(t, k)
            push!(next, t2)
            push!(vars, t2[k])
        end
        frontier = next
    end
    return vars
end

@testset "bounds" begin

    @testset "is_acyclic" begin
        @test is_acyclic(Quiver(:A, 3))
        @test is_acyclic(Seed(Quiver(:A, 3)))
        @test !is_acyclic(Quiver(_MARKOV))
        # Acyclicity is a seed property, not a class property: mutating A₃ at
        # the middle vertex produces the oriented 3-cycle.
        @test !is_acyclic(mutate(Quiver(:A, 3), 2))
        @test is_acyclic(Quiver([0 2; -2 0]))          # Kronecker
    end

    @testset "initial-seed requirement (ledger 1)" begin
        s = Seed(Quiver(:A, 2))
        @test_throws ClusterAlgebras.InvalidArgument in_upper_bound(s[1], mutate(s, 1))
        @test_throws ClusterAlgebras.InvalidArgument lower_bound_generators(mutate(s, 1))
        @test_throws ClusterAlgebras.InvalidArgument bound_certificate(mutate(s, 1))
        # The documented escape hatch works on the mutated quiver.
        @test bound_certificate(Seed(mutate(s, 1).quiver)).conclusion === :seed_acyclic
    end

    @testset "lower-bound generators" begin
        s = Seed(Quiver(:A, 2))
        g = lower_bound_generators(s)
        @test length(g) == 4
        @test g[1] == s[1] && g[2] == s[2]
        @test g[3] == mutate(s, 1)[1]
        @test g[4] == mutate(s, 2)[2]
        # Frozen vertices are units of ZP, not generators: Gr(2,5) has 2
        # mutable and 5 frozen vertices, so 2·2 = 4 generators.
        @test length(lower_bound_generators(grassmannian(2, 5))) == 4
    end

    @testset "standard monomials: count and shape" begin
        for (q, n) in ((Quiver(:A, 2), 2), (Quiver(:A, 3), 3))
            s = Seed(q)
            for d in 0:3
                sm = standard_monomials(s; max_degree = d)
                # No vertex carries both xⱼ and xⱼ′ (the defining condition).
                @test all(all(t.a .* t.b .== 0) for t in sm)
                @test all(sum(t.a) + sum(t.b) <= d for t in sm)
                # Closed form: choose exponents e ≥ 0 with Σe ≤ d, then a side
                # for each non-zero exponent - counted independently of the
                # ring code by brute-force enumeration of exponent vectors.
                expected = 0
                for ev in Iterators.product(ntuple(_ -> 0:d, n)...)
                    sum(ev) <= d || continue
                    expected += 2^count(!iszero, ev)
                end
                @test length(sm) == expected
            end
        end
        @test_throws ClusterAlgebras.InvalidArgument standard_monomials(
            Seed(Quiver(:A, 2)); max_degree = -1)
    end

    @testset "L ⊆ U on generators and standard monomials" begin
        for q in (Quiver(:A, 2), Quiver(:A, 3), Quiver(:B, 2),
                  Quiver([0 2; -2 0]), Quiver([0 1 1; -1 0 1; -1 -1 0]))
            s = Seed(q)
            @test all(in_upper_bound(g, s) for g in lower_bound_generators(s))
            @test all(in_upper_bound(t.element, s) for t in standard_monomials(s; max_degree = 3))
        end
    end

    @testset "in_upper_bound rejects non-members" begin
        s = Seed(Quiver(:A, 2))
        F = s.ring
        # 1/(1 + x₁) is not even Laurent in the initial cluster.
        @test !is_laurent(one(F) / (one(F) + s[1]), s)
        @test !in_upper_bound(one(F) / (one(F) + s[1]), s)
        # 1/x₁ is Laurent in x but not in the cluster adjacent at 1, where it
        # becomes x₁′/(1 + x₂).
        @test is_laurent(one(F) / s[1], s)
        @test !in_upper_bound(one(F) / s[1], s)
        # Laurent monomials in the *frozen* variables stay in U (ledger 3).
        gr = grassmannian(2, 5)
        @test in_upper_bound(one(gr.ring) / gr[gr.quiver.n_mutable + 1], gr)
    end

    @testset "Laurent phenomenon: every cluster variable lies in U" begin
        for (q, depth) in ((Quiver(:A, 2), 4), (Quiver(:A, 3), 4),
                           (Quiver(:B, 2), 4), (Quiver([0 2; -2 0]), 4),
                           (Quiver(_MARKOV), 3))
            s = Seed(q)
            @test all(in_upper_bound(v, s) for v in _reachable_variables(s, depth))
        end
    end

    @testset "the A₂ identity x₅ = x₁′x₂′ − 1" begin
        s  = Seed(Quiver(:A, 2))
        x5 = mutate(s, [1, 2])[2]
        e  = lower_bound_expansion(x5, s)
        @test e !== nothing
        @test sum(t.coefficient * t.element for t in e) == x5
        # Exactly two terms: the constant −1 and the standard monomial x₁′x₂′.
        @test length(e) == 2
        by = Dict((t.a, t.b) => t.coefficient for t in e)
        @test by[([0, 0], [0, 0])] == -1
        @test by[([0, 0], [1, 1])] == 1
    end

    @testset "every cluster variable expands in L (acyclic types)" begin
        for q in (Quiver(:A, 2), Quiver(:B, 2), Quiver(:A, 3))
            s = Seed(q)
            for v in _reachable_variables(s, 4)
                e = lower_bound_expansion(v, s; max_degree = 3)
                @test e !== nothing
                @test sum(t.coefficient * t.element for t in e) == v
            end
        end
    end

    @testset "lower_bound_expansion: refusals and misses" begin
        s = Seed(Quiver(:A, 2))
        F = s.ring
        # Not Laurent ⇒ certainly not in L.
        @test lower_bound_expansion(one(F) / (one(F) + s[1]), s) === nothing
        # Laurent but outside the degree budget, and outside L at any degree.
        @test lower_bound_expansion(one(F) / s[1], s; max_degree = 3) === nothing
        # Coefficients would live in ZP, not Z.
        @test_throws ClusterAlgebras.InvalidArgument lower_bound_expansion(
            grassmannian(2, 5)[1], grassmannian(2, 5))
    end

    @testset "coprimality and full rank" begin
        # A₂ coefficient-free: P₁ = 1 + x₂, P₂ = x₁ + 1, coprime.
        @test is_coprime(Seed(Quiver(:A, 2)))
        @test has_full_rank(Seed(Quiver(:A, 2)))
        # A₃ coefficient-free: P₁ = 1 + x₂ = P₃, so coprimality FAILS - the
        # coefficient-free acyclic seed is not covered by the theorem, and
        # B̃ = B is skew-symmetric of odd size, so it cannot have full rank.
        @test !is_coprime(Seed(Quiver(:A, 3)))
        @test !has_full_rank(Seed(Quiver(:A, 3)))
        # Principal coefficients repair both: P₁ = y₁ + x₂, P₃ = y₃ + x₂.
        sp = Seed(_principal_quiver([0 1 0; -1 0 1; 0 -1 0]))
        @test has_full_rank(sp)
        @test is_coprime(sp)
        # Gr(2,5): frozen variables make the exchange polynomials coprime.
        @test is_coprime(grassmannian(2, 5))
        @test has_full_rank(grassmannian(2, 5))
    end

    @testset "bound certificate" begin
        # Coprime + acyclic: the full Berenstein-Fomin-Zelevinsky conclusion.
        for q in (Quiver(:A, 2), Quiver(:B, 2), Quiver([0 2; -2 0]))
            c = bound_certificate(Seed(q))
            @test c.conclusion === :seed_acyclic
            @test c.acyclic && c.coprime && c.acyclic_class === true
        end
        @test bound_certificate(_principal_quiver([0 1 0; -1 0 1; 0 -1 0]) |> Seed).conclusion ===
              :seed_acyclic
        @test bound_certificate(grassmannian(2, 5)).conclusion === :seed_acyclic

        # Acyclic but not coprime.
        c3 = bound_certificate(Seed(Quiver(:A, 3)))
        @test c3.conclusion === :not_coprime
        @test c3.acyclic && !c3.coprime

        # Cyclic seed whose class contains an acyclic (and coprime) seed.
        sp_cyclic = Seed(mutate(_principal_quiver([0 1 0; -1 0 1; 0 -1 0]), 2))
        cc = bound_certificate(sp_cyclic)
        @test !cc.acyclic
        @test cc.acyclic_class === true
        @test cc.conclusion === :class_acyclic

        # Markov: no acyclic seed anywhere in the class.
        cm = bound_certificate(Seed(Quiver(_MARKOV)))
        @test cm.conclusion === :not_acyclic
        @test !cm.acyclic && cm.acyclic_class === false
    end

    @testset "upper_bound_generators" begin
        s = Seed(Quiver(:A, 2))
        @test upper_bound_generators(s) == lower_bound_generators(s)
        @test length(upper_bound_generators(grassmannian(2, 5))) == 4
        # Refused where the theorem does not apply, with the reason named.
        @test_throws ClusterAlgebras.InvalidArgument upper_bound_generators(Seed(Quiver(_MARKOV)))
        @test_throws ClusterAlgebras.InvalidArgument upper_bound_generators(Seed(Quiver(:A, 3)))
    end

    @testset "show" begin
        c = bound_certificate(Seed(Quiver(:A, 2)))
        @test c isa BoundCertificate
        @test occursin("seed_acyclic", repr(c))
    end
end
