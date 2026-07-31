using Test
using ClusterAlgebras

# The rank-2 quiver B = [0 b; −c 0] with its symmetrizer d = [c, b] / gcd(b, c).
function _rank2(b::Int, c::Int)
    g = gcd(b, c)
    return Quiver([0 b; -c 0], 2, [div(c, g), div(b, g)])
end

# (b, c) pairs: A₂, B₂, G₂ (finite, bc ≤ 3), Kronecker, and two wild ones.
const _BC = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))

# Cluster variables of a rank-2 seed along the alternating walk, paired with their
# denominator vectors; the two initial variables (d = −eᵢ) are dropped.
function _walk_variables(b::Int, c::Int, steps::Int)
    s0  = Seed(_rank2(b, c))
    out = Tuple{Vector{Int}, eltype(s0.cluster)}[]
    seen = Set{eltype(s0.cluster)}()
    s = s0
    for t in 1:steps
        s = mutate(s, isodd(t) ? 1 : 2)
        for k in 1:2
            v = s.cluster[k]
            v in seen && continue
            push!(seen, v)
            dv = denominator_vector(s, k)
            all(>=(0), dv) && push!(out, (dv, v))
        end
    end
    return s0, out
end

@testset "Greedy elements (rank 2)" begin

    # ─── The load-bearing oracle ─────────────────────────────────────────────
    # A cluster variable *is* the greedy element at its own denominator vector
    # (LLZ: the greedy basis contains all cluster monomials).  Exact equality in
    # Frac(ℤ[x₁,x₂]); this is what pins both conventions in the src ledger.
    @testset "cluster variables are greedy elements" begin
        for (b, c) in _BC
            s0, vars = _walk_variables(b, c, 6)
            @test !isempty(vars)
            for (dv, v) in vars
                @test greedy_element(s0, dv) == v
            end
        end
    end

    # ─── Degenerate exponents ────────────────────────────────────────────────
    @testset "initial variables and the unit" begin
        for (b, c) in _BC
            s0 = Seed(_rank2(b, c))
            @test greedy_element(s0, [-1, 0]) == s0.cluster[1]
            @test greedy_element(s0, [0, -1]) == s0.cluster[2]
            @test greedy_element(s0, [0, 0]) == one(s0.ring)
        end
    end

    # ─── Positivity (LLZ) ────────────────────────────────────────────────────
    # Regression guard on ledger point 1: with the *generalized* binomial the
    # Kronecker (1,1) coefficient is −1 and this testset fails.
    @testset "coefficients are nonnegative" begin
        for (b, c) in _BC, a1 in 0:3, a2 in 0:3
            @test all(>=(0), greedy_coefficients([a1, a2], b, c))
        end
    end

    # ─── The first imaginary element ─────────────────────────────────────────
    # Beyond finite type the greedy basis is strictly larger than the cluster
    # monomials; for the Kronecker quiver the extra element at (1,1) is
    # (1 + x₁² + x₂²)/(x₁x₂).
    @testset "Kronecker imaginary element" begin
        s0, vars = _walk_variables(2, 2, 8)
        x1, x2 = s0.cluster[1], s0.cluster[2]
        z = greedy_element(s0, [1, 1])
        @test z == (1 + x1^2 + x2^2) * inv(x1 * x2)
        @test all(v != z for (_, v) in vars)          # not a cluster variable
    end

    # ─── Cluster monomials ───────────────────────────────────────────────────
    @testset "cluster monomials sit in the basis" begin
        for (b, c) in _BC
            s0 = Seed(_rank2(b, c))
            x1, x2 = s0.cluster[1], s0.cluster[2]
            for m in 0:2, n in 0:2
                @test greedy_element(s0, [-m, -n]) == x1^m * x2^n
            end
        end
        # A finite type: an exchangeable pair's product is the greedy element at
        # the sum of the two denominator vectors.
        s0 = Seed(_rank2(1, 1))                        # A₂
        s1 = mutate(s0, 1)
        d1 = denominator_vector(s1, 1)
        s2 = mutate(s1, 2)
        d2 = denominator_vector(s2, 2)
        @test greedy_element(s0, d1 + d2) == s1.cluster[1] * s2.cluster[2]
    end

    # ─── Orientation ─────────────────────────────────────────────────────────
    # B = [0 −b; c 0] is the same algebra with the two vertices swapped, which
    # `greedy_element` absorbs internally.  Checked end to end by running the
    # cluster-variable oracle on the reversed quivers.
    @testset "reversed orientation" begin
        for (b, c) in ((1, 2), (2, 2), (2, 3))
            g   = gcd(b, c)
            rev = Quiver([0 -b; c 0], 2, [div(c, g), div(b, g)])
            s0  = Seed(rev)
            s   = s0
            seen = Set{eltype(s0.cluster)}()
            n = 0
            for t in 1:6
                s = mutate(s, isodd(t) ? 2 : 1)
                for k in 1:2
                    v = s.cluster[k]
                    v in seen && continue
                    push!(seen, v)
                    dv = denominator_vector(s, k)
                    all(>=(0), dv) || continue
                    n += 1
                    @test greedy_element(s0, dv) == v
                end
            end
            @test n > 0
        end
    end

    # ─── Coefficient tables: shape and pinned values ─────────────────────────
    @testset "coefficient table shape and pinned values" begin
        for (b, c) in _BC, a1 in -1:3, a2 in -1:3
            d = greedy_coefficients([a1, a2], b, c)
            @test size(d) == (max(a2, 0) + 1, max(a1, 0) + 1)
            @test d[1, 1] == 1
        end
        # Kronecker at (1,1): 1 + x₁² + x₂², i.e. d(0,0) = d(1,0) = d(0,1) = 1,
        # d(1,1) = 0 (ledger point 1).
        @test greedy_coefficients([1, 1], 2, 2) == BigInt[1 1; 1 0]
        # (b,c) = (3,3) at (2,2): 1 + 2x₁³ + x₁⁶ + 2x₂³ + x₂⁶ over x₁²x₂².
        @test greedy_coefficients([2, 2], 3, 3) == BigInt[1 2 1; 2 0 0; 1 0 0]
    end

    # ─── Error paths ─────────────────────────────────────────────────────────
    @testset "errors" begin
        @test_throws InvalidArgument greedy_element(Seed(Quiver(:A, 3)), [1, 0])
        @test_throws InvalidArgument greedy_element(Quiver(:A, 2), [1, 0, 0])
        @test_throws InvalidArgument greedy_coefficients([1, 0], 0, 1)
        # Frozen vertices: the greedy basis is coefficient-free.
        qf = Quiver([0 1 1; -1 0 0; -1 0 0], 2)
        @test_throws InvalidArgument greedy_element(Seed(qf), [1, 0])
    end
end
