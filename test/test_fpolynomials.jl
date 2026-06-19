# ─── Helpers ──────────────────────────────────────────────────────────────────

_coeffs(f) = Int[Int(c) for c in coefficients(f)]

# ─── Basic construction ────────────────────────────────────────────────────────

@testset "fpolynomials — initial seed" begin
    es = extend(Seed(Quiver(:A, 2)))
    F  = fpolynomials(es)
    @test length(F) == 2
    @test isone(F[1])
    @test isone(F[2])
end

# ─── A₂ oracles ───────────────────────────────────────────────────────────────
#
# B = [[0,1],[-1,0]] (arrow 1→2).  Non-initial F-polynomials:
#   μ₁   → F₁ = 1+y₁
#   μ₂   → F₂ = 1+y₂
#   μ₁μ₂ → F₂ = 1+y₁+y₁y₂   (the unique "mixed" F-poly for this quiver)
#   μ₂μ₁ → F₁ = 1+y₁+y₁y₂

@testset "fpolynomials — A₂ oracle μ₁" begin
    es = mutate(extend(Seed(Quiver(:A, 2))), 1)
    F  = fpolynomials(es)
    R  = parent(F[1])
    y1, y2 = gens(R)
    @test F[1] == 1 + y1
    @test isone(F[2])
end

@testset "fpolynomials — A₂ oracle μ₂" begin
    es = mutate(extend(Seed(Quiver(:A, 2))), 2)
    F  = fpolynomials(es)
    R  = parent(F[1])
    y1, y2 = gens(R)
    @test isone(F[1])
    @test F[2] == 1 + y2
end

@testset "fpolynomials — A₂ oracle μ₁μ₂ (mixed F-poly at position 2)" begin
    es = mutate(extend(Seed(Quiver(:A, 2))), [1, 2])
    F  = fpolynomials(es)
    R  = parent(F[1])
    y1, y2 = gens(R)
    @test F[1] == 1 + y1
    @test F[2] == 1 + y1 + y1*y2
end

@testset "fpolynomials — A₂ oracle μ₂μ₁ (mixed F-poly at position 1)" begin
    es = mutate(extend(Seed(Quiver(:A, 2))), [2, 1])
    F  = fpolynomials(es)
    R  = parent(F[1])
    y1, y2 = gens(R)
    @test F[1] == 1 + y1 + y1*y2
    @test F[2] == 1 + y2
end

@testset "fpolynomials — mixed F-poly appears at same position for symmetric paths" begin
    # μ₁μ₂ puts the mixed poly at position 2; μ₂μ₁ puts it at position 1.
    F12 = fpolynomials(mutate(extend(Seed(Quiver(:A, 2))), [1, 2]))
    F21 = fpolynomials(mutate(extend(Seed(Quiver(:A, 2))), [2, 1]))
    @test F12[2] == F21[1]   # both = 1+y₁+y₁y₂
end

# ─── Positivity and constant-term checks over A₂ exchange graph ───────────────

@testset "fpolynomials — A₂ positivity and constant-term over full graph" begin
    _key(es) = Tuple(sort([denominator_vector(es, k) for k in 1:2]))
    es0   = extend(Seed(Quiver(:A, 2)))
    seen  = Set{Any}()
    stack = [es0]
    while !isempty(stack)
        e = pop!(stack)
        key = _key(e)
        key ∈ seen && continue
        push!(seen, key)
        for f in fpolynomials(e)
            cs = _coeffs(f)
            @test cs[1] == 1          # constant term = 1
            @test all(>=(0), cs)      # positivity
        end
        for k in 1:2; push!(stack, mutate(e, k)); end
    end
    @test length(seen) == 5
end

# ─── Positivity and constant-term checks over A₃ exchange graph ───────────────

@testset "fpolynomials — A₃ positivity and constant-term over full graph" begin
    _key3(es) = Tuple(sort([denominator_vector(es, k) for k in 1:3]))
    es0   = extend(Seed(Quiver(:A, 3)))
    seen  = Set{Any}()
    stack = [es0]
    while !isempty(stack)
        e = pop!(stack)
        key = _key3(e)
        key ∈ seen && continue
        push!(seen, key)
        for f in fpolynomials(e)
            cs = _coeffs(f)
            @test cs[1] == 1
            @test all(>=(0), cs)
        end
        for k in 1:3; push!(stack, mutate(e, k)); end
    end
    @test length(seen) == 14
end

# ─── Path independence: same seed → same unordered F-poly set ─────────────────
#
# Paths [1,2] and [2,1,2] both reach the seed t₂ of the A₂ exchange graph
# (different orderings of the same cluster), so their F-polynomials are the
# same multiset.

@testset "fpolynomials — path independence (A₂, same seed via two paths)" begin
    es_a = mutate(extend(Seed(Quiver(:A, 2))), [1, 2])
    es_b = mutate(extend(Seed(Quiver(:A, 2))), [2, 1, 2])
    Fa = sort(string.(fpolynomials(es_a)))
    Fb = sort(string.(fpolynomials(es_b)))
    @test Fa == Fb
end
