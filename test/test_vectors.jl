@testset "principal-coefficient seed — construction" begin
    q  = Quiver(:A, 2)
    s  = Seed(q)
    es = extend(s)

    # extend returns a Seed{PrincipalCoefficients}; cluster variables match original
    @test collect(es) == collect(s)
    @test cmatrix(es) == [1 0; 0 1]
    @test gmatrix(es) == [1 0; 0 1]
    @test cvectors(es) == [[1,0], [0,1]]
    @test gvectors(es) == [[1,0], [0,1]]
    @test length(es) == 2
    @test is_sign_coherent(es)
end

@testset "principal-coefficient seed — A₂ mutation oracle" begin
    # B = [[0,1],[-1,0]], all oracle values computed by hand
    q  = Quiver(:A, 2)
    s  = Seed(q)
    es = extend(s)

    # after μ₁
    es1 = mutate(es, 1)
    @test cmatrix(es1) == [-1 1; 0 1]
    @test gmatrix(es1) == [-1 0; 1 1]
    @test is_sign_coherent(es1)

    # after μ₁μ₂
    es12 = mutate(es1, 2)
    @test cmatrix(es12) == [0 -1; 1 -1]
    @test gmatrix(es12) == [-1 -1; 1 0]
    @test is_sign_coherent(es12)

    # after μ₁μ₂μ₁
    es121 = mutate(es12, 1)
    @test cmatrix(es121) == [0 -1; -1 0]
    @test gmatrix(es121) == [0 -1; -1 0]
    @test is_sign_coherent(es121)

    # after μ₁μ₂μ₁μ₂
    es1212 = mutate(es121, 2)
    @test cmatrix(es1212) == [0 1; -1 0]
    @test gmatrix(es1212) == [0 1; -1 0]
    @test is_sign_coherent(es1212)

    # after μ₁μ₂μ₁μ₂μ₁  — A₂ period-5 sequence swaps variables
    es12121 = mutate(es1212, 1)
    @test cmatrix(es12121) == [0 1; 1 0]
    @test gmatrix(es12121) == [0 1; 1 0]
    @test is_sign_coherent(es12121)

    # sequence-form mutate should match step-by-step
    es_seq = mutate(extend(Seed(Quiver(:A, 2))), [1,2,1,2,1])
    @test cmatrix(es_seq) == cmatrix(es12121)
    @test gmatrix(es_seq) == gmatrix(es12121)
end

@testset "principal-coefficient seed — label-based mutation" begin
    # Quiver(:A,2) labels are "1", "2"
    q  = Quiver(:A, 2)
    es = extend(Seed(q))
    @test cmatrix(mutate(es, "1")) == cmatrix(mutate(es, 1))
    @test gmatrix(mutate(es, "1")) == gmatrix(mutate(es, 1))
end

# Deduplication key: sorted d-vector tuple, same as exchange_graph internals.
_es_key(es) =
    Tuple(sort([denominator_vector(es, k) for k in 1:es.quiver.n_mutable]))

@testset "principal-coefficient seed — A₂ sign coherence over full exchange graph" begin
    es    = extend(Seed(Quiver(:A, 2)))
    seen  = Set{Any}()
    stack = [es]
    count = 0
    while !isempty(stack)
        e = pop!(stack)
        key = _es_key(e)
        key ∈ seen && continue
        push!(seen, key)
        count += 1
        @test is_sign_coherent(e)
        for k in 1:e.quiver.n_mutable
            push!(stack, mutate(e, k))
        end
    end
    @test count == 5   # A₂ has 5 distinct seeds
end

@testset "principal-coefficient seed — A₃ sign coherence over full exchange graph" begin
    es    = extend(Seed(Quiver(:A, 3)))
    seen  = Set{Any}()
    stack = [es]
    count = 0
    while !isempty(stack)
        e = pop!(stack)
        key = _es_key(e)
        key ∈ seen && continue
        push!(seen, key)
        count += 1
        @test is_sign_coherent(e)
        for k in 1:e.quiver.n_mutable
            push!(stack, mutate(e, k))
        end
    end
    @test count == 14   # A₃ has 14 distinct seeds
end

@testset "principal-coefficient seed — cluster variables unchanged by extension" begin
    q  = Quiver(:A, 2)
    s  = Seed(q)
    es = mutate(extend(s), [1, 2, 1])
    s2 = mutate(s, [1, 2, 1])
    for k in 1:2
        @test es[k] == s2[k]
    end
end

@testset "is_mutation_finite" begin
    # Finite Dynkin types — small mutation classes, well within default cutoff
    @test is_mutation_finite(Quiver(:A, 2))
    @test is_mutation_finite(Quiver(:A, 4))
    @test is_mutation_finite(Quiver(:D, 4))

    # E_6 mutation class is large (> 10 000 labeled quivers), but is_finite_type
    # short-circuits before the BFS, so the default cutoff is no longer needed.
    @test is_mutation_finite(Quiver(:E, 6))

    # Oriented 3-cycle (mutation-equivalent to A₃) — mutation-finite
    B_cycle = [0 1 -1; -1 0 1; 1 -1 0]
    @test is_mutation_finite(Quiver(B_cycle))

    # Genuinely affine Ã₂ — acyclic orientation of the triangle, mutation-finite
    B_affine = [0 1 1; -1 0 1; -1 -1 0]
    @test is_affine_type(Quiver(B_affine))
    @test is_mutation_finite(Quiver(B_affine))

    # Mutation-infinite: rank-3 quiver with large arrow multiplicities generates
    # > 200 distinct exchange matrices; truncated well before any cutoff ≥ 500.
    B_inf = [0 3 -1; -3 0 2; 1 -2 0]
    @test !is_mutation_finite(Quiver(B_inf); max_quivers = 500)
end
