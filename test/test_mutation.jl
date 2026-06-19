@testset "Cluster mutation - A2" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B))
    F  = s.ring
    x1, x2 = s.cluster

    # Mutate at 1: x1' = (1 + x2) / x1
    s1 = mutate(s, 1)
    @test s1.cluster[2] == x2
    @test s1.cluster[1] * x1 == 1 + x2

    # Mutate at 2: x2' = (1 + x1) / x2
    s2 = mutate(s, 2)
    @test s2.cluster[1] == x1
    @test s2.cluster[2] * x2 == 1 + x1

    # Involution: mutating twice in the same direction recovers the original
    @test mutate(s1, 1).cluster[1] == x1
    @test mutate(s2, 2).cluster[2] == x2
end

@testset "Cluster mutation - A2 orbit" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B))
    x1, x2 = s.cluster

    # Second step: mutate s1 at 2 → x2'' = (1 + x1 + x2) / (x1 * x2)
    s12 = mutate(s, [1, 2])
    @test s12.cluster[2] * x1 * x2 == 1 + x1 + x2

    # Mutation sequence produces same result as chained calls
    @test s12.cluster[1] == mutate(mutate(s, 1), 2).cluster[1]
    @test s12.cluster[2] == mutate(mutate(s, 1), 2).cluster[2]
end

@testset "Laurent positivity - A2" begin
    # All cluster variables in the A2 exchange graph must be Laurent polynomials
    # with non-negative integer coefficients when written in the initial cluster.
    # We verify by clearing denominators and checking that all polynomial
    # coefficients are non-negative.
    B = [0 1; -1 0]
    s = Seed(Quiver(B))

    is_positive_laurent(f) = all(>=(0), coefficients(numerator(f)))

    seen = Set{String}()
    queue = [s]
    all_positive = true

    while !isempty(queue)
        cur = popfirst!(queue)
        for k in 1:cur.quiver.n_mutable
            s_new = mutate(cur, k)
            key = string(s_new.cluster[k])
            if key ∉ seen
                push!(seen, key)
                all_positive &= is_positive_laurent(s_new.cluster[k])
                length(seen) < 10 && push!(queue, s_new)  # bound the BFS for tests
            end
        end
    end
    @test all_positive
end

@testset "Label-based Seed mutation" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B, 2, [1, 1], ["p", "q"]))

    sp = mutate(s, "p")
    @test sp.cluster == mutate(s, 1).cluster
    @test sp.mutation_path == [1]

    @test_throws InvalidArgument mutate(s, "z")
end

@testset "Frozen vertices - mutation does not change frozen cluster" begin
    # 2 mutable vertices, 1 frozen
    B = [0 1 1; -1 0 0; -1 0 0]
    q = Quiver(B, 2)
    s = Seed(q)
    x3_frozen = s.cluster[3]

    s1 = mutate(s, 1)
    @test s1.cluster[3] == x3_frozen

    s2 = mutate(s, 2)
    @test s2.cluster[3] == x3_frozen
end
