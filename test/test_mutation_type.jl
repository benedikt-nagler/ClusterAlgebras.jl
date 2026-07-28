using Test
using ClusterAlgebras
using Random

# Star with `k` leaves attached to vertex 1, all arrows pointing outwards.
function star_quiver(k::Int)
    B = zeros(Int, k + 1, k + 1)
    for j in 2:(k + 1)
        B[1, j], B[j, 1] = 1, -1
    end
    return Quiver(B)
end

# Tree on `n` vertices from an edge list, all arrows oriented i -> j.
function tree_quiver(n::Int, edges)
    B = zeros(Int, n, n)
    for (i, j) in edges
        B[i, j], B[j, i] = 1, -1
    end
    return Quiver(B)
end

@testset "MutationType" begin

    @testset "printing" begin
        @test string(MutationType(:A, 3)) == "A3"
        @test string(MutationType(:D, 4, 1)) == "D4^(1)"
        @test string(MutationType(:A, (1, 2), 1)) == "A(1,2)^(1)"
        @test string(MutationType(:E, 6, (1, 1))) == "E6^(1,1)"
        @test string(MutationType(:X, 7)) == "X7"
        @test string(MutationType(:R2, (1, 5))) == "R2(1,5)"
    end

    @testset "equality and hashing" begin
        @test MutationType(:A, 3) == MutationType(:A, 3, nothing)
        @test MutationType(:A, 3) != MutationType(:A, 3, 1)
        @test hash(MutationType(:E, 6, (1, 1))) == hash(MutationType(:E, 6, (1, 1)))
        @test length(Set([MutationType(:A, 3), MutationType(:A, 3)])) == 1
    end
end

@testset "mutation_type: finite" begin

    @testset "agrees with cartan_type on the Dynkin quivers" begin
        for (letter, rank) in [(:A, 2), (:A, 5), (:D, 4), (:D, 6), (:E, 6), (:E, 7),
                               (:E, 8), (:B, 3), (:C, 4), (:F, 4), (:G, 2)]
            q = Quiver(letter, rank)
            @test mutation_type(q) == MutationType(cartan_type(q)...)
        end
    end

    @testset "is a mutation invariant" begin
        rng = MersenneTwister(20260728)
        for (letter, rank) in [(:A, 4), (:D, 5), (:E, 6), (:B, 3), (:G, 2)]
            q = Quiver(letter, rank)
            expected = mutation_type(q)
            p = q
            for _ in 1:12
                p = mutate(p, rand(rng, 1:rank))
                @test mutation_type(p) == expected
            end
        end
    end

    @testset "rank 1 and rank 2" begin
        @test mutation_type(Quiver(zeros(Int, 1, 1))) == MutationType(:A, 1)
        @test mutation_type(Quiver([0 1; -1 0])) == MutationType(:A, 2)
        @test mutation_type(Quiver([0 1; -2 0], 2, [2, 1])) == MutationType(:B, 2)
        @test mutation_type(Quiver([0 1; -3 0], 2, [3, 1])) == MutationType(:G, 2)
        @test mutation_type(Quiver([0 2; -2 0])) == MutationType(:A, (1, 1), 1)
        @test mutation_type(Quiver([0 1; -4 0], 2, [4, 1])) == MutationType(:BC, 1, 1)
        @test mutation_type(Quiver([0 5; -5 0])) == MutationType(:R2, (5, 5))
        @test mutation_type(Quiver([0 1; -5 0], 2, [5, 1])) == MutationType(:R2, (1, 5))
    end

    @testset "frozen vertices are ignored" begin
        B = [0 1 0 1; -1 0 1 0; 0 -1 0 0; -1 0 0 0]
        @test mutation_type(Quiver(B, 2)) == MutationType(:A, 2)
    end
end

@testset "mutation_type: affine" begin

    @testset "cycles are Ã(p, q)" begin
        # 1 -> 2 -> 3, 1 -> 3: two arrows with the traversal, one against.
        @test mutation_type(tree_quiver(3, [(1, 2), (2, 3), (1, 3)])) ==
              MutationType(:A, (1, 2), 1)
        # A cycle on four vertices with two arrows each way.
        @test mutation_type(tree_quiver(4, [(1, 2), (3, 2), (3, 4), (1, 4)])) ==
              MutationType(:A, (2, 2), 1)
        @test is_affine_type(tree_quiver(4, [(1, 2), (3, 2), (3, 4), (1, 4)]))
    end

    @testset "extended Dynkin trees" begin
        @test mutation_type(star_quiver(4)) == MutationType(:D, 4, 1)
        # D5^(1): two leaves at each end of the edge 3-4.
        d5 = tree_quiver(6, [(1, 3), (2, 3), (3, 4), (4, 5), (4, 6)])
        @test mutation_type(d5) == MutationType(:D, 5, 1)
        d6 = tree_quiver(7, [(1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)])
        @test mutation_type(d6) == MutationType(:D, 6, 1)
        # E6^(1) = T(3,3,3), E7^(1) = T(2,4,4), E8^(1) = T(2,3,6).
        e6 = tree_quiver(7, [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 7)])
        @test mutation_type(e6) == MutationType(:E, 6, 1)
        e7 = tree_quiver(8, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (4, 8)])
        @test mutation_type(e7) == MutationType(:E, 7, 1)
        e8 = tree_quiver(9, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (3, 9)])
        @test mutation_type(e8) == MutationType(:E, 8, 1)
    end

    @testset "is a mutation invariant" begin
        rng = MersenneTwister(7)
        for q in (Quiver([0 2; -2 0]), star_quiver(4),
                  tree_quiver(7, [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (6, 7)]))
            expected = mutation_type(q)
            @test expected !== nothing
            p = q
            for _ in 1:8
                p = mutate(p, rand(rng, 1:nvertices(q)))
                @test mutation_type(p) == expected
            end
        end
    end
end

@testset "mutation_type: exceptional" begin

    @testset "the stored representatives are the classes SageMath reports" begin
        # Class sizes up to isomorphism, from SageMath's QuiverMutationType.
        # E8^(1,1) (rank 10) is left out of the suite: its class is 5739 quivers
        # and takes ~10s to enumerate, which was checked once by hand.
        for (n, size) in [(6, 5), (7, 2), (8, 49), (9, 506)]
            rep = ClusterAlgebras._exceptional_quiver(n)
            mc = mutation_class(rep; up_to_isomorphism = true, max_quivers = 1000)
            @test !is_truncated(mc)
            @test length(mc) == size
        end
    end

    @testset "elliptic E7" begin
        e7 = ClusterAlgebras._exceptional_quiver(9)
        @test mutation_type(e7) == MutationType(:E, 7, (1, 1))
        @test mutation_type(mutate(mutate(e7, 3), 6)) == MutationType(:E, 7, (1, 1))
    end

    @testset "X6 and X7" begin
        x6 = ClusterAlgebras._exceptional_quiver(6)
        x7 = ClusterAlgebras._exceptional_quiver(7)
        @test mutation_type(x6) == MutationType(:X, 6)
        @test mutation_type(x7) == MutationType(:X, 7)
        @test is_mutation_finite(x6)
        @test !is_finite_type(x6)
        for k in 1:6
            @test mutation_type(mutate(x6, k)) == MutationType(:X, 6)
        end
        for k in 1:7
            @test mutation_type(mutate(x7, k)) == MutationType(:X, 7)
        end
    end

    @testset "elliptic E6" begin
        e6 = ClusterAlgebras._exceptional_quiver(8)
        @test mutation_type(e6) == MutationType(:E, 6, (1, 1))
        rng = MersenneTwister(11)
        p = e6
        for _ in 1:6
            p = mutate(p, rand(rng, 1:8))
            @test mutation_type(p) == MutationType(:E, 6, (1, 1))
        end
    end

    @testset "a surface-type quiver is left unnamed, not mistaken for one" begin
        # Markov (the once-punctured torus) is mutation-finite with no acyclic
        # representative, exactly like the exceptionals, and must not match one.
        markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
        @test mutation_type(markov) === nothing
        @test is_mutation_finite(markov)
    end
end

@testset "mutation_type: unidentified" begin

    @testset "mutation-infinite quivers give nothing" begin
        wild = Quiver([0 3 -3; -3 0 3; 3 -3 0])
        @test mutation_type(wild) === nothing
        @test mutation_type(Quiver([0 1 -1; -1 0 3; 1 -3 0])) === nothing
    end

    @testset "errors" begin
        @test_throws InvalidArgument mutation_type(Quiver(zeros(Int, 0, 0)))
        disconnected = Quiver([0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0])
        @test_throws InvalidArgument mutation_type(disconnected)
    end
end

@testset "mutation_types" begin
    disconnected = Quiver([0 1 0 0; -1 0 0 0; 0 0 0 1; 0 0 -1 0])
    @test mutation_types(disconnected) == [MutationType(:A, 2), MutationType(:A, 2)]

    mixed = Quiver([0 1 0 0 0; -1 0 0 0 0; 0 0 0 1 0; 0 0 -1 0 1; 0 0 0 -1 0])
    @test mutation_types(mixed) == [MutationType(:A, 2), MutationType(:A, 3)]

    @test mutation_types(Quiver(zeros(Int, 0, 0))) == []
end
