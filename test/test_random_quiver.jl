using Test
using ClusterAlgebras
using Random

@testset "random_quiver" begin
    @testset "reproducible from an RNG seed" begin
        # The dataset-manifest convention depends on this: seed + parameters must
        # be enough to regenerate a dataset exactly.
        a = random_quiver(MersenneTwister(7), 5)
        b = random_quiver(MersenneTwister(7), 5)
        @test a == b
        @test random_quiver(MersenneTwister(8), 5).B != a.B
    end

    @testset "shape" begin
        rng = MersenneTwister(1)
        q = random_quiver(rng, 6)
        @test nvertices(q) == 6
        @test q.n_mutable == 6
        @test q.n_frozen == 0

        qf = random_quiver(rng, 4; n_frozen = 2)
        @test nvertices(qf) == 6
        @test qf.n_mutable == 4
        @test qf.n_frozen == 2
        @test all(qf.B[5:6, 5:6] .== 0)   # frozen-frozen block carries no data
    end

    @testset "skew-symmetric, within max_arrows, and mutable" begin
        rng = MersenneTwister(2)
        for _ in 1:50
            q = random_quiver(rng, 5; max_arrows = 2)
            @test q.B == -transpose(q.B)
            @test all(abs.(q.B) .<= 2)
            @test all(q.B[i, i] == 0 for i in 1:5)
            @test mutate(q, 1) isa Quiver   # the constructor's invariants hold
        end
    end

    @testset "density controls arrow count" begin
        rng = MersenneTwister(3)
        @test all(random_quiver(rng, 5; density = 0.0).B .== 0)
        dense = random_quiver(rng, 5; density = 1.0)
        @test all(dense.B[i, j] != 0 for i in 1:5, j in 1:5 if i != j)
    end

    @testset ":acyclic model orients every arrow along the vertex order" begin
        rng = MersenneTwister(4)
        for _ in 1:20
            q = random_quiver(rng, 5; model = :acyclic, density = 0.8)
            # Every arrow runs i → j with i < j, so the vertex order is a
            # topological order and no directed cycle can close.
            @test all(q.B[i, j] >= 0 for i in 1:5, j in 1:5 if i < j)
            @test is_finite_type(q) isa Bool   # acyclic quivers are well-formed
        end
    end

    @testset "rank 0 and rank 1" begin
        @test nvertices(random_quiver(MersenneTwister(5), 0)) == 0
        @test random_quiver(MersenneTwister(5), 1).B == zeros(Int, 1, 1)
    end

    @testset "invalid arguments throw" begin
        rng = MersenneTwister(6)
        @test_throws InvalidArgument random_quiver(rng, -1)
        @test_throws InvalidArgument random_quiver(rng, 3; n_frozen = -1)
        @test_throws InvalidArgument random_quiver(rng, 3; max_arrows = 0)
        @test_throws InvalidArgument random_quiver(rng, 3; density = 1.5)
        @test_throws InvalidArgument random_quiver(rng, 3; model = :nonsense)
    end
end

@testset "random_mutate" begin
    @testset "reproducible from an RNG seed" begin
        a = random_mutate(MersenneTwister(7), Quiver(:A, 4), 12)
        b = random_mutate(MersenneTwister(7), Quiver(:A, 4), 12)
        @test a.B == b.B
    end

    @testset "stays in the mutation class" begin
        # Mutation class is a mutation invariant, so a random walk from a
        # finite-type quiver can only ever produce finite-type quivers.
        rng = MersenneTwister(8)
        for _ in 1:20
            @test is_finite_type(random_mutate(rng, Quiver(:A, 4), 15))
            @test is_finite_type(random_mutate(rng, Quiver(:D, 4), 15))
        end
    end

    @testset "zero steps is the identity" begin
        q = Quiver(:A, 3)
        @test random_mutate(MersenneTwister(9), q, 0).B == q.B
    end

    @testset "walks on seeds too" begin
        s = random_mutate(MersenneTwister(10), Seed(Quiver(:A, 3)), 8)
        @test s isa Seed
        @test is_finite_type(s.quiver)
    end

    @testset "actually moves" begin
        # A long walk on A₄ should leave the initial quiver at least sometimes;
        # a walk that never moved would silently pass every test above.
        rng = MersenneTwister(11)
        q = Quiver(:A, 4)
        @test any(random_mutate(rng, q, 5).B != q.B for _ in 1:20)
    end

    @testset "invalid arguments throw" begin
        rng = MersenneTwister(12)
        @test_throws InvalidArgument random_mutate(rng, Quiver(:A, 3), -1)
        # No mutable vertices: nothing to walk on.
        @test_throws InvalidArgument random_mutate(rng, Quiver(zeros(Int, 2, 2), 0), 3)
    end
end
