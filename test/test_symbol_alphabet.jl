using Test
using ClusterAlgebras

@testset "Symbol alphabets & cluster adjacency" begin

    @testset "symbol_alphabet - oracle sizes" begin
        # finite-type oracle: alphabet size == n_cluster_variables
        for (k, n) in ((2, 5), (2, 6), (3, 6))
            s = grassmannian(k, n)
            alpha = symbol_alphabet(s)
            @test length(alpha) == n_cluster_variables(s.quiver)
        end
    end

    @testset "symbol_alphabet - all elements distinct" begin
        for (k, n) in ((2, 5), (2, 6), (3, 6))
            alpha = symbol_alphabet(grassmannian(k, n))
            @test length(unique(alpha)) == length(alpha)
        end
    end

    @testset "symbol_alphabet - initial variables included" begin
        s = grassmannian(2, 5)
        alpha = symbol_alphabet(s)
        for k in 1:s.quiver.n_mutable
            @test s[k] in alpha
        end
    end

    @testset "cluster_adjacency_matrix - structural properties" begin
        s = grassmannian(2, 5)
        alpha, adj = cluster_adjacency_matrix(s)

        # alphabet matches symbol_alphabet
        @test alpha == symbol_alphabet(s)

        # symmetric
        @test adj == adj'

        # diagonal: every variable is in a cluster with itself
        @test all(adj[i, i] for i in axes(adj, 1))

        # correct size
        @test size(adj) == (length(alpha), length(alpha))
    end

    @testset "cluster_adjacency_matrix - initial cluster mutually adjacent" begin
        s = grassmannian(2, 6)   # A_3, rank 3
        alpha, adj = cluster_adjacency_matrix(s)
        idx = Dict(v => i for (i, v) in enumerate(alpha))
        # all pairs from the initial seed share the initial cluster
        for j in 1:s.quiver.n_mutable, k in 1:s.quiver.n_mutable
            i1, i2 = idx[s[j]], idx[s[k]]
            @test adj[i1, i2]
        end
    end

    @testset "cluster_adjacent - consistent with adjacency matrix" begin
        s = grassmannian(2, 5)
        alpha, adj = cluster_adjacency_matrix(s)
        for i in eachindex(alpha), j in eachindex(alpha)
            @test cluster_adjacent(s, alpha[i], alpha[j]) == adj[i, j]
        end
    end

    @testset "cluster_adjacent - initial variables are adjacent" begin
        s = grassmannian(2, 5)
        @test cluster_adjacent(s, s[1], s[2])
    end

end
