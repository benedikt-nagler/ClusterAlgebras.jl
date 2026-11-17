@testset "cluster fan & cluster complex" begin

    @testset "counts against the enumerative closed forms" begin
        for (t, r) in [(:A, 2), (:A, 3), (:B, 2), (:B, 3), (:C, 3), (:D, 4), (:G, 2)]
            q   = Quiver(t, r)
            fan = g_vector_fan(q)
            @test !is_truncated(fan)
            # One ray per cluster variable, one maximal cone per cluster.
            @test length(rays(fan))          == n_cluster_variables(q)
            @test length(maximal_cones(fan)) == n_clusters(q)
            # Distinct cluster variables have distinct g-vectors.
            @test length(unique(rays(fan)))  == length(rays(fan))
            @test all(!iszero, rays(fan))
            # Every cone is spanned by exactly `rank` rays.
            @test all(c -> length(c) == r, maximal_cones(fan))
            # Finite type: the fan is complete.
            @test is_complete(fan)
        end
    end

    @testset "reducible type: A1 x A1" begin
        q   = Quiver([0 0; 0 0])
        fan = g_vector_fan(q)
        @test length(rays(fan))          == n_cluster_variables(q) == 4
        @test length(maximal_cones(fan)) == n_clusters(q)          == 4
        @test is_complete(fan)
    end

    @testset "A2 fan is the pentagon fan" begin
        fan = g_vector_fan(Quiver(:A, 2))
        @test length(rays(fan)) == 5
        # The initial cluster's cone is spanned by the standard basis vectors.
        @test Set([[1, 0], [0, 1]]) ⊆ Set(rays(fan))
        # Sum over the rays of a complete rank-2 cluster fan: the five g-vectors
        # of A2 are e1, e2, -e1, -e2, e1-e2 up to the initial-seed basis.
        @test sum(sum(r) for r in rays(fan)) == 0
    end

    @testset "rays are a basis-dependent labelling" begin
        # Mutating the initial seed changes the g-vector basis but not the
        # combinatorics of the fan (ledger, src/cluster_fan.jl).
        f1 = g_vector_fan(Seed(Quiver(:A, 3)))
        f2 = g_vector_fan(mutate(Seed(Quiver(:A, 3)), 2))
        @test length(rays(f1))          == length(rays(f2))
        @test length(maximal_cones(f1)) == length(maximal_cones(f2))
        @test is_complete(f2)
        # Isomorphic fans: each ray lies in the same multiset of cone counts.
        degrees(f) = sort([count(c -> i in c, maximal_cones(f)) for i in 1:length(rays(f))])
        @test degrees(f1) == degrees(f2)
        @test rays(f1) != rays(f2)                          # but not in the same basis
    end

    @testset "truncation flags an incomplete fan" begin
        kronecker = Quiver([0 2; -2 0])                     # affine, infinite type
        fan       = g_vector_fan(kronecker; max_seeds = 20)
        @test is_truncated(fan)
        @test !is_complete(fan)
        @test length(maximal_cones(fan)) == 20
    end

    @testset "cluster complex agrees with the exchange graph" begin
        for (t, r) in [(:A, 2), (:A, 3), (:B, 3), (:D, 4)]
            q  = Quiver(t, r)
            cc = cluster_complex(q)
            eg = exchange_graph(Seed(q))
            @test nvertices(cc)     == n_cluster_variables(q)
            @test length(facets(cc)) == n_clusters(q)
            @test f_vector(cc)      == f_vector(eg)
            @test f_vector(cc)[1]   == 1
            @test f_vector(cc)[2]   == n_cluster_variables(q)
            @test f_vector(cc)[end] == n_clusters(q)
            @test length(vertex_labels(cc)) == nvertices(cc)
        end
    end

    @testset "the fan and the complex have the same face numbers" begin
        # Same cluster combinatorics, two independent labellings (g-vectors vs
        # d-vectors), so the facet index sets must agree up to relabelling.
        q   = Quiver(:A, 3)
        fan = g_vector_fan(q)
        cc  = cluster_complex(q)
        @test length(maximal_cones(fan)) == length(facets(cc))
        @test sort(map(length, maximal_cones(fan))) == sort(map(length, facets(cc)))
    end

    @testset "display" begin
        @test occursin("ClusterFan: 5 rays", sprint(show, g_vector_fan(Quiver(:A, 2))))
        @test occursin("(truncated)",
                       sprint(show, g_vector_fan(Quiver([0 2; -2 0]); max_seeds = 5)))
        @test occursin("ClusterComplex: 9 vertices, 14 facets",
                       sprint(show, cluster_complex(Quiver(:A, 3))))
    end
end
