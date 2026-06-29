using Test
using ClusterAlgebras

@testset "Grassmannian cluster algebras" begin

    @testset "Plücker label helpers" begin
        @test plucker_label([1, 3, 5]) == "p_{135}"
        @test plucker_label([3, 1, 5]) == "p_{135}"   # sorting
        @test plucker_label([1, 2])    == "p_{12}"
        @test plucker_subset("p_{135}") == [1, 3, 5]
        @test plucker_subset("p_{12}")  == [1, 2]
        @test is_plucker_label("p_{135}")
        @test is_plucker_label("p_{12}")
        @test !is_plucker_label("x_1")
        @test !is_plucker_label("p135")
        @test_throws ClusterAlgebraError plucker_subset("x_1")
    end

    @testset "Quiver structure" begin
        # Gr(2,5) ≅ A_2: 2 mutable, 5 frozen
        q25 = Quiver(:Grassmannian, 2, 5)
        @test q25.n_mutable == 2
        @test q25.n_frozen  == 5
        @test all(is_frozen(q25, i) for i in 3:7)

        # Gr(2,6) ≅ A_3: 3 mutable, 6 frozen
        q26 = Quiver(:Grassmannian, 2, 6)
        @test q26.n_mutable == 3
        @test q26.n_frozen  == 6

        # Gr(3,6) ≅ D_4: 4 mutable, 6 frozen
        q36 = Quiver(:Grassmannian, 3, 6)
        @test q36.n_mutable == 4
        @test q36.n_frozen  == 6

        # Gr(4,7) ≅ E_6: 6 mutable, 7 frozen
        q47 = Quiver(:Grassmannian, 4, 7)
        @test q47.n_mutable == 6
        @test q47.n_frozen  == 7

        # B-matrix is antisymmetric
        for q in (q25, q26, q36, q47)
            @test q.B == -q.B'
        end
    end

    @testset "Plücker labels on initial quiver" begin
        q25 = Quiver(:Grassmannian, 2, 5)
        @test q25.labels[1] == "p_{13}"
        @test q25.labels[2] == "p_{14}"
        @test "p_{12}" in q25.labels

        q36 = Quiver(:Grassmannian, 3, 6)
        @test q36.labels[1] == "p_{136}"
        @test q36.labels[2] == "p_{145}"
        @test q36.labels[3] == "p_{236}"
        @test q36.labels[4] == "p_{245}"
    end

    @testset "grassmannian() equals Quiver(:Grassmannian, k, n)" begin
        @test grassmannian(2, 5).quiver == Quiver(:Grassmannian, 2, 5)
        @test grassmannian(3, 6).quiver == Quiver(:Grassmannian, 3, 6)
    end

    @testset "Type recognition — oracle table" begin
        @test is_finite_type(grassmannian(2, 5).quiver)
        @test cartan_type(grassmannian(2, 5).quiver) == (:A, 2)
        @test cartan_type(grassmannian(2, 6).quiver) == (:A, 3)
        @test cartan_type(grassmannian(3, 6).quiver) == (:D, 4)
        @test cartan_type(grassmannian(4, 6).quiver) == (:A, 3)  # duality k ↔ n-k
        @test cartan_type(grassmannian(4, 7).quiver) == (:E, 6)
    end

    @testset "Enumerative invariants — oracle table" begin
        @test n_cluster_variables(grassmannian(2, 5).quiver) == 5
        @test n_clusters(grassmannian(2, 5).quiver)          == 5
        @test n_cluster_variables(grassmannian(2, 6).quiver) == 9
        @test n_clusters(grassmannian(2, 6).quiver)          == 14
        @test n_cluster_variables(grassmannian(3, 6).quiver) == 16
        @test n_clusters(grassmannian(3, 6).quiver)          == 50
    end

    @testset "grassmannian() Seed construction" begin
        s = grassmannian(2, 5)
        @test length(s) == 2 + 5   # mutable + frozen variables
        @test all(isone(denominator(s[i])) for i in 1:length(s))
    end

    @testset "x_coordinates" begin
        xc = x_coordinates(extend(grassmannian(2, 5)))
        @test length(xc) == 2   # rank of A_2
    end

    @testset "Invalid arguments" begin
        @test_throws ClusterAlgebraError grassmannian(1, 5)
        @test_throws ClusterAlgebraError grassmannian(4, 5)
    end

end
