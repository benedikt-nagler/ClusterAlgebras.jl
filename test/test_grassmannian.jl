using Test
using ClusterAlgebras
using AbstractAlgebra

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

        # Rectangles seed: mutable labels row-major over interior rectangles,
        # I(a,b) = {1,…,k−a} ∪ {k−a+b+1,…,k+b}
        q36 = Quiver(:Grassmannian, 3, 6)
        @test q36.labels[1] == "p_{124}"
        @test q36.labels[2] == "p_{125}"
        @test q36.labels[3] == "p_{134}"
        @test q36.labels[4] == "p_{145}"
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

    # ─── Oracle helpers: numeric Plücker minors of a totally positive matrix ──

    # Recursive cofactor determinant (exact, Rational{BigInt}).
    function _det(M)
        m = size(M, 1)
        m == 1 && return M[1, 1]
        d = zero(M[1, 1])
        for j in 1:m
            d += (isodd(j) ? 1 : -1) * M[1, j] * _det(M[2:end, [1:j-1; j+1:m]])
        end
        return d
    end

    # k×n Vandermonde with nodes 1,…,n: totally positive, so ALL Plücker
    # minors are strictly positive (no accidental vanishing, unambiguous signs).
    _gr_matrix(k, n) = [Rational{BigInt}(j)^(i - 1) for i in 1:k, j in 1:n]

    _minor(A, cols) = _det(A[:, collect(cols)])

    # Evaluate a cluster variable (element of Frac(ZZ[p...])) at the numeric
    # Plücker values of the seed's ring variables.
    _eval_var(v, vals) = evaluate(numerator(v), vals) // evaluate(denominator(v), vals)

    _plucker_values(s, A) =
        [_minor(A, plucker_subset(lbl)) for lbl in s.quiver.labels]

    # All k-element subsets of the vector v.
    function _ksubsets(v::Vector{Int}, k::Int)
        k == 0 && return [Int[]]
        length(v) < k && return Vector{Int}[]
        head, rest = v[1], v[2:end]
        return vcat([vcat([head], s) for s in _ksubsets(rest, k - 1)],
                    _ksubsets(rest, k))
    end

    # S, T weakly separated: no a < b < c < d with a,c ∈ S∖T, b,d ∈ T∖S or vice versa.
    function _weakly_separated(S, T)
        Sd, Td = setdiff(S, T), setdiff(T, S)
        for (X, Y) in ((Sd, Td), (Td, Sd))
            for a in X, c in X, b in Y, d in Y
                a < b < c < d && return false
            end
        end
        return true
    end

    @testset "initial cluster is weakly separated (Oh–Postnikov–Speyer)" begin
        for (k, n) in ((2, 5), (2, 6), (3, 6), (3, 7), (4, 7), (4, 8))
            q = Quiver(:Grassmannian, k, n)
            subsets = [plucker_subset(lbl) for lbl in q.labels]
            @test all(_weakly_separated(subsets[i], subsets[j])
                      for i in 1:length(subsets) for j in i+1:length(subsets))
        end
    end

    @testset "mutation yields Plücker minors — numeric oracle" begin
        # Gr(2,5): μ at p13 must give p24 (Plücker relation
        # p13·p24 = p12·p34 + p14·p23); μ at p14 must give p35.
        s  = grassmannian(2, 5)
        A  = _gr_matrix(2, 5)
        vv = _plucker_values(s, A)
        @test _eval_var(mutate(s, 1)[1], vv) == _minor(A, [2, 4])
        @test _eval_var(mutate(s, 2)[2], vv) == _minor(A, [3, 5])

        # Gr(3,6): μ at the three rectangles whose exchange is a short Plücker
        # relation.  p124→p135, p125→p146, p134→p245 (hand-verified relations).
        s  = grassmannian(3, 6)
        A  = _gr_matrix(3, 6)
        vv = _plucker_values(s, A)
        @test _eval_var(mutate(s, 1)[1], vv) == _minor(A, [1, 3, 5])
        @test _eval_var(mutate(s, 2)[2], vv) == _minor(A, [1, 4, 6])
        @test _eval_var(mutate(s, 3)[3], vv) == _minor(A, [2, 4, 5])
    end

    @testset "symbol alphabet as Plücker functions — numeric oracle" begin
        # Gr(4,6) ≅ Gr(2,6) is type A₃: every one of the 9 letters must be a
        # Plücker minor of the 4×6 matrix (this exercises k = 4 frozen arrows).
        s  = grassmannian(4, 6)
        A  = _gr_matrix(4, 6)
        vv = _plucker_values(s, A)
        minors46 = Set(_minor(A, cols) for cols in _ksubsets(collect(1:6), 4))
        alpha = symbol_alphabet(s)
        @test length(alpha) == 9
        @test all(_eval_var(v, vv) in minors46 for v in alpha)

        # Gr(3,6) has 16 letters: 14 mutable Plückers + 2 quadratics.
        s  = grassmannian(3, 6)
        A  = _gr_matrix(3, 6)
        vv = _plucker_values(s, A)
        minors36 = Set(_minor(A, cols) for cols in _ksubsets(collect(1:6), 3))
        alpha = symbol_alphabet(s)
        @test length(alpha) == 16
        @test count(_eval_var(v, vv) in minors36 for v in alpha) == 14
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
