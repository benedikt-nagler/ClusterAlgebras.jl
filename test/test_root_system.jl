@testset "cartan_companion" begin
    @test cartan_companion(Quiver(:A, 2)) == [2 -1; -1 2]
    @test cartan_companion(Quiver(:B, 2)) == [2 -1; -2 2]
    @test cartan_companion(Quiver(:G, 2)) == [2 -1; -3 2]
    @test cartan_companion(Quiver(:A, 3)) == [2 -1 0; -1 2 -1; 0 -1 2]
end

@testset "RootSystem positive roots" begin
    # A_2: 3 positive roots, h=3, exponents [1,2], |W|=6
    rs = RootSystem(:A, 2)
    @test length(rs.positive_roots) == 3
    @test rs.coxeter_number == 3
    @test rs.exponents == [1, 2]
    @test rs.weyl_group_order == 6
    # check: n·h/2 = 2·3/2 = 3 ✓
    @test length(rs.positive_roots) == rs.n * rs.coxeter_number ÷ 2

    # A_3: 6 positive roots, h=4, |W|=24
    rs = RootSystem(:A, 3)
    @test length(rs.positive_roots) == 6
    @test rs.coxeter_number == 4
    @test rs.exponents == [1, 2, 3]
    @test rs.weyl_group_order == 24

    # B_2 = C_2: 4 positive roots, h=4, exponents [1,3], |W|=8
    rs = RootSystem(:B, 2)
    @test length(rs.positive_roots) == 4
    @test rs.coxeter_number == 4
    @test rs.exponents == [1, 3]
    @test rs.weyl_group_order == 8

    # G_2: 6 positive roots, h=6, exponents [1,5], |W|=12
    rs = RootSystem(:G, 2)
    @test length(rs.positive_roots) == 6
    @test rs.coxeter_number == 6
    @test rs.exponents == [1, 5]
    @test rs.weyl_group_order == 12

    # D_4: 12 positive roots, h=6
    rs = RootSystem(:D, 4)
    @test length(rs.positive_roots) == 12
    @test rs.coxeter_number == 6
    @test rs.exponents == [1, 3, 3, 5]
    @test rs.weyl_group_order == 192   # 2^3 * 4! = 8*24

    # E_6: 36 positive roots, h=12, |W|=51840
    rs = RootSystem(:E, 6)
    @test length(rs.positive_roots) == 36
    @test rs.coxeter_number == 12
    @test rs.weyl_group_order == 51840

    # F_4: 24 positive roots, h=12
    rs = RootSystem(:F, 4)
    @test length(rs.positive_roots) == 24
    @test rs.coxeter_number == 12
    @test rs.exponents == [1, 5, 7, 11]
    @test rs.weyl_group_order == 1152

    # n·h/2 invariant for all finite types
    for (t, n) in [(:A,2),(:A,3),(:B,2),(:C,3),(:D,4),(:F,4),(:G,2)]
        rs = RootSystem(t, n)
        @test length(rs.positive_roots) == rs.n * rs.coxeter_number ÷ 2
    end
end

@testset "almost_positive_roots" begin
    rs = RootSystem(:A, 2)
    aprs = almost_positive_roots(rs)
    # 2 negative simples + 3 positive roots = 5
    @test length(aprs) == 5
    @test [-1,  0] ∈ aprs
    @test [ 0, -1] ∈ aprs
    # all positive roots also present
    for r in rs.positive_roots
        @test r ∈ aprs
    end

    # G_2: 2 + 6 = 8 almost-positive roots
    rs = RootSystem(:G, 2)
    @test length(almost_positive_roots(rs)) == 8
end

@testset "is_finite_type" begin
    # All named finite Dynkin types must pass
    for (t, n) in [(:A,2),(:A,3),(:B,2),(:C,2),(:D,4),(:E,6),(:E,7),(:E,8),(:F,4),(:G,2)]
        @test is_finite_type(Quiver(t, n))
    end

    # Kronecker quiver (affine Ã₁): mutation-finite but NOT finite type (infinite clusters)
    @test !is_finite_type(Quiver([0 2; -2 0]))

    # Markov quiver: mutation-infinite, definitely not finite type
    @test !is_finite_type(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
end

@testset "is_affine_type" begin
    # Kronecker quiver (affine Ã₁): symmetrized Cartan companion is [2 -2; -2 2], det=0
    @test is_affine_type(Quiver([0 2; -2 0]))

    # Markov quiver: indefinite - not affine
    @test !is_affine_type(Quiver([0 2 -2; -2 0 2; 2 -2 0]))

    # Finite types are not affine
    @test !is_affine_type(Quiver(:A, 2))
    @test !is_affine_type(Quiver(:D, 4))
end

@testset "cartan_type" begin
    @test cartan_type(Quiver(:A, 2)) == (:A, 2)
    @test cartan_type(Quiver(:A, 3)) == (:A, 3)
    @test cartan_type(Quiver(:B, 2)) == (:B, 2)
    @test cartan_type(Quiver(:C, 2)) == (:C, 2)
    @test cartan_type(Quiver(:C, 3)) == (:C, 3)
    @test cartan_type(Quiver(:D, 4)) == (:D, 4)
    @test cartan_type(Quiver(:E, 6)) == (:E, 6)
    @test cartan_type(Quiver(:E, 7)) == (:E, 7)
    @test cartan_type(Quiver(:E, 8)) == (:E, 8)
    @test cartan_type(Quiver(:F, 4)) == (:F, 4)
    @test cartan_type(Quiver(:G, 2)) == (:G, 2)

    # Non-finite type throws
    @test_throws ClusterAlgebraError cartan_type(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
end

# ─── Reducible (disconnected) finite types ────────────────────────────────────
#
# Regression for a silent-mislabel bug: cartan_type read the type from rank alone
# in the simply-laced fall-through, so every rank-6 finite-type quiver that was
# not A₆/D₆ (including the disconnected A1⊔A5, 16 positive roots) was returned as
# E₆ (36 positive roots).  A single (type, rank) pair cannot name a direct sum;
# cartan_types decomposes, and cartan_type now throws on reducible input.
@testset "cartan_types - reducible decompositions" begin
    # Block-diagonal direct sum of named Dynkin quivers.
    direct_sum(types) = begin
        n = sum(r for (_, r) in types)
        B = zeros(Int, n, n); o = 0
        for (t, r) in types
            B[o+1:o+r, o+1:o+r] = Quiver(t, r).B; o += r
        end
        Quiver(B)
    end

    @test cartan_types(Quiver(:A, 5)) == [(:A, 5)]                  # connected ⇒ singleton
    @test cartan_types(direct_sum([(:A, 1), (:A, 5)])) == [(:A, 1), (:A, 5)]
    @test cartan_types(direct_sum([(:A, 3), (:A, 3)])) == [(:A, 3), (:A, 3)]
    @test cartan_types(direct_sum([(:D, 4), (:A, 2)])) == [(:A, 2), (:D, 4)]  # sorted
    @test cartan_types(direct_sum([(:A, 1) for _ in 1:6])) == [(:A, 1) for _ in 1:6]
    @test cartan_types(Quiver(zeros(Int, 0, 0))) == Tuple{Symbol, Int}[]      # empty

    # The bug: A1⊔A5 must NOT be E₆, and cartan_type must refuse a reducible type.
    @test_throws ClusterAlgebraError cartan_type(direct_sum([(:A, 1), (:A, 5)]))

    # Non-finite-type components still throw.
    @test_throws ClusterAlgebraError cartan_types(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
end

@testset "n_clusters / n_cluster_variables - reducible types" begin
    direct_sum(types) = begin
        n = sum(r for (_, r) in types)
        B = zeros(Int, n, n); o = 0
        for (t, r) in types
            B[o+1:o+r, o+1:o+r] = Quiver(t, r).B; o += r
        end
        Quiver(B)
    end
    # n_clusters multiplies over components, n_cluster_variables adds.
    @test n_clusters(direct_sum([(:A, 1), (:A, 5)])) == 2 * 132          # = 264, not E₆'s 833
    @test n_clusters(direct_sum([(:A, 3), (:A, 3)])) == 14 * 14
    @test n_cluster_variables(direct_sum([(:A, 1), (:A, 5)])) ==
          n_cluster_variables(Quiver(:A, 1)) + n_cluster_variables(Quiver(:A, 5))
    @test n_clusters(direct_sum([(:A, 1) for _ in 1:6])) == 2^6
end

# ─── Soundness on non-acyclic quivers (the headline bug fix) ──────────────────
#
# Before this fix, is_finite_type / is_affine_type tested only the Cartan
# companion of the *given* seed - valid only for acyclic quivers.  The oriented
# 3-cycle below is mutation-equivalent to A₃ (finite type) but has an affine-
# looking Cartan companion, so the old code returned is_finite_type=false and
# is_affine_type=true - both wrong.

@testset "is_finite_type - non-acyclic input (regression)" begin
    # Oriented 3-cycle: mutate linear A₃ at vertex 2 → still mutation-class A₃
    q_cyc = mutate(Quiver([0 1 0; -1 0 1; 0 -1 0]), 2)
    @test q_cyc.B == [0 -1 1; 1 0 -1; -1 1 0]   # confirm the quiver
    @test is_finite_type(q_cyc)                   # must be true (was: false - BUG)
    @test !is_affine_type(q_cyc)                  # must be false (was: true - BUG)
    @test cartan_type(q_cyc) == (:A, 3)           # must round-trip correctly
end

@testset "is_finite_type - all acyclic named types still pass" begin
    for (t, n) in [(:A,2),(:A,3),(:B,2),(:C,2),(:D,4),(:E,6),(:F,4),(:G,2)]
        @test is_finite_type(Quiver(t, n))
    end
end

@testset "is_affine_type - non-acyclic input stays correct" begin
    # Ã₁: Kronecker quiver, acyclic - fast path.
    @test is_affine_type(Quiver([0 2; -2 0]))
    # Genuinely affine rank 3: an ACYCLIC orientation of the triangle is Ã₂
    # (the cyclically oriented triangle is instead mutation-equivalent to A₃).
    q_A2tilde = Quiver([0 1 1; -1 0 1; -1 -1 0])
    @test is_affine_type(q_A2tilde)
    @test !is_finite_type(q_A2tilde)
    # Cyclic mutation of Ã₂ must still be recognized via the BFS path.
    @test is_affine_type(mutate(q_A2tilde, 2))
    # Confirm finite types are not affine via the non-acyclic path too.
    q_cyc = mutate(Quiver([0 1 0; -1 0 1; 0 -1 0]), 2)
    @test !is_affine_type(q_cyc)
end

@testset "cartan_type - B/C tie-break is permutation-invariant" begin
    # Baseline orientation/labeling from the named constructors.
    @test cartan_type(Quiver(:B, 3)) == (:B, 3)
    @test cartan_type(Quiver(:C, 3)) == (:C, 3)

    # Vertex order reversed (long/short root relabeled).  A positional
    # tie-break like d[1] ≥ d[end] misclassifies both of these.
    B_C3_rev = [0 -1 0; 2 0 -1; 0 1 0]   # reversed C₃, d = (2,1,1)
    @test cartan_type(Quiver(B_C3_rev, 3, [2, 1, 1])) == (:C, 3)
    B_B3_rev = [0 -2 0; 1 0 -1; 0 1 0]   # reversed B₃, d = (1,2,2)
    @test cartan_type(Quiver(B_B3_rev, 3, [1, 2, 2])) == (:B, 3)
end

@testset "n_clusters / n_cluster_variables - non-acyclic finite-type input" begin
    # These delegate to cartan_type, so they inherit the soundness fix.
    q_cyc = mutate(Quiver([0 1 0; -1 0 1; 0 -1 0]), 2)
    @test n_cluster_variables(q_cyc) == 9    # same as A₃
    @test n_clusters(q_cyc)          == 14   # same as A₃
end
