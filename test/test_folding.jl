# ─── test_folding.jl ──────────────────────────────────────────────────────────
# Folding a symmetric quiver by an admissible automorphism → skew-symmetrizable
# (non-simply-laced) quiver. Classic A/D/E → B/C/F/G foldings as exact oracles.

@testset "folding - admissible automorphisms" begin
    # A₃ → rank-2 double bond. Orientation 1→2←3 is (1 3)-invariant.
    A3 = Quiver([0 1 0; -1 0 -1; 0 1 0])
    σ  = [3, 2, 1]                    # swap 1 ↔ 3, fix 2
    @test is_admissible_folding(A3, σ)
    f  = fold(A3, σ)
    @test f.B == [0 1; -2 0]
    @test f.d == [2, 1]
    @test f.labels == ["1+3", "2"]
    @test is_finite_type(f)
    @test cartan_type(f) == (:B, 2)   # rank-2 double bond (B₂ ≅ C₂)

    # A₅ → rank-3. Orientation 1→2→3←4←5 is (1 5)(2 4)-invariant.
    A5 = Quiver([ 0  1  0  0  0;
                 -1  0  1  0  0;
                  0 -1  0 -1  0;
                  0  0  1  0 -1;
                  0  0  0  1  0])
    σ5 = [5, 4, 3, 2, 1]
    @test is_admissible_folding(A5, σ5)
    f5 = fold(A5, σ5)
    @test f5.B == [0 1 0; -1 0 1; 0 -2 0]
    @test f5.d == [2, 2, 1]
    @test is_finite_type(f5)
    @test cartan_type(f5) == (:B, 3)

    # D₄ → G₂. Center 1, legs 2,3,4 all oriented outward; fold by the 3-cycle.
    D4 = Quiver([0 1 1 1; -1 0 0 0; -1 0 0 0; -1 0 0 0])
    σg = [1, 3, 4, 2]                 # fix 1, cycle 2→3→4
    @test is_admissible_folding(D4, σg)
    fg = fold(D4, σg)
    @test fg.B == [0 3; -1 0]
    @test fg.d == [1, 3]
    @test is_finite_type(fg)
    @test cartan_type(fg) == (:G, 2)

    # Identity fold is a no-op on B and d.
    id = fold(A3, [1, 2, 3])
    @test id.B == A3.B && id.d == A3.d

    # Guards.
    A2 = Quiver(:A, 2)
    @test_throws InvalidArgument fold(A2, [2, 1])          # (1 2) is not a B-automorphism
    @test_throws InvalidArgument fold(A3, [2, 2, 1])       # not a permutation
    # Oriented 3-cycle: (1 2 3) is an automorphism but its single orbit is
    # connected ⇒ not admissible.
    cyc = Quiver([0 1 -1; -1 0 1; 1 -1 0])
    @test !is_admissible_folding(cyc, [2, 3, 1])
    @test_throws InvalidArgument fold(cyc, [2, 3, 1])
    # Frozen vertices are rejected.
    withfrozen = Quiver([0 1 0; -1 0 0; 0 0 0], 2)         # vertex 3 frozen
    @test !is_admissible_folding(withfrozen, [1, 2])
end
