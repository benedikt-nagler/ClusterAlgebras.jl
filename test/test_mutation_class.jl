@testset "mutation_class — quiver-level BFS" begin

    # ── A₂: two distinct exchange matrices ─────────────────────────────────────
    mc = mutation_class(Quiver(:A, 2))
    @test length(mc)         == 2
    @test !is_truncated(mc)
    # first quiver is the initial one
    @test mc[1] == Quiver(:A, 2)
    # each quiver has n_mutable = 2 adjacency entries
    @test length(mc.adj[1])  == 2
    @test length(mc.adj[2])  == 2
    # adjacency indices are valid
    @test all(all(1 <= j <= length(mc) for j in mc.adj[i]) for i in eachindex(mc.quivers))

    # ── A₁: trivial quiver, one exchange matrix ─────────────────────────────────
    mc1 = mutation_class(Quiver([0;;]))
    @test length(mc1) == 1

    # ── Kronecker quiver [0 2; -2 0]: mutation-finite (two quivers) ────────────
    mc_kron = mutation_class(Quiver([0 2; -2 0]))
    @test length(mc_kron) == 2
    @test !is_truncated(mc_kron)

    # ── (1,1,3) cyclic triangle: mutation-infinite (entries grow) → truncation ──
    # B_{12}=1, B_{23}=1, B_{31}=3 and antisymmetric; abc=3 ≥ 4 fails for integers
    # but mutation produces entries ≥5, confirming infinite mutation class
    B_wild = [0 1 -3; -1 0 1; 3 -1 0]
    mc_inf  = mutation_class(Quiver(B_wild); max_quivers = 20)
    @test is_truncated(mc_inf)
    @test length(mc_inf) == 20

end

@testset "exchange_graph — seed-level BFS" begin

    # ── A₂: 5 seeds (pentagon) ─────────────────────────────────────────────────
    eg2 = exchange_graph(Seed(Quiver(:A, 2)))
    @test length(eg2)         == 5
    @test !is_truncated(eg2)
    # each vertex has exactly n_mutable = 2 adjacency entries
    @test all(length(eg2.adj[i]) == 2 for i in 1:5)
    # each vertex has n_mutable = 2 adjacency entries, all pointing to valid indices
    @test all(all(1 <= j <= 5 for j in eg2.adj[i]) for i in 1:5)
    # all 2 neighbors of each vertex are distinct (exchange graph is simple)
    @test all(length(unique(eg2.adj[i])) == 2 for i in 1:5)

    # ── A₃: 14 seeds ───────────────────────────────────────────────────────────
    eg3 = exchange_graph(Seed(Quiver(:A, 3)))
    @test length(eg3) == 14
    @test !is_truncated(eg3)

    # ── B₂: 6 seeds ────────────────────────────────────────────────────────────
    eg_B2 = exchange_graph(Seed(Quiver(:B, 2)))
    @test length(eg_B2) == 6
    @test !is_truncated(eg_B2)

    # ── G₂: 8 seeds ────────────────────────────────────────────────────────────
    eg_G2 = exchange_graph(Seed(Quiver(:G, 2)))
    @test length(eg_G2) == 8
    @test !is_truncated(eg_G2)

    # ── D₄: 50 seeds ───────────────────────────────────────────────────────────
    eg_D4 = exchange_graph(Seed(Quiver(:D, 4)))
    @test length(eg_D4) == 50
    @test !is_truncated(eg_D4)

    # ── getindex returns a Seed with the right quiver ──────────────────────────
    @test eg2[1] isa Seed
    @test eg2[1].quiver == Quiver(:A, 2)

    # ── Markov quiver: infinite type → truncation ───────────────────────────────
    B_markov = [0 2 -2; -2 0 2; 2 -2 0]
    eg_inf   = exchange_graph(Seed(Quiver(B_markov)); max_seeds = 15)
    @test is_truncated(eg_inf)
    @test length(eg_inf) == 15

end
