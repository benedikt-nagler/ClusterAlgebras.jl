# Tests for enumerative invariants (Track A: A1).
#
# Three independent routes must agree on every oracle:
#   Route 1 — closed form via RootSystem (Coxeter h, degrees dᵢ)
#   Route 2 — BFS exchange-graph count (length of ExchangeGraph)
#   Route 3 — f_vector (f₀ = n_cluster_variables, f_{n-1} = n_clusters)

# ─── Oracle table ─────────────────────────────────────────────────────────────
#
#  Type  | h  | n_cluster_variables | n_clusters
#  A₂    | 3  | 5                   | 5
#  A₃    | 4  | 9                   | 14
#  B₂=C₂ | 4  | 6                   | 6
#  G₂    | 6  | 8                   | 8
#  D₄    | 6  | 16                  | 50

const _ENUM_ORACLES = [
    # (type, n, h, ncv, nc)
    (:A, 2, 3,  5,  5),
    (:A, 3, 4,  9, 14),
    (:B, 2, 4,  6,  6),
    (:G, 2, 6,  8,  8),
    (:D, 4, 6, 16, 50),
]

# ─── Route 1 — closed forms via RootSystem ────────────────────────────────────

@testset "enumerative — closed forms (RootSystem)" begin
    for (type, n, h, ncv_oracle, nc_oracle) in _ENUM_ORACLES
        rs = RootSystem(type, n)
        @test rs.coxeter_number == h
        @test n_cluster_variables(rs) == ncv_oracle
        @test n_clusters(rs)          == nc_oracle
    end
end

# ─── Almost-positive roots count matches n_cluster_variables ─────────────────

@testset "enumerative — n_cluster_variables == length(almost_positive_roots)" begin
    for (type, n, _, ncv_oracle, _) in _ENUM_ORACLES
        rs = RootSystem(type, n)
        @test length(almost_positive_roots(rs)) == ncv_oracle
    end
end

# ─── Quiver overloads delegate correctly ──────────────────────────────────────

@testset "enumerative — Quiver overloads" begin
    @test n_cluster_variables(Quiver(:A, 2)) == 5
    @test n_clusters(Quiver(:A, 2))          == 5
    @test n_cluster_variables(Quiver(:A, 3)) == 9
    @test n_clusters(Quiver(:A, 3))          == 14
    @test n_cluster_variables(Quiver(:D, 4)) == 16
    @test n_clusters(Quiver(:D, 4))          == 50
end

# ─── Route 2 — BFS exchange-graph counts ─────────────────────────────────────

@testset "enumerative — exchange graph counts" begin
    for (type, n, _, ncv_oracle, nc_oracle) in _ENUM_ORACLES
        q  = Quiver(type, n)
        eg = exchange_graph(Seed(q); max_seeds = 200)
        @test !is_truncated(eg)
        # number of seeds in the exchange graph = number of clusters
        @test length(eg) == nc_oracle
    end
end

# ─── Route 3 — f-vector cross-checks ─────────────────────────────────────────

@testset "enumerative — f_vector structure" begin
    for (type, n, _, ncv_oracle, nc_oracle) in _ENUM_ORACLES
        s  = Seed(Quiver(type, n))
        fv = f_vector(s; max_seeds = 200)
        @test length(fv) == n + 1          # [f_{-1}, f_0, …, f_{n-1}]
        @test fv[1]   == 1                 # f_{-1} = 1 (empty face)
        @test fv[2]   == ncv_oracle        # f_0 = number of cluster variables
        @test fv[end] == nc_oracle         # f_{n-1} = number of clusters
        @test all(>=(0), fv)
    end
end

# ─── h-vector: non-negative entries, sum = n_clusters ────────────────────────

@testset "enumerative — h_vector non-negative and sums to Cat(W)" begin
    for (type, n, _, _, nc_oracle) in _ENUM_ORACLES
        s  = Seed(Quiver(type, n))
        hv = h_vector(s; max_seeds = 200)
        @test length(hv) == n + 1
        @test all(>=(0), hv)
        @test sum(hv) == nc_oracle
    end
end

# ─── A₂ Narayana numbers: h = [1, 1, …] wait — A₂ gives [1,3,1] for A₃ ─────
# Classical Narayana N(m+1, k+1) for type A_m (0-indexed k=0..m):
#   A₂ (m=2): h = [1, 1] → sum=2 ✗ …
#   Actually for A_m: n_clusters = Cat_{m+1}, h_k = N(m+1, k) for k=1..m,
#   and h_0=1.  A₂ has 5 clusters: [1,3,1]. A₃ has 14: [1,6,6,1].

@testset "enumerative — A₃ Narayana h-vector oracle" begin
    # Type A₃: h_vector = [1, 6, 6, 1]  (Narayana N(4,1)..N(4,4) = 1,6,6,1)
    hv = h_vector(Seed(Quiver(:A, 3)); max_seeds = 50)
    @test hv == [1, 6, 6, 1]
end

@testset "enumerative — A₂ Narayana h-vector oracle" begin
    # Type A₂: h_vector = [1, 3, 1]  (Narayana N(3,1)..N(3,3) = 1,3,1)
    hv = h_vector(Seed(Quiver(:A, 2)); max_seeds = 20)
    @test hv == [1, 3, 1]
end

# ─── Three-route agreement (closed form, BFS count, f-vector) ─────────────────

@testset "enumerative — three routes agree" begin
    for (type, n, _, _, _) in _ENUM_ORACLES
        rs       = RootSystem(type, n)
        q        = Quiver(type, n)
        s        = Seed(q)
        eg       = exchange_graph(s; max_seeds = 200)
        fv       = f_vector(eg)

        nc_closed = n_clusters(rs)
        nc_bfs    = length(eg)
        nc_fvec   = fv[end]

        @test nc_closed == nc_bfs
        @test nc_closed == nc_fvec

        ncv_closed = n_cluster_variables(rs)
        ncv_fvec   = fv[2]
        @test ncv_closed == ncv_fvec
    end
end

# ─── Error path: non-finite-type quiver ──────────────────────────────────────

@testset "enumerative — throws on non-finite-type quiver" begin
    # Markov quiver: 3 vertices, all double arrows in a cycle — indefinite
    markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
    @test_throws ClusterAlgebraError n_cluster_variables(markov)
    @test_throws ClusterAlgebraError n_clusters(markov)
end

@testset "enumerative — throws on truncated exchange graph" begin
    # An infinite-type quiver (Markov) will truncate; f_vector should throw
    markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
    @test_throws ClusterAlgebraError f_vector(Seed(markov); max_seeds = 5)
end
