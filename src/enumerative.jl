# Enumerative invariants of finite-type cluster algebras (Track A: A1).
#
# Two independent derivations are provided for each invariant:
#   Route 1 - closed forms from the root system (Coxeter number + degrees).
#   Route 2 - combinatorial counts from the exchange-graph BFS.
# They must agree on every finite type; use both in tests as mutual cross-checks.

# ─── Closed forms from the root system ───────────────────────────────────────

"""
    n_cluster_variables(rs::RootSystem) -> Int

Return the number of cluster variables of the finite-type cluster algebra whose
Cartan type equals `rs`.  Closed form: `n(h+2)/2`, where `n` is the rank and
`h` is the Coxeter number.  This also equals `length(almost_positive_roots(rs))`.
"""
function n_cluster_variables(rs::RootSystem)
    n, h = rs.n, rs.coxeter_number
    return n * (h + 2) ÷ 2
end

"""
    n_cluster_variables(q::Quiver) -> Int

Return the number of cluster variables of the finite-type cluster algebra
defined by the quiver `q`.

Requires `is_finite_type(q)`; throws `InvalidArgument` otherwise.  Reducible
(disconnected) types are handled: the cluster variables of a direct sum are the
disjoint union of the summands', so the count is **additive** over components.
"""
function n_cluster_variables(q::Quiver)
    # Additive over irreducible components (almost-positive roots are a disjoint
    # union).  cartan_types throws unless q is finite type.
    return sum(n_cluster_variables(RootSystem(t, r)) for (t, r) in cartan_types(q); init = 0)
end

"""
    n_clusters(rs::RootSystem) -> Int

Return the number of clusters (W-Catalan number) of the finite-type cluster
algebra whose Cartan type equals `rs`.

Closed form: `Cat(W) = ∏ᵢ (h + dᵢ) / ∏ᵢ dᵢ`, where `dᵢ = eᵢ + 1` are the
degrees and `h` is the Coxeter number.  The product is always an integer.
Computed exactly in `Int128` to handle the large values of types E₇ and E₈.
"""
function n_clusters(rs::RootSystem)
    h  = rs.coxeter_number
    ds = rs.exponents .+ 1          # degrees dᵢ = eᵢ + 1
    num = prod(Int128(h + d) for d in ds)
    den = prod(Int128(d)     for d in ds)
    return Int(num ÷ den)
end

"""
    n_clusters(q::Quiver) -> Int

Return the number of clusters of the finite-type cluster algebra defined by the
quiver `q`.

Requires `is_finite_type(q)`; throws `InvalidArgument` otherwise.  Reducible
(disconnected) types are handled: the exchange graph of a direct sum is the
Cartesian product of the components', so the cluster count is **multiplicative**
over components (`A1 ⊔ A5` has 2·132 = 264 clusters, not the 833 of E₆).
"""
function n_clusters(q::Quiver)
    # Multiplicative over irreducible components (Cartesian product of exchange
    # graphs).  cartan_types throws unless q is finite type.
    return prod(n_clusters(RootSystem(t, r)) for (t, r) in cartan_types(q); init = 1)
end

# ─── Combinatorial counts from the exchange graph ────────────────────────────

# Collect all distinct d-vectors across the full exchange graph, and for each
# seed the set of d-vectors it contains, as integer vectors.
function _cluster_complex_data(eg::ExchangeGraph)
    all_dvecs = Vector{Int}[]
    dvec_index = Dict{Vector{Int}, Int}()
    clusters   = Vector{Int}[]          # each entry = sorted indices into all_dvecs

    for i in 1:length(eg)
        s     = eg[i]::Seed
        n     = length(s.cluster)
        cidxs = Int[]
        for k in 1:n
            dv = denominator_vector(s, k)
            if !haskey(dvec_index, dv)
                push!(all_dvecs, dv)
                dvec_index[dv] = length(all_dvecs)
            end
            push!(cidxs, dvec_index[dv])
        end
        push!(clusters, sort!(cidxs))
    end

    return all_dvecs, clusters
end

"""
    f_vector(eg::ExchangeGraph) -> Vector{Int}

Return the f-vector `[f_{-1}, f_0, f_1, …, f_{n-1}]` of the cluster complex
associated with the exchange graph `eg`, where `n` is the rank.

`f_{k-1}` is the number of (k-1)-dimensional faces (k-element subsets of
cluster variables that appear together in some cluster).  In particular:
- `f_{-1} = 1` (the empty face),
- `f_0` = number of cluster variables,
- `f_{n-1}` = number of clusters (maximal faces).

For a finite-type cluster algebra `eg` must not be truncated; the result is
exact only when `is_truncated(eg) == false`.
"""
function f_vector(eg::ExchangeGraph)
    _, clusters = _cluster_complex_data(eg)
    isempty(clusters) && return Int[1]
    n = length(clusters[1])          # rank = size of each cluster

    # faces[k] = set of sorted k-subsets (as tuples) that appear in some cluster
    faces = [Set{NTuple{k, Int}}() for k in 0:n]

    for cl in clusters
        # all subsets of cl of every size
        for k in 0:n
            for combo in _combinations(cl, k)
                push!(faces[k+1], NTuple{k, Int}(combo))
            end
        end
    end

    # f_{-1}=1 (empty), then f_{k-1} = |faces[k]| for k=1..n  (1-indexed: f[k+1])
    result = Vector{Int}(undef, n + 1)
    result[1] = 1                        # f_{-1}
    for k in 1:n
        result[k+1] = length(faces[k+1])
    end
    return result
end

"""
    f_vector(s::Seed; max_seeds::Int = 5000) -> Vector{Int}

Compute the f-vector via a full BFS exchange graph starting from `s`.
See `f_vector(::ExchangeGraph)` for details.
"""
function f_vector(s::Seed; max_seeds::Int = 5000)
    eg = exchange_graph(s; max_seeds)
    is_truncated(eg) && throw(InvalidArgument(
        "exchange graph was truncated at $max_seeds seeds; increase max_seeds or " *
        "use a finite-type seed"))
    return f_vector(eg)
end

"""
    h_vector(eg::ExchangeGraph) -> Vector{Int}

Return the h-vector `[h_0, h_1, …, h_n]` of the cluster complex, computed
from the f-vector via the standard transform:

    h_k = Σ_{i=0}^{k} (-1)^{k-i} C(n-i, k-i) · f_{i-1}

The entries are the W-Narayana numbers; they are non-negative and sum to the
Catalan number `Cat(W)`.  In type `A_m` they are the classical Narayana
numbers `N(m+1, k)`.
"""
function h_vector(eg::ExchangeGraph)
    fv = f_vector(eg)
    n  = length(fv) - 1              # rank
    hv = Vector{Int}(undef, n + 1)
    for k in 0:n
        s = 0
        for i in 0:k
            s += (-1)^(k - i) * binomial(n - i, k - i) * fv[i + 1]
        end
        hv[k + 1] = s
    end
    return hv
end

"""
    h_vector(s::Seed; max_seeds::Int = 5000) -> Vector{Int}

Compute the h-vector via a full BFS exchange graph starting from `s`.
See `h_vector(::ExchangeGraph)` for details.
"""
function h_vector(s::Seed; max_seeds::Int = 5000)
    eg = exchange_graph(s; max_seeds)
    is_truncated(eg) && throw(InvalidArgument(
        "exchange graph was truncated at $max_seeds seeds; increase max_seeds or " *
        "use a finite-type seed"))
    return h_vector(eg)
end

# ─── Internal: k-combinations ─────────────────────────────────────────────────

# Generate all k-element subsets of `items` as Vector{T}.
function _combinations(items::Vector{T}, k::Int) where {T}
    n = length(items)
    k == 0 && return [T[]]
    k > n  && return Vector{T}[]
    result = Vector{T}[]
    combo  = Vector{T}(undef, k)
    _combinations_rec!(result, combo, items, 1, 1, k)
    return result
end

function _combinations_rec!(result, combo, items, start, depth, k)
    if depth > k
        push!(result, copy(combo))
        return
    end
    for i in start:length(items)-(k-depth)
        combo[depth] = items[i]
        _combinations_rec!(result, combo, items, i + 1, depth + 1, k)
    end
end
