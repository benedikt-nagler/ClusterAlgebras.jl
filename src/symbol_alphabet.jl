"""
    symbol_alphabet(s::Seed; max_seeds::Int = 1000) → Vector

Return all distinct mutable cluster variables appearing in the exchange
graph rooted at `s`.  For a finite-type cluster algebra this is the full
symbol alphabet used in scattering-amplitude computations.

Issues a warning when `max_seeds` is reached before the graph is exhausted
(infinite-type or undersampled finite-type case).
"""
function symbol_alphabet(s::Seed; max_seeds::Int = 1000)
    eg = exchange_graph(s; max_seeds = max_seeds)
    eg.truncated && @warn "exchange_graph truncated at $max_seeds seeds; alphabet may be incomplete"
    n_mut = s.quiver.n_mutable
    T = eltype(s.cluster)
    seen   = Set{T}()
    result = T[]
    for seed in eg.seeds
        for k in 1:n_mut
            v = seed[k]
            if v ∉ seen
                push!(seen, v)
                push!(result, v)
            end
        end
    end
    return result
end

"""
    cluster_adjacency_matrix(s::Seed; max_seeds::Int = 1000) → (Vector, BitMatrix)

Return `(alphabet, adj)` where `alphabet` is `symbol_alphabet(s)` and
`adj[i,j]` is `true` iff the `i`-th and `j`-th alphabet letters appear
together in some cluster.  The diagonal is always `true`.
"""
function cluster_adjacency_matrix(s::Seed; max_seeds::Int = 1000)
    eg = exchange_graph(s; max_seeds = max_seeds)
    eg.truncated && @warn "exchange_graph truncated at $max_seeds seeds; adjacency matrix may be incomplete"
    n_mut = s.quiver.n_mutable
    T = eltype(s.cluster)

    idx  = Dict{T, Int}()
    alpha = T[]
    for seed in eg.seeds
        for k in 1:n_mut
            v = seed[k]
            if !haskey(idx, v)
                idx[v] = length(alpha) + 1
                push!(alpha, v)
            end
        end
    end

    N   = length(alpha)
    adj = falses(N, N)
    for seed in eg.seeds
        is = [idx[seed[k]] for k in 1:n_mut]
        for i in is, j in is
            adj[i, j] = true
        end
    end

    return alpha, adj
end

"""
    cluster_adjacent(s::Seed, v, w; max_seeds::Int = 1000) → Bool

Return `true` if cluster variables `v` and `w` (ring elements from the
same cluster algebra) appear together in some cluster reachable from `s`.
"""
function cluster_adjacent(s::Seed, v, w; max_seeds::Int = 1000)
    eg    = exchange_graph(s; max_seeds = max_seeds)
    n_mut = s.quiver.n_mutable
    for seed in eg.seeds
        has_v = false
        has_w = false
        for k in 1:n_mut
            vk = seed[k]
            has_v |= (vk == v)
            has_w |= (vk == w)
        end
        has_v && has_w && return true
    end
    return false
end
