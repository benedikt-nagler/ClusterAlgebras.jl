# ─── MutationClass ─────────────────────────────────────────────────────────────
# Quiver-level BFS: tracks distinct exchange matrices reachable by mutation.
#
# NOTE(design): MutationClass (quiver-level) and ExchangeGraph (seed-level) are
# kept as separate types because they have different deduplication semantics and
# return different element types.  Unifying them into a single parameterised
# struct is worth revisiting once the API stabilises.

"""
    MutationClass

Result of a quiver-level BFS starting from a `Quiver`.  Vertices are distinct
exchange matrices reachable by mutation; `adj[i][k]` is the index of the quiver
obtained by mutating `quivers[i]` at mutable vertex `k`.

When produced with `up_to_isomorphism = true`, vertices are quiver *isomorphism
classes* rather than labeled exchange matrices: each entry of `quivers` is the
first-seen representative of its class, and `adj` is the *unlabeled* mutation
graph built from those stored representatives (mutating a representative at `k`
may land on a relabeling of an already-seen class).

`truncated` is `true` when the BFS was stopped early because the number of
distinct vertices reached `max_quivers`.
"""
struct MutationClass
    quivers   :: Vector{Quiver}
    adj       :: Vector{Vector{Int}}  # adj[i][k] = index of μ_k(quivers[i])
    truncated :: Bool
end

Base.length(mc::MutationClass)           = length(mc.quivers)
Base.getindex(mc::MutationClass, i::Int) = mc.quivers[i]
is_truncated(mc::MutationClass)          = mc.truncated

"""
    mutation_class(q::Quiver; max_quivers::Int = 1000,
                   up_to_isomorphism::Bool = false) → MutationClass

BFS over all quivers mutation-equivalent to `q`.  Stops early
(`truncated = true`) once `max_quivers` distinct vertices have been found.

By default vertices are deduplicated by exact exchange-matrix equality, so
vertex-relabeled copies of the same quiver count as distinct vertices.  With
`up_to_isomorphism = true` the dedup key is the isomorphism invariant
`(canonical_form(qi).B, qi.d)` (matching [`is_isomorphic`](@ref)), so each
isomorphism class is stored once - this gives literature-consistent class sizes.
The stored representative is the first-seen labeled quiver, so `mc[1] == q`, and
`adj` is then the unlabeled mutation graph (see [`MutationClass`](@ref)).

For mutation-finite types (e.g. all finite Dynkin types) the BFS terminates
naturally.  For mutation-infinite types the `max_quivers` cutoff prevents
the search from hanging.
"""
function mutation_class(q::Quiver; max_quivers::Int = 1000,
                        up_to_isomorphism::Bool = false)
    _key(qi) = up_to_isomorphism ? (canonical_form(qi).B, qi.d) : qi.B

    seen    = Dict{Any, Int}()
    quivers = Quiver[]
    adj     = Vector{Vector{Int}}()

    function _add!(qi)
        push!(quivers, qi)
        push!(adj, Int[])
        j = length(quivers)
        seen[_key(qi)] = j
        j
    end

    _add!(q)
    queue = [1]
    head  = 1

    while head <= length(queue)
        i  = queue[head]; head += 1
        qi = quivers[i]

        for k in 1:qi.n_mutable
            qk  = mutate(qi, k)
            key = _key(qk)
            if !haskey(seen, key)
                length(quivers) >= max_quivers &&
                    return MutationClass(quivers, adj, true)
                j = _add!(qk)
                push!(queue, j)
            end
            push!(adj[i], seen[key])
        end
    end

    MutationClass(quivers, adj, false)
end

# ─── ExchangeGraph ─────────────────────────────────────────────────────────────

"""
    ExchangeGraph

Result of a seed-level BFS starting from a `Seed`.  Vertices are distinct seeds
(identified by their ordered d-vector tuple); `adj[i][k]` is the index of the
seed obtained by mutating `seeds[i]` at mutable vertex `k`.

`truncated` is `true` when the BFS was stopped early because the number of
distinct seeds reached `max_seeds`.

For finite-type cluster algebras the BFS terminates and `length(eg)` equals the
number of clusters.  Deduplication uses ordered d-vector tuples, which correctly
separates all seeds for finite types (d-vector conjecture) and is safe up to the
cutoff for infinite types.
"""
struct ExchangeGraph{S <: AbstractSeed}
    seeds     :: Vector{S}
    adj       :: Vector{Vector{Int}}  # adj[i][k] = index of μ_k(seeds[i])
    truncated :: Bool
end

Base.length(eg::ExchangeGraph)           = length(eg.seeds)
Base.getindex(eg::ExchangeGraph, i::Int) = eg.seeds[i]
is_truncated(eg::ExchangeGraph)          = eg.truncated

# Identify a seed by its ordered tuple of d-vectors.
# Using integer vectors rather than the cluster variable expressions themselves
# avoids expensive rational-function comparisons while remaining correct for all
# finite types.
function _seed_key(s::Seed)
    dvecs = [denominator_vector(s, k) for k in 1:length(s.cluster)]
    Tuple(sort(dvecs))   # unordered: identifies clusters regardless of variable order
end

"""
    exchange_graph(s::Seed; max_seeds::Int = 1000) → ExchangeGraph

BFS over all seeds reachable from `s` by single mutations, deduplicating by
ordered d-vector tuple.  Stops early (`truncated = true`) once `max_seeds`
distinct seeds have been found.

For finite-type cluster algebras `length(exchange_graph(s))` equals the number
of clusters (e.g. 5 for A₂, 14 for A₃, 50 for D₄).
"""
function exchange_graph(s::S; max_seeds::Int = 1000) where {S <: Seed}
    seen  = Dict{Any, Int}()
    seeds = S[]
    adj   = Vector{Vector{Int}}()

    function _add!(si)
        push!(seeds, si)
        push!(adj, Int[])
        j = length(seeds)
        seen[_seed_key(si)] = j
        j
    end

    _add!(s)
    queue = [1]
    head  = 1

    while head <= length(queue)
        i  = queue[head]; head += 1
        si = seeds[i]

        for k in 1:si.quiver.n_mutable
            sk  = mutate(si, k)
            key = _seed_key(sk)
            if !haskey(seen, key)
                length(seeds) >= max_seeds &&
                    return ExchangeGraph(seeds, adj, true)
                j = _add!(sk)
                push!(queue, j)
            end
            push!(adj[i], seen[key])
        end
    end

    ExchangeGraph(seeds, adj, false)
end

# ─── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, mc::MutationClass)
    n = mc.quivers[1].n_mutable
    suffix = mc.truncated ? " (truncated)" : ""
    print(io, "MutationClass: $(length(mc)) quivers, rank $n$suffix")
end

function Base.show(io::IO, eg::ExchangeGraph)
    n = eg.seeds[1].quiver.n_mutable
    suffix = eg.truncated ? " (truncated)" : ""
    print(io, "ExchangeGraph: $(length(eg)) seeds, rank $n$suffix")
end
