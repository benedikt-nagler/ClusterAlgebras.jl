# ─── Relabeling detection ─────────────────────────────────────────────────────

# Find the permutation σ with labels2[σ[i]] == labels1[i] and
# B2[σ[i], σ[j]] == B1[i, j] for all i, j, or return `nothing` if none exists.
# Assumes the labels within each vector are pairwise distinct (true for cluster
# variables and for y-variables along a mutation orbit).
function _find_relabeling(B1::Matrix{Int}, labels1::AbstractVector,
                          B2::Matrix{Int}, labels2::AbstractVector)
    n = length(labels1)
    length(labels2) == n || return nothing
    σ = Vector{Int}(undef, n)
    for i in 1:n
        j = findfirst(==(labels1[i]), labels2)
        j === nothing && return nothing
        σ[i] = j
    end
    allunique(σ) || return nothing
    for i in 1:n, j in 1:n
        B2[σ[i], σ[j]] == B1[i, j] || return nothing
    end
    return σ
end

# ─── mutation_period ──────────────────────────────────────────────────────────

"""
    mutation_period(s::Seed, seq::AbstractVector{Int}; max_period=100) → Union{Int, Nothing}

Smallest `m ≥ 1` such that applying the mutation sequence `seq` `m` times to `s`
returns to the initial seed **up to relabeling** — i.e. there is a permutation of
the vertices carrying the cluster to the initial cluster and the exchange matrix
to the initial exchange matrix. Returns `nothing` if no return is detected
within `max_period` applications (cluster variables of non-recurrent sequences
grow quickly, so keep the cutoff modest).

In finite type, periodic mutation sequences govern Y- and T-systems; see
[`y_system`](@ref) for the bipartite Zamolodchikov dynamics with period `h + 2`.

# Example
```julia
julia> mutation_period(Seed(Quiver(:A, 2)), [1, 2])   # pentagon recurrence
5
```
"""
function mutation_period(s::Seed, seq::AbstractVector{Int}; max_period::Int=100)
    max_period >= 1 ||
        throw(InvalidArgument("max_period must be ≥ 1, got $max_period"))
    isempty(seq) &&
        throw(InvalidArgument("mutation sequence must be non-empty"))
    B0, cluster0 = s.quiver.B, s.cluster
    current = s
    for m in 1:max_period
        current = mutate(current, seq)
        if _find_relabeling(B0, cluster0, current.quiver.B, current.cluster) !== nothing
            return m
        end
    end
    return nothing
end

# ─── Bipartite structure ──────────────────────────────────────────────────────

# Mutable vertices with no outgoing (resp. no incoming) arrows within the
# mutable block; arrows to frozen vertices are ignored. A vertex with no
# mutable neighbors is reported in both.
_sinks(q::Quiver)   = [k for k in 1:q.n_mutable if all(q.B[k, j] <= 0 for j in 1:q.n_mutable)]
_sources(q::Quiver) = [k for k in 1:q.n_mutable if all(q.B[j, k] <= 0 for j in 1:q.n_mutable)]

"""
    is_bipartite(q::Quiver) → Bool

Return `true` if every mutable vertex of `q` is a source or a sink within the
mutable block (arrows to frozen vertices are ignored). Bipartite quivers carry
the "all sinks, then all sources" Y-system dynamics; see [`y_system`](@ref).
"""
function is_bipartite(q::Quiver)
    sinks = Set(_sinks(q))
    sources = Set(_sources(q))
    return all(k in sinks || k in sources for k in 1:q.n_mutable)
end

# ─── y_system ─────────────────────────────────────────────────────────────────

"""
    y_system(s::Seed{PrincipalCoefficients}; max_period=100)
    y_system(q::Quiver; max_period=100)

Iterate the bipartite Y-system dynamics on `s`: alternately mutate at **all
sinks**, then at **all sources** (vertices within each part are pairwise
non-adjacent, so each half-step is order-independent). Requires
`is_bipartite(s.quiver)`; throws `InvalidArgument` otherwise.

Returns a NamedTuple `(period, y_trajectory)`:
- `period::Union{Int, Nothing}` — the number of **half-steps** after which the
  Y-seed `(B, y)` first returns to the initial one up to relabeling, or
  `nothing` if not detected within `max_period` half-steps.
- `y_trajectory` — the rational y-variables at each half-step, starting with
  the initial ones (`period + 1` entries when the period is found).

By Zamolodchikov periodicity, a bipartite quiver of finite type has
`period == h + 2` where `h` is the Coxeter number: `5` for `A₂` (the pentagon
recurrence), `6` for `A₃`, `8` for `D₄`.
"""
function y_system(s::Seed{PrincipalCoefficients}; max_period::Int=100)
    q = s.quiver
    is_bipartite(q) ||
        throw(InvalidArgument("y_system requires a bipartite quiver (every mutable " *
                              "vertex a source or a sink); got a non-bipartite quiver"))
    max_period >= 1 ||
        throw(InvalidArgument("max_period must be ≥ 1, got $max_period"))

    sinks = _sinks(q)
    parts = (sinks, setdiff(1:q.n_mutable, sinks))   # sinks first, then sources
    B0, y0 = q.B, y_variables(s)
    trajectory = [y0]
    current = s
    for t in 1:max_period
        part = parts[mod1(t, 2)]
        # After a half-step every mutated part swaps orientation, so the sink
        # set of the CURRENT quiver alternates between the two initial parts.
        current = mutate(current, part)
        y = y_variables(current)
        push!(trajectory, y)
        if _find_relabeling(B0, y0, current.quiver.B, y) !== nothing
            return (period = t, y_trajectory = trajectory)
        end
    end
    return (period = nothing, y_trajectory = trajectory)
end

y_system(q::Quiver; kwargs...) = y_system(extend(Seed(q)); kwargs...)
