# ─── Forks, pre-forks, and the forkless part ───────────────────────────────────
#
# Warkentin's forks [War14] and the mutation-graph structure around them.  A fork
# is a local, purely combinatorial pattern, yet it controls the global shape of a
# mutation class: past a fork lie only forks (Lemma 2.8), and a connected quiver
# is mutation-infinite exactly when its class contains one (Thm 3.2).  Deleting
# those subtrees leaves the *forkless part*, which is finite for many
# mutation-infinite quivers and is a hereditary, mutation-invariant property
# [Erv23].
#
# Convention ledger (choices that are ours, not forced by the mathematics):
#   * The forkless part is computed in the LABELED mutation graph, as in
#     [Erv23] Def 3.1 - vertex-relabeled copies are distinct members.  Pinned by
#     the mutation-invariance test in test_forks.jl (the set is identical from
#     every starting point in it).
#   * Forks and pre-forks are boundary, not members: they are traversed (at the
#     point of return only) but never returned.  Pinned by Def 3.1, which counts
#     non-forks.
#   * A budget exhaustion is reported as `nothing`, never as "infinite".  No
#     finite search can certify an infinite forkless part, the same asymmetry
#     `mgs_search` carries for maximal green sequences.
#
# [War14] M. Warkentin. *Exchange graphs via quiver mutation.* Dissertation,
#         TU Chemnitz 2014.  (forks, Def 2.1; pre-forks and stopovers, Def 3.7;
#         the Tree Lemma 2.8; mutation-infinite iff a fork, Thm 3.2)
# [Erv23] J. Ervin. *New hereditary and mutation-invariant properties arising
#         from forks.* arXiv:2306.07502.  (forkless part, Def 3.1; hereditarity,
#         Thm 3.5 / Cor 3.6; pre-forkless part, Def 5.1; Cor 5.5)

"""
    ForklessPart

Result of a [`forkless_part`](@ref) or [`pre_forkless_part`](@ref) search.

`quivers` holds the members - the non-forks (resp. non-forks and non-pre-forks)
of the labeled mutation class, each a distinct exchange matrix.  `adj[i][k]` is
the index of the quiver obtained by mutating `quivers[i]` at vertex `k`, or `0`
when that mutation leaves the part.  `kind` is `:fork` or `:pre_fork`, recording
which boundary was used, and `truncated` is `true` when the search stopped at
`max_quivers` instead of closing.
"""
struct ForklessPart
    quivers   :: Vector{Quiver}
    adj       :: Vector{Vector{Int}}
    kind      :: Symbol
    truncated :: Bool
end

Base.length(fp::ForklessPart)           = length(fp.quivers)
Base.getindex(fp::ForklessPart, i::Int) = fp.quivers[i]
is_truncated(fp::ForklessPart)          = fp.truncated

function Base.show(io::IO, fp::ForklessPart)
    name = fp.kind === :fork ? "forkless part" : "pre-forkless part"
    state = fp.truncated ? "truncated at" : "closed with"
    print(io, "$name ($state $(length(fp.quivers)) quivers)")
end

# Fork theory is stated for quivers: skew-symmetric, no frozen vertices.  Arrow
# counts between a frozen and a mutable vertex are not mutation data, and the
# valued (skew-symmetrizable) case has no published fork definition.
function _require_plain_quiver(q::Quiver, fname::String)
    q.n_frozen == 0 ||
        throw(InvalidArgument("$fname is defined for quivers without frozen vertices, got $(q.n_frozen)"))
    all(isone, q.d) ||
        throw(InvalidArgument("$fname is defined for skew-symmetric quivers, got symmetrizers $(q.d)"))
    return nothing
end

"""
    is_abundant(q::Quiver) → Bool

Whether every pair of vertices of `q` is joined by at least two arrows
(Warkentin's *abundant*, also called 2-complete).

# Example
```jldoctest
julia> using ClusterAlgebras

julia> is_abundant(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
true

julia> is_abundant(Quiver(:A, 3))
false
```
"""
function is_abundant(q::Quiver)
    n = q.n_mutable
    for i in 1:n, j in (i + 1):n
        abs(q.B[i, j]) >= 2 || return false
    end
    return true
end

# The induced subquiver on vertex set `S`, as a plain quiver.
_induced(q::Quiver, S::AbstractVector{Int}) = Quiver(q.B[S, S], length(S))

"""
    point_of_return(q::Quiver) → Union{Int, Nothing}

The point of return of `q` when `q` is a fork, and `nothing` otherwise.

A fork (Warkentin 2014, Definition 2.1) is an abundant quiver that is not
acyclic and has a vertex `r` such that every arrow between the in-neighbours
and out-neighbours of `r` runs from an out-neighbour `j` back to an in-neighbour
`i`, with `B[j,i]` strictly larger than both `B[i,r]` and `B[r,j]`, and such
that the in-neighbours and out-neighbours of `r` each induce an acyclic
subquiver.  The point of return is unique.

Mutating a fork anywhere except at `r` gives another fork with strictly more
arrows, so `r` is the only direction in which the mutation class does not grow.

# Example
```jldoctest
julia> using ClusterAlgebras

julia> point_of_return(Quiver([0 2 -3; -2 0 2; 3 -2 0]))
2

julia> point_of_return(Quiver([0 2 -2; -2 0 2; 2 -2 0]))  # the Markov quiver
```
"""
function point_of_return(q::Quiver)
    _require_plain_quiver(q, "point_of_return")
    n = q.n_mutable
    n >= 3 || return nothing
    is_abundant(q) || return nothing
    _is_acyclic(q) && return nothing

    B = q.B
    for r in 1:n
        into = [i for i in 1:n if i != r && B[i, r] > 0]
        outof = [j for j in 1:n if j != r && B[r, j] > 0]
        ok = true
        for i in into, j in outof
            if !(B[j, i] > B[i, r] && B[j, i] > B[r, j])
                ok = false
                break
            end
        end
        ok || continue
        _is_acyclic(_induced(q, into)) || continue
        _is_acyclic(_induced(q, outof)) || continue
        return r
    end
    return nothing
end

"""
    is_fork(q::Quiver) → Bool

Whether `q` is a fork in the sense of Warkentin 2014 - see
[`point_of_return`](@ref), which also names the vertex.

# Example
```jldoctest
julia> using ClusterAlgebras

julia> is_fork(Quiver([0 2 -3; -2 0 2; 3 -2 0]))
true

julia> is_fork(Quiver(:A, 3))
false
```
"""
is_fork(q::Quiver) = point_of_return(q) !== nothing

"""
    pre_fork_vertices(q::Quiver) → Union{NamedTuple, Nothing}

The distinguished vertices `(r, k, k′)` when `q` is a pre-fork, and `nothing`
otherwise.

A pre-fork (Warkentin 2014, Definition 3.7) is a quiver with two vertices `k`
and `k′` such that deleting either one leaves a fork, both with the same point
of return `r`, and such that no third vertex is a *stopover*: no `i` carries
arrows `k′ → i → k` or `k → i → k′`.  Pre-forks generalise forks and bound the
mutation class the same way, which is what makes the pre-forkless part a
coarser, and more often finite, invariant than the forkless part.

`k < k′` in the returned tuple; a quiver can be a pre-fork on several pairs and
the first in that order is reported.

# Example
```jldoctest
julia> using ClusterAlgebras

julia> pre_fork_vertices(Quiver(:A, 3))
```
"""
function pre_fork_vertices(q::Quiver)
    _require_plain_quiver(q, "pre_fork_vertices")
    n = q.n_mutable
    n >= 4 || return nothing
    B = q.B
    rest = Vector{Int}(undef, n - 1)

    for k in 1:n, kp in (k + 1):n
        # Both deletions must be forks with a common point of return.  Vertices
        # are renumbered by deletion, so map the returned index back.
        Sk = [v for v in 1:n if v != k]
        Skp = [v for v in 1:n if v != kp]
        rk = point_of_return(_induced(q, Sk))
        rk === nothing && continue
        rkp = point_of_return(_induced(q, Skp))
        rkp === nothing && continue
        Sk[rk] == Skp[rkp] || continue
        r = Sk[rk]

        # Stopovers: k′ → i → k or k → i → k′ for some third vertex i.
        stopover = false
        for i in 1:n
            (i == k || i == kp) && continue
            if (B[kp, i] > 0 && B[i, k] > 0) || (B[i, kp] > 0 && B[k, i] > 0)
                stopover = true
                break
            end
        end
        stopover && continue
        return (r = r, k = k, k′ = kp)
    end
    return nothing
end

"""
    is_pre_fork(q::Quiver) → Bool

Whether `q` is a pre-fork - see [`pre_fork_vertices`](@ref), which also names
the distinguished vertices.
"""
is_pre_fork(q::Quiver) = pre_fork_vertices(q) !== nothing

# The shared BFS.  Non-boundary quivers are expanded at every vertex; a boundary
# quiver is expanded only at its point of return, which by Warkentin's Tree
# Lemma 2.8 (and Ervin Cor 4.25 for pre-forks) is the single edge leading out of
# the subtree it heads - everything else behind it is boundary again, so no
# member of the part is lost by stopping there.
function _part_search(q::Quiver, kind::Symbol, max_quivers::Int)
    _require_plain_quiver(q, kind === :fork ? "forkless_part" : "pre_forkless_part")
    max_quivers >= 1 ||
        throw(InvalidArgument("max_quivers must be at least 1, got $max_quivers"))

    # boundary(qi) → the point of return to escape through, or nothing.
    function boundary(qi::Quiver)
        r = point_of_return(qi)
        r === nothing || return r
        kind === :fork && return nothing
        pf = pre_fork_vertices(qi)
        return pf === nothing ? nothing : pf.r
    end

    index = Dict{Matrix{Int},Int}()      # members only
    members = Quiver[]
    adj = Vector{Vector{Int}}()
    queue = Tuple{Quiver,Union{Int,Nothing}}[]

    function visit!(qi, esc)
        if esc === nothing
            haskey(index, qi.B) && return index[qi.B]
            push!(members, qi)
            push!(adj, zeros(Int, qi.n_mutable))
            index[qi.B] = length(members)
            push!(queue, (qi, nothing))
            return length(members)
        end
        push!(queue, (qi, esc))
        return 0
    end

    seen_boundary = Set{Matrix{Int}}()
    visit!(q, boundary(q))
    head = 1
    while head <= length(queue)
        qi, esc = queue[head]; head += 1
        if esc !== nothing
            qi.B in seen_boundary && continue
            push!(seen_boundary, qi.B)
        end
        ks = esc === nothing ? (1:qi.n_mutable) : (esc:esc)
        for k in ks
            qk = try
                mutate(qi,k)
            catch e
                e isa OverflowError || rethrow()
                @warn "forkless search reached integer storage; returning an incomplete result"
                return ForklessPart(members,adj,kind,true)
            end
            e = boundary(qk)
            if e === nothing && !haskey(index, qk.B) && length(members) >= max_quivers
                return ForklessPart(members, adj, kind, true)
            end
            j = visit!(qk, e)
            esc === nothing && (adj[index[qi.B]][k] = j)
        end
    end
    return ForklessPart(members, adj, kind, false)
end

"""
    forkless_part(q::Quiver; max_quivers::Int = 10_000) → ForklessPart

The forkless part of the labeled mutation class of `q`: every quiver
mutation-equivalent to `q` that is not a fork (Ervin 2023, Definition 3.1).

Forks are boundary rather than members.  Warkentin's Tree Lemma makes the search
exact: everything behind a fork, except through its point of return, is again a
fork, so the returned part is complete whenever the search closes.

The search closes for a quiver with a finite forkless part - which includes every
mutation-finite quiver, since a mutation-finite connected quiver is
mutation-equivalent to no fork at all - and runs to `max_quivers` otherwise, when
the result is a truncated piece and `is_truncated` is `true`.
The same holds when an exact mutation exceeds the stored integer type.  Use
[`has_finite_forkless_part`](@ref) for the yes/undecided verdict.

# Example
```jldoctest
julia> using ClusterAlgebras

julia> fp = forkless_part(Quiver(:A, 3));

julia> length(fp), is_truncated(fp)
(14, false)
```
"""
forkless_part(q::Quiver; max_quivers::Int = 10_000) =
    _part_search(q, :fork, max_quivers)

"""
    pre_forkless_part(q::Quiver; max_quivers::Int = 10_000) → ForklessPart

The pre-forkless part of the labeled mutation class of `q`: every quiver
mutation-equivalent to `q` that is neither a fork nor a pre-fork (Ervin 2023,
Definition 5.1).

Pre-forks bound larger regions than forks, so this part is contained in the
forkless part and closes strictly more often - Warkentin's `(a,b)` box quivers
with `a, b ≥ 2` have a finite pre-forkless part and an infinite forkless part.
See [`forkless_part`](@ref) for the search's semantics.
"""
pre_forkless_part(q::Quiver; max_quivers::Int = 10_000) =
    _part_search(q, :pre_fork, max_quivers)

"""
    has_finite_forkless_part(q::Quiver; max_quivers::Int = 10_000)
        → Union{Bool, Nothing}

`true` when the forkless part of `q` closes within `max_quivers`, `nothing`
when the budget runs out first.

Never `false`: no finite search certifies an infinite forkless part, so the
undecided verdict is reported as `nothing` rather than folded into `false`.
Raise `max_quivers` to push the boundary.

Having a finite forkless part is hereditary and invariant under mutation
(Ervin 2023, Corollary 3.6), so a `true` here also settles every full subquiver,
and a `nothing` says nothing about them.

# Example
```jldoctest
julia> using ClusterAlgebras

julia> has_finite_forkless_part(Quiver(:A, 3))
true
```
"""
function has_finite_forkless_part(q::Quiver; max_quivers::Int = 10_000)
    fp = forkless_part(q; max_quivers)
    return fp.truncated ? nothing : true
end

"""
    has_finite_pre_forkless_part(q::Quiver; max_quivers::Int = 10_000)
        → Union{Bool, Nothing}

`true` when the pre-forkless part of `q` closes within `max_quivers`, `nothing`
when the budget runs out first, with the same asymmetry as
[`has_finite_forkless_part`](@ref).

Having a finite pre-forkless part is hereditary and invariant under mutation
(Ervin 2023, Corollary 5.5), and is implied by a finite forkless part.
"""
function has_finite_pre_forkless_part(q::Quiver; max_quivers::Int = 10_000)
    fp = pre_forkless_part(q; max_quivers)
    return fp.truncated ? nothing : true
end
