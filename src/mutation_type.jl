# Identification of the mutation class a quiver belongs to.
#
# `cartan_type` already names the finite types; this file names the rest of what
# is nameable: affine (simply-laced), rank 2, and the exceptional mutation-finite
# classes that come from neither a Dynkin diagram nor a surface.
#
# The naming follows SageMath's `QuiverMutationType`, a triple
# `(letter, rank, twist)`, so the answer round-trips through the interop layer.
# Note that `rank` is the rank of the *Dynkin datum*, not the number of vertices:
# an affine type has one vertex more, an elliptic type two more.

"""
    MutationType(letter, rank, twist = nothing)

The mutation class of a quiver, named as SageMath's `QuiverMutationType` names
it: a `letter`, a `rank`, and a `twist`.

- `twist === nothing` is a finite type (`MutationType(:A, 3)` is `A3`), the
  Derksen-Owen exceptionals `X6`/`X7`, or the wild rank-2 class `R2` whose
  `rank` is the pair of edge labels.
- `twist == 1` is an affine type (`D4^(1)`); for `Ã(p, q)` the `rank` is the
  pair `(p, q)` of arrow counts around the cycle, `p <= q`.
- `twist == (1, 1)` is an elliptic type (`E6^(1,1)`).

`rank` counts the Dynkin datum, so the quiver has `rank` vertices in the finite
case, `rank + 1` in the affine case (`sum(rank)` for `Ã(p, q)`), and `rank + 2`
in the elliptic case.

See [`mutation_type`](@ref).
"""
struct MutationType
    letter :: Symbol
    rank   :: Union{Int, Tuple{Int, Int}}
    twist  :: Union{Nothing, Int, Tuple{Int, Int}}
end

MutationType(letter::Symbol, rank) = MutationType(letter, rank, nothing)

Base.:(==)(a::MutationType, b::MutationType) =
    a.letter == b.letter && a.rank == b.rank && a.twist == b.twist

Base.hash(mt::MutationType, h::UInt) =
    hash(mt.twist, hash(mt.rank, hash(mt.letter, hash(:MutationType, h))))

function _mutation_type_string(mt::MutationType)
    body = mt.rank isa Tuple ? "$(mt.letter)($(mt.rank[1]),$(mt.rank[2]))" :
                               "$(mt.letter)$(mt.rank)"
    mt.twist === nothing && return body
    mt.twist isa Tuple && return body * "^($(mt.twist[1]),$(mt.twist[2]))"
    return body * "^($(mt.twist))"
end

Base.show(io::IO, mt::MutationType) = print(io, _mutation_type_string(mt))

Base.isless(a::MutationType, b::MutationType) =
    _mutation_type_string(a) < _mutation_type_string(b)

# The undirected neighbour lists of the mutable block.
function _undirected_neighbours(B::AbstractMatrix, n::Int)
    return [[j for j in 1:n if j != i && (B[i, j] != 0 || B[j, i] != 0)] for i in 1:n]
end

_is_simply_laced(q::Quiver) =
    all(isone, q.d) && all(abs(q.B[i, j]) <= 1 for i in 1:q.n_mutable, j in 1:q.n_mutable)

_is_skew_symmetric(q::Quiver) = all(isone, q.d) && q.B == -transpose(q.B)

"""
    _affine_type(rep::Quiver) → Union{MutationType, Nothing}

Name the affine type of an **acyclic** representative already known to be of
affine type. Simply-laced only: the twisted and non-simply-laced affine diagrams
(SageMath's `BB`/`CC`/`BC`/`BD`/`CD` family) are not identified and give
`nothing`.
"""
function _affine_type(rep::Quiver)
    n = rep.n_mutable
    _is_simply_laced(rep) || return nothing
    adj = _undirected_neighbours(rep.B, n)
    degrees = length.(adj)

    # A cycle is affine type Ã(p, q), with p and q the arrow counts in the two
    # directions of a traversal. The acyclic representative is never cyclically
    # oriented, so both are nonzero.
    if all(==(2), degrees)
        order = [1, adj[1][1]]
        while length(order) < n
            i = order[end]
            push!(order, adj[i][1] == order[end - 1] ? adj[i][2] : adj[i][1])
        end
        forward = 0
        for k in 1:n
            i, j = order[k], order[mod1(k + 1, n)]
            rep.B[i, j] > 0 && (forward += 1)
        end
        p, q = minmax(forward, n - forward)
        return MutationType(:A, (p, q), 1)
    end

    branch = count(>(2), degrees)
    if branch == 1
        v = findfirst(>(2), degrees)
        degrees[v] == 4 && n == 5 && return MutationType(:D, 4, 1)
        degrees[v] == 3 && 7 <= n <= 9 && return MutationType(:E, n - 1, 1)
        return nothing
    elseif branch == 2
        all(d -> d == 3, degrees[degrees .> 2]) || return nothing
        return MutationType(:D, n - 1, 1)
    end
    return nothing
end

# The exceptional mutation-finite skew-symmetric classes that are neither
# Dynkin nor from a surface: the Derksen-Owen quivers X6, X7 and the elliptic
# types E6, E7, E8 with twist (1,1). Edge lists (source, target, arrows) are
# SageMath's, translated from its 0-based vertices.
const _EXCEPTIONAL_EDGES = Dict{Int, Tuple{MutationType, Vector{NTuple{3, Int}}}}(
    6  => (MutationType(:X, 6),
           [(0, 1, 2), (1, 2, 1), (2, 0, 1), (2, 3, 1), (3, 4, 2), (4, 2, 1), (2, 5, 1)]),
    7  => (MutationType(:X, 7),
           [(0, 1, 2), (1, 2, 1), (2, 0, 1), (2, 3, 1), (3, 4, 2), (4, 2, 1), (2, 5, 1),
            (5, 6, 2), (6, 2, 1)]),
    8  => (MutationType(:E, 6, (1, 1)),
           [(0, 1, 1), (1, 2, 1), (3, 2, 1), (3, 4, 1), (5, 6, 1), (6, 7, 1), (5, 1, 1),
            (2, 5, 2), (5, 3, 1), (6, 2, 1)]),
    9  => (MutationType(:E, 7, (1, 1)),
           [(1, 0, 1), (1, 2, 1), (2, 3, 1), (4, 3, 1), (4, 5, 1), (6, 5, 1), (7, 8, 1),
            (3, 7, 2), (7, 2, 1), (7, 4, 1), (8, 3, 1)]),
    10 => (MutationType(:E, 8, (1, 1)),
           [(0, 1, 1), (1, 9, 1), (3, 9, 1), (3, 4, 1), (2, 8, 1), (2, 1, 1), (9, 2, 2),
            (2, 3, 1), (8, 9, 1), (5, 4, 1), (5, 6, 1), (7, 6, 1)]),
)

"""
    _exceptional_quiver(n) → Union{Quiver, Nothing}

The standard representative of the exceptional mutation-finite class on `n`
vertices; there is at most one for each `n`.
"""
function _exceptional_quiver(n::Int)
    haskey(_EXCEPTIONAL_EDGES, n) || return nothing
    B = zeros(Int, n, n)
    for (i, j, w) in _EXCEPTIONAL_EDGES[n][2]
        B[i + 1, j + 1], B[j + 1, i + 1] = w, -w
    end
    return Quiver(B)
end

# Canonical forms of each exceptional class, computed once on demand.
const _EXCEPTIONAL_CLASS_CACHE = Dict{Int, Union{Set{Matrix{Int}}, Nothing}}()

function _exceptional_class(n::Int, max_quivers::Int)
    get!(_EXCEPTIONAL_CLASS_CACHE, n) do
        rep = _exceptional_quiver(n)
        rep === nothing && return nothing
        mc = mutation_class(rep; max_quivers = max_quivers, up_to_isomorphism = true)
        is_truncated(mc) && return nothing
        return Set(canonical_form(qi).B for qi in mc.quivers)
    end
end

"""
    mutation_type(q::Quiver; max_quivers::Int = 20_000) → Union{MutationType, Nothing}

Identify the mutation class of the **connected** quiver `q` (its frozen vertices
are ignored), or `nothing` when the class has no name here.

Named classes are the finite types (via [`cartan_type`](@ref)), the simply-laced
affine types, the rank-2 classes, and the exceptional mutation-finite classes
`X6`, `X7` and elliptic `E6`, `E7`, `E8` of twist `(1, 1)`.

`nothing` means *not identified*, which covers three genuinely different
situations: `q` is mutation-infinite (no type exists), `q` is mutation-finite of
surface type (a name exists, but block decomposition is not implemented yet), or
`q` is an affine type that is not simply laced. Use
[`is_mutation_finite`](@ref) and [`is_affine_type`](@ref) to tell them apart.

Recognizing an exceptional class enumerates its mutation class once, which is
cached; `max_quivers` bounds that enumeration, and exceeding it gives `nothing`.

```jldoctest
julia> mutation_type(Quiver(:A, 3))
A3

julia> mutation_type(Quiver([0 2; -2 0]))
A(1,1)^(1)
```
"""
function mutation_type(q::Quiver; max_quivers::Int = 20_000)
    n = q.n_mutable
    n >= 1 || throw(InvalidArgument("mutation_type requires n_mutable ≥ 1"))
    length(_undirected_components(q.B, n)) == 1 || throw(InvalidArgument(
        "mutation_type is only defined for connected quivers; this one has " *
        "$(length(_undirected_components(q.B, n))) components. Use mutation_types(q)."))

    qm = Quiver(q.B[1:n, 1:n], n, q.d)

    n == 1 && return MutationType(:A, 1)

    if n == 2
        p, r = abs(qm.B[1, 2]), abs(qm.B[2, 1])
        lo, hi = minmax(p, r)
        product = p * r
        product == 1 && return MutationType(:A, 2)
        product == 2 && return MutationType(:B, 2)
        product == 3 && return MutationType(:G, 2)
        product == 4 && return lo == 2 ? MutationType(:A, (1, 1), 1) :
                                         MutationType(:BC, 1, 1)
        return MutationType(:R2, (lo, hi))
    end

    # A skew-symmetric quiver of rank ≥ 3 with an entry beyond 2 is mutation
    # infinite, which rules out every named class before any search.
    if _is_skew_symmetric(qm) && maximum(abs, qm.B) > 2
        return nothing
    end

    rep = _acyclic_representative(qm)
    if rep isa Quiver
        S = _symmetrized_cartan(rep)
        minors = [_det_int(S[1:k, 1:k]) for k in 1:n]
        all(>(0), minors) && return MutationType(_connected_cartan_type(rep)...)
        all(>(0), minors[1:(n - 1)]) && minors[n] == 0 && return _affine_type(rep)
        return nothing
    end

    # No acyclic representative: the only named classes left are exceptional.
    _is_skew_symmetric(qm) || return nothing
    class = _exceptional_class(n, max_quivers)
    class === nothing && return nothing
    return canonical_form(qm).B in class ? _EXCEPTIONAL_EDGES[n][1] : nothing
end

"""
    mutation_types(q::Quiver; max_quivers::Int = 20_000) → Vector

Identify each connected component of `q` with [`mutation_type`](@ref), sorted.
Components that cannot be identified contribute `nothing`.
"""
function mutation_types(q::Quiver; max_quivers::Int = 20_000)
    n = q.n_mutable
    n == 0 && return Union{MutationType, Nothing}[]
    comps = _undirected_components(q.B, n)
    types = Union{MutationType, Nothing}[
        mutation_type(Quiver(q.B[c, c], length(c), q.d[c]); max_quivers = max_quivers)
        for c in comps]
    return sort(types; by = t -> t === nothing ? "" : _mutation_type_string(t))
end
