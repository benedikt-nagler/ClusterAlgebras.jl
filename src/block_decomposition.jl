# Block decomposition of skew-symmetric quivers.
#
# [FeSTu12] Felikson-Shapiro-Tumarkin: a connected skew-symmetric quiver is the adjacency
# quiver of an ideal (tagged) triangulation of a bordered marked surface exactly
# when it can be assembled from the seven blocks below by identifying outlets of
# *different* blocks along a partial matching, arrows between identified vertices
# adding up (so opposite arrows cancel and equal ones give weight 2).
#
# Convention ledger (the block table is a *transcription*, so it is pinned twice):
#   • The blocks come from Figure 2.1 of [FeSTu12] (A. Felikson, M. Shapiro,
#     P. Tumarkin, "Skew-symmetric cluster algebras of finite mutation type",
#     JEMS 14 (2012) 1135-1180, arXiv:0811.1703); outlets are the white vertices
#     there, dead ends the black ones. The block criterion itself is [FST08] §13.
#   • Two independent checks on the transcription: block IV is itself of type D4 and
#     block V of type D4^(1), as the literature records (`test_block_decomposition.jl`).
#
# This is the missing half of the surface story. `mutation_type` names a quiver
# only when its class carries a Dynkin-style label, so it answers `nothing` both
# for a mutation-infinite quiver and for a surface quiver with no such label (a
# once-punctured torus, say). `is_block_decomposable` separates the two.

"""
    Block

One of the seven Felikson-Shapiro-Tumarkin blocks: a small quiver whose vertices
are split into `outlets`, which may be identified with an outlet of one other
block, and dead ends, which may not.

`label` is the name the classification gives the block (`:I`, `:II`, `:IIIa`,
`:IIIb`, `:IV`, `:V`, and `:point` for the single vertex).

See [`FST_BLOCKS`](@ref) and [`block_decomposition`](@ref).
"""
struct Block
    label   :: Symbol
    B       :: Matrix{Int}
    outlets :: Vector{Int}
end

nvertices(b::Block) = size(b.B, 1)
is_outlet(b::Block, i::Int) = i in b.outlets

Base.show(io::IO, b::Block) = print(io, "Block ", b.label, " (", nvertices(b), " vertices)")

# Build a block from an arrow list, so the table below reads like the figure.
function _block(label::Symbol, n::Int, arrows::Vector{Tuple{Int, Int}}, outlets::Vector{Int})
    B = zeros(Int, n, n)
    for (i, j) in arrows
        B[i, j] += 1
        B[j, i] -= 1
    end
    return Block(label, B, sort(outlets))
end

"""
    FST_BLOCKS

The seven blocks of the Felikson-Shapiro-Tumarkin classification, in the order
`:point`, `:I`, `:II`, `:IIIa`, `:IIIb`, `:IV`, `:V`.

Geometrically each block is a triangulated piece of a surface: `:point`, `:I`
and `:II` are a triangle with two, one and no sides on the boundary, and the
remaining four are the pieces carrying a puncture.
"""
const FST_BLOCKS = [
    _block(:point, 1, Tuple{Int, Int}[],                         [1]),
    _block(:I,     2, [(1, 2)],                                  [1, 2]),
    _block(:II,    3, [(1, 2), (2, 3), (3, 1)],                  [1, 2, 3]),
    # outlet 1, dead ends 2 and 3, both arrows pointing in
    _block(:IIIa,  3, [(2, 1), (3, 1)],                          [1]),
    # the same with both arrows pointing out
    _block(:IIIb,  3, [(1, 2), (1, 3)],                          [1]),
    # outlets 1, 2; two oriented triangles glued along the arrow 2 -> 1
    _block(:IV,    4, [(2, 1), (1, 3), (3, 2), (1, 4), (4, 2)],  [1, 2]),
    # outlet 1 only; dead ends 2, 3, 4, 5
    _block(:V,     5, [(2, 1), (1, 3), (3, 2), (5, 1), (1, 4),
                       (4, 5), (4, 2), (3, 5)],                  [1]),
]

"""
    BlockDecomposition

A witness that a quiver is block-decomposable: the `blocks` it is assembled
from, and for each one the `placements` vector sending block vertex `i` to the
quiver vertex it is identified with.

Reassembling is exact - summing each block's arrows over its placement returns
the quiver's B-matrix, cancellations included. See [`block_decomposition`](@ref).
"""
struct BlockDecomposition
    blocks     :: Vector{Block}
    placements :: Vector{Vector{Int}}
end

Base.length(d::BlockDecomposition) = length(d.blocks)

function Base.show(io::IO, d::BlockDecomposition)
    print(io, "BlockDecomposition(")
    print(io, join((string(b.label) for b in d.blocks), " + "))
    print(io, ")")
end

"""
    reassemble(d::BlockDecomposition, n::Int) → Matrix{Int}

The B-matrix obtained by gluing the blocks of `d` on `n` vertices. Arrows
between identified vertices add up, so opposite ones cancel and equal ones give
weight 2. This is the inverse of [`block_decomposition`](@ref) and the roundtrip
oracle for it.
"""
function reassemble(d::BlockDecomposition, n::Int)
    B = zeros(Int, n, n)
    for (blk, φ) in zip(d.blocks, d.placements)
        m = nvertices(blk)
        for i in 1:m, j in 1:m
            B[φ[i], φ[j]] += blk.B[i, j]
        end
    end
    return B
end

# The search assigns blocks vertex by vertex. Once vertex v has been passed no
# later block may touch it, so every block is discovered at its smallest vertex
# and each decomposition is reached exactly once.
struct _SearchState
    B     :: Matrix{Int}
    n     :: Int
    usage :: Vector{Int}          # how many blocks each quiver vertex sits in
    R     :: Matrix{Int}          # B minus what the placed blocks already give
end

# Arrows at a dead end can never cancel or double up: the dead end lies in one
# block only, so no other block contributes an arrow there. That makes the
# neighbourhood of a dead end an exact match against the quiver, which is what
# keeps the search small.
function _deadend_ok(st::_SearchState, blk::Block, φ::Vector{Int}, i::Int)
    v = φ[i]
    st.usage[v] == 0 || return false
    m = nvertices(blk)
    for j in 1:m
        blk.B[i, j] == 0 && continue
        st.B[v, φ[j]] == blk.B[i, j] || return false
    end
    # and no arrows in the quiver beyond the ones the block accounts for
    return count(!=(0), @view st.B[v, :]) == count(!=(0), @view blk.B[i, :])
end

# Can block vertex `i` of `blk` sit on quiver vertex `v`? Checks that are cheap
# and depend only on `v` itself; the arrow conditions come later.
function _slot_ok(st::_SearchState, blk::Block, v::Int, i::Int, vmin::Int)
    v >= vmin || return false
    if is_outlet(blk, i)
        return st.usage[v] <= 1
    else
        return st.usage[v] == 0
    end
end

# All placements of `blk` that use quiver vertex `v` and no vertex below `vmin`.
function _placements(st::_SearchState, blk::Block, v::Int, vmin::Int)
    m = nvertices(blk)
    out = Vector{Vector{Int}}()
    for anchor in 1:m
        _slot_ok(st, blk, v, anchor, vmin) || continue
        φ = zeros(Int, m)
        φ[anchor] = v
        _extend!(out, st, blk, φ, vmin, anchor)
    end
    return out
end

function _extend!(out, st::_SearchState, blk::Block, φ::Vector{Int}, vmin::Int, anchor::Int)
    m = nvertices(blk)
    i = findfirst(==(0), φ)
    if i === nothing
        # dead-end neighbourhoods must match the quiver exactly
        for k in 1:m
            is_outlet(blk, k) || _deadend_ok(st, blk, φ, k) || return
        end
        push!(out, copy(φ))
        return
    end
    for v in vmin:st.n
        v in φ && continue
        _slot_ok(st, blk, v, i, vmin) || continue
        φ[i] = v
        # An arrow of the block touching a dead end must be present in the
        # quiver; between two outlets it may cancel against another block.
        ok = true
        for j in 1:m
            (φ[j] == 0 || blk.B[i, j] == 0) && continue
            if !is_outlet(blk, i) || !is_outlet(blk, j)
                if st.B[φ[i], φ[j]] != blk.B[i, j]
                    ok = false
                    break
                end
            end
        end
        ok && _extend!(out, st, blk, φ, vmin, anchor)
        φ[i] = 0
    end
    return
end

function _place!(st::_SearchState, blk::Block, φ::Vector{Int})
    m = nvertices(blk)
    for i in 1:m
        st.usage[φ[i]] += 1
        for j in 1:m
            st.R[φ[i], φ[j]] -= blk.B[i, j]
        end
    end
    return
end

function _unplace!(st::_SearchState, blk::Block, φ::Vector{Int})
    m = nvertices(blk)
    for i in 1:m
        st.usage[φ[i]] -= 1
        for j in 1:m
            st.R[φ[i], φ[j]] += blk.B[i, j]
        end
    end
    return
end

# Close out vertex v: it must sit in one or two blocks, and once we move on its
# row of the residual has to be zero, because no later block may touch it.
function _search(st::_SearchState, v::Int, blocks::Vector{Block},
                 places::Vector{Vector{Int}})
    if v > st.n
        return all(iszero, st.R) ? BlockDecomposition(copy(blocks), copy(places)) : nothing
    end
    if st.usage[v] >= 1 && all(iszero, @view st.R[v, :])
        found = _search(st, v + 1, blocks, places)
        found === nothing || return found
    end
    st.usage[v] < 2 || return nothing
    for blk in FST_BLOCKS
        for φ in _placements(st, blk, v, v)
            _place!(st, blk, φ)
            push!(blocks, blk)
            push!(places, φ)
            found = _search(st, v, blocks, places)
            pop!(blocks)
            pop!(places)
            _unplace!(st, blk, φ)
            found === nothing || return found
        end
    end
    return nothing
end

"""
    block_decomposition(q::Quiver) → Union{BlockDecomposition, Nothing}

Decompose the skew-symmetric quiver `q` into Felikson-Shapiro-Tumarkin blocks,
or return `nothing` when no decomposition exists. Frozen vertices are ignored.

A decomposition exists exactly when `q` is the adjacency quiver of an ideal
tagged triangulation of a bordered surface with marked points, so it is a
certificate of surface type - and in particular of mutation finiteness. It is
not unique: an oriented triangle decomposes both as a single block `:II` and as
three blocks `:I` glued in a cycle, and this returns whichever the search meets
first.

Throws for a quiver that is not skew-symmetric, where blocks say nothing.

```jldoctest
julia> d = block_decomposition(Quiver(:A, 3));

julia> reassemble(d, 3) == Quiver(:A, 3).B
true

julia> block_decomposition(Quiver(:E, 6)) === nothing
true
```

See [`is_block_decomposable`](@ref), [`is_surface_type`](@ref).
"""
function block_decomposition(q::Quiver)
    n = q.n_mutable
    n >= 1 || throw(InvalidArgument("block_decomposition requires n_mutable ≥ 1"))
    qm = Quiver(q.B[1:n, 1:n], n, q.d)
    _is_skew_symmetric(qm) || throw(InvalidArgument(
        "block_decomposition is only defined for skew-symmetric quivers"))

    B = Matrix{Int}(qm.B)
    # Blocks contribute at most weight 2 between any pair, so anything larger
    # rules out a decomposition before the search starts.
    maximum(abs, B; init = 0) <= 2 || return nothing

    st = _SearchState(B, n, zeros(Int, n), copy(B))
    return _search(st, 1, Block[], Vector{Int}[])
end

"""
    is_block_decomposable(q::Quiver) → Bool

Whether the skew-symmetric quiver `q` decomposes into Felikson-Shapiro-Tumarkin
blocks, i.e. whether it comes from a triangulated surface.

See [`block_decomposition`](@ref) for what the decomposition records.
"""
is_block_decomposable(q::Quiver) = block_decomposition(q) !== nothing

"""
    is_surface_type(q::Quiver) → Bool

Whether the skew-symmetric quiver `q` is the adjacency quiver of an ideal tagged
triangulation of a bordered surface with marked points.

This is [`is_block_decomposable`](@ref) under the name the classification gives
it. Together with [`mutation_type`](@ref) it separates the three reasons that
function answers `nothing`: a surface type with no Dynkin-style label is of
surface type and mutation finite, an exceptional class is neither, and a
mutation-infinite quiver is neither.

Every mutation-finite skew-symmetric quiver of rank at least 3 is either of
surface type or mutation equivalent to one of the eleven exceptional quivers
`X6`, `X7`, `E6`, `E7`, `E8`, `E6^(1)`, `E7^(1)`, `E8^(1)`, `E6^(1,1)`,
`E7^(1,1)`, `E8^(1,1)`.

```jldoctest
julia> is_surface_type(Quiver([0 2 -2; -2 0 2; 2 -2 0]))   # once-punctured torus
true

julia> mutation_type(Quiver([0 2 -2; -2 0 2; 2 -2 0])) === nothing
true

julia> is_surface_type(Quiver(:E, 6))
false
```
"""
is_surface_type(q::Quiver) = is_block_decomposable(q)
