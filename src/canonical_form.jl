# ─── Vertex relabeling ────────────────────────────────────────────────────────

"""
    permute_vertices(q::Quiver, σ::AbstractVector{Int}) → Quiver

Relabel the vertices of `q` by the permutation `σ`, so that vertex `i` of `q`
becomes vertex `σ[i]` of the result: `B′[σ[i], σ[j]] == B[i, j]`.

`σ` must permute the mutable block `1:n_mutable` and the frozen block onto
themselves - the mutable/frozen split is structural, not a labeling choice.
Symmetrizers and labels are carried along with their vertices, so the result is
the same quiver presented in a different vertex order.

# Example
```julia
julia> q = Quiver(:A, 3);

julia> is_isomorphic(permute_vertices(q, [3, 1, 2]), q)
true
```
"""
function permute_vertices(q::Quiver, σ::AbstractVector{Int})
    n = nvertices(q)
    length(σ) == n ||
        throw(InvalidArgument("permutation must have length $n, got $(length(σ))"))
    isperm(σ) ||
        throw(InvalidArgument("σ is not a permutation of 1:$n"))
    for i in 1:n
        is_frozen(q, i) == is_frozen(q, σ[i]) ||
            throw(InvalidArgument(
                "σ maps vertex $i across the mutable/frozen split (to $(σ[i]))"))
    end

    B = zeros(Int, n, n)
    for i in 1:n, j in 1:n
        B[σ[i], σ[j]] = q.B[i, j]
    end
    d = Vector{Int}(undef, q.n_mutable)
    for i in 1:q.n_mutable
        d[σ[i]] = q.d[i]
    end
    labels = Vector{String}(undef, n)
    for i in 1:n
        labels[σ[i]] = q.labels[i]
    end
    return Quiver(B, q.n_mutable, d, labels)
end

# ─── Colour refinement ────────────────────────────────────────────────────────

# Initial vertex colours: everything an isomorphism must preserve but that costs
# nothing to read off - the mutable/frozen split and the symmetrizer.  Frozen
# vertices have no symmetrizer, so they take 0.
function _initial_colours(q::Quiver)
    n = nvertices(q)
    keys = [(is_frozen(q, i), i <= q.n_mutable ? q.d[i] : 0) for i in 1:n]
    return _renumber(keys)
end

# Map an arbitrary vector of keys to colours 1:k, ordered by sorted key, so the
# numbering depends only on the multiset of keys and never on vertex order.
function _renumber(keys::AbstractVector)
    order = sort(unique(keys))
    index = Dict(k => c for (c, k) in enumerate(order))
    return [index[k] for k in keys]
end

# One round of 1-WL refinement: recolour each vertex by its own colour together
# with the multiset of (colour, arrow weight) pairs on its outgoing and incoming
# arrows.  Iterated to a stable partition.
function _refine_colours(B::Matrix{Int}, colours::Vector{Int})
    n = length(colours)
    while true
        keys = Vector{Tuple{Int, Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}}(undef, n)
        for i in 1:n
            out = sort([(colours[j], B[i, j]) for j in 1:n if B[i, j] != 0])
            inc = sort([(colours[j], B[j, i]) for j in 1:n if B[j, i] != 0])
            keys[i] = (colours[i], out, inc)
        end
        new_colours = _renumber(keys)
        new_colours == colours && return colours
        colours = new_colours
    end
end

# ─── Canonical form ───────────────────────────────────────────────────────────

# The exchange matrix read out under a candidate ordering, shell by shell: for
# each vertex k in order, its row against all earlier vertices, then its column.
#
# Reading in shells rather than column-major is what makes the search prunable:
# the key of an m-vertex prefix is a genuine *prefix* of the key of any ordering
# extending it, so a prefix that already loses lexicographically can never win.
# A column-major flattening does not have this property - its first column
# interleaves entries from vertices outside the prefix.
function _shell_key(B::Matrix{Int}, order::AbstractVector{Int})
    key = Int[]
    for k in eachindex(order)
        for j in 1:k
            push!(key, B[order[k], order[j]])
        end
        for j in 1:(k - 1)
            push!(key, B[order[j], order[k]])
        end
    end
    return key
end

"""
    canonical_permutation(q::Quiver) → Vector{Int}

The permutation `σ` sending `q` to its canonical form: `permute_vertices(q, σ)`
is [`canonical_form(q)`](@ref).

Vertices are first partitioned by 1-WL colour refinement (seeded by the
mutable/frozen split and the symmetrizers), then the orderings compatible with
that partition are searched, and the one minimising the flattened exchange
matrix lexicographically is returned. The refinement makes the search cheap on
the ranks this package is used at; it is not a nauty-grade implementation and is
not intended for large graphs.
"""
function canonical_permutation(q::Quiver)
    n = nvertices(q)
    n == 0 && return Int[]

    colours = _refine_colours(q.B, _initial_colours(q))

    # Candidate orderings visit colour classes in colour order; within a class,
    # every arrangement is a candidate.  Frozen vertices get a strictly larger
    # initial colour than mutable ones, so they always land in the frozen block.
    classes = [findall(==(c), colours) for c in 1:maximum(colours)]

    best_order = Int[]
    best_key   = Int[]

    # Depth-first over the classes, pruning any prefix already losing on the
    # shells it fixes.
    function extend!(prefix::Vector{Int}, ci::Int)
        if ci > length(classes)
            key = _shell_key(q.B, prefix)
            if isempty(best_key) || key < best_key
                best_key   = key
                best_order = copy(prefix)
            end
            return
        end
        for perm in _permutations(classes[ci])
            append!(prefix, perm)
            key = _shell_key(q.B, prefix)
            if isempty(best_key) || key <= best_key[1:length(key)]
                extend!(prefix, ci + 1)
            end
            resize!(prefix, length(prefix) - length(perm))
        end
    end

    extend!(Int[], 1)

    # best_order[k] is the vertex placed at position k; invert to get σ.
    σ = Vector{Int}(undef, n)
    for (pos, v) in enumerate(best_order)
        σ[v] = pos
    end
    return σ
end

# All arrangements of `v`.  Recursive and allocating, which is fine: it is only
# ever called on a single colour class, and refinement keeps those small.
function _permutations(v::Vector{Int})
    length(v) <= 1 && return [copy(v)]
    out = Vector{Vector{Int}}()
    for (i, x) in enumerate(v)
        rest = deleteat!(copy(v), i)
        for p in _permutations(rest)
            push!(out, pushfirst!(p, x))
        end
    end
    return out
end

"""
    canonical_form(q::Quiver) → Quiver

A canonical representative of `q`'s isomorphism class: two quivers are isomorphic
(equal up to relabeling vertices within the mutable and frozen blocks) exactly
when their canonical forms have equal exchange matrices and symmetrizers.

Labels are carried along with their vertices rather than reset, so the canonical
form is the same quiver in a canonical vertex order. Note that this means
`canonical_form(q1) == canonical_form(q2)` - which compares labels - is *stricter*
than isomorphism; use [`is_isomorphic`](@ref) to test isomorphism, and
`canonical_form(q).B` as a hash key for isomorphism-class deduplication.

Canonical forms give normal forms for dataset deduplication and for splitting
quivers into isomorphism-disjoint groups, neither of which
[`mutation_class`](@ref) provides - it deduplicates by raw exchange matrix, hence
not up to relabeling.

# Example
```julia
julia> q = Quiver(:A, 3);

julia> canonical_form(permute_vertices(q, [2, 3, 1])).B == canonical_form(q).B
true
```
"""
canonical_form(q::Quiver) = permute_vertices(q, canonical_permutation(q))

"""
    is_isomorphic(q1::Quiver, q2::Quiver) → Bool

Whether `q1` and `q2` agree up to relabeling of vertices, within the mutable and
frozen blocks separately. Labels are ignored; symmetrizers are not.

This is isomorphism of quivers, *not* mutation equivalence: isomorphic quivers
lie in the same mutation class, but mutation-equivalent quivers are generally not
isomorphic.

# Example
```julia
julia> is_isomorphic(Quiver(:A, 3), permute_vertices(Quiver(:A, 3), [3, 2, 1]))
true
```
"""
function AbstractAlgebra.is_isomorphic(q1::Quiver, q2::Quiver)
    nvertices(q1) == nvertices(q2)   || return false
    q1.n_mutable == q2.n_mutable     || return false
    c1 = canonical_form(q1)
    c2 = canonical_form(q2)
    return c1.B == c2.B && c1.d == c2.d
end
