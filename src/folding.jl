# Folding a symmetric quiver by an admissible automorphism to obtain a
# skew-symmetrizable (non-simply-laced) quiver - the standard route from A/D/E to
# B/C/F/G. Given a skew-symmetric B and a vertex permutation σ that preserves B
# (an automorphism) whose orbits are independent sets, the folded exchange matrix
# on the σ-orbits is
#
#     C[a,b] = Σ_{j ∈ orbit(b)} B[i, j]          (i any representative of orbit a)
#
# with symmetrizer d[a] = |orbit(a)|.  One checks d[a]·C[a,b] = −d[b]·C[b,a], so C
# is skew-symmetrizable (the Quiver constructor re-verifies this).

# σ-orbits (cycles of the permutation), each sorted, ordered by least element.
function _orbits(σ::Vector{Int})
    n = length(σ)
    seen = falses(n)
    orbs = Vector{Vector{Int}}()
    for i in 1:n
        seen[i] && continue
        orb = Int[]
        j = i
        while !seen[j]
            seen[j] = true
            push!(orb, j)
            j = σ[j]
        end
        push!(orbs, sort!(orb))
    end
    return sort!(orbs; by = first)
end

_is_permutation(σ::Vector{Int}, n::Int) = length(σ) == n && sort(σ) == collect(1:n)

# σ preserves the exchange matrix: B[σ(i), σ(j)] == B[i, j] for all i, j.
function _is_automorphism(B::Matrix{Int}, σ::Vector{Int})
    n = length(σ)
    for i in 1:n, j in 1:n
        B[σ[i], σ[j]] == B[i, j] || return false
    end
    return true
end

# Every orbit is an independent set: no arrows between vertices of one orbit.
function _orbits_independent(B::Matrix{Int}, orbs::Vector{Vector{Int}})
    for orb in orbs, i in orb, j in orb
        i == j && continue
        B[i, j] == 0 || return false
    end
    return true
end

"""
    is_admissible_folding(q::Quiver, σ::Vector{Int}) -> Bool

Whether `σ` (a permutation given as `σ[i] = image of i`) is an **admissible**
automorphism of the skew-symmetric quiver `q` for folding: `q` is all-mutable and
skew-symmetric, `σ` is a permutation of the vertices that preserves the exchange
matrix, and each `σ`-orbit is an independent set (no arrows inside an orbit).
See [`fold`](@ref).
"""
function is_admissible_folding(q::Quiver, σ::Vector{Int})
    n = q.n_mutable
    q.n_frozen == 0            || return false
    all(==(1), q.d)            || return false
    _is_permutation(σ, n)      || return false
    _is_automorphism(q.B, σ)   || return false
    return _orbits_independent(q.B, _orbits(σ))
end

"""
    fold(q::Quiver, σ::Vector{Int}) -> Quiver

Fold the skew-symmetric quiver `q` by the admissible automorphism `σ` (a
permutation, `σ[i] = image of i`), returning the skew-symmetrizable quiver on the
`σ`-orbits.  Classic foldings: `A_{2n-1} → C_n`, `D_{n+1} → B_n`, `E_6 → F_4`,
`D_4 → G_2`.

The orbits become the vertices (ordered by least element), the folded exchange
matrix is `C[a,b] = Σ_{j ∈ orbit(b)} B[rep(a), j]`, and the symmetrizer is
`d[a] = |orbit(a)|`.  Throws `InvalidArgument` if `q` has frozen vertices, is not
skew-symmetric, or `σ` is not an admissible automorphism (see
[`is_admissible_folding`](@ref)).
"""
function fold(q::Quiver, σ::Vector{Int})
    n = q.n_mutable
    q.n_frozen == 0 || throw(InvalidArgument(
        "fold requires an all-mutable quiver (n_frozen = $(q.n_frozen))"))
    all(==(1), q.d) || throw(InvalidArgument(
        "fold requires a skew-symmetric quiver (symmetrizer d = $(q.d) ≠ 1)"))
    _is_permutation(σ, n) || throw(InvalidArgument(
        "σ = $σ is not a permutation of 1:$n"))
    _is_automorphism(q.B, σ) || throw(InvalidArgument(
        "σ = $σ is not an automorphism of the exchange matrix"))

    orbs = _orbits(σ)
    _orbits_independent(q.B, orbs) || throw(InvalidArgument(
        "σ = $σ is not admissible: some orbit has an internal arrow"))

    m = length(orbs)
    C = zeros(Int, m, m)
    for a in 1:m, b in 1:m
        i = first(orbs[a])                       # representative of orbit a
        C[a, b] = sum(q.B[i, j] for j in orbs[b])
    end
    d = [length(orb) for orb in orbs]
    labels = [join(orb, "+") for orb in orbs]
    return Quiver(C, m, d, labels)
end
