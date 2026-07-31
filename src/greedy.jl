# ─── Greedy elements / the greedy basis in rank 2 ──────────────────────────────
#
# Lee-Li-Zelevinsky, "Greedy elements in rank 2 cluster algebras" (arXiv:1208.2391,
# Selecta Math 2014).  For B = [0 b; −c 0] with b, c ≥ 1 and any (a₁, a₂) ∈ ℤ²,
#
#   x[a₁,a₂] = x₁^{−a₁} x₂^{−a₂} · Σ_{p,q} d(p,q) x₁^{b p} x₂^{c q},
#
# with the coefficients given by LLZ's recursion (1.5).  The greedy elements form a
# positive ℤ-basis of the rank-2 cluster algebra containing every cluster monomial;
# for b·c ≤ 3 the basis *is* the set of cluster monomials, and beyond that it adds
# "imaginary" elements (the first one being (1 + x₁² + x₂²)/(x₁x₂) for b = c = 2).
#
# Convention ledger (both pinned by the cluster-variable oracle in
# test/test_greedy.jl - a cluster variable equals the greedy element at its own
# denominator vector, exactly, so a wrong choice shows up immediately):
#
#   1. C(m,k) is the ORDINARY binomial, zero for m < 0.  With the generalized
#      binomial instead, Kronecker (a₁,a₂) = (1,1) gives d(1,1) = −1 plus an
#      unterminated alternating tail, contradicting LLZ positivity.
#   2. The support bounds CROSS: 0 ≤ p ≤ [a₂]₊ and 0 ≤ q ≤ [a₁]₊.  This is LLZ
#      Theorem 1.11: d(p,q) counts compatible pairs (S₁,S₂) in the maximal Dyck
#      path D^{[a₁]₊ × [a₂]₊} with |S₂| = p vertical and |S₁| = q horizontal edges.
#      Using p ≤ [a₁]₊, q ≤ [a₂]₊ instead breaks even A₂.

# Ordinary binomial extended by zero to negative top arguments (ledger point 1).
_greedy_binomial(m::Int, k::Int) =
    (m < 0 || k > m) ? big(0) : binomial(big(m), big(k))

"""
    greedy_coefficients(a::AbstractVector{Int}, b::Int, c::Int) → Matrix{BigInt}

The greedy coefficients `d(p,q)` of the greedy element `x[a₁,a₂]` of the rank-2
cluster algebra with exchange matrix `B = [0 b; −c 0]`, as a table indexed
`M[p+1, q+1]` of size `([a₂]₊ + 1) × ([a₁]₊ + 1)` (`[·]₊ = max(·, 0)`).

`d(p,q)` is the coefficient of `x₁^{b·p} x₂^{c·q}` in
`x₁^{a₁} x₂^{a₂} · x[a₁,a₂]`, computed by the Lee-Li-Zelevinsky recursion

```
d(0,0) = 1,
d(p,q) = max( Σ_{k=1}^{p} (−1)^{k−1} d(p−k, q) C(a₂ − c·q + k − 1, k),
              Σ_{k=1}^{q} (−1)^{k−1} d(p, q−k) C(a₁ − b·p + k − 1, k) )
```

with the ordinary binomial (zero for a negative top argument).  Coefficients are
nonnegative (LLZ positivity) and returned as `BigInt`, since they grow quickly.

# Example
```julia
julia> greedy_coefficients([1, 1], 2, 2)      # Kronecker: 1 + x₁² + x₂²
2×2 Matrix{BigInt}:
 1  1
 1  0
```
"""
function greedy_coefficients(a::AbstractVector{Int}, b::Int, c::Int)
    length(a) == 2 || throw(InvalidArgument(
        "greedy_coefficients needs a rank-2 exponent vector, got length $(length(a))"))
    b >= 1 && c >= 1 || throw(InvalidArgument(
        "greedy_coefficients needs b, c ≥ 1, got b = $b, c = $c"))
    a1, a2 = a[1], a[2]
    P, Q = max(a2, 0), max(a1, 0)              # ledger point 2: the bounds cross
    d = zeros(BigInt, P + 1, Q + 1)
    d[1, 1] = big(1)
    for p in 0:P, q in 0:Q
        (p, q) == (0, 0) && continue
        s1 = sum((-1)^(k - 1) * d[p - k + 1, q + 1] *
                 _greedy_binomial(a2 - c * q + k - 1, k) for k in 1:p; init = big(0))
        s2 = sum((-1)^(k - 1) * d[p + 1, q - k + 1] *
                 _greedy_binomial(a1 - b * p + k - 1, k) for k in 1:q; init = big(0))
        d[p + 1, q + 1] = max(s1, s2)
    end
    return d
end

"""
    greedy_element(s::Seed, a::AbstractVector{Int})
    greedy_element(q::Quiver, a)

The greedy element `x[a₁,a₂]` of the rank-2 cluster algebra of `s`, as an element
of the seed's ambient field `Frac(ℤ[x₁,x₂])` - directly comparable to the entries
of `s.cluster`.

`(b, c)` are read off the exchange matrix; both orientations are accepted (for
`B = [0 −b; c 0]` the two vertices, and correspondingly `a₁` and `a₂`, swap
roles).  The greedy elements for `a ∈ ℤ²` form a positive ℤ-basis containing every
cluster monomial: `x[−m,−n] = x₁^m x₂^n`, and a cluster variable is the greedy
element at its own [`denominator_vector`](@ref).  See
[`greedy_coefficients`](@ref) for the expansion.

Throws `InvalidArgument` unless the seed has exactly two mutable vertices, no
frozen vertices (the greedy basis is coefficient-free), and a connected quiver.

# Example
```julia
julia> s = Seed(Quiver(:A, 2));

julia> greedy_element(s, [1, 0]) == mutate(s, 1).cluster[1]
true
```
"""
function greedy_element(s::Seed, a::AbstractVector{Int})
    length(a) == 2 || throw(InvalidArgument(
        "greedy_element needs a rank-2 exponent vector, got length $(length(a))"))
    s.quiver.n_mutable == 2 || throw(InvalidArgument(
        "greedy elements are defined for rank 2, got n_mutable = $(s.quiver.n_mutable)"))
    s.quiver.n_frozen == 0 || throw(InvalidArgument(
        "greedy elements are defined for a coefficient-free cluster algebra, but " *
        "the seed has $(s.quiver.n_frozen) frozen vertices"))

    B = s.quiver.B
    # LLZ orientation is B = [0 b; −c 0].  The reversed one is the same algebra
    # with the two vertices swapped, so swap (b,c), the exponents and the variables.
    forward = B[1, 2] >= 0
    b, c    = forward ? (B[1, 2], -B[2, 1]) : (B[2, 1], -B[1, 2])
    a1, a2  = forward ? (a[1], a[2]) : (a[2], a[1])
    x1, x2  = forward ? (s.cluster[1], s.cluster[2]) : (s.cluster[2], s.cluster[1])
    b >= 1 && c >= 1 || throw(InvalidArgument(
        "greedy elements need a connected rank-2 quiver, got B = $B"))

    d   = greedy_coefficients([a1, a2], b, c)
    acc = zero(s.ring)
    for p in axes(d, 1), q in axes(d, 2)
        iszero(d[p, q]) && continue
        acc += s.ring(d[p, q]) * x1^(b * (p - 1)) * x2^(c * (q - 1))
    end
    return acc * x1^(-a1) * x2^(-a2)
end

greedy_element(q::Quiver, a::AbstractVector{Int}) = greedy_element(Seed(q), a)
