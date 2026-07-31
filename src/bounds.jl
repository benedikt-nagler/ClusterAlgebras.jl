# ─── Upper and lower bounds of a cluster algebra ───────────────────────────────
#
# Berenstein-Fomin-Zelevinsky, "Cluster algebras III: upper bounds and double
# Bruhat cells" ([BFZ05]).  For a seed Σ = (x, B̃) of geometric type, with
# ZP = Z[x_{n+1}^{±1}, …, x_m^{±1}] the coefficient ring,
#
#   upper bound   U(Σ) = ZP[x^{±1}] ∩ ZP[x_1^{±1}] ∩ ⋯ ∩ ZP[x_n^{±1}]   ([BFZ05] Def. 1.1)
#   lower bound   L(Σ) = ZP[x_1, x_1', …, x_n, x_n']                    ([BFZ05] Def. 1.10)
#
# where x_k is the cluster adjacent to x at k.  Always L(Σ) ⊆ A ⊆ U(Σ); the
# content of [BFZ05] is that a coprime acyclic seed collapses the sandwich:
# L(Σ) = U(Σ) (Thm. 1.18), so A equals the upper cluster algebra (Cor. 1.19)
# and the standard monomials are a ZP-basis (Thm. 1.16, an *iff*).
#
# What this file deliberately does NOT do: compute generators of U(Σ) in the
# non-acyclic case.  That is Sage's `find_upper_bound`, which needs a Groebner
# basis; it belongs behind the OSCAR firewall in a downstream package, not in
# this dependency-light core.  Everything here is exact and Groebner-free.
#
# Convention ledger (each pinned by an oracle in test/test_bounds.jl):
#
#   1. All bound functions take an INITIAL seed (empty `mutation_path`).  The
#      cluster of a mutated `Seed` is a vector of rational functions of the
#      *initial* variables, so "Laurent in this seed's cluster" would silently
#      mean the wrong thing; we refuse rather than guess.  Use `Seed(s.quiver)`
#      to get the initial seed on a mutated quiver.  Pinned by the
#      round-trip oracle: bounds of `Seed(mutate(q,k))` vs `mutate(Seed(q),k)`.
#   2. `in_upper_bound` re-expresses f in the adjacent cluster by SUBSTITUTING
#      x_k ↦ P_k/x_k into f - no new ring, no renaming.  This is legitimate
#      because μ_k is an involution: P_k does not involve x_k, so the
#      substitution is its own inverse and sends the cluster x to the adjacent
#      cluster under the same variable names.  Pinned by the Laurent-phenomenon
#      oracle (every cluster variable in the class must pass).
#   3. Units of ZP are ± monomials in the FROZEN variables.  So a Laurent
#      denominator is any single term with coefficient ±1, and two exchange
#      polynomials count as coprime when their gcd is such a unit.  Pinned by
#      the geometric-seed oracle on `grassmannian(2, 5)`.
#   4. `lower_bound_generators` returns the 2n mutable generators only.  The
#      frozen variables are units of ZP ([BFZ05] Def. 1.10 takes L over ZP),
#      not algebra generators.

# ─── Acyclicity ───────────────────────────────────────────────────────────────

"""
    is_acyclic(q::Quiver) → Bool
    is_acyclic(s::Seed)   → Bool

Whether the mutable part of the quiver has no oriented cycle.

Acyclicity is a property of a single seed, not of its mutation class: the
oriented 3-cycle is not acyclic but the `A₃` quiver it is mutation equivalent
to is.  Use [`bound_certificate`](@ref) for the class-level question, which is
what Berenstein-Fomin-Zelevinsky's `𝒜 = 𝒰` theorem actually needs.
"""
is_acyclic(q::Quiver) = _is_acyclic(q)
is_acyclic(s::Seed)   = _is_acyclic(s.quiver)

# ─── Seed and ring plumbing ───────────────────────────────────────────────────

# Ledger point 1.
function _require_initial_seed(s::Seed, fn::String)
    isempty(s.mutation_path) || throw(InvalidArgument(
        "$fn requires an initial seed, but this seed carries the mutation path " *
        "$(s.mutation_path); the bounds are defined relative to a seed's own " *
        "cluster, so pass `Seed(s.quiver)` to work with the mutated quiver"))
    return nothing
end

_bounds_poly_ring(s::Seed) = base_ring(s.ring)

# A unit of ZP: ±1 times a monomial in the frozen variables only (ledger 3).
function _is_zp_unit(p, n_mutable::Int)
    length(p) == 1 || return false
    isone(abs(coeff(p, 1))) || return false
    ev = first(exponent_vectors(p))
    return all(iszero, @view ev[1:n_mutable])
end

# The exchange polynomial P_k = x_k x_k' as an element of the polynomial ring.
function _exchange_polynomial(s::Seed, k::Int)
    pos, neg = _exchange_monomials(s, k)
    return numerator(pos + neg)
end

# ─── Laurentness and upper-bound membership ───────────────────────────────────

"""
    is_laurent(f, s::Seed) → Bool

Whether `f` is a Laurent polynomial in the cluster of the initial seed `s`,
i.e. whether it lies in `ZP[x^{±1}]`.  A reduced fraction qualifies exactly
when its denominator is a single monomial with coefficient `±1`.
"""
function is_laurent(f, s::Seed)
    _require_initial_seed(s, "is_laurent")
    den = denominator(f)
    return length(den) == 1 && isone(abs(coeff(den, 1)))
end

# f re-expressed in the cluster adjacent to x at k (ledger point 2).
function _substitute_adjacent(f, s::Seed, k::Int)
    pos, neg = _exchange_monomials(s, k)
    vals     = copy(s.cluster)
    vals[k]  = (pos + neg) / s.cluster[k]
    return evaluate(numerator(f), vals) / evaluate(denominator(f), vals)
end

"""
    in_upper_bound(f, s::Seed) → Bool

Whether `f` lies in the upper bound `𝒰(Σ)` of the initial seed `s`: Berenstein-
Fomin-Zelevinsky define `𝒰(Σ)` as the intersection of the Laurent ring of the
cluster of `s` with the Laurent rings of the `n` adjacent clusters, so this is
`n + 1` exact Laurentness checks and needs no Groebner basis.

Every cluster variable of the algebra lies in `𝒰(Σ)` (the Laurent phenomenon),
so this is a decision procedure for the *upper* bound whatever the quiver;
what needs acyclicity is only the identification of `𝒰(Σ)` with the cluster
algebra itself - see [`bound_certificate`](@ref).

```jldoctest
julia> s = Seed(Quiver(:A, 2));

julia> x5 = mutate(s, [1, 2])[2];

julia> in_upper_bound(x5, s)
true
```
"""
function in_upper_bound(f, s::Seed)
    _require_initial_seed(s, "in_upper_bound")
    is_laurent(f, s) || return false
    for k in 1:s.quiver.n_mutable
        is_laurent(_substitute_adjacent(f, s, k), s) || return false
    end
    return true
end

# ─── The lower bound ──────────────────────────────────────────────────────────

"""
    lower_bound_generators(s::Seed) → Vector

The `2n` generators `x₁, …, xₙ, x₁′, …, xₙ′` of the lower bound
`ℒ(Σ) = ZP[x₁, x₁′, …, xₙ, xₙ′]` of the initial seed `s`, where `xₖ′` is the
cluster variable obtained by mutating at `k`.

Only mutable vertices contribute: the frozen variables are units of the
coefficient ring `ZP`, over which `ℒ(Σ)` is taken, not algebra generators.
"""
function lower_bound_generators(s::Seed)
    _require_initial_seed(s, "lower_bound_generators")
    n = s.quiver.n_mutable
    gens_x  = s.cluster[1:n]
    gens_xp = [mutate(s, k)[k] for k in 1:n]
    return vcat(gens_x, gens_xp)
end

"""
    standard_monomials(s::Seed; max_degree::Int = 4)

The standard monomials of the initial seed `s` of total degree at most
`max_degree`: monomials in `x₁, x₁′, …, xₙ, xₙ′` containing no product
`xⱼ xⱼ′` (Berenstein-Fomin-Zelevinsky's definition).  Each entry is a named
tuple `(a, b, element)` for `∏ xⱼ^{aⱼ} ∏ (xⱼ′)^{bⱼ}`, so `aⱼ · bⱼ = 0` for
every `j`.

These form a `ZP`-basis of `ℒ(Σ)` **if and only if** the seed is acyclic, so
for a cyclic seed the list still enumerates the right monomials but they are
no longer independent.
"""
function standard_monomials(s::Seed; max_degree::Int = 4)
    _require_initial_seed(s, "standard_monomials")
    max_degree >= 0 || throw(InvalidArgument(
        "standard_monomials: max_degree must be non-negative, got $max_degree"))

    n    = s.quiver.n_mutable
    gs   = lower_bound_generators(s)
    F    = s.ring
    out  = Vector{@NamedTuple{a::Vector{Int}, b::Vector{Int}, element::eltype(gs)}}()

    # Per vertex j the exponent is (aⱼ, 0) or (0, bⱼ) - never both non-zero.
    function rec(j::Int, budget::Int, a::Vector{Int}, b::Vector{Int})
        if j > n
            elt = one(F)
            for i in 1:n
                a[i] > 0 && (elt *= gs[i]^a[i])
                b[i] > 0 && (elt *= gs[n + i]^b[i])
            end
            push!(out, (a = copy(a), b = copy(b), element = elt))
            return
        end
        rec(j + 1, budget, a, b)                       # exponent 0 at j
        for e in 1:budget
            a[j] = e; rec(j + 1, budget - e, a, b); a[j] = 0
            b[j] = e; rec(j + 1, budget - e, a, b); b[j] = 0
        end
    end
    rec(1, max_degree, zeros(Int, n), zeros(Int, n))
    return out
end

"""
    lower_bound_expansion(f, s::Seed; max_degree::Int = 4)

Expand `f` in the standard monomials of the initial seed `s` of total degree at
most `max_degree`, returning a vector of named tuples
`(a, b, coefficient, element)` with non-zero integer coefficients, or `nothing`
if no such expansion exists.

`nothing` does **not** prove that `f ∉ ℒ(Σ)` - only that no expansion fits
inside the degree budget.  When the seed is acyclic the standard monomials are
a basis, so an expansion that is found is unique.

Requires a coefficient-free seed (no frozen vertices): with coefficients the
scalars live in `ZP = Z[x_{n+1}^{±1}, …]` rather than in `Z`, which is a linear
algebra problem over a Laurent ring rather than over `Q`.

```jldoctest
julia> s = Seed(Quiver(:A, 2));

julia> e = lower_bound_expansion(mutate(s, [1, 2])[2], s);

julia> [(t.a, t.b, t.coefficient) for t in e]
2-element Vector{Tuple{Vector{Int64}, Vector{Int64}, BigInt}}:
 ([0, 0], [0, 0], -1)
 ([0, 0], [1, 1], 1)
```
"""
function lower_bound_expansion(f, s::Seed; max_degree::Int = 4)
    _require_initial_seed(s, "lower_bound_expansion")
    iszero(s.quiver.n_frozen) || throw(InvalidArgument(
        "lower_bound_expansion requires a coefficient-free seed; this one has " *
        "$(s.quiver.n_frozen) frozen vertices, so the coefficients of an " *
        "expansion would lie in the Laurent ring ZP rather than in Z"))
    is_laurent(f, s) || return nothing

    n     = s.quiver.n_mutable
    R     = _bounds_poly_ring(s)
    cands = standard_monomials(s; max_degree = max_degree)

    # Clear denominators: every standard monomial has denominator dividing
    # ∏ xⱼ^{max_degree}, and f's is a monomial by the Laurentness check above.
    fden  = denominator(f)
    fev   = first(exponent_vectors(fden))
    shift = [max(max_degree, fev[j]) for j in 1:n]
    N     = prod(gens(R)[j]^shift[j] for j in 1:n; init = one(R))

    polys  = [numerator(t.element * N) for t in cands]
    target = numerator(f * N)

    # Row index: every exponent vector occurring anywhere.
    rows = Dict{Vector{Int}, Int}()
    for p in polys, ev in exponent_vectors(p)
        get!(rows, ev, length(rows) + 1)
    end
    for ev in exponent_vectors(target)
        haskey(rows, ev) || return nothing   # a monomial no candidate can reach
    end

    M = zero_matrix(QQ, length(rows), length(polys))
    for (j, p) in enumerate(polys)
        for (c, ev) in zip(coefficients(p), exponent_vectors(p))
            M[rows[ev], j] = QQ(c)
        end
    end
    v = zero_matrix(QQ, length(rows), 1)
    for (c, ev) in zip(coefficients(target), exponent_vectors(target))
        v[rows[ev], 1] = QQ(c)
    end

    ok, sol = can_solve_with_solution(M, v; side = :right)
    ok || return nothing

    out = Vector{@NamedTuple{a::Vector{Int}, b::Vector{Int},
                             coefficient::BigInt, element::eltype(s.cluster)}}()
    for j in eachindex(cands)
        c = sol[j, 1]
        iszero(c) && continue
        isone(denominator(c)) || return nothing   # not a Z-combination
        push!(out, (a = cands[j].a, b = cands[j].b,
                    coefficient = BigInt(numerator(c)), element = cands[j].element))
    end
    return out
end

# ─── Coprimality, full rank, and the BFZ certificate ──────────────────────────

"""
    is_coprime(s::Seed) → Bool

Whether the exchange polynomials `P₁, …, Pₙ` of the initial seed `s` are
pairwise coprime in `ZP[x]`, i.e. whether every common divisor of `Pᵢ` and `Pⱼ`
is a unit of the coefficient ring.  This is the hypothesis, alongside
acyclicity, of Berenstein-Fomin-Zelevinsky's `ℒ(Σ) = 𝒰(Σ)` theorem.

It is not automatic: `Pᵢ` is a binomial, and binomials do factor (`x₁³ + x₂³`
has the factor `x₁ + x₂`), so a gcd is genuinely computed.
"""
function is_coprime(s::Seed)
    _require_initial_seed(s, "is_coprime")
    n = s.quiver.n_mutable
    P = [_exchange_polynomial(s, k) for k in 1:n]
    for i in 1:n, j in (i + 1):n
        _is_zp_unit(gcd(P[i], P[j]), n) || return false
    end
    return true
end

"""
    has_full_rank(s::Seed) → Bool
    has_full_rank(q::Quiver) → Bool

Whether the extended exchange matrix `B̃` (all rows, mutable columns) has full
rank `n`.  By Berenstein-Fomin-Zelevinsky this is a sufficient condition for
*every* seed in the mutation class to be coprime, so it upgrades
[`is_coprime`](@ref) from a statement about one seed to one about the class.

A coefficient-free seed has `B̃ = B` skew-symmetrizable, hence of even rank, so
full rank fails for every odd rank - which is why the certificate tests
coprimality directly instead of relying on this.
"""
function has_full_rank(q::Quiver)
    n = q.n_mutable
    n == 0 && return true
    B_tilde = matrix(QQ, [QQ(q.B[i, j]) for i in 1:size(q.B, 1), j in 1:n])
    return rank(B_tilde) == n
end
has_full_rank(s::Seed) = has_full_rank(s.quiver)

"""
    BoundCertificate

The outcome of [`bound_certificate`](@ref): which hypotheses of Berenstein-
Fomin-Zelevinsky's `𝒜 = 𝒰` theorem hold for a seed, and what therefore follows.

Fields: `acyclic`, `coprime` and `full_rank` for the seed itself,
`acyclic_class` (`true`/`false`/`missing` - whether the mutation class contains
an acyclic seed, `missing` when the search budget was exhausted), and
`conclusion`, one of

- `:seed_acyclic` - the seed is coprime and acyclic, so
  `ℒ(Σ) = 𝒜 = 𝒰(Σ) = 𝒰` and the standard monomials are a basis of all of them;
- `:class_acyclic` - this seed is not acyclic but the class contains a coprime
  acyclic seed, so `𝒜 = 𝒰` still holds as algebras, though nothing is claimed
  about `ℒ(Σ)` at *this* seed;
- `:not_coprime` - coprimality fails, so the theorem does not apply here;
- `:not_acyclic` - no seed in the mutation class is acyclic (the Markov quiver
  is the standard example), so nothing is concluded;
- `:inconclusive` - the search for an acyclic representative ran out of budget.
"""
struct BoundCertificate
    acyclic::Bool
    coprime::Bool
    full_rank::Bool
    acyclic_class::Union{Bool, Missing}
    conclusion::Symbol
end

"""
    bound_certificate(s::Seed) → BoundCertificate

Decide which of Berenstein-Fomin-Zelevinsky's hypotheses the initial seed `s`
satisfies, and hence whether the cluster algebra equals its upper cluster
algebra.  See [`BoundCertificate`](@ref) for the possible conclusions.

The class-level search reuses the same acyclic-representative BFS as the
finite-type classification, so a cyclic seed whose class contains an acyclic
one (an oriented 3-cycle of `A₃`, say) is still certified.

```jldoctest
julia> bound_certificate(Seed(Quiver(:A, 2))).conclusion
:seed_acyclic

julia> bound_certificate(Seed(Quiver(:A, 3))).conclusion   # P₁ = 1 + x₂ = P₃
:not_coprime
```
"""
function bound_certificate(s::Seed)
    _require_initial_seed(s, "bound_certificate")
    acyclic   = is_acyclic(s)
    coprime   = is_coprime(s)
    full_rank = has_full_rank(s)

    if acyclic && coprime
        return BoundCertificate(true, true, full_rank, true, :seed_acyclic)
    elseif !coprime && acyclic
        return BoundCertificate(acyclic, false, full_rank, true, :not_coprime)
    end

    rep = _acyclic_representative(s.quiver)
    if rep === missing
        @warn "bound_certificate: acyclic-representative search budget exhausted " *
              "before the mutation class was covered; reporting :inconclusive " *
              "rather than guessing that no acyclic seed exists"
        return BoundCertificate(acyclic, coprime, full_rank, missing, :inconclusive)
    elseif rep === nothing
        return BoundCertificate(acyclic, coprime, full_rank, false, :not_acyclic)
    end

    # A coprime acyclic seed elsewhere in the class still gives 𝒜 = 𝒰.
    is_coprime(Seed(rep)) ||
        return BoundCertificate(acyclic, coprime, full_rank, true, :not_coprime)
    return BoundCertificate(acyclic, coprime, full_rank, true, :class_acyclic)
end

"""
    upper_bound_generators(s::Seed) → Vector

Generators of the upper bound `𝒰(Σ)` of the initial seed `s`.

For a coprime acyclic seed the upper and lower bounds coincide, so the answer
is [`lower_bound_generators`](@ref): the `2n` variables `xₖ`, `xₖ′`.  For any
other seed this throws - computing generators of `𝒰(Σ)` in the cyclic case
needs a Groebner basis over the defining ideal and is out of scope for this
package.  Call [`bound_certificate`](@ref) to see which hypothesis failed.
"""
function upper_bound_generators(s::Seed)
    cert = bound_certificate(s)
    cert.conclusion === :seed_acyclic || throw(InvalidArgument(
        "upper_bound_generators: generators of the upper bound are known only " *
        "for a coprime acyclic seed, and `bound_certificate` reports " *
        ":$(cert.conclusion) (acyclic = $(cert.acyclic), coprime = $(cert.coprime)); " *
        "the cyclic case needs a Groebner basis and is out of scope here"))
    return lower_bound_generators(s)
end
