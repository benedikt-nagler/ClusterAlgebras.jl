# ─── Ordered c-vectors along a maximal green sequence ─────────────────────────

"""
    ordered_c_vectors(s::Seed, seq::AbstractVector{Int}) → Vector{Vector{Int}}
    ordered_c_vectors(q::Quiver, seq)

Replay the maximal green sequence `seq` on `s` and return the c-vector of each
mutated vertex, recorded immediately *before* its mutation.

Physically these are the BPS charges `γ₁, …, γ_ℓ` of the finite chamber
described by the sequence, listed in phase order; each occurrence carries
DT invariant `Ω(γ) = 1` (see [`omega`](@ref)).  A seed with
`TrivialCoefficients` is automatically extended to principal coefficients.

Throws `InvalidArgument` if `seq` mutates a non-green vertex, or is green but
not maximal (green vertices remain at the end) - the charge list is only a
chamber's BPS spectrum for a *maximal* green sequence.

# Example
```julia
julia> ordered_c_vectors(Quiver(:A, 2), [1, 2, 1])
3-element Vector{Vector{Int64}}:
 [1, 0]
 [1, 1]
 [0, 1]
```
"""
function ordered_c_vectors(s::Seed, seq::AbstractVector{Int})
    s isa Seed{PrincipalCoefficients} || (s = extend(s))
    return _replay_green(s, seq).charges
end

ordered_c_vectors(q::Quiver, seq::AbstractVector{Int}) =
    ordered_c_vectors(extend(Seed(q)), seq)

# ─── The quantum dilogarithm word ─────────────────────────────────────────────

"""
    QuantumDilogWord

The ordered quantum-dilogarithm product `𝔼_{q_1}(ŷ^γ₁) ⋯ 𝔼_{q_ℓ}(ŷ^γ_ℓ)`
attached to a maximal green sequence, stored symbolically.  Fields:

- `charges::Vector{Vector{Int}}` - the ordered c-vectors `γ₁, …, γ_ℓ`
  (see [`ordered_c_vectors`](@ref)); `charges[t]` belongs to the `t`-th factor,
  leftmost first.
- `skew::Matrix{Int}` - the skew form `Λ = −D·B` on the charge lattice (`B` the
  mutable principal part of the initial exchange matrix, `D = diag(d)` the
  symmetrizer), defining the quantum torus `ŷ^α ŷ^β = q^{Λ(α,β)} ŷ^β ŷ^α`.
  `D·B` is skew for every skew-symmetrizable `B`, and the overall sign is the
  one that makes the sequence-ordered product the DT invariant in this package's
  B-matrix/c-vector conventions - pinned by the pentagon test.  For a
  skew-symmetric quiver (`d ≡ 1`) it reduces to `Bᵀ`.
- `weights::Vector{Int}` - aligned with `charges`: factor `t` is the
  `q^{d_k}`-dilogarithm, `k` the vertex mutated at step `t`.  All `1` for a
  skew-symmetric quiver; for a merely skew-symmetrizable one the weighting is
  what makes wall-crossing invariance hold (the unweighted product fails it
  already on `B₂`).
- `sequence::Vector{Int}` - the maximal green sequence the word came from.
- `shifts::Vector{Int}` - aligned with `charges`: factor `t` has argument
  `v^{shifts[t]} ŷ^γ` (`v = q^{1/2}`).  All `0` for a word built from a green
  sequence; a nonzero shift is how a state of nontrivial *refined* index
  `Ω(γ, y)` splits into dilogarithm factors (a vector multiplet has
  `Ω = −y − y^{-1}`, i.e. shifts `±1` at exponent `−1`).
- `exponents::Vector{Int}` - aligned with `charges`: the integer power the factor
  is raised to.  All `1` for a green sequence; `−1` is what a vector multiplet
  needs.
- `signs::Vector{Int}` - aligned with `charges`: the sign `σ = ±1` in the argument
  `σ v^{shift} ŷ^γ`, i.e. the quadratic refinement.  All `+1` for a green sequence;
  `−1` is what the imaginary root of an affine quiver carries (the `SU(2)`
  W-boson - see the `ExactWKB.jl` bridge).

By the Reineke/Keller theorem the evaluated product ([`ks_dilog_product`](@ref))
is the refined DT invariant of the quiver - independent of which maximal green
sequence is chosen.  Constructed by [`quantum_dilog_word`](@ref).
"""
struct QuantumDilogWord
    charges::Vector{Vector{Int}}
    skew::Matrix{Int}
    sequence::Vector{Int}
    weights::Vector{Int}
    shifts::Vector{Int}
    exponents::Vector{Int}
    signs::Vector{Int}
end

QuantumDilogWord(charges::Vector{Vector{Int}}, skew::Matrix{Int},
                 sequence::Vector{Int}, weights::Vector{Int},
                 shifts::Vector{Int}, exponents::Vector{Int}) =
    QuantumDilogWord(charges, skew, sequence, weights, shifts, exponents,
                     ones(Int, length(charges)))

# Unshifted, unit-exponent default: every factor is a plain q-dilogarithm.
QuantumDilogWord(charges::Vector{Vector{Int}}, skew::Matrix{Int},
                 sequence::Vector{Int}, weights::Vector{Int}) =
    QuantumDilogWord(charges, skew, sequence, weights,
                     zeros(Int, length(charges)), ones(Int, length(charges)))

# Skew-symmetric default: every factor is the plain q-dilogarithm.
QuantumDilogWord(charges::Vector{Vector{Int}}, skew::Matrix{Int},
                 sequence::Vector{Int}) =
    QuantumDilogWord(charges, skew, sequence, ones(Int, length(charges)))

function Base.show(io::IO, w::QuantumDilogWord)
    factors = join(("E$(d == 1 ? "" : "_q^$d")($(σ == 1 ? "" : "−")$(s == 0 ? "" : "v^$s ")ŷ^$(γ))" *
                    (e == 1 ? "" : "^$e")
                    for (γ, d, s, e, σ) in
                        zip(w.charges, w.weights, w.shifts, w.exponents, w.signs)), "·")
    print(io, "QuantumDilogWord ", factors)
end

"""
    quantum_dilog_word(s::Seed, seq::AbstractVector{Int}) → QuantumDilogWord
    quantum_dilog_word(q::Quiver, seq)

Build the [`QuantumDilogWord`](@ref) of the maximal green sequence `seq`: the
ordered charges from [`ordered_c_vectors`](@ref) together with the skew form
`Λ = −D·B` of the initial exchange matrix and the per-factor dilogarithm weights
`d_k`.  Validation is as in `ordered_c_vectors`.

Skew-symmetrizable quivers (`B₂`, `G₂`, `F₄`, …) are supported: `D·B` is skew for
all of them, and factor `t` carries the `q^{d_k}`-dilogarithm for the vertex `k`
mutated at that step.  For a skew-symmetric quiver this is exactly the older
`Λ = Bᵀ` with unit weights.
"""
function quantum_dilog_word(s::Seed, seq::AbstractVector{Int})
    s isa Seed{PrincipalCoefficients} || (s = extend(s))
    n = s.quiver.n_mutable
    Bm = s.quiver.B[1:n, 1:n]
    d = s.quiver.d
    Λ = [-d[i] * Bm[i, j] for i in 1:n, j in 1:n]   # −D·B, skew by skew-symmetrizability
    @assert Λ == -permutedims(Λ)
    charges = _replay_green(s, seq).charges
    return QuantumDilogWord(charges, Λ, collect(seq), [d[k] for k in seq])
end

quantum_dilog_word(q::Quiver, seq::AbstractVector{Int}) =
    quantum_dilog_word(extend(Seed(q)), seq)

"""
    omega(w::QuantumDilogWord) → Dict{Vector{Int}, Int}

The DT invariants `Ω(γ)` read off the word: the multiplicity of each charge
among `w.charges`.  Every factor of a maximal green sequence contributes `1`;
in finite type each charge occurs exactly once per chamber, so all values
are `1`.
"""
function omega(w::QuantumDilogWord)
    Ω = Dict{Vector{Int}, Int}()
    for γ in w.charges
        Ω[γ] = get(Ω, γ, 0) + 1
    end
    return Ω
end

# ─── Truncated quantum-torus evaluation ───────────────────────────────────────

# Skew form Λ(α, β) = αᵀ (w.skew) β of the word's quantum torus.
_skew(B::Matrix{Int}, α::Vector{Int}, β::Vector{Int}) =
    sum(α[i] * B[i, j] * β[j] for i in eachindex(α), j in eachindex(β))

# Multiply two truncated quantum-torus elements (Dict: exponent vector →
# coefficient in Frac(ZZ[v]), v = q^{1/2}), normal-ordering monomials via
# ŷ^α · ŷ^β = v^⟨α,β⟩ ŷ^(α+β) and dropping total y-degree > D.
function _qt_mul(a::Dict{Vector{Int}, T}, b::Dict{Vector{Int}, T},
                 B::Matrix{Int}, v, D::Int) where {T}
    c = Dict{Vector{Int}, T}()
    for (α, ca) in a, (β, cb) in b
        sum(α) + sum(β) <= D || continue
        coeff = ca * cb * v^_skew(B, α, β)
        γ = α + β
        acc = get(c, γ, nothing)
        newc = acc === nothing ? coeff : acc + coeff
        iszero(newc) ? delete!(c, γ) : (c[γ] = newc)
    end
    return c
end

# Inverse of a truncated quantum-torus element with constant term 1: the positive-degree
# part is nilpotent modulo degree D, so the geometric series terminates.
function _qt_inv(a::Dict{Vector{Int}, T}, B::Matrix{Int}, v, D::Int) where {T}
    n = size(B, 1)
    unit = zeros(Int, n)
    isone(get(a, unit, zero(v))) || throw(InvalidArgument(
        "only a quantum-torus element with constant term 1 can be inverted here"))
    negrest = Dict{Vector{Int}, T}(α => -c for (α, c) in a if α != unit)   # a = 1 − (−r)
    inv = Dict{Vector{Int}, T}(unit => one(v))
    term = Dict{Vector{Int}, T}(unit => one(v))
    for _ in 1:D                                    # (1 + r)⁻¹ = Σ_{k≥0} (−r)ᵏ, deg r ≥ 1
        term = _qt_mul(term, negrest, B, v, D)
        isempty(term) && break
        for (α, c) in term
            acc = get(inv, α, nothing)
            newc = (acc === nothing ? c : acc + c)
            iszero(newc) ? delete!(inv, α) : (inv[α] = newc)
        end
    end
    return inv
end

# `a^e` for any integer `e`, by repeated multiplication (and one inversion when e < 0).
function _qt_pow(a::Dict{Vector{Int}, T}, e::Int, B::Matrix{Int}, v, D::Int) where {T}
    e == 1 && return a
    base = e < 0 ? _qt_inv(a, B, v, D) : a
    out = Dict{Vector{Int}, T}(zeros(Int, size(B, 1)) => one(v))
    for _ in 1:abs(e)
        out = _qt_mul(out, base, B, v, D)
    end
    return out
end

# One factor 𝔼_{q_k}(v^s ŷ^γ) truncated at total y-degree D, in the convention
#   𝔼_{q_k}(y) = Σ_{n≥0} (−q_k^{1/2})ⁿ yⁿ / ∏_{i=1}^n (1 − q_kⁱ)
#              = Σ_{n≥0} q_k^{n/2} yⁿ / ∏_{i=1}^n (q_kⁱ − 1),   q_k = q^{d_k} = v^{2d_k}.
# The q^{n/2} (not q^{n²/2}) power is what the pentagon test pins for the
# normal-ordered monomials used here, and the weight d_k (which only rescales the
# series parameter, never the normal ordering) is what the B₂/G₂ tests pin.
# The shift `s` and sign `σ` are the substitution y → σ v^s y: term n picks up σⁿ v^{ns}.
function _dilog_factor(γ::Vector{Int}, F, v, D::Int, dk::Int, shift::Int = 0,
                       sign::Int = 1)
    e = Dict{Vector{Int}, elem_type(F)}(zero(γ) => one(F))
    vk = v^dk
    d = sum(γ)
    denom = one(F)
    for n in 1:(d > 0 ? D ÷ d : 0)
        denom *= vk^(2n) - 1
        e[n * γ] = divexact(vk^n, denom) * v^(n * shift) * sign^n
    end
    return e
end

# The coefficient field Frac(ℤ[v]) shared by every truncated evaluation here.
function _ks_ring(truncation_degree::Int)
    truncation_degree >= 1 ||
        throw(InvalidArgument("truncation_degree must be ≥ 1"))
    R, v0 = polynomial_ring(ZZ, "v")
    F = fraction_field(R)
    return R, F(v0)
end

"""
    ks_dilog_product(w::QuantumDilogWord; truncation_degree::Int=6)
        → Dict{Vector{Int}, coefficient}

Evaluate the ordered product `𝔼_{q_1}(ŷ^γ₁) ⋯ 𝔼_{q_ℓ}(ŷ^γ_ℓ)` in the quantum
affine space truncated at total `ŷ`-degree `truncation_degree`.  Monomials are
normal-ordered `ŷ^α` with `ŷ^α ŷ^β = v^{Λ(α,β)} ŷ^{α+β}`, where `Λ` is the skew
form stored in `w` and `v = q^{1/2}`; factor `t` is the quantum dilogarithm
series `𝔼_{q_t}(y) = Σ_{n≥0} (−q_t^{1/2})ⁿ yⁿ / ∏_{i=1}^n (1 − q_tⁱ)` at
`q_t = q^{d_k}` (`w.weights[t]`, trivial in the skew-symmetric case), and
factors multiply in sequence order (first mutation leftmost).

A factor may additionally carry an argument shift and an integer exponent,
`𝔼_{q_t}(v^{s_t} ŷ^{γ_t})^{e_t}` (`w.shifts`, `w.exponents`; `0`/`1` for a word
built from a green sequence).  These are what a state of nontrivial **refined**
index needs: a vector multiplet, `Ω(γ, y) = −y − y^{-1}`, is the pair of factors
`𝔼(−v ŷ^γ)^{-1} 𝔼(−v^{-1} ŷ^γ)^{-1}` (`w.signs` carries the `−`, the quadratic
refinement).  A negative exponent is evaluated by the terminating geometric series
(the factor's constant term is `1`).

Returns a `Dict` mapping exponent vectors `α` to coefficients in
`Frac(ℤ[v])` (zero coefficients removed).  Kontsevich–Soibelman wall-crossing
(Reineke/Keller): all maximal green sequences of one quiver give the **same**
product - for `A₂` this equality is the pentagon identity, for `B₂` the hexagon
and for `G₂` the octagon one.  These conventions (product order, sign of the skew
form, `q^{1/2}` powers, the `q^{d_k}` weighting) are pinned by the corresponding
tests in `test/test_ks_dilog.jl`.
"""
function ks_dilog_product(w::QuantumDilogWord; truncation_degree::Int = 6)
    R, v0 = _ks_ring(truncation_degree)
    F = parent(v0)
    v = v0
    n = size(w.skew, 1)
    result = Dict{Vector{Int}, elem_type(F)}(zeros(Int, n) => one(F))
    for (t, γ) in enumerate(w.charges)
        factor = _dilog_factor(γ, F, v, truncation_degree, w.weights[t],
                               w.shifts[t], w.signs[t])
        factor = _qt_pow(factor, w.exponents[t], w.skew, v, truncation_degree)
        result = _qt_mul(result, factor, w.skew, v, truncation_degree)
    end
    return result
end

"""
    ks_dilog_adjoint(w::QuantumDilogWord, μ::AbstractVector{Int};
                     truncation_degree::Int=6) → Dict{Vector{Int}, coefficient}

The image of the monomial `ŷ^μ` under conjugation by the word's product
`P =` [`ks_dilog_product`](@ref)`(w)`, i.e. `P ŷ^μ P⁻¹`, truncated at total
`ŷ`-degree `truncation_degree` (degrees are counted from `μ`, so the answer is
`ŷ^μ` times a series of degree `≤ truncation_degree`).

This is the Kontsevich–Soibelman *transformation*, and unlike the product itself
it has a classical limit: every coefficient is a Laurent polynomial in `v` and
[`ks_classical_limit`](@ref) evaluates it at `v = 1`, where a single factor acts
by `X_μ ↦ X_μ (1 + X_γ)^{Λ(γ,μ)}` - the Poisson automorphism the DDP jump uses.
The product's own coefficients diverge as `v → 1`, so the adjoint is the only
place the two pictures can be compared.
"""
function ks_dilog_adjoint(w::QuantumDilogWord, μ::AbstractVector{Int};
                          truncation_degree::Int = 6)
    _, v = _ks_ring(truncation_degree)
    F = parent(v)
    D = truncation_degree + sum(μ)
    P = ks_dilog_product(w; truncation_degree = D)
    mono = Dict{Vector{Int}, elem_type(F)}(collect(Int, μ) => one(F))
    # P ŷ^μ P⁻¹, not the other order: this is the direction in which a single
    # factor acts by (1 + X_γ)^{+Λ(γ,μ)}, matching the DDP jump (pinned in
    # test/test_ks_dilog.jl - the opposite order gives the inverse exponent).
    out = _qt_mul(P, mono, w.skew, v, D)
    return _qt_mul(out, _qt_inv(P, w.skew, v, D), w.skew, v, D)
end

"""
    ks_classical_limit(series::Dict) → Dict{Vector{Int}, Rational{BigInt}}

Specialize a truncated quantum-torus element at `v = 1` (`q → 1`), returning the
commutative Poisson-torus series.  Throws `InvalidArgument` on a coefficient with
a pole there - which is the honest outcome for [`ks_dilog_product`](@ref) itself;
use it on [`ks_dilog_adjoint`](@ref), whose coefficients are regular at `v = 1`.
"""
function ks_classical_limit(series::Dict{Vector{Int}, T}) where {T}
    out = Dict{Vector{Int}, Rational{BigInt}}()
    for (α, c) in series
        den = evaluate(denominator(c), one(ZZ))
        iszero(den) && throw(InvalidArgument(
            "coefficient of ŷ^$α has a pole at v = 1: this series has no classical limit"))
        val = Rational{BigInt}(BigInt(evaluate(numerator(c), one(ZZ))), BigInt(den))
        iszero(val) || (out[α] = val)
    end
    return out
end
