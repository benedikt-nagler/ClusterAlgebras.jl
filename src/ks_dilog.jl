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
not maximal (green vertices remain at the end) — the charge list is only a
chamber's BPS spectrum for a *maximal* green sequence.

# Example
```julia
julia> ordered_c_vectors(Quiver(:A, 2), [2, 1, 2])
3-element Vector{Vector{Int64}}:
 [0, 1]
 [1, 1]
 [1, 0]
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

The ordered quantum-dilogarithm product `𝔼(ŷ^γ₁) ⋯ 𝔼(ŷ^γ_ℓ)` attached to a
maximal green sequence, stored symbolically.  Fields:

- `charges::Vector{Vector{Int}}` — the ordered c-vectors `γ₁, …, γ_ℓ`
  (see [`ordered_c_vectors`](@ref)); `charges[t]` belongs to the `t`-th factor,
  leftmost first.
- `skew::Matrix{Int}` — the skew form `Λ = Bᵀ` on the charge lattice (`B` the
  mutable principal part of the initial exchange matrix), defining the quantum
  torus `ŷ^α ŷ^β = q^{Λ(α,β)} ŷ^β ŷ^α`.  The sign (`Bᵀ`, not `B`) is the one
  that makes the sequence-ordered product the DT invariant in this package's
  B-matrix/c-vector conventions — pinned by the pentagon test.
- `sequence::Vector{Int}` — the maximal green sequence the word came from.

By the Reineke/Keller theorem the evaluated product ([`ks_dilog_product`](@ref))
is the refined DT invariant of the quiver — independent of which maximal green
sequence is chosen.  Constructed by [`quantum_dilog_word`](@ref).
"""
struct QuantumDilogWord
    charges::Vector{Vector{Int}}
    skew::Matrix{Int}
    sequence::Vector{Int}
end

function Base.show(io::IO, w::QuantumDilogWord)
    factors = join(("E(ŷ^$(γ))" for γ in w.charges), "·")
    print(io, "QuantumDilogWord ", factors)
end

"""
    quantum_dilog_word(s::Seed, seq::AbstractVector{Int}) → QuantumDilogWord
    quantum_dilog_word(q::Quiver, seq)

Build the [`QuantumDilogWord`](@ref) of the maximal green sequence `seq`: the
ordered charges from [`ordered_c_vectors`](@ref) together with the skew form of
the initial exchange matrix.  Validation is as in `ordered_c_vectors`.
"""
function quantum_dilog_word(s::Seed, seq::AbstractVector{Int})
    s isa Seed{PrincipalCoefficients} || (s = extend(s))
    n = s.quiver.n_mutable
    charges = _replay_green(s, seq).charges
    return QuantumDilogWord(charges, permutedims(s.quiver.B[1:n, 1:n]), collect(seq))
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

# One factor 𝔼(ŷ^γ) truncated at total y-degree D, in the convention
#   𝔼(y) = Σ_{n≥0} (−q^{1/2})ⁿ yⁿ / ∏_{i=1}^n (1 − qⁱ)
#        = Σ_{n≥0} q^{n/2} yⁿ / ∏_{i=1}^n (qⁱ − 1),   q = v².
# The q^{n/2} (not q^{n²/2}) power is what the pentagon test pins for the
# normal-ordered monomials used here.
function _dilog_factor(γ::Vector{Int}, F, v, D::Int)
    e = Dict{Vector{Int}, elem_type(F)}(zero(γ) => one(F))
    d = sum(γ)
    denom = one(F)
    for n in 1:(d > 0 ? D ÷ d : 0)
        denom *= v^(2n) - 1
        e[n * γ] = divexact(v^n, denom)
    end
    return e
end

"""
    ks_dilog_product(w::QuantumDilogWord; truncation_degree::Int=6)
        → Dict{Vector{Int}, coefficient}

Evaluate the ordered product `𝔼(ŷ^γ₁) ⋯ 𝔼(ŷ^γ_ℓ)` in the quantum affine space
truncated at total `ŷ`-degree `truncation_degree`.  Monomials are normal-ordered
`ŷ^α` with `ŷ^α ŷ^β = v^{Λ(α,β)} ŷ^{α+β}`, where `Λ` is the skew form stored in
`w` and `v = q^{1/2}`; each factor is the quantum dilogarithm series
`𝔼(y) = Σ_{n≥0} (−q^{1/2})ⁿ yⁿ / ∏_{i=1}^n (1 − qⁱ)`, and factors multiply in
sequence order (first mutation leftmost).

Returns a `Dict` mapping exponent vectors `α` to coefficients in
`Frac(ℤ[v])` (zero coefficients removed).  Kontsevich–Soibelman wall-crossing
(Reineke/Keller): all maximal green sequences of one quiver give the **same**
product — for `A₂` this equality is the pentagon identity.  These conventions
(product order, sign of the skew form, `q^{1/2}` powers) are pinned by the
pentagon test in `test/test_ks_dilog.jl`.
"""
function ks_dilog_product(w::QuantumDilogWord; truncation_degree::Int = 6)
    truncation_degree >= 1 ||
        throw(InvalidArgument("truncation_degree must be ≥ 1"))
    R, v0 = polynomial_ring(ZZ, "v")
    F = fraction_field(R)
    v = F(v0)
    n = size(w.skew, 1)
    result = Dict{Vector{Int}, elem_type(F)}(zeros(Int, n) => one(F))
    for γ in w.charges
        result = _qt_mul(result, _dilog_factor(γ, F, v, truncation_degree),
                         w.skew, v, truncation_degree)
    end
    return result
end
