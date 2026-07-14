# ─── σ extraction ─────────────────────────────────────────────────────────────

# Read the permutation σ off the final C-matrix of a maximal green sequence:
# C = −P_σ with (P_σ)[i, j] = δ_{i, σ[j]}, i.e. column j has a single nonzero
# entry, equal to −1, at row σ[j]. Return `nothing` if C is not minus a
# permutation matrix.
function _signed_permutation(C::Matrix{Int})
    n = size(C, 1)
    size(C, 2) == n || return nothing
    σ = Vector{Int}(undef, n)
    for j in 1:n
        rows = findall(!=(0), view(C, :, j))
        length(rows) == 1 || return nothing
        C[rows[1], j] == -1 || return nothing
        σ[j] = rows[1]
    end
    allunique(σ) || return nothing
    return σ
end

# ─── Green replay ─────────────────────────────────────────────────────────────

# Replay `seq` on a principal-coefficient seed, validating that it is a maximal
# green sequence: every step mutates a green vertex, and no green vertex remains
# at the end. Returns `(final, charges)` where `charges[t]` is the c-vector of
# the vertex mutated at step `t`, recorded *before* that mutation - the ordered
# BPS charges consumed by the quantum-dilogarithm export in `ks_dilog.jl`.
function _replay_green(s::Seed{PrincipalCoefficients}, seq::AbstractVector{Int})
    n = s.quiver.n_mutable
    current = s
    charges = Vector{Vector{Int}}()
    for (t, k) in enumerate(seq)
        1 <= k <= n || throw(InvalidVertex(k, n))
        is_green(current, k) ||
            throw(InvalidArgument("not a green sequence: vertex $k is not green " *
                                  "at step $t"))
        push!(charges, c_vector(current, k))
        current = mutate(current, k)
    end
    is_all_red(current) ||
        throw(InvalidArgument("the sequence is green but not maximal: green " *
                              "vertices remain after replaying it"))
    return (final = current, charges = charges)
end

# ─── DT transformation ────────────────────────────────────────────────────────

"""
    dt_transformation(s::Seed; seq=nothing, max_length=100) → NamedTuple
    dt_transformation(q::Quiver; kwargs...)

Compute the Donaldson–Thomas (DT) transformation of the initial seed `s` as a
birational map on the y-variables, by composing a maximal green sequence (MGS)
and un-permuting the result.

For any MGS `i = (i₁, …, iₗ)` the final C-matrix equals `−P_σ` for a
permutation `σ`, and the DT transformation is `σ⁻¹ ∘ μ_{iₗ} ∘ ⋯ ∘ μ_{i₁}`.
By Keller's theorem (*On cluster theory and quantum dilogarithm identities*,
2011) the result is independent of which MGS is chosen - this invariance is
covered by tests.

A seed with `TrivialCoefficients` is automatically extended to principal
coefficients. A `Seed{PrincipalCoefficients}` must be initial (C-matrix = Iₙ);
otherwise `InvalidArgument` is thrown, since the DT transformation is attached
to the initial seed. If `seq` is given it is validated to be an MGS (each step
green, terminal all-red); if `seq === nothing` the first MGS found by
[`maximal_green_sequences`](@ref) within `max_length` is used (which one is
irrelevant, by Keller invariance), and `InvalidArgument` is thrown when none
exists - e.g. for the Markov quiver, which admits no MGS.

Returns a NamedTuple `(sequence, sigma, seed, y_images)`:
- `sequence::Vector{Int}` - the MGS that was composed.
- `sigma::Vector{Int}` - the permutation, pinned by `C_final[σ[j], j] == -1`
  (the final c-vector of vertex `j` is `−e_{σ[j]}`); equivalently the final
  quiver satisfies `B_final[i, j] == B0[σ[i], σ[j]]`.
- `seed::Seed{PrincipalCoefficients}` - the final, un-relabeled seed; its
  C-/G-matrix, F-polynomials and tropical y-variables are available through
  the usual accessors.
- `y_images` - the DT map on the initial y-torus, `y_j ↦ y_images[j]`, as
  elements of `Frac(ZZ[y₁, …, yₙ])`: `y_images[j] = y_variables(seed)[σ⁻¹[j]]`,
  the final y-variable at the vertex whose c-vector is `−e_j`. Its
  tropicalization is `y_j⁻¹`.

Iterating the DT transformation means substituting `y_images` into itself -
**not** replaying the raw mutation sequence twice, which composes differently
whenever `σ ≠ id`.

# Example
```julia
julia> r = dt_transformation(Quiver(:A, 2); seq = [2, 1]);

julia> r.sigma
2-element Vector{Int64}:
 1
 2

julia> r.y_images   # DT: y₁ ↦ 1/(y₁(1+y₂)),  y₂ ↦ (1+y₁+y₁y₂)/y₂
2-element Vector{AbstractAlgebra.Generic.FracFieldElem{...}}:
 1//(y1*y2 + y1)
 (y1*y2 + y1 + 1)//y2
```
"""
function dt_transformation(s::Seed;
                           seq::Union{Nothing, AbstractVector{Int}} = nothing,
                           max_length::Int = 100)
    s isa Seed{PrincipalCoefficients} || (s = extend(s))
    n = s.quiver.n_mutable
    all(cmatrix(s)[i, j] == (i == j ? 1 : 0) for i in 1:n, j in 1:n) ||
        throw(InvalidArgument("dt_transformation requires the initial principal " *
                              "seed (C-matrix = identity); mutate after, not before"))

    if seq === nothing
        found = maximal_green_sequences(s; max_length, max_count = 1)
        isempty(found) &&
            throw(InvalidArgument("no maximal green sequence found within " *
                                  "max_length = $max_length; none may exist " *
                                  "(e.g. the Markov quiver has no MGS)"))
        seq = found[1]
    end

    final = _replay_green(s, seq).final
    σ = _signed_permutation(cmatrix(final))
    σ === nothing &&
        error("final C-matrix of a maximal green sequence is not −P_σ - " *
              "this contradicts the theorem and indicates a bug in the " *
              "C-matrix recurrence")

    # Relabel by σ⁻¹ so the result is a self-map of the initial y-torus: the
    # image at index j is the final y-variable at the vertex whose c-vector is
    # −e_j, hence tropicalizes to y_j⁻¹ (DT sends the C-matrix to −Iₙ).
    y = y_variables(final)
    σinv = invperm(σ)
    return (sequence = collect(seq), sigma = σ, seed = final,
            y_images = [y[σinv[j]] for j in 1:n])
end

dt_transformation(q::Quiver; kwargs...) = dt_transformation(extend(Seed(q)); kwargs...)
