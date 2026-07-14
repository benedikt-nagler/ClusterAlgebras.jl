# ─── Green / red vertex predicates ───────────────────────────────────────────

"""
    is_green(s::Seed{PrincipalCoefficients}, k::Int) → Bool

Return `true` if vertex `k` is green in `s`, meaning its c-vector is
entirely non-negative.  By sign coherence, every c-vector is either all ≥ 0
(green) or all ≤ 0 (red).
"""
is_green(s::Seed{PrincipalCoefficients}, k::Int) = all(>=(0), c_vector(s, k))

"""
    is_red(s::Seed{PrincipalCoefficients}, k::Int) → Bool

Return `true` if vertex `k` is red in `s`, meaning its c-vector is entirely
non-positive.
"""
is_red(s::Seed{PrincipalCoefficients}, k::Int) = all(<=(0), c_vector(s, k))

"""
    is_all_red(s::Seed{PrincipalCoefficients}) → Bool

Return `true` if every mutable vertex of `s` is red.  This is the terminal
condition for a maximal green sequence.
"""
function is_all_red(s::Seed{PrincipalCoefficients})
    all(k -> is_red(s, k), 1:s.quiver.n_mutable)
end

_green_vertices(s::Seed{PrincipalCoefficients}) =
    [k for k in 1:s.quiver.n_mutable if is_green(s, k)]

# ─── Maximal green sequences ──────────────────────────────────────────────────

"""
    maximal_green_sequences(s::Seed; max_length::Int=100, max_count::Int=1000)
                            → Vector{Vector{Int}}

Enumerate all maximal green sequences (MGS) starting from `s` via
depth-first search.

A *maximal green sequence* is a mutation sequence where every mutation is at a
green vertex (c-vector all ≥ 0) and the sequence is maximal: it terminates
when no green vertex remains (equivalently, all c-vectors are red).

If `s` has `TrivialCoefficients` it is automatically extended to principal
coefficients before the search.

The search prunes any branch whose length exceeds `max_length`, and stops
globally once `max_count` complete sequences have been collected.  For
mutation-finite quivers of moderate rank the defaults are sufficient; raise
`max_length` or lower `max_count` as needed.

Returns a `Vector{Vector{Int}}`; each inner vector lists the mutable vertex
indices mutated in order.
"""
function maximal_green_sequences(s::Seed; max_length::Int = 100, max_count::Int = 1000)
    n = s.quiver.n_mutable
    B = s.quiver.B[1:n, 1:n]
    if s isa Seed{PrincipalCoefficients}
        C = copy(cmatrix(s))
    else
        C = zeros(Int, n, n); for i in 1:n; C[i, i] = 1; end
    end
    results = Vector{Vector{Int}}()

    # DFS on (B, C) integer matrices only — greenness and the mutation
    # recurrences depend solely on these, so no symbolic cluster arithmetic
    # is needed.  Stack entries: (B, C, mutation sequence leading here).
    stack = [(B, C, Int[])]

    while !isempty(stack) && length(results) < max_count
        B, C, seq = pop!(stack)
        greens = [k for k in 1:n if all(>=(0), view(C, :, k))]

        if isempty(greens)
            all(k -> all(<=(0), view(C, :, k)), 1:n) && push!(results, seq)
            continue
        end

        length(seq) >= max_length && continue

        for k in greens
            push!(stack, (_mutate_matrix(B, k), _mutate_C(C, B, k), push!(copy(seq), k)))
        end
    end

    return results
end

# Sign type of one c-vector column: c-vectors are sign-coherent (never mixed for
# a valid seed), but classify defensively so externally-built seeds render too.
function _cvector_sign(col)
    all(>=(0), col) && return :green
    all(<=(0), col) && return :red
    return :mixed
end

"""
    green_sequence_signs(s::Seed, seq::AbstractVector{<:Integer}) → Matrix{Symbol}

Track the sign type of every c-vector as the mutation sequence `seq` is applied
starting from `s`.  Returns a `(length(seq)+1) × n` matrix (`n = s.quiver.n_mutable`)
of `:green` (c-vector all ≥ 0), `:red` (all ≤ 0), or `:mixed`.  Row `t` is the
state *before* the `t`-th mutation; the final row is the state after the whole
sequence.

If `s` has `TrivialCoefficients` it is first extended to principal coefficients
(where c-vectors live).  For a maximal green sequence the first row is all
`:green` and the last all `:red` — the "all green → all red" theorem made into a
matrix.  Rendered by `plot_green_sequence` in the Makie extension.
"""
function green_sequence_signs(s::Seed, seq::AbstractVector{<:Integer})
    sp = s isa Seed{PrincipalCoefficients} ? s : extend(s)
    n = sp.quiver.n_mutable
    signs = Matrix{Symbol}(undef, length(seq) + 1, n)
    for k in 1:n
        signs[1, k] = _cvector_sign(c_vector(sp, k))
    end
    for (t, j) in enumerate(seq)
        sp = mutate(sp, j)
        for k in 1:n
            signs[t + 1, k] = _cvector_sign(c_vector(sp, k))
        end
    end
    return signs
end
