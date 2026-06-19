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
    sp = s isa Seed{PrincipalCoefficients} ? s : extend(s)
    results = Vector{Vector{Int}}()

    # Stack entries: (current seed, mutation sequence leading here)
    stack = [(sp, Int[])]

    while !isempty(stack) && length(results) < max_count
        seed, seq = pop!(stack)
        greens = _green_vertices(seed)

        if isempty(greens)
            is_all_red(seed) && push!(results, seq)
            continue
        end

        length(seq) >= max_length && continue

        for k in greens
            push!(stack, (mutate(seed, k), push!(copy(seq), k)))
        end
    end

    return results
end
