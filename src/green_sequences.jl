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

# Entries above this bound get their branch pruned: 2^62 keeps every product
# of two in-bound entries inside Int128, so the mutation arithmetic in the
# green-sequence searches below stays exact.
const _MGS_ENTRY_LIMIT = Int128(2)^62

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
    B = Int128.(s.quiver.B[1:n, 1:n])
    if s isa Seed{PrincipalCoefficients}
        C = Int128.(cmatrix(s))
    else
        C = zeros(Int128, n, n); for i in 1:n; C[i, i] = 1; end
    end
    results = Vector{Vector{Int}}()

    # DFS on (B, C) integer matrices only — greenness and the mutation
    # recurrences depend solely on these, so no symbolic cluster arithmetic
    # is needed.  Int128 with the entry guard below, for the same reason as
    # `mgs_search`: silent wraparound can fabricate an all-red state, i.e. a
    # sequence that is not actually green.  Guard-pruned branches can only
    # lose sequences, never invent them.
    # Stack entries: (B, C, mutation sequence leading here).
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
            state = (_mutate_matrix(B, k), _mutate_C(C, B, k))
            max(maximum(abs, state[1]), maximum(abs, state[2])) > _MGS_ENTRY_LIMIT && continue
            push!(stack, (state[1], state[2], push!(copy(seq), k)))
        end
    end

    return results
end

"""
    mgs_search(s::Seed; max_length::Int=100, max_nodes::Int=10^6)
        → (status = ::Symbol, min_length = ::Union{Int,Nothing},
           sequence = ::Union{Vector{Int},Nothing})

Decide whether `s` admits a maximal green sequence, by breadth-first search on
`(B, C)` states with visited-state deduplication.

Because the search is breadth-first, a found sequence has **minimal length**
among all MGS of `s`. `status` is one of:

- `:found` — an MGS exists; `sequence` is a minimal-length one (the sequence
  itself is the certificate).
- `:none` — the reachable state space was exhausted below both budgets and
  contains no all-red state: certified non-existence. Note that by sign
  coherence a state with no green vertex is all-red — every dead end of the
  search *is* an MGS — so `:none` requires the deduplicated state graph to
  be exhausted without any dead end, which for well-behaved quivers does not
  occur: non-existence (e.g. the Markov quiver) typically manifests as an
  *infinite* search, i.e. a budget hit at every finite budget.
- `:none_within_length` — `max_nodes` was never hit and the search was
  exhausted except for branches pruned at `max_length`: certified
  non-existence of any MGS of length ≤ `max_length` (BFS reaches every state
  at its minimal depth, so the pruning cannot hide a shorter sequence).
  Longer MGS may still exist. This is the exact label for the bounded
  question "does an MGS of length ≤ L exist?".
- `:unknown` — more than `max_nodes` states were expanded, **or** a branch was
  pruned because its matrix entries exceeded the overflow guard (see below);
  no conclusion.

`min_length` and `sequence` are `nothing` unless `status == :found`.

The search state is `Int128` with branches pruned once any `|entry|` exceeds
`2^62` (so all mutation arithmetic stays exact). C-matrix entries grow
doubly-exponentially on wild quivers, and silent `Int64` wraparound can turn
an overflowed column "red" and fabricate an all-red state — a real bug this
guard fixed: `[0 2 0 2 0; -2 0 2 -2 -2; 0 -2 0 2 0; -2 2 -2 0 -2; 0 2 0 2 0]`
produced a fake length-10 MGS under `Int64` (BigInt ground truth: no MGS of
length ≤ 10; entries reach ~2^66). Overflow-pruning can only lose sequences,
never invent them, so a `:found` certificate is always exact; when pruning
occurred and nothing was found the result is honestly `:unknown`.

Unlike [`maximal_green_sequences`](@ref), which enumerates sequences by
depth-first search, this function only decides existence and minimality —
use it when the answer, not the list, is wanted.
"""
function mgs_search(s::Seed; max_length::Int = 100, max_nodes::Int = 10^6)
    n = s.quiver.n_mutable
    B0 = Int128.(s.quiver.B[1:n, 1:n])
    if s isa Seed{PrincipalCoefficients}
        C0 = Int128.(cmatrix(s))
    else
        C0 = zeros(Int128, n, n); for i in 1:n; C0[i, i] = 1; end
    end

    notfound = (status = :unknown, min_length = nothing, sequence = nothing)

    # BFS queue with a head index (O(1) dequeue); states deduplicated by
    # (B, C).  BFS reaches every state first at its minimal depth, so the
    # first all-red state yields a minimal-length MGS.
    queue = [(B0, C0, Int[])]
    visited = Set{Tuple{Matrix{Int128},Matrix{Int128}}}([(B0, C0)])
    head = 1
    truncated = false
    overflowed = false

    while head <= length(queue)
        length(queue) > max_nodes && return notfound
        B, C, seq = queue[head]
        head += 1

        greens = [k for k in 1:n if all(>=(0), view(C, :, k))]
        if isempty(greens)
            all(k -> all(<=(0), view(C, :, k)), 1:n) &&
                return (status = :found, min_length = length(seq), sequence = seq)
            continue    # sign-incoherent dead end (defensive; not an MGS)
        end

        if length(seq) >= max_length
            truncated = true
            continue
        end

        for k in greens
            state = (_mutate_matrix(B, k), _mutate_C(C, B, k))
            if max(maximum(abs, state[1]), maximum(abs, state[2])) > _MGS_ENTRY_LIMIT
                overflowed = true    # prune: can only lose sequences, never invent them
            elseif state ∉ visited
                push!(visited, state)
                push!(queue, (state[1], state[2], vcat(seq, k)))
            end
        end
    end

    overflowed && return notfound
    truncated && return (status = :none_within_length, min_length = nothing, sequence = nothing)
    return (status = :none, min_length = nothing, sequence = nothing)
end

"""
    verify_mutation_sequence(s::Seed, word::AbstractVector{<:Integer};
                             require_maximal::Bool = true)
        → (valid = ::Bool, maximal = ::Bool,
           charges = ::Union{Vector{Vector{BigInt}}, Nothing},
           reason = ::Union{String, Nothing})

Certify a candidate green sequence by exact replay: check that `word` mutates
only green vertices starting from `s`, and (with `require_maximal = true`)
that the sequence is maximal, i.e. ends in the all-red state.  When `valid`,
`charges` holds the c-vector of each mutated vertex recorded immediately
before its mutation — the ordered BPS charges of the chamber, identical in
meaning to [`ordered_c_vectors`](@ref).

This is the certification entry point for externally proposed sequences
(search heuristics, reinforcement-learning agents, hand computations): a
found sequence is a proof only if the verifier is trusted code, so the
replay runs on `Int128` `(B, C)` matrices with the same `2^62` entry guard
as [`mgs_search`](@ref) — and when the guard trips, the single word is
replayed again in `BigInt` (exact at any size, cheap for one sequence), so
the verifier has no ceiling and never refuses a genuine certificate.
Failure reasons name the first offending step.

Unlike [`ordered_c_vectors`](@ref), which throws on a bad sequence and
replays through the full symbolic seed, this function reports verdicts as
data and touches only integer matrices.

# Example
```julia
julia> verify_mutation_sequence(Seed(Quiver(:A, 2)), [2, 1]).valid
true

julia> verify_mutation_sequence(Seed(Quiver(:A, 2)), [1, 2]).reason
"sequence is green but not maximal: green vertices remain after step 2"
```
"""
function verify_mutation_sequence(s::Seed, word::AbstractVector{<:Integer};
                                  require_maximal::Bool = true)
    r = _verify_replay(Int128, s, word, require_maximal)
    # Guard trip on the fast path → exact BigInt replay of this one word.
    r === :overflow && return _verify_replay(BigInt, s, word, require_maximal)
    return r
end

function _verify_replay(::Type{T}, s::Seed, word::AbstractVector{<:Integer},
                        require_maximal::Bool) where {T <: Integer}
    n = s.quiver.n_mutable
    B = T.(s.quiver.B[1:n, 1:n])
    if s isa Seed{PrincipalCoefficients}
        C = T.(cmatrix(s))
    else
        C = zeros(T, n, n); for i in 1:n; C[i, i] = 1; end
    end

    reject(reason) = (valid = false, maximal = false, charges = nothing, reason = reason)

    charges = Vector{Vector{BigInt}}()
    for (t, k) in enumerate(word)
        1 <= k <= n || return reject("vertex $k out of range 1:$n at step $t")
        col = view(C, :, k)
        all(>=(0), col) || return reject("vertex $k not green at step $t")
        push!(charges, BigInt.(col))
        state = (_mutate_matrix(B, k), _mutate_C(C, B, k))
        if T !== BigInt &&
           max(maximum(abs, state[1]), maximum(abs, state[2])) > _MGS_ENTRY_LIMIT
            return :overflow
        end
        B, C = state
    end

    maximal = all(k -> all(<=(0), view(C, :, k)), 1:n)
    if require_maximal && !maximal
        return reject("sequence is green but not maximal: green vertices remain " *
                      "after step $(length(word))")
    end
    return (valid = true, maximal = maximal, charges = charges, reason = nothing)
end

verify_mutation_sequence(q::Quiver, word::AbstractVector{<:Integer}; kwargs...) =
    verify_mutation_sequence(Seed(q), word; kwargs...)

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
