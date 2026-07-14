# Cartan companion, root systems, and finite/affine type detection (Track A: P2 + A2).

# ─── Cartan companion ─────────────────────────────────────────────────────────

"""
    cartan_companion(q::Quiver) -> Matrix{Int}

Return the Cartan companion matrix of the mutable block of `q`.
Diagonal entries are 2; off-diagonal entries are `A[i,j] = -|B[i,j]|`.
"""
function cartan_companion(q::Quiver)
    n = q.n_mutable
    A = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n
        A[i, j] = i == j ? 2 : -abs(q.B[i, j])
    end
    return A
end

# ─── BFS for positive roots ───────────────────────────────────────────────────
#
# Roots are integer vectors in the simple-root basis.
# Criterion: β + e_k is a root when p − c ≥ 1, where
#   p = max{r ≥ 0 : β − r·e_k is already a discovered root}
#   c = A[k,:] · β  (the k-th row of A dotted with β)
#
# The simple "c < 0 ⟹ extend" rule is insufficient for non-simply-laced types
# (e.g. G₂); the p-correction is required.

function _compute_positive_roots(A::Matrix{Int})
    n = size(A, 1)
    root_set = Set{Vector{Int}}()
    queue    = Vector{Vector{Int}}()
    for i in 1:n
        ei = zeros(Int, n); ei[i] = 1
        push!(root_set, ei)
        push!(queue,    ei)
    end
    idx = 1
    while idx <= length(queue)
        β = queue[idx]; idx += 1
        for k in 1:n
            # compute p
            p = 0
            while β[k] - p - 1 >= 0
                γ = copy(β); γ[k] -= (p + 1)
                γ ∉ root_set && break
                p += 1
            end
            # c = A[k,:] · β
            c = 0; for j in 1:n; c += A[k, j] * β[j]; end
            # extend if q = p − c ≥ 1
            if p - c >= 1
                γ = copy(β); γ[k] += 1
                if γ ∉ root_set
                    push!(root_set, γ)
                    push!(queue,    γ)
                end
            end
        end
    end
    return sort!(collect(root_set))
end

# ─── Hardcoded Lie-theory data ────────────────────────────────────────────────

function _coxeter_number(type::Symbol, n::Int)
    type === :A && return n + 1
    type === :B && return 2n
    type === :C && return 2n
    type === :D && return 2(n - 1)
    type === :E && n == 6 && return 12
    type === :E && n == 7 && return 18
    type === :E && n == 8 && return 30
    type === :F && return 12
    type === :G && return 6
    throw(InvalidArgument("unknown finite Dynkin type :$type (n=$n)"))
end

function _exponents(type::Symbol, n::Int)
    type === :A && return collect(1:n)
    type === :B && return collect(1:2:2n-1)
    type === :C && return collect(1:2:2n-1)
    type === :D && return sort!(vcat(collect(1:2:2n-3), [n - 1]))
    type === :E && n == 6 && return [1, 4, 5, 7, 8, 11]
    type === :E && n == 7 && return [1, 5, 7, 9, 11, 13, 17]
    type === :E && n == 8 && return [1, 7, 11, 13, 17, 19, 23, 29]
    type === :F && return [1, 5, 7, 11]
    type === :G && return [1, 5]
    throw(InvalidArgument("unknown finite Dynkin type :$type (n=$n)"))
end

function _weyl_group_order(type::Symbol, n::Int)
    type === :A && return factorial(n + 1)
    type === :B && return 2^n  * factorial(n)
    type === :C && return 2^n  * factorial(n)
    type === :D && return 2^(n-1) * factorial(n)
    type === :E && n == 6 && return 51840
    type === :E && n == 7 && return 2903040
    type === :E && n == 8 && return 696729600
    type === :F && return 1152
    type === :G && return 12
    throw(InvalidArgument("unknown finite Dynkin type :$type (n=$n)"))
end

# ─── RootSystem ───────────────────────────────────────────────────────────────

"""
    RootSystem

Root system of a finite Dynkin type, computed from the standard Cartan matrix.

Fields:
- `type::Symbol` — one of `:A, :B, :C, :D, :E, :F, :G`
- `n::Int` — rank
- `positive_roots::Vector{Vector{Int}}` — positive roots in simple-root coordinates (sorted)
- `coxeter_number::Int` — Coxeter number `h`
- `exponents::Vector{Int}` — fundamental exponents `e_1 ≤ … ≤ e_n` (sorted)
- `weyl_group_order::Int` — order of the Weyl group
"""
struct RootSystem
    type::Symbol
    n::Int
    positive_roots::Vector{Vector{Int}}
    coxeter_number::Int
    exponents::Vector{Int}
    weyl_group_order::Int
end

"""
    RootSystem(type::Symbol, n::Int) -> RootSystem

Construct the root system of the given finite Dynkin type and rank.
Positive roots are computed via BFS on the Cartan matrix.
Coxeter number, exponents, and Weyl group order come from classical tables.
"""
function RootSystem(type::Symbol, n::Int)
    B, _ = _dynkin_exchange_matrix(type, n)
    A = Matrix{Int}(undef, n, n)
    for i in 1:n, j in 1:n; A[i, j] = i == j ? 2 : -abs(B[i, j]); end
    return RootSystem(
        type, n,
        _compute_positive_roots(A),
        _coxeter_number(type, n),
        _exponents(type, n),
        _weyl_group_order(type, n),
    )
end

# ─── Almost-positive roots ────────────────────────────────────────────────────

"""
    almost_positive_roots(rs::RootSystem) -> Vector{Vector{Int}}

Return the almost-positive roots: the `n` negative simple roots `{−α_i}` together
with all positive roots, as vectors in the simple-root basis.
`−α_i` is encoded as the vector with `−1` at position `i` and `0` elsewhere.
"""
function almost_positive_roots(rs::RootSystem)
    n = rs.n
    neg_simples = [[-Int(i == j) for j in 1:n] for i in 1:n]
    return vcat(neg_simples, rs.positive_roots)
end

# ─── Acyclicity and acyclic-representative search ────────────────────────────
#
# Finite type is a mutation-class invariant, not a property of one seed.
# By Fomin–Zelevinsky II, A(B) is finite type iff some mutation-equivalent
# quiver is acyclic with positive-definite Cartan companion.  These helpers
# implement that search so is_finite_type / is_affine_type / cartan_type
# are correct on non-acyclic (e.g. oriented-3-cycle) input.

# True iff the mutable block of q has no oriented cycle (i.e. no directed cycle
# in the graph where i→j iff B[i,j]>0).
function _is_acyclic(q::Quiver)
    n = q.n_mutable
    # state: 0=unvisited, 1=on current DFS stack, 2=finished
    state = zeros(Int8, n)
    function dfs(v)
        state[v] = 1
        for u in 1:n
            q.B[v, u] > 0 || continue
            state[u] == 1 && return false   # back edge = directed cycle
            state[u] == 0 && !dfs(u) && return false
        end
        state[v] = 2
        return true
    end
    for v in 1:n
        state[v] == 0 && !dfs(v) && return false
    end
    return true
end

# BFS over the mutation class of q, returning the first acyclic quiver found.
# Returns `nothing` if the class is exhausted without finding one (a definitive
# answer: no acyclic representative exists), and `missing` if the search budget
# is hit first (inconclusive — callers must not silently report `false`).
# The dedup key is the LABELED mutable block, so the enumeration may revisit
# vertex-relabeled copies; the budget is sized to absorb that redundancy for
# the finite types (including E₇/E₈).
function _acyclic_representative(q::Quiver; max_quivers::Int = 30_000)
    _is_acyclic(q) && return q            # fast path: already acyclic
    seen  = Set{Matrix{Int}}()
    queue = Quiver[q]
    # Key on the mutable sub-block only: the mutable block mutates independently
    # of frozen vertices, so this correctly identifies distinct mutable states
    # while avoiding exponential growth from frozen-vertex interaction terms.
    _mut_key(qi) = qi.B[1:qi.n_mutable, 1:qi.n_mutable]
    push!(seen, _mut_key(q))
    while !isempty(queue)
        qi = popfirst!(queue)
        for k in 1:qi.n_mutable
            qk = mutate(qi, k)
            _mut_key(qk) ∈ seen && continue
            length(seen) >= max_quivers && return missing
            push!(seen, _mut_key(qk))
            _is_acyclic(qk) && return qk
            push!(queue, qk)
        end
    end
    return nothing                        # class exhausted, no acyclic rep
end

# Shared handling of an inconclusive acyclic-representative search.
function _warn_budget(fn::AbstractString)
    @warn "$fn: acyclic-representative search budget exhausted before the " *
          "mutation class was covered; returning false, but this may be a " *
          "false negative for deep non-acyclic quivers"
end

# ─── Exact integer determinant (cofactor expansion, n ≤ 8 in practice) ──────

function _det_int(M::Matrix{Int})
    n = size(M, 1)
    n == 0 && return 1
    n == 1 && return M[1, 1]
    n == 2 && return M[1,1]*M[2,2] - M[1,2]*M[2,1]
    d = 0
    for j in 1:n
        M[1, j] == 0 && continue
        cols = [c for c in 1:n if c != j]
        d += (-1)^(j + 1) * M[1, j] * _det_int(M[2:end, cols])
    end
    return d
end

# ─── Finite / affine type detection ──────────────────────────────────────────

function _symmetrized_cartan(q::Quiver)
    n = q.n_mutable
    A = cartan_companion(q)
    return [q.d[i] * A[i, j] for i in 1:n, j in 1:n]
end

"""
    is_finite_type(q::Quiver) -> Bool

Return `true` if `q` is of finite cluster type (finitely many cluster variables).

Finite type is a mutation-class invariant: by Fomin–Zelevinsky II, A(B) is
finite type iff some mutation-equivalent quiver is acyclic with positive-definite
symmetrized Cartan companion. This function searches the mutation class for such
a representative (up to 30 000 quivers) and applies Sylvester's criterion to it.

Named constructors `Quiver(:A, n)` etc. always produce acyclic quivers, so the
search terminates instantly for them.

Note: *finite type* (finitely many cluster variables) is strictly stronger than
*mutation-finite* (finitely many quivers in the class). The Markov quiver is
mutation-finite but not finite type.
"""
function is_finite_type(q::Quiver)
    n = q.n_mutable
    n == 0 && return true
    rep = _acyclic_representative(q)
    rep === missing && (_warn_budget("is_finite_type"); return false)
    rep === nothing && return false
    S = _symmetrized_cartan(rep)
    for k in 1:n
        _det_int(S[1:k, 1:k]) > 0 || return false
    end
    return true
end

"""
    is_affine_type(q::Quiver) -> Bool

Return `true` if `q` is of affine (tame) cluster type: the symmetrized Cartan
companion `D·A` of an acyclic mutation representative is positive semidefinite
of corank 1 (all leading principal minors of size 1:n−1 positive; determinant 0).

Searches the mutation class for an acyclic representative (up to 30 000 quivers).

Note: affine type differs from finite type; the Markov quiver is neither.
"""
function is_affine_type(q::Quiver)
    n = q.n_mutable
    n == 0 && return false
    rep = _acyclic_representative(q)
    rep === missing && (_warn_budget("is_affine_type"); return false)
    rep === nothing && return false
    S = _symmetrized_cartan(rep)
    for k in 1:n-1
        _det_int(S[1:k, 1:k]) > 0 || return false
    end
    return _det_int(S) == 0
end

# ─── Type recognition ─────────────────────────────────────────────────────────

"""
    cartan_type(q::Quiver) -> (Symbol, Int)

Return the Dynkin type `(type, rank)` of a finite-type quiver.
`type` is one of `:A, :B, :C, :D, :E, :F, :G` and `rank = q.n_mutable`.

Searches the mutation class for an acyclic representative (up to 30 000 quivers)
and reads the type from it.  Throws `InvalidArgument` if `q` is not of finite
type or if the type cannot be determined.
"""
function cartan_type(q::Quiver)
    n = q.n_mutable
    n >= 1 || throw(InvalidArgument("cartan_type requires n_mutable ≥ 1"))

    rep = _acyclic_representative(q)
    (rep === nothing || rep === missing) && throw(InvalidArgument(
        "cartan_type requires a finite-type quiver; no acyclic mutation " *
        "representative was found" *
        (rep === missing ? " within the search budget (inconclusive)" : "")))

    # Verify positive definiteness on the acyclic representative.
    S = _symmetrized_cartan(rep)
    for k in 1:n
        _det_int(S[1:k, 1:k]) > 0 || throw(InvalidArgument(
            "cartan_type requires a finite-type quiver; is_finite_type(q) is false"))
    end

    A        = cartan_companion(rep)
    n_pos    = length(_compute_positive_roots(A))
    max_bond = maximum(-A[i, j] for i in 1:n, j in 1:n if i != j; init=0)

    if max_bond <= 1
        n_pos == n * (n + 1) ÷ 2 && return (:A, n)
        n_pos == n * (n - 1)      && return (:D, n)
        n == 6                    && return (:E, 6)
        n == 7                    && return (:E, 7)
        n == 8                    && return (:E, 8)
    elseif max_bond == 2
        n == 4 && n_pos == 24 && return (:F, 4)
        if n_pos == n^2
            # Distinguish Bₙ from Cₙ by the symmetrizer multiset, which is
            # permutation-invariant (positional checks like d[1] vs d[end] are
            # not): Bₙ has d-multiset {big × (n−1), small × 1}, Cₙ the reverse.
            # B₂ ≅ C₂: the multiset {2,1} is symmetric, so fall back to the
            # positional convention of the named constructors (long root last
            # for B, short root last for C).  Either answer is correct.
            n == 2 && return rep.d[1] >= rep.d[2] ? (:B, 2) : (:C, 2)
            hi = maximum(rep.d)
            return count(==(hi), rep.d) == n - 1 ? (:B, n) : (:C, n)
        end
    elseif max_bond == 3
        n == 2 && n_pos == 6 && return (:G, 2)
    end

    throw(InvalidArgument(
        "unable to determine Cartan type (n_mutable=$n, n_pos=$n_pos, max_bond=$max_bond)"))
end
