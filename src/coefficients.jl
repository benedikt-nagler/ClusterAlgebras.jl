# ─── PrincipalData ────────────────────────────────────────────────────────────

"""
    PrincipalData{S, P}

Payload for `Seed{PrincipalCoefficients}`.  Tracks:
- `C`: C-matrix, columns are c-vectors (tropical y-variables).
- `yring`: the fraction field `Frac(ZZ[y1,…,yn])`.
- `yvars`: rational y-variables, elements of `yring`.
- `B0`: the INITIAL mutable exchange matrix (constant along mutation paths).
- `B_ext`: the CURRENT principal-extended 2n×2n exchange matrix.
- `fpolys`: the current F-polynomials, elements of `ZZ[y1,…,yn]`.

`B_ext` and `fpolys` are updated incrementally by one recurrence step per
mutation, so accessors never replay the whole mutation path.

G-vectors (g-matrix) are derived on demand from the cluster variables and
F-polynomials via the separation formula; they are not stored incrementally
because the update formula requires sign information that g-vectors (unlike
c-vectors) do not satisfy sign-coherence for.

Constructed via `extend(s::Seed{TrivialCoefficients})`.
"""
struct PrincipalData{S, P}
    C      :: Matrix{Int}
    yring  :: Any          # FracField{<:MPolyRing}
    yvars  :: Vector{S}    # elements of yring
    B0     :: Matrix{Int}  # initial mutable exchange matrix
    B_ext  :: Matrix{Int}  # current principal-extended matrix (2n × 2n)
    fpolys :: Vector{P}    # current F-polynomials in ZZ[y1,…,yn]
end

# ─── extend ───────────────────────────────────────────────────────────────────

"""
    extend(s::Seed{TrivialCoefficients}) → Seed{PrincipalCoefficients}

Attach principal-coefficient data (C = Iₙ, rational y-variables) to `s`.
The new seed wraps the SAME quiver, cluster, ring, and mutation path.
"""
function extend(s::Seed{TrivialCoefficients})
    n = s.quiver.n_mutable
    C = zeros(Int, n, n); for i in 1:n; C[i,i] = 1; end
    Ry, ygen = polynomial_ring(ZZ, ["y$i" for i in 1:n])
    yring    = fraction_field(Ry)
    yvars    = yring.(ygen)

    # Initial mutable exchange matrix: un-mutate the current one along the path.
    B0 = s.quiver.B[1:n, 1:n]
    for k in reverse(s.mutation_path); B0 = _mutate_matrix(B0, k); end

    # Principal-coefficient extended matrix at the initial seed, then replay the
    # path once to obtain the current B_ext and F-polynomials.  From here on,
    # `mutate` updates both incrementally (one recurrence step per mutation).
    B_ext = zeros(Int, 2n, 2n)
    B_ext[1:n, 1:n] = B0
    for i in 1:n
        B_ext[n+i, i] =  1
        B_ext[i, n+i] = -1
    end
    fpolys = fill(one(Ry), n)
    for k in s.mutation_path
        fpolys = _mutate_fpolys(fpolys, B_ext, ygen, k)
        B_ext  = _mutate_matrix(B_ext, k)
    end

    pd = PrincipalData(C, yring, yvars, B0, B_ext, fpolys)
    T  = eltype(s.cluster)
    F  = typeof(s.ring)
    S  = eltype(yvars)
    P  = eltype(fpolys)
    return Seed{PrincipalCoefficients, T, F, PrincipalData{S, P}}(
        s.quiver, s.cluster, s.ring, s.mutation_path, pd)
end

# ─── Mutation internals ────────────────────────────────────────────────────────

# ε_k ∈ {+1, -1}: whether the k-th c-vector is non-negative or non-positive.
# Sign coherence (a theorem) guarantees exactly one case holds; we default to
# +1 for the zero vector (initial identity diagonal is always positive).
# A mixed-sign c-vector means the C-matrix recurrence itself is broken, so we
# error loudly instead of silently returning −1.
function _epsilon(c::AbstractVector{<:Integer})
    all(>=(0), c) && return 1
    all(<=(0), c) && return -1
    error("sign coherence violated: c-vector $c has mixed signs - " *
          "this indicates a bug in the C-matrix recurrence")
end

# One step of the F-polynomial recurrence (FZ-IV Prop. 5.1) at vertex k,
# using the CURRENT principal-extended matrix B_ext (before mutating it).
function _mutate_fpolys(F::Vector{P}, B_ext::Matrix{Int}, ygens, k::Int) where {P}
    n  = length(F)
    R  = parent(F[1])
    y_pos = prod(ygens[i]^max( B_ext[n+i, k], 0) for i in 1:n; init=one(R))
    y_neg = prod(ygens[i]^max(-B_ext[n+i, k], 0) for i in 1:n; init=one(R))

    M_pos = y_pos
    M_neg = y_neg
    for j in 1:n
        j == k && continue
        b = B_ext[j, k]
        b > 0 && (M_pos *= F[j]^b)
        b < 0 && (M_neg *= F[j]^(-b))
    end

    F′ = copy(F)
    F′[k] = divexact(M_pos + M_neg, F[k])
    return F′
end

function _mutate_C(C::Matrix{T}, B::Matrix{T}, k::Int) where {T <: Integer}
    n  = size(C, 1)
    εC = _epsilon(C[:, k])
    C′ = copy(C)
    C′[:, k] = -C[:, k]
    for j in 1:n
        j == k && continue
        bump = max(εC * B[k, j], 0)
        iszero(bump) && continue
        @. C′[:, j] = C[:, j] + bump * C[:, k]
    end
    return C′
end

# ─── G-matrix: on-demand computation from cluster + F-polynomials ─────────────
#
# g-vectors are defined via the separation formula:
#   x_k = x^{g_k} * F_k(ŷ₁,…,ŷₙ)
# where ŷ_j = ∏_i xᵢ^{B₀[i,j]} (column j of the INITIAL exchange matrix).
# g-vectors are NOT required to be sign-coherent, so they cannot be tracked
# incrementally by an εG-based formula.  We derive them exactly by solving
# for x^{g_k} = x_k / F_k(ŷ)|_{y=1}, which is always a Laurent monomial.

function _gmatrix_and_B0(s::Seed{PrincipalCoefficients})
    # With frozen vertices, _mutate_cluster includes frozen variables in the
    # exchange relations, so x_k / F_k(ŷ) built from the mutable block alone is
    # NOT a Laurent monomial - silently taking its first exponent vector would
    # produce garbage g-vectors.  Error loudly until a coefficient-aware
    # separation formula is implemented.  (C-matrix and y-variables remain
    # well-defined and available for such seeds.)
    s.quiver.n_frozen == 0 || throw(InvalidArgument(
        "g-vectors and the separation formula are not implemented for seeds " *
        "with frozen vertices (n_frozen = $(s.quiver.n_frozen)); " *
        "c-vectors and y-variables remain available"))
    n     = s.quiver.n_mutable
    n_tot = n + s.quiver.n_frozen
    Fps   = s.coeffs.fpolys
    B0    = s.coeffs.B0

    Fx    = s.ring
    xvars = Fx.(gens(base_ring(Fx)))
    yhat  = [prod(xvars[i]^B0[i, j] for i in 1:n; init = one(Fx)) for j in 1:n]

    G = zeros(Int, n, n)
    for k in 1:n
        Fk_x = evaluate(Fps[k], yhat)
        gmon = s.cluster[k] * inv(Fk_x)
        num  = numerator(gmon)
        den  = denominator(gmon)
        ev_n = isone(num) ? zeros(Int, n_tot) : collect(exponent_vectors(num))[1]
        ev_d = isone(den) ? zeros(Int, n_tot) : collect(exponent_vectors(den))[1]
        G[:, k] = (ev_n - ev_d)[1:n]
    end
    return G, B0
end

"""
    _mutate_yvars(yvars, B, k) → new yvars

Mutate rational y-variables in the universal semifield Frac(ZZ[y₁,…,yₙ]).

Convention (Fomin–Zelevinsky IV, Prop. 3.9): using b_{kj} = B[k, j] (the
(k,j) entry of the CURRENT exchange matrix):

    y′_k = y_k⁻¹
    y′_j = y_j · y_k^[b_{kj}]₊ · (1 + y_k)^{-b_{kj}}    for j ≠ k

where [x]₊ = max(x, 0).  The tropicalization of this convention (componentwise
min of exponents in numerator minus min in denominator) equals the c-vector
maintained by `_mutate_C`, and the resulting C-matrix satisfies tropical
duality C = (Gᵀ)⁻¹ in the skew-symmetric case - both covered by tests.
"""
function _mutate_yvars(yvars::Vector{S}, B::Matrix{Int}, k::Int) where {S}
    n    = length(yvars)
    yk   = yvars[k]
    ynew = copy(yvars)

    ynew[k] = inv(yk)

    for j in 1:n
        j == k && continue
        bkj = B[k, j]          # entry (k, j) of current exchange matrix
        # y′_j = y_j * y_k^[b_{kj}]_+ * (1 + y_k)^{-b_{kj}}
        factor = yk^max(bkj, 0)
        if bkj > 0
            factor = factor // (1 + yk)^bkj
        elseif bkj < 0
            factor = factor * (1 + yk)^(-bkj)
        end
        ynew[j] = yvars[j] * factor
    end
    return ynew
end

# ─── Mutation for PrincipalCoefficients ───────────────────────────────────────

"""
    mutate(s::Seed{PrincipalCoefficients}, k::Int) → Seed{PrincipalCoefficients}

Mutate at mutable vertex `k`, updating the cluster, C-matrix (c-vectors),
G-matrix (g-vectors), and rational y-variables simultaneously.
"""
function mutate(s::Seed{PrincipalCoefficients}, k::Int)
    q_new, cluster_new, path_new = _mutate_cluster(s, k)
    pd     = s.coeffs
    C′     = _mutate_C(pd.C, s.quiver.B, k)
    ynew   = _mutate_yvars(pd.yvars, s.quiver.B, k)
    ygens  = gens(parent(pd.fpolys[1]))
    fpolys′ = _mutate_fpolys(pd.fpolys, pd.B_ext, ygens, k)
    B_ext′  = _mutate_matrix(pd.B_ext, k)
    pd_new = PrincipalData(C′, pd.yring, ynew, pd.B0, B_ext′, fpolys′)
    T  = eltype(s.cluster)
    F  = typeof(s.ring)
    S  = eltype(ynew)
    P  = eltype(fpolys′)
    return Seed{PrincipalCoefficients, T, F, PrincipalData{S, P}}(
        q_new, cluster_new, s.ring, path_new, pd_new)
end

# ─── Accessors ────────────────────────────────────────────────────────────────

"""Return the C-matrix (columns are c-vectors)."""
cmatrix(s::Seed{PrincipalCoefficients})  = s.coeffs.C

"""Return the G-matrix (columns are g-vectors), computed on demand."""
gmatrix(s::Seed{PrincipalCoefficients})  = _gmatrix_and_B0(s)[1]

"""Return the c-vectors as a `Vector` of integer column vectors."""
cvectors(s::Seed{PrincipalCoefficients}) = [s.coeffs.C[:, k] for k in 1:size(s.coeffs.C, 2)]

"""Return the g-vectors as a `Vector` of integer column vectors."""
gvectors(s::Seed{PrincipalCoefficients}) = [gmatrix(s)[:, k] for k in 1:s.quiver.n_mutable]

"""Return the k-th c-vector (column k of the C-matrix)."""
c_vector(s::Seed{PrincipalCoefficients}, k::Int) = s.coeffs.C[:, k]

"""Return the k-th g-vector (column k of the G-matrix), computed on demand."""
g_vector(s::Seed{PrincipalCoefficients}, k::Int) = _gmatrix_and_B0(s)[1][:, k]

"""
    is_sign_coherent(s::Seed{PrincipalCoefficients}) → Bool

Return `true` if every c-vector is either entirely ≥ 0 or entirely ≤ 0.
"""
function is_sign_coherent(s::Seed{PrincipalCoefficients})
    n = size(s.coeffs.C, 2)
    all(1:n) do k
        c = s.coeffs.C[:, k]
        all(>=(0), c) || all(<=(0), c)
    end
end

"""
    y_variables(s::Seed{PrincipalCoefficients}; semifield=:rational)

Return the y-variables of `s`.
- `semifield=:rational` (default): the rational y-variables in `Frac(ZZ[y₁,…,yₙ])`.
- `semifield=:tropical`: the c-vectors (columns of the C-matrix) as `Vector{Vector{Int}}`.
"""
function y_variables(s::Seed{PrincipalCoefficients}; semifield::Symbol=:rational)
    if semifield === :rational
        return s.coeffs.yvars
    elseif semifield === :tropical
        return cvectors(s)
    else
        throw(InvalidArgument("unknown semifield=$semifield; use :rational or :tropical"))
    end
end

# ─── F-polynomials ────────────────────────────────────────────────────────────

"""
    fpolynomials(s::Seed{PrincipalCoefficients}) → Vector{<:MPolyRingElem}

Return the F-polynomials of the cluster variables in `s`, as elements of
`ZZ[y₁,…,yₙ]` where `n = n_mutable`.

The k-th entry is the F-polynomial of the k-th cluster variable. Initial cluster
variables have F-polynomial 1.

The F-polynomials are maintained incrementally by `mutate` (one recurrence
step per mutation, see `_mutate_fpolys`), so this accessor is O(1).

See `fpolynomials(::Seed{TrivialCoefficients})` for algorithm documentation.
"""
fpolynomials(s::Seed{PrincipalCoefficients}) = s.coeffs.fpolys

"""
    f_polynomial(s::Seed{PrincipalCoefficients}, k::Int) → MPolyRingElem

Return the k-th F-polynomial of `s`.
"""
f_polynomial(s::Seed{PrincipalCoefficients}, k::Int) = fpolynomials(s)[k]

# ─── Separation formula ───────────────────────────────────────────────────────

"""
    separation_formula(s::Seed{PrincipalCoefficients}, k::Int) → FracElem

Return the k-th cluster variable of `s` expressed via the separation formula
in the ring `Frac(ZZ[x₁,…,xₙ, y₁,…,yₙ])`:

    x_k = (∏ᵢ xᵢ^{gᵢ}) · F_k(ŷ₁,…,ŷₙ)

where ŷⱼ = yⱼ · ∏ᵢ xᵢ^{B₀[i,j]} (column j of the INITIAL exchange matrix B₀),
gᵢ is the FZ g-vector of x_k, equal to column k of the G-matrix G[:,k], and
F_k|_trop = 1 (for principal coefficients).
"""
function separation_formula(s::Seed{PrincipalCoefficients}, k::Int)
    n      = s.quiver.n_mutable
    n_tot  = n + s.quiver.n_frozen
    G, B0  = _gmatrix_and_B0(s)
    gvec   = G[:, k]
    Fk     = fpolynomials(s)[k]

    # Build an xy-ring for the full expression.
    x_names = ["x_$i" for i in 1:n_tot]
    y_names = ["y$i" for i in 1:n]
    all_names = vcat(x_names, y_names)
    Rxy, xy_gens = polynomial_ring(ZZ, all_names)
    Fxy = fraction_field(Rxy)

    xvars = Fxy.(xy_gens[1:n_tot])
    yvars = Fxy.(xy_gens[n_tot+1:n_tot+n])

    # ŷ_j = y_j · ∏ᵢ xᵢ^{B₀[i,j]}   (column j of B₀ gives the x-exponents)
    # In code: for each j, sum over rows i.
    yhat = Vector{eltype(xvars)}(undef, n)
    for j in 1:n
        yhat[j] = yvars[j]
        for i in 1:n
            b = B0[i, j]
            if b > 0
                yhat[j] *= xvars[i]^b
            elseif b < 0
                yhat[j] *= inv(xvars[i]^(-b))
            end
        end
    end

    # Evaluate F_k at ŷ
    Fk_xy = evaluate(Fk, yhat)

    # ∏ xᵢ^{gᵢ} (g-vector, may have negative entries → fractions)
    xpow = one(Fxy)
    for i in 1:n
        gi = gvec[i]
        if gi > 0
            xpow *= xvars[i]^gi
        elseif gi < 0
            xpow *= inv(xvars[i]^(-gi))
        end
    end

    return xpow * Fk_xy  # F_k|_trop = 1 for principal coefficients
end

"""
    separation_formula_trivial(s::Seed{PrincipalCoefficients}, k::Int) → FracElem

Return the k-th cluster variable via the separation formula specialized to
trivial coefficients (y → 1, i.e. ŷⱼ = ∏ᵢ xᵢ^{B₀[i,j]}), as an element of
`Frac(ZZ[x₁,…,xₙ])`.  This ring is the SAME as `Seed(q).ring`, so the result
can be `==`-compared (as strings) to `mutate(Seed(q), path)[k]`.

The FZ g-vector equals column k of the G-matrix (see `separation_formula`).
"""
function separation_formula_trivial(s::Seed{PrincipalCoefficients}, k::Int)
    n     = s.quiver.n_mutable
    n_tot = n + s.quiver.n_frozen
    G, B0 = _gmatrix_and_B0(s)
    gvec  = G[:, k]
    Fk    = fpolynomials(s)[k]

    # Build the x-ring (same variable names as the initial trivial seed).
    x_names = ["x_$i" for i in 1:n_tot]
    Rx, x_gens = polynomial_ring(ZZ, x_names)
    Fx = fraction_field(Rx)
    xvars = Fx.(x_gens)

    # ŷ_j = ∏ᵢ xᵢ^{B₀[i,j]}   (y → 1 specialization; column j of B₀)
    yhat = Vector{eltype(xvars)}(undef, n)
    for j in 1:n
        yhat[j] = one(Fx)
        for i in 1:n
            b = B0[i, j]
            if b > 0
                yhat[j] *= xvars[i]^b
            elseif b < 0
                yhat[j] *= inv(xvars[i]^(-b))
            end
        end
    end

    # Evaluate F_k at ŷ
    Fk_x = evaluate(Fk, yhat)

    # ∏ xᵢ^{gᵢ}
    xpow = one(Fx)
    for i in 1:n
        gi = gvec[i]
        if gi > 0
            xpow *= xvars[i]^gi
        elseif gi < 0
            xpow *= inv(xvars[i]^(-gi))
        end
    end

    return xpow * Fk_x
end

# ─── Display for Seed{PrincipalCoefficients} ──────────────────────────────────

# Compact one-line form (used inside collections, arrays, etc.)
function Base.show(io::IO, s::Seed{PrincipalCoefficients})
    n = length(s.cluster)
    print(io, "Seed($n cluster variables, $(s.quiver.n_mutable) mutable, principal coefficients)")
    isempty(s.mutation_path) || print(io, " after μ$(s.mutation_path)")
end

# Verbose form for the REPL
function Base.show(io::IO, ::MIME"text/plain", s::Seed{PrincipalCoefficients})
    show(io, MIME"text/plain"(), s.quiver)
    println(io, "Cluster variables:")
    for (i, x) in enumerate(s.cluster)
        tag = i > s.quiver.n_mutable ? "  [frozen]" : ""
        println(io, "  [$i]$tag: $x")
    end
    isempty(s.mutation_path) || println(io, "Mutation path: $(s.mutation_path)")
    n = size(s.coeffs.C, 2)
    println(io, "C-matrix (c-vectors as columns):")
    for i in 1:n
        print(io, " ")
        for j in 1:n; print(io, lpad(s.coeffs.C[i, j], 4)); end
        println(io)
    end
    if s.quiver.n_frozen == 0
        G = _gmatrix_and_B0(s)[1]
        println(io, "G-matrix (g-vectors as columns):")
        for i in 1:n
            print(io, " ")
            for j in 1:n; print(io, lpad(G[i, j], 4)); end
            println(io)
        end
    end
end

function Base.show(io::IO, ::MIME"text/latex", s::Seed{PrincipalCoefficients})
    n     = length(s.cluster)
    n_mut = s.quiver.n_mutable
    m     = size(s.coeffs.C, 1)
    B     = s.quiver.B

    println(io, "\\begin{aligned}")

    # cluster row vector
    print(io, "  \\mathbf{x} &= \\begin{pmatrix} ")
    for i in 1:n
        _show_latex(io, s.cluster[i])
        i > n_mut && print(io, "^{\\ast}")
        i < n     && print(io, " & ")
    end
    println(io, " \\end{pmatrix} \\\\[6pt]")

    # exchange matrix
    print(io, "  B &= \\begin{pmatrix} ")
    for i in 1:n
        print(io, join(string.(B[i, :]), " & "))
        i < n && print(io, " \\\\ ")
    end
    println(io, " \\end{pmatrix} \\\\[6pt]")

    # C-matrix
    print(io, "  C &= \\begin{pmatrix} ")
    for i in 1:m
        print(io, join(string.(s.coeffs.C[i, :]), " & "))
        i < m && print(io, " \\\\ ")
    end
    print(io, " \\end{pmatrix}")

    # G-matrix (only defined for seeds without frozen vertices)
    if s.quiver.n_frozen == 0
        G = _gmatrix_and_B0(s)[1]
        print(io, " \\\\[6pt]\n  G &= \\begin{pmatrix} ")
        for i in 1:m
            print(io, join(string.(G[i, :]), " & "))
            i < m && print(io, " \\\\ ")
        end
        print(io, " \\end{pmatrix}")
    end

    # optional mutation path
    if !isempty(s.mutation_path)
        path_str = join(string.(s.mutation_path), ",\\,")
        print(io, " \\\\[4pt]\n  \\text{mutations} &= ($(path_str))")
    end

    print(io, "\n\\end{aligned}")
end

# ─── Mutation-finiteness ───────────────────────────────────────────────────────

"""
    is_mutation_finite(q::Quiver; max_quivers::Int = 10_000) → Bool

Return `true` if `q` is mutation-finite (its mutation class is finite).

Finite-type quivers are detected first via `is_finite_type` (which searches for
an acyclic representative) and return `true` immediately without a full BFS - this
correctly handles all finite Dynkin types including E₆, E₇, E₈ whose mutation
classes exceed the default cutoff.

For non-finite-type quivers, falls back to BFS with an early-exit cutoff;
returns `false` when more than `max_quivers` distinct exchange matrices are
found before the BFS terminates.  The default cutoff of 10 000 covers the
largest known mutation-finite, non-finite-type classes.
"""
function is_mutation_finite(q::Quiver; max_quivers::Int = 10_000)
    is_finite_type(q) && return true
    !is_truncated(mutation_class(q; max_quivers))
end

# ─── ExtendedCoefficients (geometric type) ────────────────────────────────────
#
# A cluster algebra of *geometric type*: the seed's quiver carries m ≥ 1 frozen
# vertices whose variables generate the coefficient (tropical) semifield.  The
# extended exchange matrix is the quiver's full (n+m)×(n+m) matrix; its mutable
# block `1:n` drives cluster mutation while its frozen block `n+1:n+m` tracks the
# coefficient tropical y-variables (c-vectors) automatically under matrix
# mutation.  Nothing extra is stored - every coefficient is recomputed from the
# quiver and cluster, so no incremental state can drift.

"""
    extend_geometric(s::Seed{TrivialCoefficients}) → Seed{ExtendedCoefficients}

Reinterpret a seed whose quiver has frozen vertices as a cluster algebra of
**geometric type**: the frozen variables become the coefficients, exposing the
geometric y-variables ([`y_variables`](@ref)`(s; semifield=:geometric)`), their
tropical shadow (the coefficient c-vectors), and the full ŷ cross-ratios
([`y_hat`](@ref)).  Requires `s.quiver.n_frozen ≥ 1`.  The quiver, cluster, ring,
and mutation path are shared unchanged.
"""
function extend_geometric(s::Seed{TrivialCoefficients})
    s.quiver.n_frozen >= 1 || throw(InvalidArgument(
        "extend_geometric requires a quiver with frozen (coefficient) vertices; " *
        "n_frozen = 0. Use extend for principal coefficients."))
    T = eltype(s.cluster)
    F = typeof(s.ring)
    return Seed{ExtendedCoefficients, T, F, Nothing}(
        s.quiver, s.cluster, s.ring, s.mutation_path, nothing)
end

"""
    mutate(s::Seed{ExtendedCoefficients}, k::Int) → Seed{ExtendedCoefficients}

Mutate at mutable vertex `k`.  The cluster and the extended exchange matrix
(quiver) mutate together; the coefficients follow from the mutated frozen block.
Mutating a frozen (coefficient) vertex throws `FrozenVertexMutation`.
"""
function mutate(s::Seed{ExtendedCoefficients}, k::Int)
    q_new, cluster_new, path_new = _mutate_cluster(s, k)
    T = eltype(cluster_new)
    F = typeof(s.ring)
    return Seed{ExtendedCoefficients, T, F, Nothing}(
        q_new, cluster_new, s.ring, path_new, nothing)
end

# ∏_i base[i]^{e[i]} in the fraction field, negative exponents via division.
function _laurent_monomial(base::Vector{T}, e::AbstractVector{<:Integer}) where {T}
    R = parent(base[1])
    m = one(R)
    for i in eachindex(e)
        b = e[i]
        b > 0 && (m *= base[i]^b)
        b < 0 && (m //= base[i]^(-b))
    end
    return m
end

"""
    y_variables(s::Seed{ExtendedCoefficients}; semifield=:geometric)

The coefficient y-variables of a geometric-type seed.
- `:geometric` (default): `y_j = ∏_{frozen i} x_i^{B[i,j]}` - Laurent monomials in
  the frozen (coefficient) variables, one per mutable vertex `j`.
- `:tropical`: the coefficient c-vectors (columns of the frozen block of the
  extended exchange matrix), as `Vector{Vector{Int}}`.
"""
function y_variables(s::Seed{ExtendedCoefficients}; semifield::Symbol=:geometric)
    n = s.quiver.n_mutable
    m = s.quiver.n_frozen
    B = s.quiver.B
    if semifield === :tropical
        return [B[n+1:n+m, j] for j in 1:n]
    elseif semifield === :geometric
        frozen = s.cluster[n+1:n+m]
        return [_laurent_monomial(frozen, B[n+1:n+m, j]) for j in 1:n]
    else
        throw(InvalidArgument("unknown semifield=$semifield; use :geometric or :tropical"))
    end
end

"""
    y_hat(s::Seed{ExtendedCoefficients}) → Vector

The full geometric ŷ-variables `ŷ_j = ∏_i x_i^{B[i,j]}` over **all** vertices
(mutable and frozen), one per mutable vertex `j`.  For a triangulated surface
these are the Thurston shear coordinates (cross-ratios of lambda lengths).
"""
function y_hat(s::Seed{ExtendedCoefficients})
    n = s.quiver.n_mutable
    B = s.quiver.B
    return [_laurent_monomial(s.cluster, B[:, j]) for j in 1:n]
end

"""Return the coefficient C-matrix (frozen block of the extended exchange matrix)."""
cmatrix(s::Seed{ExtendedCoefficients}) =
    s.quiver.B[s.quiver.n_mutable+1:end, 1:s.quiver.n_mutable]

"""Return the coefficient c-vectors (columns of [`cmatrix`](@ref))."""
cvectors(s::Seed{ExtendedCoefficients}) =
    [cmatrix(s)[:, j] for j in 1:s.quiver.n_mutable]

function Base.show(io::IO, s::Seed{ExtendedCoefficients})
    n = length(s.cluster)
    print(io, "Seed($n cluster variables, $(s.quiver.n_mutable) mutable, " *
              "$(s.quiver.n_frozen) frozen, extended coefficients)")
    isempty(s.mutation_path) || print(io, " after μ$(s.mutation_path)")
end

function Base.show(io::IO, ::MIME"text/plain", s::Seed{ExtendedCoefficients})
    show(io, MIME"text/plain"(), s.quiver)
    println(io, "Cluster variables:")
    for (i, x) in enumerate(s.cluster)
        tag = i > s.quiver.n_mutable ? "  [frozen]" : ""
        println(io, "  [$i]$tag: $x")
    end
    isempty(s.mutation_path) || println(io, "Mutation path: $(s.mutation_path)")
    println(io, "Coefficient c-vectors (columns): ", cvectors(s))
end
