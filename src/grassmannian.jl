# Grassmannian cluster algebras - Track C: GR.
#
# Reference: J. Scott, "Grassmannians and Cluster Algebras," Proc. London Math. Soc. 92 (2006).
#
# Initial seed: the "rectangles seed".  Vertices sit on the k × (n−k) grid of
# rectangle Young diagrams inside the k × (n−k) box, plus one extra vertex for
# the empty diagram.  The a × b rectangle corresponds to the Plücker coordinate
#
#   I(a,b) = {1, …, k−a} ∪ {k−a+b+1, …, k+b},        I(∅) = {1, …, k}.
#
# Mutable vertices: 1 ≤ a ≤ k−1, 1 ≤ b ≤ n−k−1 (interior rectangles).
# Frozen vertices: full-height (a = k) and full-width (b = n−k) rectangles and
# the empty diagram - exactly the n cyclic-interval Plückers.
#
# Quiver arrows (each unit cell of the grid carries an oriented triangle):
#   (a,b) → (a,b+1)   (right)
#   (a,b) → (a+1,b)   (down)
#   (a+1,b+1) → (a,b) (diagonal)
#   I(∅) → (1,1)
# with arrows between two frozen vertices omitted.  Mutation at an interior
# rectangle then reproduces the three-term Plücker relation
#   p_{(a,b)} · p′ = p_{(a,b−1)} p_{(a,b+1)} + p_{(a−1,b)} p_{(a+1,b)}-type
# exchanges (verified against short Plücker relations for small (k,n) and by
# the numeric minor-substitution oracle in the tests).

# ─── Plücker label helpers ────────────────────────────────────────────────────

"""
    plucker_label(subset) → String

Format a sorted integer vector as a Plücker label, e.g. `[1,3,5]` → `"p_{135}"`.
Indices must be single-digit (1–9): the label format is not parseable otherwise.
"""
function plucker_label(subset::AbstractVector{Int})
    all(i -> 1 <= i <= 9, subset) || throw(InvalidArgument(
        "plucker_label requires single-digit indices (1–9), got $subset"))
    return "p_{" * join(sort(subset), "") * "}"
end

"""
    is_plucker_label(s) → Bool

Return `true` if `s` is a Plücker label of the form `"p_{...}"`.
"""
is_plucker_label(s::AbstractString) = startswith(s, "p_{") && endswith(s, "}")

"""
    plucker_subset(s) → Vector{Int}

Parse a Plücker label back to a sorted vector of single-digit indices.
"""
function plucker_subset(s::AbstractString)
    is_plucker_label(s) || throw(InvalidArgument(
        "\"$s\" is not a Plücker label; expected format \"p_{...}\""))
    inner = s[4:end-1]   # strip "p_{" and "}"
    return sort!([parse(Int, string(c)) for c in inner])
end

# ─── Internal helpers ─────────────────────────────────────────────────────────

# Plücker k-subset of the a × b rectangle Young diagram inside the k × (n−k)
# box.  a = 0 or b = 0 is the empty diagram {1,…,k}.
function _rectangle_subset(a::Int, b::Int, k::Int)
    (a == 0 || b == 0) && return collect(1:k)
    return vcat(collect(1:k-a), collect(k-a+b+1:k+b))
end

# Frozen boundary Plücker l: the cyclic k-subset {l, l+1, …, l+k-1} (mod n).
_boundary_subset(l::Int, k::Int, n::Int) = sort!([mod1(l + m, n) for m in 0:k-1])

# Vertex index of grid position (a, b).  Mutable interior rectangles come
# first in row-major order; the n frozen vertices follow, ordered by the
# starting point l of their cyclic interval:
#   empty diagram        = {1,…,k}              → l = 1
#   full height (a = k)  = {b+1,…,b+k}          → l = b + 1
#   full width (b = n−k) = {1,…,k−a}∪{n−a+1,…,n} → l = n − a + 1
function _gr_index(a::Int, b::Int, k::Int, n::Int)
    n_mut = (k - 1) * (n - k - 1)
    (a == 0 || b == 0) && return n_mut + 1
    a == k     && return n_mut + b + 1
    b == n - k && return n_mut + n - a + 1
    return (a - 1) * (n - k - 1) + b
end

_gr_frozen(a::Int, b::Int, k::Int, n::Int) = a == 0 || b == 0 || a == k || b == n - k

# ─── Quiver constructor ───────────────────────────────────────────────────────

function _grassmannian_quiver(k::Int, n::Int)
    n <= 9 || throw(InvalidArgument(
        "Grassmannian quivers are limited to n ≤ 9 (single-digit Plücker labels), got n=$n"))
    n_mut = (k - 1) * (n - k - 1)
    N     = n_mut + n   # mutable + frozen

    # Subsets by vertex index (see _gr_index for the ordering).
    subsets = Vector{Vector{Int}}(undef, N)
    for a in 1:k-1, b in 1:n-k-1
        subsets[_gr_index(a, b, k, n)] = _rectangle_subset(a, b, k)
    end
    for l in 1:n
        subsets[n_mut + l] = _boundary_subset(l, k, n)
    end

    B = zeros(Int, N, N)
    function add_arrow!(a1, b1, a2, b2)   # (a1,b1) → (a2,b2), skip frozen–frozen
        _gr_frozen(a1, b1, k, n) && _gr_frozen(a2, b2, k, n) && return
        i, j = _gr_index(a1, b1, k, n), _gr_index(a2, b2, k, n)
        B[i, j] += 1
        B[j, i] -= 1
    end

    for a in 1:k, b in 1:n-k
        b < n - k     && add_arrow!(a, b, a, b + 1)          # right
        a < k         && add_arrow!(a, b, a + 1, b)          # down
        a > 1 && b > 1 && add_arrow!(a, b, a - 1, b - 1)     # diagonal
    end
    add_arrow!(0, 0, 1, 1)                                   # empty → (1,1)

    labels = [plucker_label(subsets[i]) for i in 1:N]
    return Quiver(B, n_mut, ones(Int, n_mut), labels)
end

# ─── Public API ───────────────────────────────────────────────────────────────

"""
    grassmannian(k, n) → Seed

Return the initial seed of the Grassmannian cluster algebra for `Gr(k,n)` (Scott 2006).

The seed is Scott's *rectangles seed*: `(k-1)*(n-k-1)` mutable cluster
variables (Plücker coordinates of interior rectangle Young diagrams inside the
`k × (n-k)` box) and `n` frozen variables (cyclic-interval boundary Plückers).
Cluster variables print as `p_S` where `S` is the corresponding `k`-subset of
`{1,…,n}`.

The quiver type corresponds to classical Dynkin diagrams for small cases:

| `Gr(k,n)` | Type |
|---|---|
| `Gr(2,5)` | `A_2` |
| `Gr(2,6)` | `A_3` |
| `Gr(3,6)` | `D_4` |
| `Gr(4,7)` | `E_6` |

# Examples
```julia
julia> s = grassmannian(2, 5)
julia> cartan_type(s.quiver)    # (:A, 2)
julia> n_clusters(s.quiver)     # 5
```
"""
function grassmannian(k::Int, n::Int)
    2 <= k <= n - 2 || throw(InvalidArgument(
        "Gr(k,n) requires 2 ≤ k ≤ n-2, got k=$k, n=$n"))
    q = _grassmannian_quiver(k, n)
    # Ring variable names: brace-free "p135" to avoid any special-character issues
    # in AbstractAlgebra's polynomial ring; the quiver labels carry the "p_{135}" form.
    ring_names = ["p" * join(plucker_subset(lbl), "") for lbl in q.labels]
    return Seed(q, ring_names)
end

"""
    x_coordinates(s::Seed{PrincipalCoefficients}) → Vector

Return the cluster X-coordinates (rational y-variables) of `s`.  These are the
arguments of the symbol in scattering-amplitude computations.
Call `extend(s)` first to attach principal coefficients.

# Example
```julia
julia> xc = x_coordinates(extend(grassmannian(2, 5)))
```
"""
x_coordinates(s::Seed{PrincipalCoefficients}) = y_variables(s; semifield=:rational)
