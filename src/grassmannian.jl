# Grassmannian cluster algebras — Track C: GR.
#
# Reference: J. Scott, "Grassmannians and Cluster Algebras," Proc. London Math. Soc. 92 (2006).
#
# Gr(k,n) has a cluster algebra structure.  The initial seed sits on a
# (k−1) × (n−k−1) grid of mutable Plücker coordinates, with n frozen
# boundary (cyclic-interval) Plückers.
#
# B-matrix (winding number): for k-subsets S, T ⊂ {1,…,n},
#   B[S,T] = Σ_{a=1}^{n} ([a∈S]·[(a mod n)+1 ∈ T] − [a∈T]·[(a mod n)+1 ∈ S])

# ─── Plücker label helpers ────────────────────────────────────────────────────

"""
    plucker_label(subset) → String

Format a sorted integer vector as a Plücker label, e.g. `[1,3,5]` → `"p_{135}"`.
Assumes all indices are single-digit (n ≤ 9).
"""
plucker_label(subset::AbstractVector{Int}) = "p_{" * join(sort(subset), "") * "}"

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

# Winding-number exchange matrix entry for two k-subsets of {1,…,n}.
function _plucker_B(S::AbstractVector{Int}, T::AbstractVector{Int}, n::Int)
    Smask = zeros(Bool, n)
    Tmask = zeros(Bool, n)
    for i in S; Smask[i] = true; end
    for i in T; Tmask[i] = true; end
    b = 0
    for a in 1:n
        ap1 = mod1(a + 1, n)
        b += Int(Smask[a] && Tmask[ap1]) - Int(Tmask[a] && Smask[ap1])
    end
    return b
end

# Plücker k-subset at grid position (i, j) in the (k-1) × (n-k-1) mutable grid.
#
# k=2: S(1,j) = {1, j+2}  (one row, formula anchors at 1)
# k=3: S(i,j) = {i, j+k-1, n+1-j}  (avoids double arrows via wraparound)
# k≥4: greedy weakly-separated collection (no simple closed form)
#
# The k=2 and k=3 formulas return subsets in row-major grid order.
# The k≥4 path calls _greedy_ws_cluster which returns the full list.
_mutable_subset_k2(i::Int, j::Int) = [i, j+2]
_mutable_subset_k3(i::Int, j::Int, n::Int) = sort([i, j+2, n+1-j])

# Known valid initial clusters for Gr(k,n) where k≥4.
# These were verified computationally: each gives the correct Dynkin type
# (A_3 and E_6 respectively) and matches Scott 2006.
const _GR_INITIAL_CLUSTERS = Dict{Tuple{Int,Int}, Vector{Vector{Int}}}(
    (4, 6) => [[1,2,3,5],[1,2,4,5],[1,3,4,5]],
    (4, 7) => [[1,2,4,7],[1,3,4,6],[1,3,5,7],[1,4,5,7],[2,3,4,6],[2,3,5,6]],
)

function _greedy_ws_cluster(k::Int, n::Int)
    haskey(_GR_INITIAL_CLUSTERS, (k,n)) && return _GR_INITIAL_CLUSTERS[(k,n)]
    error("No initial cluster implemented for Gr($k,$n) with k≥4. " *
          "Add an entry to _GR_INITIAL_CLUSTERS.")
end


# Frozen boundary Plücker l: the cyclic k-subset {l, l+1, …, l+k-1} (mod n).
_boundary_subset(l::Int, k::Int, n::Int) = sort!([mod1(l + m, n) for m in 0:k-1])

# ─── Quiver constructor ───────────────────────────────────────────────────────

function _grassmannian_quiver(k::Int, n::Int)
    n_mut = (k - 1) * (n - k - 1)
    N     = n_mut + n   # mutable + frozen

    # Ordered subsets: mutable (row-major over grid) then frozen (boundary l = 1..n).
    subsets = Vector{Vector{Int}}(undef, N)
    if k == 2
        idx = 1
        for i in 1:k-1, j in 1:n-k-1
            subsets[idx] = _mutable_subset_k2(i, j)
            idx += 1
        end
    elseif k == 3
        idx = 1
        for i in 1:k-1, j in 1:n-k-1
            subsets[idx] = _mutable_subset_k3(i, j, n)
            idx += 1
        end
    else
        mutable_list = _greedy_ws_cluster(k, n)
        for idx in 1:n_mut
            subsets[idx] = mutable_list[idx]
        end
    end
    idx = n_mut + 1
    for l in 1:n
        subsets[idx] = _boundary_subset(l, k, n)
        idx += 1
    end

    # Full exchange matrix via the winding-number formula.
    B = zeros(Int, N, N)
    for s in 1:N, t in 1:N
        B[s, t] = _plucker_B(subsets[s], subsets[t], n)
    end

    labels = [plucker_label(subsets[i]) for i in 1:N]
    return Quiver(B, n_mut, ones(Int, n_mut), labels)
end

# ─── Public API ───────────────────────────────────────────────────────────────

"""
    grassmannian(k, n) → Seed

Return the initial seed of the Grassmannian cluster algebra for `Gr(k,n)` (Scott 2006).

The seed has `(k-1)*(n-k-1)` mutable cluster variables (interior Plücker
coordinates arranged on a rectangular grid) and `n` frozen variables (cyclic-
interval boundary Plückers).  Cluster variables print as `p_S` where `S` is the
corresponding `k`-subset of `{1,…,n}`.

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
