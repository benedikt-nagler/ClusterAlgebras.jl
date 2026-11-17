# Machine-integer storage is retained. Widen intermediates before abs, products
# and negation; an unrepresentable result raises rather than wrapping.
function _mutation_entry(::Type{T}, x) where {T<:Integer}
    try
        T(x)
    catch e
        e isa InexactError || rethrow()
        throw(OverflowError("exact mutation entry exceeds the matrix integer type $T"))
    end
end

function _mutate_matrix(B::Matrix{T}, k::Int) where {T <: Integer}
    n = size(B,1)
    out = copy(B)
    W = widen(T)
    for i in 1:n, j in 1:n
        if i == k || j == k
            out[i,j] = _mutation_entry(T,-W(B[i,j]))
        else
            a,b = W(B[i,k]),W(B[k,j])
            correction = a > 0 && b > 0 ? a*b : a < 0 && b < 0 ? -a*b : zero(W)
            out[i,j] = _mutation_entry(T,W(B[i,j])+correction)
        end
    end
    out
end

"""
    mutate(x, k::Int)
    mutate(x, ks::AbstractVector{Int})
    mutate(x, label::String)

Mutate a `Quiver` or `Seed` at mutable vertex `k` (Fomin–Zelevinsky matrix
mutation; for seeds also the exchange relation on cluster variables, and for
principal-coefficient seeds additionally the C-matrix, y-variables, and
F-polynomials). Unrepresentable integer entries throw `OverflowError`.  Returns a new object; the input is never modified.

A vector of indices applies the mutations left to right
(`mutate(x, [1, 2]) == mutate(mutate(x, 1), 2)`).  A string mutates at the
vertex with that label.  Mutating a frozen vertex throws
`FrozenVertexMutation`.
"""
function mutate(q::Quiver, k::Int)
    n_total = q.n_mutable + q.n_frozen
    (1 <= k <= n_total) || throw(InvalidVertex(k, n_total))
    k <= q.n_mutable    || throw(FrozenVertexMutation(k, q.n_mutable))
    return Quiver(_mutate_matrix(q.B, k), q.n_mutable, q.d, q.labels)
end

# The two monomials of the exchange relation at vertex k,
#
#   x_k x_k' = P_k = ∏_{B[i,k]>0} x_i^{B[i,k]}  +  ∏_{B[i,k]<0} x_i^{−B[i,k]},
#
# the product running over all vertices, frozen included (so for a seed with
# coefficients each monomial carries its tropical coefficient p_k^± as the
# frozen part).  Shared by mutation and by `bounds.jl`, which needs P_k
# without performing the mutation.
function _exchange_monomials(s::Seed, k::Int)
    B       = s.quiver.B
    n_total = s.quiver.n_mutable + s.quiver.n_frozen

    pos = one(s.cluster[k])
    neg = one(s.cluster[k])
    for i in 1:n_total
        i == k && continue
        b = B[i, k]
        if b > 0
            pos *= s.cluster[i]^b
        elseif b < 0
            neg *= s.cluster[i]^(-b)
        end
    end
    return pos, neg
end

# Helper: validate k and compute the new quiver, cluster, and path for any Seed.
function _mutate_cluster(s::Seed, k::Int)
    n_total = s.quiver.n_mutable + s.quiver.n_frozen
    (1 <= k <= n_total)      || throw(InvalidVertex(k, n_total))
    k <= s.quiver.n_mutable  || throw(FrozenVertexMutation(k, s.quiver.n_mutable))

    B    = s.quiver.B

    q_new = Quiver(_mutate_matrix(B, k), s.quiver.n_mutable, s.quiver.d, s.quiver.labels)

    pos, neg = _exchange_monomials(s, k)

    cluster_new = copy(s.cluster)
    cluster_new[k] = (pos + neg) / s.cluster[k]

    path_new = push!(copy(s.mutation_path), k)
    return q_new, cluster_new, path_new
end

function mutate(s::Seed{TrivialCoefficients}, k::Int)
    q_new, cluster_new, path_new = _mutate_cluster(s, k)
    return _seed(q_new, cluster_new, s.ring, path_new)
end

# Mutation sequence - works for any Seed kind (delegates to the single-Int method)
mutate(x::Union{Quiver, Seed}, ks::AbstractVector{Int}) = foldl(mutate, ks; init=x)

# Label-based mutation - works for any Seed kind
function mutate(q::Quiver, label::String)
    k = findfirst(==(label), q.labels)
    k === nothing && throw(InvalidArgument(
        "no vertex with label \"$label\"; valid labels: $(q.labels)"))
    return mutate(q, k)
end

function mutate(s::Seed, label::String)
    k = findfirst(==(label), s.quiver.labels)
    k === nothing && throw(InvalidArgument(
        "no vertex with label \"$label\"; valid labels: $(s.quiver.labels)"))
    return mutate(s, k)
end
