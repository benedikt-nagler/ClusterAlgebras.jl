abstract type AbstractSeed end

abstract type CoefficientKind end
struct TrivialCoefficients   <: CoefficientKind end
struct PrincipalCoefficients <: CoefficientKind end
struct ExtendedCoefficients  <: CoefficientKind end   # reserved for surfaces/physics; not implemented yet

"""
    Seed(q::Quiver)
    Seed(q::Quiver, var_names::Vector{String})

A cluster-algebra seed: a `Quiver` together with a cluster of variables living
in `Frac(ZZ[x₁,…,xₙ])` (one variable per vertex, frozen vertices included).
The type parameter `K <: CoefficientKind` distinguishes trivial coefficients
(`Seed(q)`) from principal coefficients (`extend(Seed(q))`, which additionally
tracks c-/g-vectors, y-variables, and F-polynomials).

Mutate with [`mutate`](@ref); `s[k]` returns the k-th cluster variable, and
`s.mutation_path` records the mutations applied since construction.
"""
struct Seed{K<:CoefficientKind, T<:RingElem, F<:Ring, C} <: AbstractSeed
    quiver::Quiver
    cluster::Vector{T}
    ring::F
    mutation_path::Vector{Int}
    coeffs::C            # nothing for trivial; a PrincipalData for principal
end

# Internal constructor used by mutate — bypasses ring rebuilding
function _seed(q::Quiver, cluster::Vector{T}, ring::F, path::Vector{Int}) where {T <: RingElem, F <: Ring}
    Seed{TrivialCoefficients, T, F, Nothing}(q, cluster, ring, path, nothing)
end

function Seed(q::Quiver, var_names::Vector{String})
    n_total = q.n_mutable + q.n_frozen
    length(var_names) == n_total ||
        error("var_names must have length $n_total, got $(length(var_names))")
    R, _ = polynomial_ring(ZZ, var_names)
    F = fraction_field(R)
    cluster = F.(gens(R))
    return Seed{TrivialCoefficients, eltype(cluster), typeof(F), Nothing}(q, cluster, F, Int[], nothing)
end

Seed(q::Quiver) = Seed(q, ["x_$i" for i in 1:q.n_mutable + q.n_frozen])

Base.getindex(s::Seed, k::Int)       = s.cluster[k]
Base.length(s::Seed)                  = length(s.cluster)
Base.iterate(s::Seed)                 = iterate(s.cluster)
Base.iterate(s::Seed, state)          = iterate(s.cluster, state)
Base.eltype(::Type{<:Seed{K,T}}) where {K,T} = T
