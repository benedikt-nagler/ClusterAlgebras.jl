# The cluster complex and the g-vector fan (Track C, group AS - combinatorial half).
#
# The polytope *realizations* of these objects live downstream in
# ClusterPositiveGeometry.jl; nothing here needs a convex-geometry library.
#
# Convention ledger:
#   - Rays of a ClusterFan are g-vectors with respect to the initial seed handed
#     to `g_vector_fan`, i.e. the seed at the end of its `mutation_path`.  A
#     different initial seed gives the same fan in a different linear basis.
#     Pinned by: `test_cluster_fan.jl`, "rays are a basis-dependent labelling"
#     (the two initial seeds of A2 give matching cone combinatorics).
#   - Vertices of a ClusterComplex are labelled by d-vector, matching the
#     labelling `f_vector(::ExchangeGraph)` already uses, so the two index sets
#     are ordered independently of the fan's rays.

"""
    ClusterFan

The g-vector fan of a cluster algebra: one ray per cluster variable, one maximal
cone per cluster.

Fields `rays` (g-vectors), `cones` (sorted indices into `rays`, one entry per
cluster) and `truncated`. In finite type the fan is complete and is the normal
fan of the generalized associahedron, so [`is_complete`](@ref) holds; in infinite
type only the part reached within `max_seeds` seeds is present and `truncated` is
`true`. Build one with [`g_vector_fan`](@ref). Immutable.
"""
struct ClusterFan
    rays      :: Vector{Vector{Int}}
    cones     :: Vector{Vector{Int}}
    truncated :: Bool
end

"""
    ClusterComplex

The cluster complex: vertices are cluster variables, facets are clusters, and the
faces are the subsets of pairwise compatible cluster variables.

Fields `vertices` (one d-vector per cluster variable), `facets` (sorted indices
into `vertices`) and `truncated`. Build one with [`cluster_complex`](@ref); its
[`f_vector`](@ref) counts the faces by dimension. Immutable.
"""
struct ClusterComplex
    vertices  :: Vector{Vector{Int}}
    facets    :: Vector{Vector{Int}}
    truncated :: Bool
end

is_truncated(f::ClusterFan)     = f.truncated
is_truncated(cc::ClusterComplex) = cc.truncated

"""
    rays(f::ClusterFan) → Vector{Vector{Int}}

The g-vectors spanning the rays of `f`, one per cluster variable.
"""
rays(f::ClusterFan) = f.rays

"""
    maximal_cones(f::ClusterFan) → Vector{Vector{Int}}

The maximal cones of `f` as sorted index vectors into [`rays`](@ref), one per cluster.
"""
maximal_cones(f::ClusterFan) = f.cones

"""
    facets(cc::ClusterComplex) → Vector{Vector{Int}}

The facets of `cc` as sorted index vectors into `cc.vertices`, one per cluster.
"""
facets(cc::ClusterComplex) = cc.facets

"""
    vertex_labels(cc::ClusterComplex) → Vector{Vector{Int}}

The d-vector of each vertex of `cc`, which identifies its cluster variable.
"""
vertex_labels(cc::ClusterComplex) = cc.vertices

nvertices(cc::ClusterComplex) = length(cc.vertices)

# ─── Construction ─────────────────────────────────────────────────────────────

"""
    g_vector_fan(s::Seed; max_seeds::Int = 1000) → ClusterFan

The g-vector fan of `s`: the g-vectors of all cluster variables as rays, grouped
into one maximal cone per cluster.

A `TrivialCoefficients` seed is extended to principal coefficients first, since
g-vectors are defined there. Distinct cluster variables have distinct g-vectors
(Gross-Hacking-Keel-Kontsevich), so the rays are in bijection with the cluster
variables and the cones with the clusters. Rays are read off with respect to the
initial seed of `s`; see the ledger at the top of `src/cluster_fan.jl`.

For a finite-type seed this is the complete normal fan of the generalized
associahedron. Otherwise the exchange graph is explored only to `max_seeds` seeds
and the result is flagged by [`is_truncated`](@ref), which makes it a piece of the
fan and not the fan.

```jldoctest
julia> fan = g_vector_fan(Quiver(:A, 2))
ClusterFan: 5 rays, 5 maximal cones, rank 2

julia> length(rays(fan)) == n_cluster_variables(Quiver(:A, 2))
true

julia> is_complete(fan)
true
```
"""
function g_vector_fan(s::Seed{PrincipalCoefficients}; max_seeds::Int = 1000)
    eg    = exchange_graph(s; max_seeds)
    rs    = Vector{Int}[]
    index = Dict{Vector{Int}, Int}()
    cones = Vector{Int}[]

    for i in 1:length(eg)
        cone = Int[]
        for g in gvectors(eg[i])
            v = Int.(g)
            j = get(index, v, 0)
            if j == 0
                push!(rs, v)
                j = length(rs)
                index[v] = j
            end
            push!(cone, j)
        end
        push!(cones, sort!(cone))
    end

    return ClusterFan(rs, cones, is_truncated(eg))
end

g_vector_fan(s::Seed{TrivialCoefficients}; kwargs...) = g_vector_fan(extend(s); kwargs...)
g_vector_fan(q::Quiver; kwargs...)                    = g_vector_fan(Seed(q); kwargs...)

"""
    cluster_complex(s::Seed; max_seeds::Int = 1000) → ClusterComplex

The cluster complex of `s`, with cluster variables as vertices and clusters as facets.

Vertices are labelled by d-vector. For a finite-type seed the complex is the
boundary complex of the generalized associahedron; otherwise it is truncated at
`max_seeds` seeds, and [`f_vector`](@ref) then undercounts.

```jldoctest
julia> cc = cluster_complex(Quiver(:A, 3))
ClusterComplex: 9 vertices, 14 facets, rank 3

julia> f_vector(cc)
4-element Vector{Int64}:
  1
  9
 21
 14
```
"""
function cluster_complex(s::Seed; max_seeds::Int = 1000)
    eg = exchange_graph(s; max_seeds)
    return cluster_complex(eg)
end

cluster_complex(q::Quiver; kwargs...) = cluster_complex(Seed(q); kwargs...)

function cluster_complex(eg::ExchangeGraph)
    dvecs, clusters = _cluster_complex_data(eg)
    return ClusterComplex(dvecs, clusters, is_truncated(eg))
end

f_vector(cc::ClusterComplex) = _f_vector_from_facets(cc.facets)

"""
    is_complete(f::ClusterFan) → Bool

Whether every wall of `f` is shared by exactly two maximal cones, the
combinatorial criterion for a complete simplicial fan.

Cluster fans are simplicial, so for an untruncated fan this is completeness, and
it is equivalent to the exchange graph being regular of degree the rank. A
truncated fan is never complete. The test is combinatorial: it does not check
that the cones cover ``\\mathbb{R}^n`` geometrically.
"""
function is_complete(f::ClusterFan)
    f.truncated && return false
    isempty(f.cones) && return false
    n = length(f.cones[1])
    n == 0 && return false

    walls = Dict{Vector{Int}, Int}()
    for c in f.cones
        length(c) == n || return false
        for i in 1:n
            w = deleteat!(copy(c), i)
            walls[w] = get(walls, w, 0) + 1
        end
    end
    return all(==(2), values(walls))
end

# ─── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, f::ClusterFan)
    n      = isempty(f.rays) ? 0 : length(f.rays[1])
    suffix = f.truncated ? " (truncated)" : ""
    print(io, "ClusterFan: $(length(f.rays)) rays, $(length(f.cones)) maximal cones, rank $n$suffix")
end

function Base.show(io::IO, cc::ClusterComplex)
    n      = isempty(cc.facets) ? 0 : length(cc.facets[1])
    suffix = cc.truncated ? " (truncated)" : ""
    print(io, "ClusterComplex: $(length(cc.vertices)) vertices, $(length(cc.facets)) facets, rank $n$suffix")
end
