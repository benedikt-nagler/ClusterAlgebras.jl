"""
    ClusterAlgebras

Exact computation with cluster algebras: quivers, seeds and mutation, coefficients,
finite-type classification, exchange graphs, green sequences, friezes, and the cluster
structures on Grassmannians.

Cluster variables live in `Frac(ZZ[x₁,…,xₙ])`, so every result is exact. Quivers and seeds
are immutable: [`mutate`](@ref) returns a new object.

# Entry points

- [`Quiver`](@ref) - build a quiver from a `B`-matrix, a Dynkin type (`Quiver(:A, 3)`), a
  type string (`Quiver("D4")`) or an edge list; [`mutate`](@ref) mutates it.
- [`Seed`](@ref) - a quiver together with its cluster; `mutate(s, k)` or `mutate(s, [1,2,1])`.
- [`extend`](@ref) / [`extend_geometric`](@ref) - attach principal or geometric coefficients,
  giving [`cmatrix`](@ref), [`gmatrix`](@ref), [`f_polynomial`](@ref) and
  [`separation_formula`](@ref).
- [`is_finite_type`](@ref), [`cartan_type`](@ref), [`n_cluster_variables`](@ref) -
  classification and counting.
- [`mutation_class`](@ref), [`exchange_graph`](@ref) - the combinatorics of all seeds.
- [`maximal_green_sequences`](@ref), [`dt_transformation`](@ref) - green sequences and the
  Donaldson–Thomas transformation.
- [`grassmannian`](@ref), [`symbol_alphabet`](@ref) - Grassmannian seeds and their symbol
  alphabets.

Plotting (`plot_quiver`, `plot_exchange_graph`, `plot_frieze`, …) lives in package
extensions: load `Graphs`, `GraphMakie` and a Makie backend to enable it.

The manual is at <https://benedikt-nagler.github.io/ClusterAlgebras.jl/>.
"""
module ClusterAlgebras

using AbstractAlgebra
using PrecompileTools
using Random

include("errors.jl")
include("quiver.jl")
include("canonical_form.jl")
include("named_quivers.jl")
include("seed.jl")
include("mutation.jl")
include("show.jl")
include("root_system.jl")
include("folding.jl")
include("denominator_vector.jl")
include("frieze.jl")
include("mutation_class.jl")
include("coefficients.jl")
include("green_sequences.jl")
include("periodicity.jl")
include("dt_transformation.jl")
include("ks_dilog.jl")
include("enumerative.jl")
include("grassmannian.jl")
include("symbol_alphabet.jl")
include("random_quiver.jl")

export Quiver, Seed, mutate, to_dot
export ClusterAlgebraError, NotSkewSymmetrizable, FrozenVertexMutation, InvalidVertex, InvalidArgument
export nvertices, labels, is_frozen
export cartan_companion, RootSystem, almost_positive_roots
export is_finite_type, is_affine_type, cartan_type, cartan_types
export fold, is_admissible_folding
export denominator_vector
export Frieze, frieze
export MutationClass, mutation_class, ExchangeGraph, exchange_graph, is_truncated
export AbstractSeed, CoefficientKind, TrivialCoefficients, PrincipalCoefficients, ExtendedCoefficients
export extend, extend_geometric, y_hat, cmatrix, gmatrix, cvectors, gvectors
export c_vector, g_vector, f_polynomial, y_variables, separation_formula, separation_formula_trivial
export is_sign_coherent, is_mutation_finite
export is_green, is_red, is_all_red, maximal_green_sequences, mgs_search,
       green_sequence_signs, verify_mutation_sequence
export mutation_period, y_system, is_bipartite
export dt_transformation
export ordered_c_vectors, QuantumDilogWord, quantum_dilog_word, omega, ks_dilog_product
export fpolynomials
export n_cluster_variables, n_clusters, f_vector, h_vector
export grassmannian, plucker_label, is_plucker_label, plucker_subset, x_coordinates
export symbol_alphabet, cluster_adjacency_matrix, cluster_adjacent
export canonical_form, canonical_permutation, permute_vertices
# `is_isomorphic` is AbstractAlgebra's generic function, extended to Quiver in
# canonical_form.jl; re-exporting the same binding keeps `using ClusterAlgebras,
# AbstractAlgebra` free of an ambiguity.
export is_isomorphic
export random_quiver, random_mutate

function plot_quiver end
function plot_quiver! end
function plot_exchange_graph end
function plot_frieze end
function plot_green_sequence end
function plot_g_vector_fan end

plot_quiver(::Any; kwargs...) =
    throw(ArgumentError("plot_quiver requires GraphMakie: run `using GraphMakie` (and a Makie backend)."))
plot_exchange_graph(::Any; kwargs...) =
    throw(ArgumentError("plot_exchange_graph requires GraphMakie: run `using Graphs, GraphMakie` (and a Makie backend)."))
plot_frieze(::Any; kwargs...) =
    throw(ArgumentError("plot_frieze requires Makie: run `using CairoMakie` (or GLMakie)."))
plot_green_sequence(::Any, ::Any; kwargs...) =
    throw(ArgumentError("plot_green_sequence requires Makie: run `using CairoMakie` (or GLMakie)."))
plot_g_vector_fan(::Any; kwargs...) =
    throw(ArgumentError("plot_g_vector_fan requires Makie: run `using CairoMakie` (or GLMakie)."))

export plot_quiver, plot_quiver!, plot_exchange_graph, plot_frieze, plot_green_sequence
export plot_g_vector_fan

@setup_workload begin
    B = [0 1; -1 0]
    @compile_workload begin
        q = Quiver(B)
        s = Seed(q)
        mutate(q, 1)
        mutate(s, 1)
        mutate(s, [1, 2])
    end
end

end
