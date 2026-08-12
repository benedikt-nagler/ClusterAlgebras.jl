# ClusterAlgebras.jl

[![Version](https://juliahub.com/docs/General/ClusterAlgebras/stable/version.svg)](https://juliahub.com/ui/Packages/General/ClusterAlgebras)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/ClusterAlgebras.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/ClusterAlgebras.jl/dev)
[![CI](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact computation with cluster algebras in Julia.

A seed is a cluster $(x_1, \ldots, x_n)$ of rational functions together with a
skew-symmetrizable integer matrix $B$, equivalently a quiver. Mutation at index $k$ replaces
$x_k$ by

$$x_k x_k' = \prod_{b_{ik}>0} x_i^{b_{ik}} + \prod_{b_{ik}<0} x_i^{-b_{ik}},$$

and updates $B$. Cluster variables are rational functions over `Frac(ZZ[x₁,…,xₙ])`, computed
exactly with [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl). Seeds are
immutable and `mutate` returns a new seed.

## Installation

```julia
pkg> add ClusterAlgebras
```

Requires Julia 1.10 or later.

## Quick start

```julia
using ClusterAlgebras

q = Quiver(:A, 3)              # named Dynkin type; also Quiver("D4") or Quiver(B)
s = Seed(q)                    # cluster variables as rational functions
s2 = mutate(s, [1, 2, 1])      # mutate along a sequence

is_finite_type(q)              # true
cartan_type(q)                 # (:A, 3)
mutation_type(q)               # A3 - names the whole mutation class, finite or not
n_cluster_variables(q)         # 9
length(mutation_class(q))      # 14

mutation_type(Quiver([0 2; -2 0]))   # A(1,1)^(1), the Kronecker quiver
is_surface_type(Quiver(:A, 3))       # true - it triangulates a hexagon
```

Principal coefficients carry the c-/g-vector and F-polynomial data:

```julia
sp = mutate(extend(s), [1, 2, 1])    # Seed{PrincipalCoefficients}

cmatrix(sp)                    # c-vectors as columns
gmatrix(sp)
f_polynomial(sp, 1)
separation_formula(sp, 1)      # cluster variable from g-vector + F-polynomial
is_sign_coherent(sp)           # true
```

## Features

- **Quivers and seeds.** Skew-symmetrizable exchange matrices with frozen vertices, from a
  matrix, an edge list or a Dynkin name. Mutation at a vertex, a label or along a sequence.
  `verify_mutation_sequence`, `to_dot`.
- **Coefficients.** Trivial, principal and geometric systems as a type parameter on `Seed`.
  `extend`, `extend_geometric`, `cmatrix`, `gmatrix`, `f_polynomial`, `y_hat`,
  `separation_formula`, denominator vectors.
- **Classification.** `is_finite_type`, `is_affine_type`, `cartan_type`, `is_mutation_finite`.
  Each searches for an acyclic representative, so the answers hold for non-acyclic seeds.
  Cartan companion, root systems, almost-positive roots, `n_clusters`, `f_vector`, `h_vector`.
- **Mutation type.** `mutation_type` names a mutation-finite quiver by SageMath's
  `(letter, rank, twist)` triple: finite, simply-laced affine, rank 2, and the exceptional
  $X_6$, $X_7$, $E_{6,7,8}^{(1,1)}$. `mutation_types` handles the disconnected case.
- **Surface types.** `block_decomposition` decomposes a skew-symmetric quiver into the seven
  Felikson–Shapiro–Tumarkin blocks, with `reassemble` as the inverse and `is_surface_type` /
  `is_block_decomposable` as predicates.
- **Folding.** `fold(q, σ)` folds by an admissible vertex automorphism, gated on
  `is_admissible_folding`. This reaches the non-simply-laced types.
- **Mutation classes and exchange graphs.** BFS with deduplication and a truncation guard for
  mutation-infinite types.
- **Green sequences and DT theory.** Green/red predicates, `maximal_green_sequences`,
  `mgs_search`, `dt_transformation`, and Kontsevich–Soibelman dilogarithm words with their
  DT invariants $\Omega(\gamma)$.
- **Y-systems.** `mutation_period` for a sequence, `y_system` for the bipartite Zamolodchikov
  dynamics, of period $h+2$ in finite type.
- **Grassmannians.** `grassmannian(k, n)` builds the Plücker seed of $\mathrm{Gr}(k,n)$
  (Scott 2006). `symbol_alphabet` enumerates the letters, with `cluster_adjacency_matrix` for
  which may sit adjacent.
- **Bounds.** The Berenstein–Fomin–Zelevinsky bounds
  $\mathcal{L}(\Sigma) \subseteq \mathcal{A} \subseteq \mathcal{U}(\Sigma)$ via
  `lower_bound_generators`, `standard_monomials`, `lower_bound_expansion`, `is_laurent` and
  `in_upper_bound`. Upper-bound membership needs $n+1$ Laurentness checks and no Gröbner
  basis. `bound_certificate` reports which hypotheses of $\mathcal{A} = \mathcal{U}$ hold.
- **Greedy basis.** `greedy_element` and `greedy_coefficients` for the Lee–Li–Zelevinsky
  basis in rank 2.
- **Friezes.** SL₂ frieze patterns from triangulated polygons.
- **Isomorphism and sampling.** `canonical_form`, `is_isomorphic`, `random_quiver`,
  `random_mutate`. `mutation_class` itself does not deduplicate up to relabeling.
- **Interoperability.** `to_dig6` / `from_dig6` for SageMath's canonical digraph6 pair, and
  `to_qmu` / `from_qmu` / `write_qmu` / `read_qmu` for Keller's quiver-mutation applet.

## Plotting

Plotting loads on demand through package extensions. The core depends only on
AbstractAlgebra.

```julia
using Graphs, GraphMakie, GLMakie

plot_quiver(q)
plot_exchange_graph(q)
plot_green_sequence(q, seq)
plot_g_vector_fan(q)
plot_frieze(f)

mutation_explorer(s)                 # click a vertex to mutate
```

`mutation_explorer` needs an interactive backend (GLMakie or WGLMakie). It supports a history
slider, hover readout of each vertex's cluster variable and c-/g-vectors, a text box for
applying a sequence such as `1,3,2`, and `on_mutate` to capture the explored state.

Loading `Graphs` alone gives the exchange graph and quiver as `Graphs.jl` objects.

## Documentation

The [documentation](https://benedikt-nagler.github.io/ClusterAlgebras.jl/stable) has a manual
page per layer and a worked $A_2$ tutorial. `examples/` holds four notebooks.

## Related packages

Part of a family of Julia packages for exact and asymptotic methods:

- [ClusterSurfaces.jl](https://github.com/benedikt-nagler/ClusterSurfaces.jl): cluster
  algebras from marked surfaces, where an arc flip is a mutation of the seed built here.
- [QuantumClusterAlgebras.jl](https://github.com/benedikt-nagler/QuantumClusterAlgebras.jl):
  the $q$-deformation.
- [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl): divergent series,
  Borel–Padé summation, transseries. Independent of this package.
- [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl): Stokes graphs of a
  Schrödinger-type ODE, linked to this package by the Iwaki–Nakanishi dictionary.

## License

MIT
