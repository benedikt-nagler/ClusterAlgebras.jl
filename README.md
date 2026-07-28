# ClusterAlgebras.jl

[![Version](https://juliahub.com/docs/General/ClusterAlgebras/stable/version.svg)](https://juliahub.com/ui/Packages/General/ClusterAlgebras)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/ClusterAlgebras.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/ClusterAlgebras.jl/dev)
[![CI](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact computation with cluster algebras in Julia: quivers and seed mutation, coefficient
systems and c-/g-vectors, finite-type classification, mutation classes and exchange graphs,
green sequences and DT invariants, friezes, and the cluster structures on Grassmannians.

A cluster algebra is generated combinatorially. From a *seed* - a cluster
$(x_1, \ldots, x_n)$ of rational functions together with a skew-symmetrizable integer matrix
$B$, equivalently a quiver - *mutation* at index $k$ replaces $x_k$ by

$$x_k x_k' = \prod_{b_{ik}>0} x_i^{b_{ik}} + \prod_{b_{ik}<0} x_i^{-b_{ik}},$$

and updates $B$ accordingly. Iterating this generates the algebra, and every cluster variable
stays a Laurent polynomial in any initial cluster.

Everything is exact: cluster variables are rational functions over `Frac(ZZ[x₁,…,xₙ])` via
[AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl), seeds are immutable, and
`mutate` returns a new seed. Each layer - quivers, coefficients, classification, green
sequences - is usable on its own.

## Installation

From the Julia REPL:

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
n_cluster_variables(q)         # 9
length(mutation_class(q))      # 14
```

Principal coefficients carry the c-/g-vector and F-polynomial data:

```julia
sp = extend(s)                 # Seed{PrincipalCoefficients}
sp = mutate(sp, [1, 2, 1])

cmatrix(sp)                    # c-vectors as columns
gmatrix(sp)
f_polynomial(sp, 1)
separation_formula(sp, 1)      # cluster variable from g-vector + F-polynomial
is_sign_coherent(sp)           # true (Derksen–Weyman–Zelevinsky)
```

## Functionality

**Quivers and seeds.** Skew-symmetrizable exchange matrices with frozen vertices, from a
matrix, an edge list, or a Dynkin name. Mutation at a vertex, a label, or along a sequence,
plus `verify_mutation_sequence` and `to_dot`.

**Coefficients.** Trivial, principal and general extended systems as a type parameter on
`Seed`, so one `mutate` covers all three. C- and G-matrices, F-polynomials, rational and
tropical y-variables, `y_hat`, the separation formula, denominator vectors.

**Classification.** `is_finite_type` / `is_affine_type` / `cartan_type` /
`is_mutation_finite`, each searching for an acyclic representative first, so the answers hold
for non-acyclic seeds too. Cartan companion, root systems, almost-positive roots, and the
finite-type invariants `n_clusters` / `f_vector` / `h_vector`.

**Folding.** `fold(q, σ)` folds a symmetric quiver by an admissible vertex automorphism, with
`is_admissible_folding` as the precondition - the route from the simply-laced types to the
non-simply-laced ones.

**Mutation classes and exchange graphs.** BFS with deduplication and a truncation guard for
mutation-infinite types.

**Green sequences and DT theory.** Green/red predicates, maximal green sequence enumeration
(`maximal_green_sequences`, `mgs_search`), the Donaldson–Thomas transformation
(`dt_transformation`, sequence-independent by Keller's theorem), and Kontsevich–Soibelman
quantum dilogarithm words with their DT invariants $\Omega(\gamma)$.

**Y-systems and periodicity.** `mutation_period` for a sequence, `y_system` for the bipartite
Zamolodchikov dynamics - period $h+2$ in finite type.

**Grassmannians and symbol alphabets.** `grassmannian(k, n)` builds the Plücker seed of
$\mathrm{Gr}(k,n)$ (Scott 2006); `symbol_alphabet` enumerates the letters - 9 for
$\mathrm{Gr}(4,6)$, 42 for $\mathrm{Gr}(4,7)$ - with `cluster_adjacency_matrix` for which may
sit next to each other.

**Friezes.** SL₂ frieze patterns from triangulated polygons, with their integrality.

**Isomorphism and sampling.** `canonical_form` / `is_isomorphic` give normal forms up to
relabeling (`mutation_class` does not dedup); `random_quiver` / `random_mutate` sample quivers
and walks.

## Extensions

The core depends only on AbstractAlgebra; plotting is loaded on demand.

```julia
using Graphs, GraphMakie, GLMakie
plot_quiver(q)
plot_exchange_graph(q)
plot_green_sequence(q, seq)
plot_g_vector_fan(q)
plot_frieze(f)
```

`Graphs` alone (without a Makie backend) gets you the exchange graph and quiver as
`Graphs.jl` objects, for use with the rest of that ecosystem.

## Related packages

This package is self-contained. It is also the discrete foundation of a family of Julia
packages for **exact and asymptotic methods** - where a combinatorial structure computed
exactly here controls a continuous object that is only defined asymptotically:

- [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) - the continuous side:
  divergent series, Borel–Padé summation, transseries and alien calculus. Independent of this
  package.
- [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl) - the bridge: Stokes graphs
  of a Schrödinger-type ODE, connected to this package by the Iwaki–Nakanishi dictionary,
  where Voros-symbol jumps are y-mutations and BPS spectra are maximal green sequences.

## License

MIT
