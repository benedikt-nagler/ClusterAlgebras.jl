# ClusterAlgebras.jl

[![CI](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ClusterAlgebras.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact computation with cluster algebras in Julia: quivers and seed mutation, principal
coefficients, finite-type classification, exchange graphs, green sequences, and the
cluster structures on Grassmannians.

Cluster algebras were introduced by Fomin and Zelevinsky to give an algebraic framework for
canonical bases in Lie theory. One starts from a *seed*: a cluster $(x_1, \ldots, x_n)$ of
rational functions together with a skew-symmetrizable integer matrix $B$. *Mutation* at
index $k$ replaces $x_k$ using

$$x_k x_k' = \prod_{b_{ik}>0} x_i^{b_{ik}} + \prod_{b_{ik}<0} x_i^{-b_{ik}},$$

and updates $B$ accordingly. Iterating this generates the algebra. Despite the division,
every cluster variable is a Laurent polynomial in any initial cluster — the Laurent
phenomenon.

The package works over `Frac(ZZ[x₁,…,xₙ])` via
[AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl), so results are exact
rather than floating-point. Seeds are immutable and `mutate` returns a new seed.

## Installation

Not yet registered in General. From the Julia REPL:

```julia
pkg> add https://github.com/benedikt-nagler/ClusterAlgebras.jl
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

**Quivers and seeds.** Skew-symmetrizable exchange matrices with frozen vertices, built from
a matrix, an edge list, or a Dynkin name. Mutation at a vertex, a label, or along a
sequence. `to_dot` for Graphviz export.

**Coefficients.** Trivial, principal, and general extended coefficient systems; C- and
G-matrices, F-polynomials, rational and tropical y-variables, the separation formula, and
denominator vectors.

**Classification.** Cartan companion, root systems, almost-positive roots, and
`is_finite_type` / `is_affine_type` / `cartan_type`, each of which searches for an acyclic
representative first, so the answers hold for non-acyclic seeds too. Enumerative invariants
(`n_clusters`, `f_vector`, `h_vector`) for finite type.

**Mutation classes and exchange graphs.** BFS with deduplication and a truncation guard for
mutation-infinite types.

**Green sequences and DT theory.** Green/red vertex predicates, maximal green sequence
enumeration, the Donaldson–Thomas transformation (`dt_transformation`, independent of the
chosen sequence by Keller's theorem), and Kontsevich–Soibelman quantum dilogarithm words
with their DT invariants $\Omega(\gamma)$ (`quantum_dilog_word`, `omega`,
`ks_dilog_product`).

**Y-systems and periodicity.** `mutation_period` for a sequence, and `y_system` for the
bipartite Zamolodchikov dynamics — period $h+2$ in finite type. See the
[Y-systems and TBA notebook](examples/y_systems_and_tba.ipynb).

**Grassmannians and amplitudes.** `grassmannian(k, n)` builds the Plücker-coordinate seed of
$\mathrm{Gr}(k,n)$ (Scott 2006), and `symbol_alphabet` enumerates the cluster variables as
functions of Plücker coordinates — the 9-letter $\mathrm{Gr}(4,6)$ and 42-letter
$\mathrm{Gr}(4,7)$ alphabets that appear in $\mathcal{N}=4$ SYM amplitudes.

**Friezes.** SL₂ frieze patterns from triangulated polygons.

**Isomorphism and sampling.** `canonical_form` / `is_isomorphic` give normal forms up to
relabeling (useful for deduplicating datasets, which `mutation_class` does not do), and
`random_quiver` / `random_mutate` sample quivers and walks.

## Extensions

The core depends only on AbstractAlgebra; plotting is loaded on demand.

```julia
using Graphs, GraphMakie, GLMakie
plot_quiver(q)
plot_exchange_graph(q)
plot_g_vector_fan(q)
```

## Related packages

Part of a family of packages around cluster algebras and exact WKB.
[Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) handles Borel–Padé
summation and transseries; [ExactWKB.jl](https://github.com/benedikt-nagler/ExactWKB.jl)
builds Stokes graphs and bridges to this package via the Iwaki–Nakanishi dictionary, where
Voros jumps are y-mutations and BPS spectra are maximal green sequences.

## License

MIT
