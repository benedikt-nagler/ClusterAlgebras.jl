```@meta
CurrentModule = ClusterAlgebras
```

# ClusterAlgebras.jl

Exact computation with cluster algebras in Julia: quivers and seed mutation, principal
coefficients, finite-type classification, exchange graphs, green sequences, and the cluster
structures on Grassmannians.

Cluster algebras were introduced by Fomin and Zelevinsky to give an algebraic framework for
canonical bases in Lie theory. One starts from a *seed*: a cluster ``(x_1, \ldots, x_n)`` of
rational functions together with a skew-symmetrizable integer matrix ``B``. *Mutation* at
index ``k`` replaces ``x_k`` using

```math
x_k x_k' = \prod_{b_{ik}>0} x_i^{b_{ik}} + \prod_{b_{ik}<0} x_i^{-b_{ik}},
```

and updates ``B`` accordingly. Iterating this generates the algebra. Despite the division,
every cluster variable is a Laurent polynomial in any initial cluster — the Laurent
phenomenon.

The package works over ``\operatorname{Frac}(\mathbb{Z}[x_1,\ldots,x_n])`` via
[AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl), so results are exact
rather than floating-point. Seeds are immutable and [`mutate`](@ref) returns a new seed.

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

Plotting lives in package extensions: load `Graphs`, `GraphMakie` and a Makie backend to
enable `plot_quiver`, `plot_exchange_graph` and friends.

See the [API](@ref) page for the full list of exported names.
