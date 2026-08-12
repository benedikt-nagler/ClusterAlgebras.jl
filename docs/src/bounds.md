```@meta
CurrentModule = ClusterAlgebras
```

# Bounds and the greedy basis

Two questions about a cluster algebra as a *ring*: what does it sit between, and what is a
basis for it.

## Upper and lower bounds

For a geometric-type seed ``\Sigma = (x, \tilde{B})`` over the coefficient ring
``\mathbb{Z}P``, Berenstein–Fomin–Zelevinsky sandwich the cluster algebra between an
intersection of Laurent rings and a ring of explicit generators:

```math
L(\Sigma) = \mathbb{Z}P[x_1, x_1', \dots, x_n, x_n']
  \;\subseteq\; \mathcal{A} \;\subseteq\;
U(\Sigma) = \mathbb{Z}P[x^{\pm 1}] \cap \mathbb{Z}P[x_1^{\pm 1}] \cap \dots \cap \mathbb{Z}P[x_n^{\pm 1}],
```

where ``x_k`` is the cluster adjacent to ``x`` at ``k``. [`lower_bound_generators`](@ref) and
[`upper_bound_generators`](@ref) produce the two sides, and [`in_upper_bound`](@ref) /
[`is_laurent`](@ref) test membership.

The sandwich collapses exactly when the seed is **coprime and acyclic**
([`is_coprime`](@ref), [`is_acyclic`](@ref)): then ``L(\Sigma) = U(\Sigma)``, so the cluster
algebra equals the upper cluster algebra and the [`standard_monomials`](@ref) are a
``\mathbb{Z}P``-basis. The condition is an *iff*, which is what makes
[`bound_certificate`](@ref) a certificate and not a heuristic.

Everything here is exact and Gröbner-free. Generators of ``U(\Sigma)`` in the **non-acyclic**
case are deliberately out of scope: that computation needs a Gröbner basis, so it belongs
behind the OSCAR firewall in a downstream package and not in this dependency-light core.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/bounds.jl"]
```

## The greedy basis (rank 2)

In rank 2, Lee–Li–Zelevinsky's greedy elements form a positive ``\mathbb{Z}``-basis
containing every cluster monomial. [`greedy_element`](@ref) builds the element at a given
denominator vector and [`greedy_coefficients`](@ref) exposes the coefficients of the
recursion.

For ``bc \le 3`` the greedy basis *is* the set of cluster monomials; beyond that it acquires
"imaginary" elements that no cluster monomial supplies, the first being
``(1 + x_1^2 + x_2^2)/(x_1x_2)`` at ``b = c = 2``. A cluster variable equals the greedy
element at its own denominator vector exactly, which is the oracle that pins the
implementation.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/greedy.jl"]
```
