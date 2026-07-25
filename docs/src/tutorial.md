```@meta
CurrentModule = ClusterAlgebras
```

# Tutorial: the ``A_2`` pentagon

The rank-2 cluster algebra of type ``A_2`` is the smallest example that shows every feature
of the theory: finitely many cluster variables, the Laurent phenomenon, and a periodicity
that is invisible from the exchange relation alone. Everything below is small enough to
check by hand.

All output on this page is verified when the documentation is built, so it is exactly what
you will see in the REPL.

## A quiver and its seed

```jldoctest a2
julia> using ClusterAlgebras

julia> q = Quiver(:A, 2)
Quiver with 2 vertices (2 mutable, 0 frozen)
Exchange matrix B:
   0  1   (1)
  -1  0   (2)
```

A [`Seed`](@ref) attaches a cluster to the quiver. The variables live in
``\operatorname{Frac}(\mathbb{Z}[x_1,x_2])``:

```jldoctest a2
julia> s = Seed(q)
Quiver with 2 vertices (2 mutable, 0 frozen)
Exchange matrix B:
   0  1   (1)
  -1  0   (2)
Cluster variables:
  [1]: x_1
  [2]: x_2
```

## Mutation

Mutating at vertex 1 replaces ``x_1`` by ``(x_2+1)/x_1`` and reverses the arrow. Seeds are
immutable - [`mutate`](@ref) returns a new seed and records the path taken:

```jldoctest a2
julia> s1 = mutate(s, 1)
Quiver with 2 vertices (2 mutable, 0 frozen)
Exchange matrix B:
   0 -1   (1)
   1  0   (2)
Cluster variables:
  [1]: (x_2 + 1)//x_1
  [2]: x_2
Mutation path: [1]
```

Individual cluster variables are reachable by indexing the seed. Mutating again at vertex 2
produces the variable in which the Laurent phenomenon is easiest to see - the division by
``x_1 x_2`` leaves no denominator behind beyond the initial cluster:

```jldoctest a2
julia> s2 = mutate(s1, 2);

julia> s2[2]
(x_1 + x_2 + 1)//(x_1*x_2)
```

## The pentagon

Alternating mutations at 1 and 2 returns the initial cluster after five steps, with the two
variables swapped. This is the pentagon recurrence:

```jldoctest a2
julia> mutate(s, [1, 2, 1, 2, 1])
Quiver with 2 vertices (2 mutable, 0 frozen)
Exchange matrix B:
   0 -1   (1)
   1  0   (2)
Cluster variables:
  [1]: x_2
  [2]: x_1
Mutation path: [1, 2, 1, 2, 1]
```

[`mutation_period`](@ref) finds that period without being told what to expect:

```jldoctest a2
julia> mutation_period(Seed(q), [1, 2])
5
```

Five mutations, five distinct cluster variables - ``x_1``, ``x_2`` and the three produced
above. The counting functions agree, and the classification confirms the type:

```jldoctest a2
julia> cartan_type(q)
(:A, 2)

julia> is_finite_type(q)
true

julia> n_cluster_variables(q)
5

julia> n_clusters(q)
5
```

## Coefficients

[`extend`](@ref) attaches principal coefficients, which carry the c-vectors, g-vectors and
F-polynomials:

```jldoctest a2
julia> sp = mutate(extend(Seed(q)), [1, 2]);

julia> cmatrix(sp)
2×2 Matrix{Int64}:
 0  -1
 1  -1

julia> gmatrix(sp)
2×2 Matrix{Int64}:
 -1  -1
  1   0

julia> f_polynomial(sp, 1)
y1 + 1
```

The separation formula of Fomin–Zelevinsky reconstructs a cluster variable from its
g-vector and F-polynomial alone:

```jldoctest a2
julia> separation_formula(sp, 1)
(x_2 + y1)//x_1
```

Sign coherence - every c-vector is nonzero and has entries of a single sign - is a theorem
of Derksen–Weyman–Zelevinsky, and holds here:

```jldoctest a2
julia> is_sign_coherent(sp)
true
```

## Green sequences

A maximal green sequence mutates every c-vector from green to red. Type ``A_2`` has two,
of lengths 2 and 3 - the two sides of the pentagon:

```jldoctest a2
julia> maximal_green_sequences(extend(Seed(q)))
2-element Vector{Vector{Int64}}:
 [2, 1]
 [1, 2, 1]
```

The Donaldson–Thomas transformation is built from one of them:

```jldoctest a2
julia> dt = dt_transformation(q);

julia> dt.sequence
2-element Vector{Int64}:
 2
 1

julia> dt.sigma
2-element Vector{Int64}:
 1
 2
```

## Friezes

Friezes are the numerical shadow of a type ``A`` cluster algebra: every entry is a positive
integer, and each diamond satisfies ``ad - bc = 1``.

```jldoctest a2
julia> frieze(4)
Frieze(n = 4, width = 1)
Quiddity: [2, 1, 2, 1]
 0 0 0 0
 1 1 1 1
 2 1 2 1
 1 1 1 1
 0 0 0 0
```

## Where to go next

- [Quivers](quivers.md) and [Seeds and mutation](seeds.md) for the core types.
- [Coefficients](coefficients.md) for c-/g-vectors, F-polynomials and the separation formula.
- [Green sequences and DT](green.md) for maximal green sequences, Y-system periodicity and
  the Donaldson–Thomas transformation.
- [Grassmannians](grassmannian.md) for Plücker seeds and symbol alphabets.
