# ClusterAlgebras.jl

A Julia package for exact cluster algebra computation.

## Features

- **Quivers** — skew-symmetrizable exchange matrices with frozen vertices, edge-list and Dynkin-type constructors (`Quiver(:A, 4)`, `Quiver("D5")`, etc.), DOT export
- **Seeds** — cluster variables as elements of `Frac(ZZ[x₁,…,xₙ])`, mutation along single vertices or sequences
- **Coefficients** — trivial and principal coefficient systems; C-matrix, G-matrix, F-polynomials, rational and tropical y-variables, separation formula
- **Green sequences** — green/red vertex predicates, maximal green sequence enumeration
- **Finite and affine type detection** — Cartan companion, root systems, almost-positive roots
- **Mutation class and exchange graph** — BFS with deduplication; truncation guard for mutation-infinite types
- **Enumerative invariants** — cluster variable counts, f-vectors, h-vectors for finite-type algebras
- **Friezes** — Conway–Coxeter SL₂ frieze patterns from triangulated polygons
- **Optional extensions** — graph layout via Graphs.jl + GraphMakie

Requires Julia 1.10 or later.

## Quick start

```julia
using ClusterAlgebras

# Build a quiver from an exchange matrix
B = [0 1 0; -1 0 1; 0 -1 0]
q = Quiver(B)

# Or use a named Dynkin type
q = Quiver(:A, 3)       # A₃ quiver
q = Quiver("D4")        # D₄ quiver

# Mutate the quiver at vertex 2
q2 = mutate(q, 2)

# Create a seed (cluster variables as rational functions)
s = Seed(q)
s2 = mutate(s, 2)
s3 = mutate(s, [1, 2, 3])   # sequence mutation

# Attach principal coefficients
sp = extend(s)
sp = mutate(sp, [1, 2, 1])

# Inspect c-vectors, g-vectors, F-polynomials
cmatrix(sp)
gmatrix(sp)
fpolynomials(sp)
f_polynomial(sp, 1)

# Finite-type detection
is_finite_type(q)         # true for A₃
cartan_type(q)            # (:A, 3)

# Mutation class
mc = mutation_class(q)
length(mc)                # number of distinct quivers in the class

# Green sequences
is_green(sp, 1)
maximal_green_sequences(sp)
```

## Optional extensions

Load additional packages to unlock extensions:

```julia
using Graphs, GraphMakie, GLMakie
plot_quiver(q)            # interactive quiver plot
```
