# Changelog

All notable changes to this project are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.1] - 2026-07-28

Documentation-only release. No changes to the public API.

### Added

- A [Documenter](https://documenter.juliadocs.org) site: tutorial, a manual page per topic
  (quivers, seeds and mutation, coefficients, classification, mutation classes, green
  sequences, friezes, Grassmannians), an errors page and an API index.
- A module-level docstring for `ClusterAlgebras` listing the entry points.
- This changelog.

### Changed

- The README and the documentation home page now give `pkg> add ClusterAlgebras` as the
  installation instruction, since the package is in the General Registry.
- The README carries a registry version badge.

### Fixed

- The `CompatHelper` workflow.

## [0.1.0] - 2026-07-25

First release, registered on the Julia General Registry.

### Added

- **Quivers and seeds.** Skew-symmetrizable exchange matrices with frozen vertices, built
  from a matrix, an edge list or a Dynkin name. Mutation at a vertex, at a label, or along a
  sequence, plus `verify_mutation_sequence` and `to_dot`. Cluster variables are exact
  elements of `Frac(ZZ[x₁,…,xₙ])` via `AbstractAlgebra`; quivers and seeds are immutable and
  `mutate` returns a new object.
- **Coefficients.** Trivial, principal and general extended coefficient systems as a type
  parameter on `Seed`, so one `mutate` covers all three. C- and G-matrices, F-polynomials,
  rational and tropical y-variables, `y_hat`, the separation formula, denominator vectors.
- **Classification.** `is_finite_type`, `is_affine_type`, `cartan_type` and
  `is_mutation_finite`, each searching for an acyclic representative so the answers hold for
  non-acyclic seeds too. Cartan companion, root systems, almost-positive roots, and the
  finite-type invariants `n_clusters`, `f_vector`, `h_vector`.
- **Folding.** `fold(q, sigma)` for an admissible vertex automorphism, with
  `is_admissible_folding` as the precondition, giving the non-simply-laced types.
- **Mutation classes and exchange graphs.** BFS with deduplication and a truncation guard for
  mutation-infinite types.
- **Green sequences and DT theory.** Green and red predicates, maximal green sequence
  enumeration (`maximal_green_sequences`, `mgs_search`), the Donaldson-Thomas transformation
  `dt_transformation` (sequence-independent by Keller's theorem), and Kontsevich-Soibelman
  quantum dilogarithm words with their DT invariants.
- **Y-systems and periodicity.** `mutation_period` for a sequence and `y_system` for the
  bipartite Zamolodchikov dynamics, of period `h+2` in finite type.
- **Grassmannians and symbol alphabets.** `grassmannian(k, n)` for the Plücker seed of
  Gr(k,n) (Scott 2006), `symbol_alphabet` for the letters, and `cluster_adjacency_matrix`.
- **Friezes.** SL2 frieze patterns from triangulated polygons, with their integrality.
- **Isomorphism and sampling.** `canonical_form` and `is_isomorphic` for normal forms up to
  relabeling, `random_quiver` and `random_mutate` for sampling quivers and walks.
- **Plotting** through package extensions on `Graphs`, `GraphMakie` and `Makie`:
  `plot_quiver`, `plot_exchange_graph`, `plot_green_sequence`, `plot_g_vector_fan`,
  `plot_frieze`. `Graphs` alone exposes the quiver and exchange graph as `Graphs.jl` objects.

[0.1.1]: https://github.com/benedikt-nagler/ClusterAlgebras.jl/releases/tag/v0.1.1
[0.1.0]: https://github.com/benedikt-nagler/ClusterAlgebras.jl/releases/tag/v0.1.0
