```@meta
CurrentModule = ClusterAlgebras
```

# Types and classification

Cartan companions, root systems, finite- and affine-type recognition, folding to non-simply-laced types, and enumerative invariants.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/root_system.jl", "src/folding.jl", "src/enumerative.jl"]
```

## Naming a mutation type

[`mutation_type`](@ref) identifies the mutation class a quiver belongs to and returns a
[`MutationType`](@ref): the finite types, the simply-laced affine ones, the rank-2 classes,
the exceptional ``X_6``/``X_7``, and the elliptic ``E_6^{(1,1)}``, ``E_7^{(1,1)}``,
``E_8^{(1,1)}``. [`mutation_types`](@ref) lists every match when the answer is not unique.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/mutation_type.jl"]
```

## Surface types and block decomposition

The mutation-finite classes that `mutation_type` does *not* name are the ones coming from
marked surfaces, and those are recognised structurally instead. Fomin–Shapiro–Thurston show
that a quiver arises from a triangulated surface exactly when it decomposes into the blocks of
a fixed list ([`FST_BLOCKS`](@ref)), glued at their outlets. [`block_decomposition`](@ref)
finds such a decomposition, [`is_block_decomposable`](@ref) and [`is_surface_type`](@ref) are
the predicates, and [`reassemble`](@ref) rebuilds the quiver from a
[`BlockDecomposition`](@ref), so a decomposition can always be checked against the quiver it
claims to explain.

This is the reverse of the bridge that `ClusterSurfaces.jl` walks forwards: given the quiver
alone, recover the surface it could have come from.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/block_decomposition.jl"]
```
