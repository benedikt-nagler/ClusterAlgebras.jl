```@meta
CurrentModule = ClusterAlgebras
```

# Interchange formats

Quivers can be exchanged with other cluster algebra software through two formats.

`dig6` is the pair `(string, edges)` that SageMath stores mutation classes in: nauty's
`digraph6` encoding of the underlying digraph, together with the labels of the arcs that
are not simply laced. It records the digraph only, so the frozen/mutable split and the
symmetrizers are supplied on the way back in.

`.qmu` is the save format of Bernhard Keller's
[quiver mutation applet](https://webusers.imj-prg.fr/~bernhard.keller/quivermutation/),
which SageMath also writes. It carries frozen vertices and vertex positions.

`to_dot` (see [Quivers](@ref)) covers Graphviz.

```@autodocs
Modules = [ClusterAlgebras]
Private = false
Pages = ["src/interop.jl"]
```
