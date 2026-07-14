module ClusterAlgebrasGraphMakieExt

using ClusterAlgebras, Graphs, GraphMakie, Makie

function _quiver_graph_args(q::Quiver)
    n = q.n_mutable + q.n_frozen
    g = SimpleDiGraph(n)
    edge_label_map = Dict{Tuple{Int,Int},String}()
    for i in 1:n, j in 1:n
        w = q.B[i, j]
        if w > 0
            add_edge!(g, i, j)
            w > 1 && (edge_label_map[(i, j)] = string(w))
        end
    end
    node_color  = [i <= q.n_mutable ? :steelblue : :lightgray for i in 1:n]
    node_marker = [i <= q.n_mutable ? :circle : :rect for i in 1:n]
    elabels     = [get(edge_label_map, (src(e), dst(e)), "") for e in edges(g)]
    g, node_color, node_marker, elabels
end

function ClusterAlgebras.plot_quiver(q::Quiver; kwargs...)
    g, node_color, node_marker, elabels = _quiver_graph_args(q)
    p = graphplot(g; node_color, node_marker, nlabels=q.labels, elabels, kwargs...)
    hidedecorations!(p.axis)
    hidespines!(p.axis)
    p
end

function ClusterAlgebras.plot_quiver!(ax, q::Quiver; kwargs...)
    g, node_color, node_marker, elabels = _quiver_graph_args(q)
    graphplot!(ax, g; node_color, node_marker, nlabels=q.labels, elabels, kwargs...)
end

# ─── Exchange / mutation graph ────────────────────────────────────────────────

function _exchange_graph_plot(obj; nlabels, kwargs...)
    g = SimpleGraph(obj)                       # from ClusterAlgebrasGraphsExt
    n = nv(g)
    labs = nlabels === true  ? string.(1:n) :
           nlabels === false ? nothing       : nlabels
    p = graphplot(g; node_color = :seagreen, nlabels = labs, kwargs...)
    hidedecorations!(p.axis)
    hidespines!(p.axis)
    p
end

"""
    plot_exchange_graph(mc::MutationClass; nlabels = true, kwargs...)
    plot_exchange_graph(eg::ExchangeGraph; nlabels = true, kwargs...)

Render the mutation graph of `mc` (vertices = quivers) or the exchange graph of
`eg` (vertices = seeds) as an undirected graph: two vertices are joined when a
single mutation relates them. For a finite-type `ExchangeGraph` this is the
1-skeleton of the generalized associahedron.

`nlabels = true` (default) labels each vertex with its 1-based index in the
class; pass `false` for none, or a vector of custom labels. Extra keywords pass
through to GraphMakie's `graphplot`. Requires `using Graphs, GraphMakie` and a
Makie backend.
"""
ClusterAlgebras.plot_exchange_graph(mc::ClusterAlgebras.MutationClass;
                                    nlabels = true, kwargs...) =
    _exchange_graph_plot(mc; nlabels, kwargs...)

ClusterAlgebras.plot_exchange_graph(eg::ClusterAlgebras.ExchangeGraph;
                                    nlabels = true, kwargs...) =
    _exchange_graph_plot(eg; nlabels, kwargs...)

end
