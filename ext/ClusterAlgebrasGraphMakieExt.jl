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

end
