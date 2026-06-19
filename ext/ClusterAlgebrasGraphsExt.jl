module ClusterAlgebrasGraphsExt

using ClusterAlgebras, Graphs

"""
    Graphs.SimpleGraph(mc::MutationClass) → SimpleGraph

Construct an undirected `SimpleGraph` whose vertices are the quivers in `mc` and
whose edges connect quivers that are related by a single mutation.
"""
function Graphs.SimpleGraph(mc::ClusterAlgebras.MutationClass)
    n = length(mc)
    g = SimpleGraph(n)
    for i in 1:n
        for j in mc.adj[i]
            i < j && add_edge!(g, i, j)
        end
    end
    return g
end

"""
    Graphs.SimpleGraph(eg::ExchangeGraph) → SimpleGraph

Construct an undirected `SimpleGraph` whose vertices are the seeds in `eg` and
whose edges connect seeds that are related by a single mutation.  For finite-type
cluster algebras this is the 1-skeleton of the generalized associahedron.
"""
function Graphs.SimpleGraph(eg::ClusterAlgebras.ExchangeGraph)
    n = length(eg)
    g = SimpleGraph(n)
    for i in 1:n
        for j in eg.adj[i]
            i < j && add_edge!(g, i, j)
        end
    end
    return g
end

end
