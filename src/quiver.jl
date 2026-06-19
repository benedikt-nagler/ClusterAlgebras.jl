struct Quiver
    B::Matrix{Int}
    n_mutable::Int
    n_frozen::Int
    d::Vector{Int}         # symmetrizers: d[i]*B[i,j] == -d[j]*B[j,i] for mutable i,j
    labels::Vector{String} # vertex labels, length n_mutable + n_frozen

    function Quiver(B::Matrix{Int}, n_mutable::Int, d::Vector{Int}, labels::Vector{String})
        n_total = size(B, 1)
        size(B, 2) == n_total ||
            throw(InvalidArgument("exchange matrix must be square, got $(size(B))"))
        0 <= n_mutable <= n_total ||
            throw(InvalidArgument("n_mutable=$n_mutable out of range 0:$n_total"))
        length(d) == n_mutable ||
            throw(InvalidArgument("symmetrizer d must have length n_mutable=$n_mutable, got $(length(d))"))
        length(labels) == n_total ||
            throw(InvalidArgument("labels must have length $n_total, got $(length(labels))"))
        all(>(0), d) ||
            throw(InvalidArgument("symmetrizers must be positive integers"))
        for i in 1:n_mutable, j in 1:n_mutable
            lhs = d[i] * B[i, j]
            rhs = -d[j] * B[j, i]
            lhs == rhs || throw(NotSkewSymmetrizable(i, j, lhs, rhs))
        end
        new(B, n_mutable, n_total - n_mutable, d, labels)
    end
end

# All-mutable, skew-symmetric (d = ones)
Quiver(B::Matrix{Int}) =
    Quiver(B, size(B, 1), ones(Int, size(B, 1)), string.(1:size(B, 1)))

# With frozen vertices, skew-symmetric mutable block
Quiver(B::Matrix{Int}, n_mutable::Int) =
    Quiver(B, n_mutable, ones(Int, n_mutable), string.(1:size(B, 1)))

# All-mutable, skew-symmetrizable
Quiver(B::Matrix{Int}, n_mutable::Int, d::Vector{Int}) =
    Quiver(B, n_mutable, d, string.(1:size(B, 1)))

Base.:(==)(a::Quiver, b::Quiver) =
    a.n_mutable == b.n_mutable && a.B == b.B && a.d == b.d && a.labels == b.labels

nvertices(q::Quiver) = q.n_mutable + q.n_frozen
labels(q::Quiver)    = q.labels

function is_frozen(q::Quiver, k::Int)
    n_total = q.n_mutable + q.n_frozen
    1 <= k <= n_total || throw(InvalidVertex(k, n_total))
    return k > q.n_mutable
end

function to_dot(q::Quiver)::String
    n = q.n_mutable + q.n_frozen
    lines = ["digraph quiver {", "  rankdir=LR;"]
    for i in 1:n
        shape = i <= q.n_mutable ? "circle" : "box"
        push!(lines, "  $i [label=\"$(q.labels[i])\", shape=$shape];")
    end
    for i in 1:n, j in 1:n
        w = q.B[i, j]
        if w > 0
            label = w == 1 ? "" : " [label=\"$w\"]"
            push!(lines, "  $i -> $j$label;")
        end
    end
    push!(lines, "}")
    join(lines, "\n")
end
