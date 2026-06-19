abstract type ClusterAlgebraError <: Exception end

struct NotSkewSymmetrizable <: ClusterAlgebraError
    i::Int
    j::Int
    lhs::Int  # d[i]*B[i,j]
    rhs::Int  # -d[j]*B[j,i]
end

struct FrozenVertexMutation <: ClusterAlgebraError
    k::Int
    n_mutable::Int
end

struct InvalidVertex <: ClusterAlgebraError
    k::Int
    n_total::Int
end

struct InvalidArgument <: ClusterAlgebraError
    msg::String
end

function Base.showerror(io::IO, e::NotSkewSymmetrizable)
    print(io, "NotSkewSymmetrizable: exchange matrix is not skew-symmetrizable with the given " *
              "symmetrizer d at (i=$(e.i), j=$(e.j)): " *
              "d[$(e.i)]·B[$(e.i),$(e.j)] = $(e.lhs) ≠ $(e.rhs) = -d[$(e.j)]·B[$(e.j),$(e.i)]")
end

function Base.showerror(io::IO, e::FrozenVertexMutation)
    print(io, "FrozenVertexMutation: vertex $(e.k) is frozen and cannot be mutated " *
              "(mutable vertices: 1:$(e.n_mutable))")
end

function Base.showerror(io::IO, e::InvalidVertex)
    print(io, "InvalidVertex: vertex $(e.k) is out of range " *
              "(valid vertices: 1:$(e.n_total))")
end

function Base.showerror(io::IO, e::InvalidArgument)
    print(io, "InvalidArgument: $(e.msg)")
end
