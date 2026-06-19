# Convenience constructors: named Dynkin types, string parsing, edge lists, and
# permissive matrix-type coercion.

# ─── Permissive matrix coercion ─────────────────────────────────────────────

Quiver(B::AbstractMatrix{<:Integer}) =
    Quiver(Matrix{Int}(B))
Quiver(B::AbstractMatrix{<:Integer}, n_mutable::Int) =
    Quiver(Matrix{Int}(B), n_mutable)
Quiver(B::AbstractMatrix{<:Integer}, n_mutable::Int, d::Vector{Int}) =
    Quiver(Matrix{Int}(B), n_mutable, d)
Quiver(B::AbstractMatrix{<:Integer}, n_mutable::Int, d::Vector{Int}, labels::Vector{String}) =
    Quiver(Matrix{Int}(B), n_mutable, d, labels)

# ─── Edge-list constructor ───────────────────────────────────────────────────
#
# Each edge is a 2-tuple (src, dst) or 3-tuple (src, dst, multiplicity).
# Multiple edges in the same direction accumulate; opposite-direction edges
# cancel (the exchange matrix is antisymmetric by construction).

function Quiver(edges::AbstractVector{T}) where {T <: Tuple}
    isempty(edges) && throw(InvalidArgument("edge list must be non-empty"))
    n = 0
    for e in edges
        length(e) >= 2 || throw(InvalidArgument(
            "each edge must be a 2- or 3-tuple (src, dst[, weight]), got $e"))
        n = max(n, e[1], e[2])
    end
    B = zeros(Int, n, n)
    for e in edges
        i, j = Int(e[1]), Int(e[2])
        w = length(e) >= 3 ? Int(e[3]) : 1
        B[i, j] += w
        B[j, i] -= w
    end
    Quiver(B)
end

# ─── Dynkin exchange matrices ─────────────────────────────────────────────────
#
# Returns (B, d) for the standard acyclic orientation of each finite Dynkin type.
# All arrows point in the direction of increasing index; for non-simply-laced
# types the multiple bond is at the position specified below.
#
# Symmetrizer conventions (d[i]*B[i,j] == -d[j]*B[j,i]):
#   A_n : d = [1,...,1]
#   B_n : double bond B[n-1,n]=1, B[n,n-1]=-2  →  d = [2,...,2,1]
#   C_n : double bond B[n-1,n]=2, B[n,n-1]=-1  →  d = [1,...,1,2]
#   D_n : all simply laced, fork at vertex n-2  →  d = [1,...,1]
#   E_n : all simply laced, branch from 3 to n →  d = [1,...,1]
#   F_4 : double bond B[2,3]=1, B[3,2]=-2       →  d = [2,2,1,1]
#   G_2 : triple bond B[1,2]=1, B[2,1]=-3       →  d = [3,1]

function _dynkin_exchange_matrix(type::Symbol, n::Int)
    if type === :A
        n >= 1 || throw(InvalidArgument("A_n requires n ≥ 1, got n = $n"))
        B = zeros(Int, n, n)
        for i in 1:n-1; B[i,i+1] = 1; B[i+1,i] = -1; end
        return B, ones(Int, n)

    elseif type === :B
        n >= 2 || throw(InvalidArgument("B_n requires n ≥ 2, got n = $n"))
        B = zeros(Int, n, n)
        for i in 1:n-2; B[i,i+1] = 1; B[i+1,i] = -1; end
        B[n-1,n] = 1; B[n,n-1] = -2
        d = fill(2, n); d[n] = 1
        return B, d

    elseif type === :C
        n >= 2 || throw(InvalidArgument("C_n requires n ≥ 2, got n = $n"))
        B = zeros(Int, n, n)
        for i in 1:n-2; B[i,i+1] = 1; B[i+1,i] = -1; end
        B[n-1,n] = 2; B[n,n-1] = -1
        d = ones(Int, n); d[n] = 2
        return B, d

    elseif type === :D
        n >= 4 || throw(InvalidArgument("D_n requires n ≥ 4, got n = $n"))
        B = zeros(Int, n, n)
        for i in 1:n-3; B[i,i+1] = 1; B[i+1,i] = -1; end
        B[n-2,n-1] = 1; B[n-1,n-2] = -1
        B[n-2,n]   = 1; B[n,n-2]   = -1
        return B, ones(Int, n)

    elseif type === :E
        n in (6,7,8) || throw(InvalidArgument("E_n requires n ∈ {6,7,8}, got n = $n"))
        B = zeros(Int, n, n)
        for i in 1:n-2; B[i,i+1] = 1; B[i+1,i] = -1; end
        B[3,n] = 1; B[n,3] = -1
        return B, ones(Int, n)

    elseif type === :F
        n == 4 || throw(InvalidArgument("F_4 has rank 4, got n = $n"))
        B = [0 1 0 0; -1 0 1 0; 0 -2 0 1; 0 0 -1 0]
        return B, [2, 2, 1, 1]

    elseif type === :G
        n == 2 || throw(InvalidArgument("G_2 has rank 2, got n = $n"))
        return [0 1; -3 0], [3, 1]

    else
        throw(InvalidArgument(
            "unknown Dynkin type :$type; expected one of :A, :B, :C, :D, :E, :F, :G"))
    end
end

# ─── Symbol constructor ───────────────────────────────────────────────────────

function Quiver(type::Symbol, n::Int)
    B, d = _dynkin_exchange_matrix(type, n)
    Quiver(B, n, d, string.(1:n))
end

# ─── String constructor ───────────────────────────────────────────────────────

const _DYNKIN_CHAR_MAP = Dict(
    'A' => :A, 'B' => :B, 'C' => :C, 'D' => :D,
    'E' => :E, 'F' => :F, 'G' => :G,
)

function Quiver(s::AbstractString)
    isempty(s) && throw(InvalidArgument(
        "empty string is not a valid Dynkin type specification"))
    tc = uppercase(first(s))
    haskey(_DYNKIN_CHAR_MAP, tc) || throw(InvalidArgument(
        "unknown Dynkin type '$tc' in \"$s\"; expected one of A,B,C,D,E,F,G"))
    rest = s[nextind(s, firstindex(s)):end]
    n = tryparse(Int, rest)
    n !== nothing || throw(InvalidArgument(
        "cannot parse rank from \"$s\"; expected format like \"A3\" or \"D4\""))
    Quiver(_DYNKIN_CHAR_MAP[tc], n)
end
