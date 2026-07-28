# Interchange formats for quivers.
#
# Two external formats are supported, both matching what SageMath's cluster
# algebra package reads and writes, so seeds can travel between the two systems:
#
#   * `dig6` - a pair `(string, edges)`: nauty's `digraph6` encoding of the
#     underlying digraph together with the labels of its non-simply-laced arcs.
#     This is the canonical form Sage stores mutation classes in.
#   * `.qmu` - the plain-text save format of Bernhard Keller's quiver mutation
#     applet, <https://webusers.imj-prg.fr/~bernhard.keller/quivermutation/>.
#
# Vertices are 1-based here but 0-based in both formats; the conversion happens
# at the boundary.

const _DIG6_HEADER = ">>digraph6<<"

# nauty's R(x): pack a bit vector into printable characters, six bits each.
function _dig6_pack_bits(bits::Vector{Bool})::String
    padded = copy(bits)
    while length(padded) % 6 != 0
        push!(padded, false)
    end
    chars = Char[]
    for k in 1:6:length(padded)
        value = 0
        for offset in 0:5
            value = 2 * value + (padded[k + offset] ? 1 : 0)
        end
        push!(chars, Char(value + 63))
    end
    return String(chars)
end

function _dig6_unpack_bits(s::AbstractString, nbits::Int)::Vector{Bool}
    bits = Bool[]
    for c in s
        value = Int(c) - 63
        0 <= value <= 63 ||
            throw(InvalidArgument("invalid digraph6 character '$c' (code $(Int(c)))"))
        for offset in 5:-1:0
            push!(bits, (value >> offset) & 1 == 1)
        end
    end
    length(bits) >= nbits ||
        throw(InvalidArgument("digraph6 string is too short: need $nbits bits, got $(length(bits))"))
    return bits[1:nbits]
end

# nauty's N(n): the vertex count, one character up to 62, four beyond that.
function _dig6_encode_order(n::Int)::String
    n >= 0 || throw(InvalidArgument("vertex count must be nonnegative, got $n"))
    n <= 62 && return string(Char(n + 63))
    n <= 258047 || throw(InvalidArgument(
        "digraph6 supports at most 258047 vertices, got $n"))
    bits = [(n >> k) & 1 == 1 for k in 17:-1:0]
    return string(Char(126)) * _dig6_pack_bits(bits)
end

# Returns (n, rest) where `rest` is the packed adjacency data.
function _dig6_decode_order(s::AbstractString)
    isempty(s) && throw(InvalidArgument("empty digraph6 string"))
    chars = collect(s)
    if chars[1] == Char(126)
        length(chars) >= 4 ||
            throw(InvalidArgument("truncated digraph6 vertex count"))
        bits = _dig6_unpack_bits(String(chars[2:4]), 18)
        n = 0
        for b in bits
            n = 2 * n + (b ? 1 : 0)
        end
        return n, String(chars[5:end])
    end
    value = Int(chars[1]) - 63
    0 <= value <= 62 ||
        throw(InvalidArgument("invalid digraph6 vertex count character '$(chars[1])'"))
    return value, String(chars[2:end])
end

"""
    to_dig6(q::Quiver; canonical::Bool = false) → (String, Vector{Pair{Tuple{Int,Int},Tuple{Int,Int}}})

Encode `q` in SageMath's `dig6` form: the `digraph6` string of the underlying
digraph (an arc `i → j` for every `B[i,j] > 0`) paired with the labels of the
arcs that are **not** simply laced.

Each entry of the companion list is `(i, j) => (B[i,j], B[j,i])` with **0-based**
vertex numbers, sorted, exactly as Sage stores it; simply-laced arcs (label
`(1, -1)`) are omitted. Pass `canonical = true` to encode
[`canonical_form(q)`](@ref) instead, which makes the pair an isomorphism
invariant.

The format carries neither the frozen/mutable split nor the symmetrizers - see
[`from_dig6`](@ref) for restoring them. Use [`to_qmu`](@ref) when frozen
vertices matter.

```jldoctest
julia> to_dig6(Quiver(:A, 2))
("&AO", Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}}[])
```
"""
function to_dig6(q::Quiver; canonical::Bool = false)
    q = canonical ? canonical_form(q) : q
    n = nvertices(q)
    bits = Bool[]
    edges = Pair{Tuple{Int, Int}, Tuple{Int, Int}}[]
    for i in 1:n, j in 1:n
        arc = q.B[i, j] > 0
        push!(bits, arc)
        if arc && (q.B[i, j], q.B[j, i]) != (1, -1)
            push!(edges, (i - 1, j - 1) => (q.B[i, j], q.B[j, i]))
        end
    end
    return "&" * _dig6_encode_order(n) * _dig6_pack_bits(bits), edges
end

"""
    from_dig6(data; n_mutable = nothing, d = nothing) → Quiver
    from_dig6(dig6::AbstractString, edges = (); n_mutable = nothing, d = nothing) → Quiver

Decode SageMath's `dig6` form back into a `Quiver`; `data` is the
`(dig6, edges)` pair produced by [`to_dig6`](@ref). A leading `>>digraph6<<`
header is accepted. `edges` may be any iterable of `(i, j) => (b_ij, b_ji)`
pairs (0-based, as Sage writes them) or the equivalent `Dict`; arcs missing from
it are simply laced.

Since the format records no frozen block, all vertices are mutable unless
`n_mutable` says otherwise. Symmetrizers are inferred from `B` unless given as
`d`; a `B` that admits none is rejected by the `Quiver` constructor.
"""
function from_dig6(dig6::AbstractString, edges = ();
                   n_mutable::Union{Nothing, Int} = nothing,
                   d::Union{Nothing, Vector{Int}} = nothing)
    s = String(dig6)
    startswith(s, _DIG6_HEADER) && (s = s[(length(_DIG6_HEADER) + 1):end])
    startswith(s, "&") ||
        throw(InvalidArgument("a digraph6 string must start with '&', got \"$s\""))
    n, packed = _dig6_decode_order(s[2:end])
    bits = _dig6_unpack_bits(packed, n * n)

    labels = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    for (key, value) in pairs(Dict(edges))
        labels[(Int(key[1]), Int(key[2]))] = (Int(value[1]), Int(value[2]))
    end

    B = zeros(Int, n, n)
    for i in 1:n, j in 1:n
        bits[(i - 1) * n + j] || continue
        a, b = get(labels, (i - 1, j - 1), (1, -1))
        B[i, j], B[j, i] = a, b
    end

    nm = n_mutable === nothing ? n : n_mutable
    dd = d === nothing ? _infer_symmetrizer(B, nm) : d
    return Quiver(B, nm, dd)
end

from_dig6(data::Tuple; kwargs...) = from_dig6(data[1], data[2]; kwargs...)

"""
    _infer_symmetrizer(B, n_mutable) → Vector{Int}

Positive integer `d` with `d[i]*B[i,j] == -d[j]*B[j,i]` on the mutable block,
normalized per connected component to be coprime. Falls back to all ones when
the constraints are inconsistent, leaving the diagnosis to the `Quiver`
constructor's `NotSkewSymmetrizable` check.
"""
function _infer_symmetrizer(B::Matrix{Int}, n_mutable::Int)
    n_mutable == 0 && return Int[]
    d = zeros(Rational{Int}, n_mutable)
    for root in 1:n_mutable
        d[root] == 0 || continue
        d[root] = 1
        queue = [root]
        component = [root]
        while !isempty(queue)
            i = popfirst!(queue)
            for j in 1:n_mutable
                (B[i, j] == 0 || d[j] != 0) && continue
                B[j, i] == 0 && return ones(Int, n_mutable)
                d[j] = -d[i] * B[i, j] // B[j, i]
                d[j] > 0 || return ones(Int, n_mutable)
                push!(queue, j)
                push!(component, j)
            end
        end
        scale = lcm([denominator(d[i]) for i in component])
        for i in component
            d[i] *= scale
        end
        divisor = gcd([numerator(d[i]) for i in component])
        for i in component
            d[i] //= divisor
        end
    end
    return [numerator(x) for x in d]
end

"""
    to_qmu(q::Quiver) → String

Serialize `q` in the plain-text `.qmu` format of Bernhard Keller's quiver
mutation applet. Frozen vertices are written the way Sage writes them: the
exchange matrix is completed to the square block matrix
`[B_mut  -B_froz'; B_froz  0]` and the frozen points carry a trailing `1` flag.

Vertex positions are laid out on a circle (the format stores coordinates; the
applet lets the user rearrange them). See [`write_qmu`](@ref) to save directly
to a file and [`from_qmu`](@ref) to read one back.
"""
function to_qmu(q::Quiver)::String
    n, m = q.n_mutable, nvertices(q)
    M = copy(q.B)
    if q.n_frozen > 0
        M[1:n, (n + 1):m] = -transpose(q.B[(n + 1):m, 1:n])
        M[(n + 1):m, (n + 1):m] .= 0
    end

    lines = ["//Number of points", string(m), "//Vertex radius", "9",
             "//Labels shown", "1", "//Matrix", "$m $m"]
    for i in 1:m
        push!(lines, join((string(M[i, j]) for j in 1:m), " "))
    end

    push!(lines, "//Points")
    for i in 1:m
        angle = 2 * pi * (i - 1) / max(m, 1)
        x = round(Int, 100 * cos(angle))
        y = round(Int, 100 * sin(angle))
        push!(lines, i <= n ? "9 $x $y" : "9 $x $y 1")
    end

    append!(lines, ["//Historycounter", "-1", "//History", "", "//Cluster is null"])
    return join(lines, "\n")
end

"""
    write_qmu(path::AbstractString, q::Quiver) → String

Write [`to_qmu(q)`](@ref) to `path`, appending the `.qmu` extension when it is
missing, and return the file name actually used.
"""
function write_qmu(path::AbstractString, q::Quiver)::String
    filename = endswith(path, ".qmu") ? String(path) : path * ".qmu"
    open(io -> write(io, to_qmu(q)), filename, "w")
    return filename
end

"""
    from_qmu(text::AbstractString; d = nothing) → Quiver

Parse the `.qmu` text produced by [`to_qmu`](@ref), Sage's `qmu_save`, or
Keller's applet. Points flagged frozen become the frozen block, which the format
requires to come last. Symmetrizers are inferred from the matrix unless passed
as `d`.
"""
function from_qmu(text::AbstractString; d::Union{Nothing, Vector{Int}} = nothing)
    lines = strip.(split(text, '\n'))
    sections = Dict{String, Vector{String}}()
    current = ""
    for line in lines
        if startswith(line, "//")
            current = line
            sections[current] = String[]
        elseif !isempty(current) && !isempty(line)
            push!(sections[current], String(line))
        end
    end

    haskey(sections, "//Matrix") ||
        throw(InvalidArgument("no //Matrix section in this .qmu text"))
    body = sections["//Matrix"]
    isempty(body) && throw(InvalidArgument("empty //Matrix section"))
    dims = parse.(Int, split(body[1]))
    length(dims) == 2 && dims[1] == dims[2] ||
        throw(InvalidArgument("//Matrix must declare a square size, got \"$(body[1])\""))
    m = dims[1]
    length(body) >= m + 1 ||
        throw(InvalidArgument("//Matrix declares $m rows but only $(length(body) - 1) follow"))
    B = zeros(Int, m, m)
    for i in 1:m
        row = parse.(Int, split(body[i + 1]))
        length(row) == m ||
            throw(InvalidArgument("row $i of //Matrix has $(length(row)) entries, expected $m"))
        B[i, :] = row
    end

    n_mutable = m
    if haskey(sections, "//Points")
        points = sections["//Points"]
        length(points) == m ||
            throw(InvalidArgument("//Points has $(length(points)) entries, expected $m"))
        frozen = [length(split(p)) >= 4 && split(p)[4] == "1" for p in points]
        n_mutable = something(findfirst(frozen), m + 1) - 1
        all(frozen[(n_mutable + 1):m]) ||
            throw(InvalidArgument("frozen points must come last in a .qmu file"))
    end

    dd = d === nothing ? _infer_symmetrizer(B, n_mutable) : d
    return Quiver(B, n_mutable, dd)
end

"""
    read_qmu(path::AbstractString; d = nothing) → Quiver

Read a `.qmu` file from `path` and parse it with [`from_qmu`](@ref).
"""
read_qmu(path::AbstractString; kwargs...) = from_qmu(read(path, String); kwargs...)
