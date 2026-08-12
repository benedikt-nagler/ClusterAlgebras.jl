# Interactive mutation explorer - backend-free core.
#
# `mutation_explorer` itself lives in the GraphMakie extension; everything that
# can be decided without a plotting backend lives here, so the explorer's
# behaviour (what a click does, what a vertex is called, what is drawn in the
# matrix panel) is unit-tested headlessly and the extension is only wiring.

const EXPLORER_LABEL_MODES  = (:index, :name, :variable, :cvector, :gvector)
const EXPLORER_COLOR_MODES  = (:green_red, :mutable_frozen, :source_sink)
const EXPLORER_MATRIX_KINDS = (:C, :B, :G)

_explorer_quiver(q::Quiver)      = q
_explorer_quiver(s::AbstractSeed) = s.quiver

_has_principal(x) = x isa Seed{PrincipalCoefficients}

"""Modes offered by the explorer's colour menu for `x`, best first."""
_explorer_color_modes(x) =
    _has_principal(x) ? (:green_red, :mutable_frozen, :source_sink) :
                        (:mutable_frozen, :source_sink)

"""Modes offered by the explorer's vertex-label menu for `x`, best first."""
function _explorer_label_modes(x)
    modes = [:index, :name]
    x isa AbstractSeed && push!(modes, :variable)
    _has_principal(x) && append!(modes, [:cvector, :gvector])
    Tuple(modes)
end

"""Matrices the explorer's matrix panel can show for `x`, best first."""
_explorer_matrix_kinds(x) = _has_principal(x) ? (:C, :B, :G) : (:B,)

"""
One explorer click: mutate at vertex `k`, or return `x` unchanged when `k` is
frozen or out of range.  A stray click is not an error.
"""
function _explorer_step(x, k::Int)
    q = _explorer_quiver(x)
    1 <= k <= q.n_mutable ? mutate(x, k) : x
end

"""
Vertex kinds under colour `mode`, one symbol per vertex, drawn by the extension
as a (colour, marker) pair - both differ per kind, so the encoding survives
colour-blindness.
"""
function _explorer_kinds(x, mode::Symbol)
    q = _explorer_quiver(x)
    n = nvertices(q)
    if mode === :mutable_frozen
        return [i <= q.n_mutable ? :mutable : :frozen for i in 1:n]
    elseif mode === :green_red
        _has_principal(x) || throw(InvalidArgument(
            "colour mode :green_red needs principal coefficients; call extend(seed) first"))
        return [i <= q.n_mutable ? (is_green(x, i) ? :green : :red) : :frozen for i in 1:n]
    elseif mode === :source_sink
        kinds = Symbol[]
        for i in 1:n
            if i > q.n_mutable
                push!(kinds, :frozen)
                continue
            end
            row = @view q.B[i, :]
            push!(kinds, all(iszero, row) ? :mutable :
                         all(<=(0), row)  ? :sink    :
                         all(>=(0), row)  ? :source  : :mutable)
        end
        return kinds
    end
    throw(InvalidArgument("unknown colour mode :$mode; expected one of $EXPLORER_COLOR_MODES"))
end

_vector_string(v) = "(" * join(v, ",") * ")"

function _truncate_string(s::AbstractString, maxlen::Int)
    length(s) <= maxlen ? String(s) : String(first(s, max(maxlen - 1, 1))) * "…"
end

"""Vertex labels under label `mode`, one string per vertex."""
function _explorer_labels(x, mode::Symbol; maxlen::Int = 22)
    q = _explorer_quiver(x)
    n = nvertices(q)
    if mode === :index
        return string.(1:n)
    elseif mode === :name
        return copy(q.labels)
    elseif mode === :variable
        x isa AbstractSeed || throw(InvalidArgument("label mode :variable needs a Seed"))
        return [_truncate_string(string(x.cluster[i]), maxlen) for i in 1:n]
    elseif mode === :cvector
        _has_principal(x) || throw(InvalidArgument(
            "label mode :cvector needs principal coefficients; call extend(seed) first"))
        return [i <= q.n_mutable ? _vector_string(c_vector(x, i)) : "" for i in 1:n]
    elseif mode === :gvector
        _has_principal(x) || throw(InvalidArgument(
            "label mode :gvector needs principal coefficients; call extend(seed) first"))
        G = gmatrix(x)
        return [i <= q.n_mutable ? _vector_string(G[:, i]) : "" for i in 1:n]
    end
    throw(InvalidArgument("unknown label mode :$mode; expected one of $EXPLORER_LABEL_MODES"))
end

"""The matrix panel's content: `(matrix, title)` for matrix kind `kind`."""
function _explorer_matrix(x, kind::Symbol)
    if kind === :B
        return Matrix{Int}(_explorer_quiver(x).B), "B (exchange matrix)"
    elseif kind === :C
        _has_principal(x) || throw(InvalidArgument(
            "matrix :C needs principal coefficients; call extend(seed) first"))
        return Matrix{Int}(cmatrix(x)), "C (c-vectors in columns)"
    elseif kind === :G
        _has_principal(x) || throw(InvalidArgument(
            "matrix :G needs principal coefficients; call extend(seed) first"))
        return Matrix{Int}(gmatrix(x)), "G (g-vectors in columns)"
    end
    throw(InvalidArgument("unknown matrix kind :$kind; expected one of $EXPLORER_MATRIX_KINDS"))
end

"""
Arrow multiplicities along `edge_list`, a fixed list of ordered vertex pairs.
Zero means "no arrow", which the extension draws as a zero-width invisible edge:
the drawn graph keeps a fixed topology so every mutation is an attribute update
rather than a re-plot, and vertices never jump.
"""
_explorer_edge_weights(x, edge_list) =
    [max(_explorer_quiver(x).B[i, j], 0) for (i, j) in edge_list]

"""Ordered pairs of the arrows actually present in `x`, for layout algorithms."""
function _explorer_visible_edges(x)
    q = _explorer_quiver(x)
    n = nvertices(q)
    pairs = Tuple{Int,Int}[]
    for i in 1:n, j in 1:n
        q.B[i, j] > 0 && push!(pairs, (i, j))
    end
    pairs
end

"""Hover text for vertex `i`: everything known about it, one fact per line."""
function _explorer_tooltip(x, i::Int; maxlen::Int = 60)
    q = _explorer_quiver(x)
    lines = ["vertex $i" * (q.labels[i] == string(i) ? "" : " ($(q.labels[i]))")]
    push!(lines, i <= q.n_mutable ? "mutable" : "frozen")
    if x isa AbstractSeed
        push!(lines, "x = " * _truncate_string(string(x.cluster[i]), maxlen))
    end
    if _has_principal(x) && i <= q.n_mutable
        push!(lines, "c = " * _vector_string(c_vector(x, i)) *
                     (is_green(x, i) ? "  (green)" : "  (red)"))
        push!(lines, "g = " * _vector_string(g_vector(x, i)))
    end
    join(lines, "\n")
end

"""
Parse a user-typed mutation sequence such as `"1,3,2"` or `"1 3 2"`.
Throws [`InvalidArgument`](@ref) on anything that is not a mutable vertex.
"""
function _parse_mutation_sequence(str::AbstractString, n_mutable::Int)
    ks = Int[]
    for token in split(str, r"[,\s]+"; keepempty = false)
        k = tryparse(Int, token)
        k === nothing && throw(InvalidArgument("not a vertex index: \"$token\""))
        1 <= k <= n_mutable ||
            throw(InvalidArgument("vertex $k is not mutable (expected 1:$n_mutable)"))
        push!(ks, k)
    end
    ks
end

"""The explorer's status line: where we are, and whether the walk has finished."""
function _explorer_status(x, path)
    q = _explorer_quiver(x)
    parts = ["path: " * (isempty(path) ? "(none)" : join(path, ", "))]
    if _has_principal(x)
        n_green = count(k -> is_green(x, k), 1:q.n_mutable)
        push!(parts, "$n_green/$(q.n_mutable) green")
        is_all_red(x) && push!(parts, "all red - maximal green sequence complete")
    end
    join(parts, "    |    ")
end
