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

# ─── Interactive mutation explorer ────────────────────────────────────────────
#
# Vertex kinds are encoded by colour *and* marker shape, never colour alone; the
# palette is Okabe-Ito.  The drawn graph is the complete digraph on the vertex
# set with absent arrows given zero width and zero alpha, so a mutation is a pure
# attribute update: vertices never jump, and the registered interactions survive.

const _EXPLORER_STYLE = (
    mutable = ("#0072B2", :circle),
    frozen  = ("#BBBBBB", :rect),
    green   = ("#009E73", :circle),
    red     = ("#D55E00", :diamond),
    source  = ("#E69F00", :utriangle),
    sink    = ("#56B4E9", :dtriangle),
)

_explorer_colors(kinds)  = [Makie.to_color(getproperty(_EXPLORER_STYLE, k)[1]) for k in kinds]
_explorer_markers(kinds) = [getproperty(_EXPLORER_STYLE, k)[2] for k in kinds]

const _EXPLORER_LAYOUTS = ["shell", "spring", "stress", "grid"]

function _explorer_positions(name::AbstractString, n::Int, present)
    NL = GraphMakie.NetworkLayout
    g = SimpleDiGraph(n)
    for (i, j) in present
        add_edge!(g, i, j)
    end
    alg = name == "spring" ? NL.Spring() :
          name == "stress" ? NL.Stress() :
          name == "grid"   ? NL.SquareGrid() : NL.Shell()
    Makie.Point2f.(alg(g))
end

const _EXPLORER_DIGITS = (Makie.Keyboard._1, Makie.Keyboard._2, Makie.Keyboard._3,
                          Makie.Keyboard._4, Makie.Keyboard._5, Makie.Keyboard._6,
                          Makie.Keyboard._7, Makie.Keyboard._8, Makie.Keyboard._9)

"""
    mutation_explorer(x; kwargs...) → Figure

Interactive explorer for a `Quiver` or `Seed`: **click a mutable vertex to
mutate there**, and watch the quiver, the C-/B-/G-matrix and (optionally) the
position in the mutation class update together.

Requires `using Graphs, GraphMakie` and an interactive backend (GLMakie or
WGLMakie); under CairoMakie the figure is built but nothing responds.

# Interaction
- **click** a mutable vertex to mutate; **drag** a vertex to rearrange; **hover**
  for that vertex's cluster variable, c-vector and g-vector.
- keys `1`-`9` mutate by index, `u` undoes, `r` resets.
- the **history slider** scrubs back and forth through the whole walk; mutating
  from an earlier point discards the states after it.
- the text box applies a whole sequence, e.g. `1,3,2`.
- menus choose the vertex colouring (green/red needs principal coefficients),
  the vertex labels, the matrix shown, and the layout algorithm.

# Keywords
- `layout = "shell"`: one of `$(_EXPLORER_LAYOUTS)`, or a vector of positions.
- `exchange_graph = false`: draw a mutation-class minimap with the current quiver
  highlighted. Off by default because the class must be enumerated first; it is
  skipped, with a note, when the class exceeds `max_class = 400` quivers.
- `on_mutate = nothing`: called with the new object after every mutation, which
  is how the explored state is captured (`res = Ref(s); mutation_explorer(s;
  on_mutate = x -> res[] = x)`).
- `node_size`, `figure_size`; anything else is forwarded to `graphplot`.
"""
function ClusterAlgebras.mutation_explorer(
        x0::Union{ClusterAlgebras.Quiver, ClusterAlgebras.AbstractSeed};
        layout = "shell",
        exchange_graph::Bool = false,
        max_class::Int = 400,
        on_mutate = nothing,
        node_size = 32,
        figure_size = (1240, 820),
        kwargs...)

    CA = ClusterAlgebras
    q0 = CA._explorer_quiver(x0)
    n, nm = nvertices(q0), q0.n_mutable
    nm >= 1 || throw(CA.InvalidArgument(
        "mutation_explorer needs at least one mutable vertex"))

    # ── state: the walk, and where we are on it ──────────────────────────────
    states  = Any[x0]
    pathvec = Int[]
    pos     = Ref(1)
    current = Observable{Any}(x0)

    # ── display modes ────────────────────────────────────────────────────────
    color_modes  = CA._explorer_color_modes(x0)
    label_modes  = CA._explorer_label_modes(x0)
    matrix_kinds = CA._explorer_matrix_kinds(x0)
    color_mode   = Observable(first(color_modes))
    label_mode   = Observable(:index)
    matrix_kind  = Observable(first(matrix_kinds))
    show_mult    = Observable(true)

    # ── the drawn graph: fixed topology, arrows switched by width/alpha ──────
    g = SimpleDiGraph(n)
    for i in 1:n, j in 1:n
        i != j && add_edge!(g, i, j)
    end
    edge_list = [(src(e), dst(e)) for e in edges(g)]
    positions = Observable(layout isa AbstractString ?
        _explorer_positions(layout, n, CA._explorer_visible_edges(x0)) :
        Makie.Point2f.(layout))

    kinds       = lift((x, m) -> CA._explorer_kinds(x, m), current, color_mode)
    node_color  = lift(_explorer_colors, kinds)
    node_marker = lift(_explorer_markers, kinds)
    nlabels     = lift((x, m) -> CA._explorer_labels(x, m), current, label_mode)
    weights     = lift(x -> CA._explorer_edge_weights(x, edge_list), current)
    edge_width  = lift(w -> Float32[wi == 0 ? 0 : 2.5 for wi in w], weights)
    edge_color  = lift(w -> [Makie.RGBAf(0.22, 0.22, 0.28, wi == 0 ? 0 : 1) for wi in w], weights)
    arrow_size  = lift(w -> Float32[wi == 0 ? 0 : 20 for wi in w], weights)
    elabels     = lift((w, s) -> [(s && wi > 1) ? string(wi) : "" for wi in w], weights, show_mult)

    fig = Figure(size = figure_size)
    ax  = Axis(fig[1, 1]; title = "click a mutable vertex to mutate")
    hidedecorations!(ax)
    hidespines!(ax)
    node_sizes = Observable(fill(Float32(node_size), n))   # a vector: NodeHoverHighlight needs one
    p = graphplot!(ax, g;
        layout = positions,
        node_color, node_marker, node_size = node_sizes,
        nlabels, nlabels_distance = 12, nlabels_fontsize = 15,
        edge_width, edge_color, arrow_size, arrow_shift = 0.55,
        elabels, elabels_fontsize = 14,
        force_straight_edges = true,
        node_attr = (; inspector_label = (_, i, _) -> CA._explorer_tooltip(current[], i)),
        kwargs...)

    # ── matrix panel: padded to a fixed n×n grid so no update can race ───────
    right   = GridLayout(fig[1, 2])
    matdata = lift((x, k) -> CA._explorer_matrix(x, k), current, matrix_kind)
    mat_ax  = Axis(right[1, 1]; title = lift(d -> d[2], matdata),
                   yreversed = true, aspect = Makie.DataAspect())
    hidedecorations!(mat_ax)
    hidespines!(mat_ax)
    padded = lift(matdata) do d
        M = d[1]
        Z = fill(NaN32, n, n)
        Z[1:size(M, 1), 1:size(M, 2)] .= Float32.(M)
        Z
    end
    heatmap!(mat_ax, 1:n, 1:n, lift(permutedims, padded);
             colormap = :RdBu,
             colorrange = lift(Z -> begin
                 m = maximum(abs, Iterators.filter(!isnan, Z); init = 1.0f0)
                 (-m, m)
             end, padded),
             nan_color = :transparent)
    text!(mat_ax, [Makie.Point2f(c, r) for r in 1:n for c in 1:n];
          text = lift(Z -> [isnan(Z[r, c]) ? "" : string(round(Int, Z[r, c]))
                            for r in 1:n for c in 1:n], padded),
          align = (:center, :center), fontsize = 13)

    # ── optional mutation-class minimap ──────────────────────────────────────
    class     = nothing
    class_idx = Int[1]
    if exchange_graph
        mc = mutation_class(q0; max_quivers = max_class)
        i0 = findfirst(==(q0), mc.quivers)
        if is_truncated(mc) || i0 === nothing
            Label(right[2, 1], "mutation class larger than $max_class - minimap skipped";
                  tellwidth = false)
        else
            class_idx[1] = i0
            highlight = Observable(i0)
            mm_ax = Axis(right[2, 1]; title = "mutation class ($(length(mc)) quivers)")
            hidedecorations!(mm_ax)
            hidespines!(mm_ax)
            graphplot!(mm_ax, SimpleGraph(mc);
                node_size = 11,
                node_color = lift(i -> [j == i ? Makie.to_color("#D55E00") :
                                        Makie.to_color("#BBBBBB") for j in 1:length(mc)],
                                  highlight))
            class = (mc = mc, highlight = highlight)
        end
    end

    # ── controls ─────────────────────────────────────────────────────────────
    ctrl = GridLayout(fig[2, 1:2]; tellheight = true)
    col  = 0
    next!() = (col += 1)
    b_undo  = Button(ctrl[1, next!()]; label = "undo (u)")
    b_reset = Button(ctrl[1, next!()]; label = "reset (r)")
    b_rand  = Button(ctrl[1, next!()]; label = "random")
    b_mgs   = CA._has_principal(x0) ? Button(ctrl[1, next!()]; label = "green sequence") : nothing
    Label(ctrl[1, next!()], "colour:"; halign = :right)
    m_color = Menu(ctrl[1, next!()]; options = string.(collect(color_modes)),
                   default = string(color_mode[]), width = 145)
    Label(ctrl[1, next!()], "labels:"; halign = :right)
    m_label = Menu(ctrl[1, next!()]; options = string.(collect(label_modes)),
                   default = string(label_mode[]), width = 110)
    Label(ctrl[1, next!()], "matrix:"; halign = :right)
    m_matrix = Menu(ctrl[1, next!()]; options = string.(collect(matrix_kinds)),
                    default = string(matrix_kind[]), width = 70)
    Label(ctrl[1, next!()], "layout:"; halign = :right)
    m_layout = Menu(ctrl[1, next!()]; options = _EXPLORER_LAYOUTS,
                    default = layout isa AbstractString ? layout : "shell", width = 100)
    Label(ctrl[1, next!()], "arrow ×n:"; halign = :right)
    t_mult = Toggle(ctrl[1, next!()]; active = true)

    row2 = GridLayout(fig[3, 1:2])
    Label(row2[1, 1], "sequence:"; halign = :right)
    tb = Textbox(row2[1, 2]; placeholder = "e.g. 1,3,2 then press enter", width = Makie.Auto())
    Label(row2[1, 3], "history:"; halign = :right)
    sl = Slider(row2[1, 4]; range = 1:1, startvalue = 1)
    colsize!(row2, 2, Relative(0.3))
    colsize!(row2, 4, Relative(0.45))

    status = Label(fig[4, 1:2], CA._explorer_status(x0, Int[]);
                   halign = :left, tellwidth = false)
    colsize!(fig.layout, 1, Relative(0.6))

    # ── the state machine ────────────────────────────────────────────────────
    syncing = Ref(false)

    function goto!(i::Int)
        i = clamp(i, 1, length(states))
        pos[] = i
        current[] = states[i]
        status.text[] = CA._explorer_status(states[i], view(pathvec, 1:i-1))
        class === nothing || (class.highlight[] = class_idx[i])
        syncing[] = true
        sl.range[] = 1:length(states)
        Makie.set_close_to!(sl, i)
        syncing[] = false
        return i
    end

    function reset!()
        resize!(states, 1)
        resize!(pathvec, 0)
        resize!(class_idx, 1)
        goto!(1)
    end

    function do_mutate!(k::Int)
        i = pos[]
        x = states[i]
        1 <= k <= CA._explorer_quiver(x).n_mutable || return
        resize!(states, i)                    # a mutation from the past drops the future
        resize!(pathvec, i - 1)
        resize!(class_idx, min(i, length(class_idx)))
        xn = CA._explorer_step(x, k)
        push!(states, xn)
        push!(pathvec, k)
        class === nothing || push!(class_idx, class.mc.adj[class_idx[i]][k])
        goto!(i + 1)
        on_mutate === nothing || on_mutate(xn)
        return
    end

    apply_sequence!(ks) = foreach(do_mutate!, ks)

    on(sl.value) do i
        syncing[] || goto!(i)
    end
    on(b_undo.clicks)  do _; pos[] > 1 && goto!(pos[] - 1); end
    on(b_reset.clicks) do _; reset!(); end
    on(b_rand.clicks) do _
        # non-backtracking, as `random_mutate` is: repeating a vertex is the
        # identity, and an explorer that appears to do nothing reads as broken
        nmut = CA._explorer_quiver(states[pos[]]).n_mutable
        last = pos[] > 1 ? pathvec[pos[]-1] : 0
        k = rand(1:nmut)
        while nmut > 1 && k == last
            k = rand(1:nmut)
        end
        do_mutate!(k)
    end
    if b_mgs !== nothing
        on(b_mgs.clicks) do _
            seqs = maximal_green_sequences(states[pos[]]; max_length = 4nm, max_count = 1)
            isempty(seqs) ? (status.text[] = "no maximal green sequence of length ≤ $(4nm)") :
                            apply_sequence!(first(seqs))
        end
    end
    on(m_color.selection)  do s; s === nothing || (color_mode[]  = Symbol(s)); end
    on(m_label.selection)  do s; s === nothing || (label_mode[]  = Symbol(s)); end
    on(m_matrix.selection) do s; s === nothing || (matrix_kind[] = Symbol(s)); end
    on(m_layout.selection) do s
        s === nothing && return
        positions[] = _explorer_positions(s, n, CA._explorer_visible_edges(states[pos[]]))
        autolimits!(ax)
    end
    on(t_mult.active) do a; show_mult[] = a; end
    on(tb.stored_string) do s
        isnothing(s) && return
        try
            apply_sequence!(CA._parse_mutation_sequence(s, CA._explorer_quiver(states[pos[]]).n_mutable))
        catch err
            err isa CA.ClusterAlgebraError || rethrow()
            status.text[] = sprint(showerror, err)
        end
    end

    # ── pointer and keyboard ─────────────────────────────────────────────────
    deregister_interaction!(ax, :rectanglezoom)
    register_interaction!(ax, :nodeclick, NodeClickHandler((idx, _, _) -> (do_mutate!(idx); true)))
    register_interaction!(ax, :nodedrag, NodeDragHandler((_, idx, event, _) -> begin
        positions[][idx] = Makie.Point2f(event.data)
        notify(positions)
        true
    end))
    register_interaction!(ax, :nodehover, NodeHoverHighlight(p))

    on(events(fig).keyboardbutton) do e
        (e.action == Makie.Keyboard.press && !tb.focused[]) || return Consume(false)
        k = findfirst(==(e.key), _EXPLORER_DIGITS)
        if k !== nothing
            do_mutate!(k)
        elseif e.key == Makie.Keyboard.u
            pos[] > 1 && goto!(pos[] - 1)
        elseif e.key == Makie.Keyboard.r
            reset!()
        end
        return Consume(false)
    end

    DataInspector(fig)
    return fig
end

end
