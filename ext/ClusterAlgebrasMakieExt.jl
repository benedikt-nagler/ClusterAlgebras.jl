# Makie weakdep extension: classical staggered rendering of a Conway–Coxeter
# frieze pattern.  Makie-only (no Graphs/GraphMakie needed).
module ClusterAlgebrasMakieExt

using ClusterAlgebras, Makie

"""
    plot_frieze(f::Frieze; interior_color = :black, border_color = :gray,
                fontsize = 20, kwargs...)

Draw the frieze `f` as the classical staggered integer diagram: each successive
row is offset by half a step, so the SL₂ unimodular diamonds line up visually.
Border rows (the bounding 0s and 1s) are drawn in `border_color`, the interior
entries in `interior_color`.

Returns a `Makie.Figure`. Requires a Makie backend (`using CairoMakie` / GLMakie).
"""
function ClusterAlgebras.plot_frieze(f::Frieze;
                                     interior_color = :black,
                                     border_color = :gray,
                                     fontsize::Real = 20,
                                     kwargs...)
    n = f.n
    F = f.entries
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; title = "Frieze (n = $n), quiddity $(f.quiddity)")
    for r in 1:n+1                       # matrix row r = math row r-1
        is_border = r <= 2 || r >= n     # rows of 0s and 1s at top/bottom
        col = is_border ? border_color : interior_color
        y = -(r - 1)
        for i in 1:n
            x = i + (r - 1) / 2          # half-step stagger
            Makie.text!(ax, x, y; text = string(F[r, i]),
                        align = (:center, :center), color = col, fontsize)
        end
    end
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)
    ax.aspect = Makie.DataAspect()
    fig
end

# ─── Green-sequence sign chart ────────────────────────────────────────────────

_sign_color(sym, green, red, mixed) =
    sym === :green ? green : sym === :red ? red : mixed

"""
    plot_green_sequence(s::Seed, seq; green = :seagreen, red = :crimson,
                        mixed = :lightgray)

Chart the c-vector signs along the mutation sequence `seq` (see
[`green_sequence_signs`](@ref)) as a grid: one column per mutable vertex, one row
per step (top = start, bottom = final state).  Each cell is coloured green
(c-vector ≥ 0), red (≤ 0), or grey (mixed); the vertex mutated at each step is
ringed.  For a maximal green sequence the top row is all green and the bottom all
red - the "all green → all red" theorem made visible.

Returns a `Makie.Figure`. Requires a Makie backend (`using CairoMakie` / GLMakie).
"""
function ClusterAlgebras.plot_green_sequence(s::ClusterAlgebras.Seed, seq;
                                             green = :seagreen,
                                             red = :crimson,
                                             mixed = :lightgray,
                                             kwargs...)
    signs = ClusterAlgebras.green_sequence_signs(s, seq)
    nsteps, n = size(signs)
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; xlabel = "mutable vertex", ylabel = "step",
                    title = "Green sequence $(collect(seq))")
    xs = Float64[]; ys = Float64[]; cs = Any[]
    for r in 1:nsteps, k in 1:n
        push!(xs, k); push!(ys, -(r - 1))
        push!(cs, _sign_color(signs[r, k], green, red, mixed))
    end
    Makie.scatter!(ax, xs, ys; color = cs, marker = :rect, markersize = 34)
    # ring the vertex mutated at each step (row r acts before mutation r)
    for (t, j) in enumerate(seq)
        Makie.scatter!(ax, [Float64(j)], [-(t - 1)];
                       marker = :rect, markersize = 40,
                       color = (:black, 0.0), strokecolor = :black, strokewidth = 2)
    end
    ax.yticks = (Float64.(-(0:nsteps-1)), string.(0:nsteps-1))
    ax.xticks = (Float64.(1:n), string.(1:n))
    fig
end

# ─── g-/d-vector fan (rank 2) ─────────────────────────────────────────────────

# All distinct g-/d-vectors across the exchange graph, as Vector{Int} of length n.
function _fan_vectors(s, kind, max_seeds)
    n = s.quiver.n_mutable
    seed0 = (kind === :g && !(s isa ClusterAlgebras.Seed{ClusterAlgebras.PrincipalCoefficients})) ?
            ClusterAlgebras.extend(s) : s
    eg = ClusterAlgebras.exchange_graph(seed0; max_seeds)
    vecs = Set{Vector{Int}}()
    for i in 1:length(eg)
        seed = eg[i]
        if kind === :g
            for g in ClusterAlgebras.gvectors(seed)
                push!(vecs, Int.(g))
            end
        else
            for k in 1:n
                push!(vecs, Int.(ClusterAlgebras.denominator_vector(seed, k)))
            end
        end
    end
    collect(vecs)
end

"""
    plot_g_vector_fan(s::Seed; kind = :g, max_seeds = 1000, kwargs...)

Draw the **g-vector fan** (`kind = :g`) or **d-vector fan** (`kind = :d`): every
g-/d-vector of every cluster variable across the exchange graph of `s`, plotted
as a ray from the origin.  For a finite-type seed this is the complete cluster
fan; for infinite type it is the portion reached within `max_seeds` seeds
(necessarily partial).

Rank 2 renders in a 2D axis, rank 3 in a 3D axis (`Axis3`); higher ranks are not
supported (no faithful low-dim picture).  `kind = :g` auto-extends a
`TrivialCoefficients` seed to principal coefficients (where g-vectors live).
Returns a `Makie.Figure`.
"""
function ClusterAlgebras.plot_g_vector_fan(s::ClusterAlgebras.Seed;
                                           kind::Symbol = :g,
                                           max_seeds::Int = 1000,
                                           kwargs...)
    kind in (:g, :d) ||
        throw(ArgumentError("plot_g_vector_fan: kind must be :g or :d, got :$kind"))
    n = s.quiver.n_mutable
    n in (2, 3) ||
        throw(ClusterAlgebras.InvalidArgument(
            "plot_g_vector_fan supports rank 2 or 3; got rank $n"))
    vecs = _fan_vectors(s, kind, max_seeds)
    title = "$(kind == :g ? "g" : "d")-vector fan (rank $n)"
    fig = Makie.Figure()
    if n == 2
        ax = Makie.Axis(fig[1, 1]; xlabel = "e₁", ylabel = "e₂", title)
        for v in vecs
            Makie.lines!(ax, [0.0, Float64(v[1])], [0.0, Float64(v[2])];
                         color = :steelblue)
        end
        isempty(vecs) ||
            Makie.scatter!(ax, [Float64(v[1]) for v in vecs],
                           [Float64(v[2]) for v in vecs]; color = :crimson)
        ax.aspect = Makie.DataAspect()
    else
        ax = Makie.Axis3(fig[1, 1]; xlabel = "e₁", ylabel = "e₂", zlabel = "e₃",
                         title, aspect = :data)
        for v in vecs
            Makie.lines!(ax, [0.0, Float64(v[1])], [0.0, Float64(v[2])],
                         [0.0, Float64(v[3])]; color = :steelblue)
        end
        isempty(vecs) ||
            Makie.scatter!(ax, [Float64(v[1]) for v in vecs],
                           [Float64(v[2]) for v in vecs],
                           [Float64(v[3]) for v in vecs]; color = :crimson)
    end
    fig
end

end # module ClusterAlgebrasMakieExt
