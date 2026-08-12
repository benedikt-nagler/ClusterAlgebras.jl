# ═══════════════════════════════════════════════════════════════════════════
#  ClusterAlgebras.jl - the interactive mutation explorer
#
#  Six guided demos of `mutation_explorer`: click a mutable vertex and the seed
#  mutates, with the quiver, the C-/B-/G-matrix and the position in the mutation
#  class all following along.
#
#  Interactivity needs GLMakie (or WGLMakie):
#
#      ] add GLMakie GraphMakie Graphs
#      julia --project=. examples/explore_mutations.jl        # all six in turn
#      julia --project=. examples/explore_mutations.jl 2      # just demo 2
#
#  Without GLMakie the script still runs: it falls back to CairoMakie, walks
#  each demo a few mutations by itself, and writes a still to examples/figures/.
#  The test environment already has CairoMakie, so this works out of the box:
#
#      julia --project=test examples/explore_mutations.jl
#
#  Set EXPLORER_STILLS=1 to force that path even where GLMakie is installed -
#  each demo then runs to completion instead of waiting on a window, which is
#  what you want from a script, from CI, or over ssh:
#
#      EXPLORER_STILLS=1 julia --project=test examples/explore_mutations.jl
# ═══════════════════════════════════════════════════════════════════════════

using ClusterAlgebras
using Graphs, GraphMakie

const HAS_GL = get(ENV, "EXPLORER_STILLS", "0") == "0" &&
               Base.find_package("GLMakie") !== nothing
@eval using $(HAS_GL ? :GLMakie : :CairoMakie)
const MK  = HAS_GL ? GLMakie.Makie : CairoMakie.Makie
const OUT = joinpath(@__DIR__, "figures")
mkpath(OUT)

"""
Show `fig` and block until the window is closed. With no interactive backend,
type `walk` into the explorer's own sequence box instead - the still then shows
a walk in progress rather than the initial seed - and save the result.
"""
function present(fig, name, instructions; walk = "")
    println(instructions)
    if HAS_GL
        println("  (close the window to continue)\n")
        wait(display(fig))
    else
        isempty(walk) || (only(filter(b -> b isa MK.Textbox, fig.content)).stored_string[] = walk)
        path = joinpath(OUT, name)
        save(path, fig)
        println("  still written to $path\n")
    end
    return fig
end

# ── 1 ────────────────────────────────────────────────────────────────────────
# The whole subject in two vertices: mutation is an involution, and the cluster
# variables are Laurent polynomials in the initial ones.
function demo_first_mutations()
    q = Quiver([0 1; -1 0], 2, [1, 1], ["x", "y"])
    present(mutation_explorer(Seed(q)), "explorer_01_first_mutations.png", """
    [1] First mutations - A₂

      Click vertex 1. The arrow reverses, and the cluster variable at 1 becomes
      (1 + y)/x. Click it again: you are back where you started, because μ_k is
      an involution.

      Set  labels: variable  to watch the cluster itself rather than the indices,
      and hover a vertex for its variable, c-vector and g-vector at once.
      A₂ is finite type, so you can only ever reach five distinct seeds - try to
      find them all with the history slider.
    """; walk = "1,2,1")
end

# ── 2 ────────────────────────────────────────────────────────────────────────
# Green sequences are the reason the explorer colours vertices at all.
function demo_green_sequences()
    s = extend(Seed(Quiver(:A, 4)))
    present(mutation_explorer(s), "explorer_02_green_sequences.png", """
    [2] Green and red - A₄ with principal coefficients

      Every vertex starts green (its c-vector is ≥ 0). Mutating at a green vertex
      keeps you on a green sequence; the status line counts how many are left.
      Colour is backed by marker shape - green circles, red diamonds - so the
      picture survives being printed in grey or read by a colour-blind viewer.

      Try mutating greedily at green vertices until none remain: that is a
      *maximal* green sequence, and the status line says so. The 'green sequence'
      button finds one for you and walks it. Then drag the history slider back
      and watch the c-vectors flip sign one at a time in the C panel.
    """; walk = "1,2,3,4")
end

# ── 3 ────────────────────────────────────────────────────────────────────────
# The minimap turns "which seed am I on?" into a picture.
function demo_exchange_graph()
    s = extend(Seed(Quiver(:A, 3)))
    present(mutation_explorer(s; exchange_graph = true), "explorer_03_exchange_graph.png", """
    [3] Where am I? - A₃ with the mutation-class minimap

      The lower right panel is the whole mutation class, with the current quiver
      in orange. Every click moves the orange dot along one edge, so the walk you
      are taking becomes a path on the exchange graph.

      This is opt-in (`exchange_graph = true`) because the class has to be
      enumerated first; it is skipped, with a note, past `max_class = 400`.
    """; walk = "2,1,3,2")
end

# ── 4 ────────────────────────────────────────────────────────────────────────
# Frozen vertices are drawn, hoverable, and never mutate - a click on one is
# deliberately a no-op rather than an error.
function demo_frozen_vertices()
    q = Quiver([0 1 1; -1 0 1; -1 -1 0], 2, ones(Int, 2), ["x₁", "x₂", "c"])
    present(mutation_explorer(extend_geometric(Seed(q))), "explorer_04_frozen.png", """
    [4] Coefficients - a quiver with a frozen vertex

      Grey squares are frozen: they carry the coefficients and are never mutated.
      Clicking one does nothing at all (a stray click is not an error).

      With geometric coefficients the green/red colouring is not available - that
      needs principal coefficients - so the colour menu offers mutable/frozen and
      source/sink instead. Source/sink is worth a look on any quiver: it marks the
      vertices where mutation is a sink-source reflection.
    """; walk = "1,2")
end

# ── 5 ────────────────────────────────────────────────────────────────────────
# Multiple arrows, and a mutation-infinite quiver.
function demo_multiple_arrows()
    q = Quiver([0 3 -3; -3 0 3; 3 -3 0])          # mutation-infinite
    present(mutation_explorer(q; layout = "shell"), "explorer_05_growth.png", """
    [5] Multiple arrows, and how fast they grow

      A multiple arrow is one edge labelled ×n; the 'arrow ×n' toggle hides the
      labels, and  matrix: B  shows the same thing as numbers.

      This quiver is mutation-infinite, and the growth is doubly exponential:
      3, 6, 15, 87, 1299, 112998, ... Nothing here would fit in a minimap, which
      is why the minimap is off by default.

      Contrast the Kronecker quiver `Quiver([0 2; -2 0])`: in rank 2 mutation only
      reverses every arrow (B ↦ −B), so a multiplicity there can never change.

      Caveat worth knowing: `B` holds `Int`, so around the tenth mutation the
      entries pass 2^63 and wrap silently. This is a quiver to look at, not to
      click twenty times.
    """; walk = "1,2,3,1")
end

# ── 6 ────────────────────────────────────────────────────────────────────────
# The explorer is not a dead end: `on_mutate` hands the state back.
function demo_capture_state()
    s = extend(Seed(Quiver(:D, 4)))
    reached = Ref{Any}(s)
    fig = mutation_explorer(s; on_mutate = x -> reached[] = x)
    present(fig, "explorer_06_capture.png", """
    [6] Taking the state with you - D₄

      This explorer was opened with

          reached = Ref{Any}(s)
          mutation_explorer(s; on_mutate = x -> reached[] = x)

      so whatever you explore your way to is a Julia object you can keep working
      with. Mutate a few times, then close the window.
    """; walk = "1,2,3,4,2")

    x = reached[]
    println("    you stopped at:  path = ", x.mutation_path)
    println("                     C    = ", cmatrix(x))
    println("                     green vertices = ",
            [k for k in 1:x.quiver.n_mutable if is_green(x, k)])
    println("                     all red? ", is_all_red(x), "\n")
end

const DEMOS = [demo_first_mutations, demo_green_sequences, demo_exchange_graph,
               demo_frozen_vertices, demo_multiple_arrows, demo_capture_state]

function main(args)
    println("\nClusterAlgebras.jl - mutation explorer",
            HAS_GL ? " (GLMakie: interactive)" : " (CairoMakie: stills only)", "\n")
    HAS_GL || println("""
    Not interactive: each demo walks itself and is saved to $OUT.
    For the real thing, `] add GLMakie` and run without EXPLORER_STILLS set.
    """)
    chosen = isempty(args) ? eachindex(DEMOS) : parse.(Int, args)
    for i in chosen
        1 <= i <= length(DEMOS) ||
            error("no demo $i; there are $(length(DEMOS))")
        DEMOS[i]()
    end
    println("Controls, in every demo: click a mutable vertex to mutate · drag to")
    println("rearrange · hover for the algebra · keys 1-9 mutate, u undoes, r resets ·")
    println("the history slider scrubs the whole walk · the box applies a sequence.")
end

main(ARGS)
