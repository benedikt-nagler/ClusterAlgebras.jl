# ═══════════════════════════════════════════════════════════════════════════
#  ClusterAlgebras.jl - visual showcase
#
#  Renders every plotting function in the package to PNG using the headless
#  CairoMakie backend.  For live, rotatable/zoomable figures swap CairoMakie for
#  GLMakie (`using GLMakie`) and replace the `save(...)` calls with `display(...)`.
#
#  Run (from the package root), using the test environment which already has the
#  visualization dependencies resolved:
#
#      julia --project=test examples/showcase_visuals.jl
#
#  or add the deps to your own environment first:
#
#      ] add CairoMakie GraphMakie Graphs
#      julia --project=. examples/showcase_visuals.jl
#
#  Figures are written to examples/figures/.
# ═══════════════════════════════════════════════════════════════════════════

using ClusterAlgebras
using Graphs, GraphMakie, CairoMakie

const OUT = joinpath(@__DIR__, "figures")
mkpath(OUT)
save_fig(name, obj) = (path = joinpath(OUT, name); CairoMakie.save(path, obj);
                       println("  wrote $path"))

println("ClusterAlgebras.jl visual showcase → $OUT")

# ── 1. Quivers ───────────────────────────────────────────────────────────────
# Mutable vertices: steelblue circles.  Frozen: lightgray rectangles.
println("\n[1] quivers (plot_quiver)")
q_A2 = Quiver([0 1; -1 0], 2, [1, 1], ["x", "y"])
q_B2 = Quiver([0 1; -2 0], 2, [2, 1], ["a", "b"])              # skew-symmetrizable
q_fr = Quiver([0 1 1; -1 0 1; -1 -1 0], 2, ones(Int, 2), ["x₁", "x₂", "c"])
save_fig("01a_quiver_A2.png", plot_quiver(q_A2))
save_fig("01b_quiver_B2.png", plot_quiver(q_B2))
save_fig("01c_quiver_frozen.png", plot_quiver(q_fr))

# ── 2. Exchange / mutation graph ───────────────────────────────────────────────
# Vertices = seeds, edges = single mutations.  A₂ → pentagon; A₃ → 1-skeleton of
# the 3D associahedron.
println("\n[2] exchange graph (plot_exchange_graph)")
save_fig("02a_exchange_A2_pentagon.png", plot_exchange_graph(exchange_graph(Seed(Quiver(:A, 2)))))
save_fig("02b_exchange_A3_associahedron.png", plot_exchange_graph(exchange_graph(Seed(Quiver(:A, 3)))))
save_fig("02c_mutation_class_A3.png", plot_exchange_graph(mutation_class(Quiver(:A, 3))))

# ── 3. Frieze patterns ─────────────────────────────────────────────────────────
# Classical staggered Conway–Coxeter diagram; borders greyed, interior black.
println("\n[3] friezes (plot_frieze)")
save_fig("03a_frieze_n7_fan.png", plot_frieze(frieze(7)))
save_fig("03b_frieze_quiddity.png", plot_frieze(frieze([3, 1, 3, 1, 3, 1])))

# ── 4. Green-sequence sign chart ───────────────────────────────────────────────
# Grid of c-vector signs along a maximal green sequence: top row all green,
# bottom all red - the "all green → all red" theorem.  Mutated vertex ringed.
println("\n[4] green sequences (plot_green_sequence)")
s_A3 = Seed(Quiver(:A, 3))
mgs = first(maximal_green_sequences(s_A3))
println("  A₃ maximal green sequence: $mgs")
save_fig("04_green_sequence_A3.png", plot_green_sequence(s_A3, mgs))

# ── 5. g-/d-vector fans ────────────────────────────────────────────────────────
# Every g-/d-vector across the exchange graph as a ray from the origin - the
# cluster fan.  Rank 2 in 2D, rank 3 in 3D.
println("\n[5] cluster fans (plot_g_vector_fan)")
save_fig("05a_gfan_A2.png", plot_g_vector_fan(Seed(Quiver(:A, 2)); kind = :g))
save_fig("05b_dfan_A2.png", plot_g_vector_fan(Seed(Quiver(:A, 2)); kind = :d))
save_fig("05c_gfan_B2.png", plot_g_vector_fan(Seed(Quiver([0 1; -2 0], 2, [2, 1])); kind = :g))
save_fig("05d_gfan_A3_3d.png", plot_g_vector_fan(Seed(Quiver(:A, 3)); kind = :g))

println("\nDone - $(length(readdir(OUT))) figures in $OUT")
