# Visualization example for ClusterAlgebras.jl
#
# First add the optional packages to your environment:
#   ] add GLMakie GraphMakie Graphs
#
# Then run with:
#   julia --project=. examples/visualize.jl

using ClusterAlgebras
using GLMakie, GraphMakie, Graphs

# ─────────────────────────────────────────────────────────────────────────────
# 1. to_dot  --  zero-dependency DOT export
#    Paste output at graphviz.org, or: echo "..." | dot -Tpng -o quiver.png
# ─────────────────────────────────────────────────────────────────────────────

q_A2 = Quiver([0 1; -1 0], 2, [1, 1], ["x", "y"])
q_A3 = Quiver([0 1 0; -1 0 1; 0 -1 0])
q_B2 = Quiver([0 1; -2 0], 2, [2, 1], ["a", "b"])
q_fr = Quiver([0 1 1; -1 0 1; -1 -1 0], 2, ones(Int, 2), ["x₁", "x₂", "c"])

println("── A₂ DOT ──")
println(to_dot(q_A2))

println("── B₂ DOT ──")
println(to_dot(q_B2))

println("── Quiver with frozen vertex c ──")
println(to_dot(q_fr))

# ─────────────────────────────────────────────────────────────────────────────
# 2. plot_quiver  --  standalone interactive figure
#    Mutable vertices: steelblue circles    Frozen: lightgray rectangles
#    The extension loads automatically when GLMakie/GraphMakie/Graphs are loaded.
# ─────────────────────────────────────────────────────────────────────────────

f1, ax1, _ = plot_quiver(q_A2)
ax1.title = "A₂ quiver"
display(f1)

f2, ax2, _ = plot_quiver(q_B2)
ax2.title = "B₂ quiver (skew-symmetrizable)"
display(f2)

f3, ax3, _ = plot_quiver(q_fr)
ax3.title = "Quiver with frozen vertex c"
display(f3)

# ─────────────────────────────────────────────────────────────────────────────
# 3. plot_quiver!  --  compose into an existing figure (multi-panel layout)
#    A₃ before and after mutation at vertex 2, side by side.
# ─────────────────────────────────────────────────────────────────────────────

fig = Figure(size = (900, 400))
plot_quiver!(Axis(fig[1, 1], title = "A₃  (before μ₂)"), q_A3)
plot_quiver!(Axis(fig[1, 2], title = "A₃  (after μ₂)"),  mutate(q_A3, 2))
display(fig)

println("\nClose the windows to exit.")
wait(display(fig))
