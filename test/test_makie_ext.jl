@testset "makie extension" begin
    q  = Quiver([0 1; -1 0])          # A₂
    f  = frieze(5)
    mc = mutation_class(q)
    eg = exchange_graph(Seed(q))

    # before the extensions load, the stubs explain what to load
    @test_throws ArgumentError plot_quiver(q)
    @test_throws ArgumentError plot_exchange_graph(mc)
    @test_throws ArgumentError plot_frieze(f)
    @test_throws ArgumentError plot_green_sequence(Seed(q), [1, 2])
    @test_throws ArgumentError plot_g_vector_fan(Seed(q))

    @eval using Graphs, GraphMakie, CairoMakie

    # quiver (existing) — graphplot returns a FigureAxisPlot with an .axis
    pq = plot_quiver(q)
    @test pq.axis isa CairoMakie.Makie.Axis

    # exchange / mutation graph
    @test plot_exchange_graph(mc).axis isa CairoMakie.Makie.Axis
    @test plot_exchange_graph(eg; nlabels = false).axis isa CairoMakie.Makie.Axis
    @test plot_exchange_graph(eg; nlabels = ["a", "b", "c", "d", "e"]).axis isa
          CairoMakie.Makie.Axis

    # frieze staggered diagram
    @test plot_frieze(f) isa CairoMakie.Makie.Figure
    @test plot_frieze(frieze([3, 1, 3, 1, 3, 1])) isa CairoMakie.Makie.Figure

    # green-sequence sign chart
    @test plot_green_sequence(Seed(q), [1, 2]) isa CairoMakie.Makie.Figure
    @test plot_green_sequence(extend(Seed(q)), [1, 2]) isa CairoMakie.Makie.Figure

    # g-/d-vector fan (rank 2 → 2D, rank 3 → 3D)
    @test plot_g_vector_fan(Seed(q)) isa CairoMakie.Makie.Figure
    @test plot_g_vector_fan(Seed(q); kind = :d) isa CairoMakie.Makie.Figure
    @test plot_g_vector_fan(Seed(Quiver(:A, 3))) isa CairoMakie.Makie.Figure
    @test plot_g_vector_fan(Seed(Quiver(:A, 3)); kind = :d) isa CairoMakie.Makie.Figure
    @test_throws ArgumentError plot_g_vector_fan(Seed(q); kind = :bad)
    @test_throws ClusterAlgebras.InvalidArgument plot_g_vector_fan(Seed(Quiver(:A, 4)))
end
