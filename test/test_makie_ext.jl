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
    @test_throws ArgumentError mutation_explorer(q)

    @eval using Graphs, GraphMakie, CairoMakie

    # quiver (existing) - graphplot returns a FigureAxisPlot with an .axis
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

    # ── the interactive explorer ─────────────────────────────────────────────
    # CairoMakie builds the figure but nothing responds, so what is tested here
    # is that every observable, widget and interaction wires up; the behaviour
    # behind them is tested backend-free in test_explorer.jl.
    @testset "mutation explorer" begin
        @test mutation_explorer(q) isa CairoMakie.Makie.Figure
        @test mutation_explorer(Seed(q)) isa CairoMakie.Makie.Figure
        @test mutation_explorer(extend(Seed(q))) isa CairoMakie.Makie.Figure
        @test mutation_explorer(extend_geometric(Seed(Quiver([0 1 0; -1 0 1; 0 -1 0], 2)))) isa CairoMakie.Makie.Figure
        @test mutation_explorer(Quiver(:A, 4); layout = "spring") isa CairoMakie.Makie.Figure
        @test mutation_explorer(q; layout = [(0.0, 0.0), (1.0, 0.0)]) isa CairoMakie.Makie.Figure
        @test mutation_explorer(extend(Seed(Quiver(:A, 3))); exchange_graph = true) isa
              CairoMakie.Makie.Figure
        # the class is enumerated only when it fits; beyond the cap it is skipped
        @test mutation_explorer(q; exchange_graph = true, max_class = 1) isa
              CairoMakie.Makie.Figure
        @test_throws ClusterAlgebras.InvalidArgument mutation_explorer(Quiver(zeros(Int, 1, 1), 0))
    end

    # Driving the widgets is how the state machine gets tested without a mouse:
    # the text box, the buttons and the slider all funnel into the same
    # mutate / undo / scrub code path a click would take.
    @testset "explorer state machine" begin
        M    = CairoMakie.Makie
        seen = Any[]
        s3   = extend(Seed(Quiver(:A, 3)))
        fig  = mutation_explorer(s3; on_mutate = x -> push!(seen, x))
        block(T) = filter(b -> b isa T, fig.content)
        tb     = only(block(M.Textbox))
        sl     = only(block(M.Slider))
        status = last(block(M.Label))
        button(lbl) = only(filter(b -> b.label[] == lbl, block(M.Button)))

        tb.stored_string[] = "1,3"
        @test length(seen) == 2
        @test seen[end].quiver == mutate(s3.quiver, [1, 3])
        @test sl.range[] == 1:3 && sl.value[] == 3
        @test occursin("path: 1, 3", status.text[])

        # undo rewinds one step and leaves the future in place for the slider
        button("undo (u)").clicks[] += 1
        @test sl.value[] == 2
        @test startswith(status.text[], "path: 1 ")       # the 3 is gone from the path
        M.set_close_to!(sl, 3)                       # scrub forward again
        @test occursin("path: 1, 3", status.text[])

        # mutating from an earlier state discards the states after it
        M.set_close_to!(sl, 1)
        tb.stored_string[] = "2"
        @test sl.range[] == 1:2 && sl.value[] == 2
        @test occursin("path: 2", status.text[])

        button("reset (r)").clicks[] += 1
        @test sl.range[] == 1:1 && occursin("(none)", status.text[])
        @test occursin("3/3 green", status.text[])

        # a bad sequence reports on the status line instead of throwing
        tb.stored_string[] = "1,9"
        @test occursin("not mutable", status.text[])
        @test sl.range[] == 1:1

        # the green-sequence button walks to an all-red seed
        button("green sequence").clicks[] += 1
        @test occursin("all red", status.text[])
        @test is_all_red(seen[end])

        # random walks never backtrack: mutating twice at one vertex is the identity
        button("reset (r)").clicks[] += 1
        for _ in 1:12
            button("random").clicks[] += 1
        end
        path = split(split(status.text[], "|")[1], r"[:,]")[2:end]
        @test length(path) == 12
        @test all(i -> strip(path[i]) != strip(path[i+1]), 1:11)

        # the keyboard path, driven straight through the figure's event stream
        button("reset (r)").clicks[] += 1
        M.events(fig).keyboardbutton[] = M.KeyEvent(M.Keyboard._2, M.Keyboard.press)
        @test occursin("path: 2", status.text[])
        M.events(fig).keyboardbutton[] = M.KeyEvent(M.Keyboard.u, M.Keyboard.press)
        @test occursin("(none)", status.text[])
        M.events(fig).keyboardbutton[] = M.KeyEvent(M.Keyboard._9, M.Keyboard.press)
        @test occursin("(none)", status.text[])      # rank 3 has no vertex 9

        # menus retarget the derived observables without rebuilding the figure
        menus = block(M.Menu)
        for m in menus, i in 1:length(m.options[])
            m.i_selected[] = i
            @test m.selection[] isa String
        end
        @test only(block(M.Toggle)).active[] == true
    end
end
