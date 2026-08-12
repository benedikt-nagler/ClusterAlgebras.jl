@testset "mutation explorer core" begin
    CA = ClusterAlgebras
    q  = Quiver([0 1; -1 0])                 # A₂: 1 → 2
    s  = Seed(q)
    ps = extend(s)                           # principal coefficients
    qf = Quiver([0 1 0; -1 0 1; 0 -1 0], 2)  # A₂ with one frozen vertex
    gs = extend_geometric(Seed(qf))          # geometric: frozen vertices appear

    @testset "a click is a mutation, or nothing at all" begin
        @test CA._explorer_step(q, 1) == mutate(q, 1)
        @test CA._explorer_step(CA._explorer_step(q, 1), 1) == q   # involution
        @test CA._explorer_step(q, 0)  === q                       # out of range
        @test CA._explorer_step(q, 3)  === q
        @test CA._explorer_step(gs, nvertices(gs.quiver)) === gs    # frozen vertex
        @test CA._explorer_step(ps, 1).quiver == mutate(q, 1)
    end

    @testset "modes offered" begin
        @test CA._explorer_color_modes(q)   == (:mutable_frozen, :source_sink)
        @test CA._explorer_color_modes(ps)  == (:green_red, :mutable_frozen, :source_sink)
        @test CA._explorer_label_modes(q)   == (:index, :name)
        @test CA._explorer_label_modes(s)   == (:index, :name, :variable)
        @test CA._explorer_label_modes(ps)  == (:index, :name, :variable, :cvector, :gvector)
        @test CA._explorer_matrix_kinds(q)  == (:B,)
        @test CA._explorer_matrix_kinds(ps) == (:C, :B, :G)
    end

    @testset "vertex kinds" begin
        @test CA._explorer_kinds(q, :mutable_frozen) == [:mutable, :mutable]
        @test CA._explorer_kinds(q, :source_sink)    == [:source, :sink]
        @test CA._explorer_kinds(Quiver(zeros(Int, 2, 2), 2), :source_sink) ==
              [:mutable, :mutable]                       # isolated is neither
        @test CA._explorer_kinds(gs, :mutable_frozen)[end] === :frozen

        # green/red tracks the c-vectors, and A₂ runs out of green after [2,1]
        @test CA._explorer_kinds(ps, :green_red) == [:green, :green]
        @test CA._explorer_kinds(mutate(ps, 1), :green_red) == [:red, :green]
        @test CA._explorer_kinds(mutate(ps, [2, 1]), :green_red) == [:red, :red]

        @test_throws CA.InvalidArgument CA._explorer_kinds(q, :green_red)
        @test_throws CA.InvalidArgument CA._explorer_kinds(q, :nonsense)
    end

    @testset "vertex labels" begin
        @test CA._explorer_labels(q, :index)  == ["1", "2"]
        @test CA._explorer_labels(q, :name)   == q.labels
        @test CA._explorer_labels(s, :variable) == ["x_1", "x_2"]
        @test CA._explorer_labels(ps, :cvector) == ["(1,0)", "(0,1)"]
        @test CA._explorer_labels(ps, :gvector) == ["(1,0)", "(0,1)"]
        # long cluster variables are truncated so a label cannot swamp the figure
        long = CA._explorer_labels(mutate(s, [1, 2, 1]), :variable; maxlen = 6)
        @test all(l -> length(l) <= 6, long)

        @test_throws CA.InvalidArgument CA._explorer_labels(q, :variable)
        @test_throws CA.InvalidArgument CA._explorer_labels(s, :cvector)
        @test_throws CA.InvalidArgument CA._explorer_labels(q, :nonsense)
    end

    @testset "matrix panel" begin
        @test CA._explorer_matrix(q, :B)[1] == q.B
        @test CA._explorer_matrix(ps, :C)[1] == cmatrix(ps)
        @test CA._explorer_matrix(ps, :G)[1] == gmatrix(ps)
        @test CA._explorer_matrix(mutate(ps, 1), :C)[1] == cmatrix(mutate(ps, 1))
        @test CA._explorer_matrix(q, :B)[2] isa String
        @test_throws CA.InvalidArgument CA._explorer_matrix(q, :C)
        @test_throws CA.InvalidArgument CA._explorer_matrix(q, :G)
        @test_throws CA.InvalidArgument CA._explorer_matrix(q, :nonsense)
    end

    @testset "edges: fixed topology, weights switch" begin
        edge_list = [(i, j) for i in 1:2 for j in 1:2 if i != j]
        w  = CA._explorer_edge_weights(q, edge_list)
        w1 = CA._explorer_edge_weights(mutate(q, 1), edge_list)
        @test length(w) == length(w1) == length(edge_list)          # topology fixed
        @test w  == [1, 0]                                          # 1 → 2
        @test w1 == [0, 1]                                          # mutation flips it
        @test CA._explorer_visible_edges(q) == [(1, 2)]
        @test CA._explorer_visible_edges(mutate(q, 1)) == [(2, 1)]

        # a multiple arrow is a weight, not a multi-edge
        kr = Quiver([0 2; -2 0])
        @test CA._explorer_edge_weights(kr, edge_list) == [2, 0]
    end

    @testset "sequence parsing" begin
        @test CA._parse_mutation_sequence("1,2,1", 2) == [1, 2, 1]
        @test CA._parse_mutation_sequence(" 1 2  1 ", 2) == [1, 2, 1]
        @test CA._parse_mutation_sequence("", 2) == Int[]
        @test_throws CA.InvalidArgument CA._parse_mutation_sequence("1,x", 2)
        @test_throws CA.InvalidArgument CA._parse_mutation_sequence("1,3", 2)
        @test_throws CA.InvalidArgument CA._parse_mutation_sequence("0", 2)
    end

    @testset "status line and tooltip" begin
        @test occursin("(none)", CA._explorer_status(q, Int[]))
        @test occursin("1, 2", CA._explorer_status(q, [1, 2]))
        @test occursin("2/2 green", CA._explorer_status(ps, Int[]))
        @test occursin("all red", CA._explorer_status(mutate(ps, [2, 1]), [2, 1]))
        @test !occursin("green", CA._explorer_status(q, Int[]))

        tip = CA._explorer_tooltip(ps, 1)
        @test occursin("vertex 1", tip) && occursin("mutable", tip)
        @test occursin("c = (1,0)", tip) && occursin("g = (1,0)", tip)
        @test occursin("green", tip)
        @test occursin("red", CA._explorer_tooltip(mutate(ps, 1), 1))
        @test occursin("frozen", CA._explorer_tooltip(gs, nvertices(gs.quiver)))
    end
end
