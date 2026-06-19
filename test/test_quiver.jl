@testset "Quiver construction" begin
    B = [0 1; -1 0]
    q = Quiver(B)
    @test q.n_mutable == 2
    @test q.n_frozen  == 0
    @test q.d         == [1, 1]
    @test q.labels    == ["1", "2"]

    # Frozen vertex
    B3 = [0 1 0; -1 0 0; 0 0 0]
    q3 = Quiver(B3, 2)
    @test q3.n_mutable == 2
    @test q3.n_frozen  == 1

    # Skew-symmetrizable: B2 type, B = [0 1; -2 0], d = [2, 1]
    Bb2 = [0 1; -2 0]
    q_b2 = Quiver(Bb2, 2, [2, 1])
    @test q_b2.d == [2, 1]

    # Validation errors
    @test_throws NotSkewSymmetrizable Quiver([0 1; 1 0])               # not skew-symmetric
    @test_throws NotSkewSymmetrizable Quiver([0 1; -2 0], 2, [1, 1])  # wrong symmetrizer
    @test_throws InvalidArgument      Quiver([0 1; -1 0], 3)           # n_mutable > n_total
    @test_throws InvalidArgument      Quiver([0 1; -1 0], 1, [1], ["1", "2", "3"])  # wrong label count
end

@testset "Quiver mutation" begin
    B = [0 1; -1 0]
    q = Quiver(B)

    q1 = mutate(q, 1)
    @test q1.B == [0 -1; 1 0]

    q2 = mutate(q, 2)
    @test q2.B == [0 -1; 1 0]

    # Mutation is an involution
    @test mutate(q1, 1).B == q.B
    @test mutate(q2, 2).B == q.B

    # Labels and symmetrizer are preserved
    @test q1.labels    == q.labels
    @test q1.d         == q.d
    @test q1.n_mutable == q.n_mutable

    # Cannot mutate a frozen vertex
    B3 = [0 1 0; -1 0 0; 0 0 0]
    q3 = Quiver(B3, 2)
    @test_throws FrozenVertexMutation mutate(q3, 3)
    @test_throws InvalidVertex        mutate(q3, 0)
end

@testset "Quiver mutation sequence" begin
    B = [0 1; -1 0]
    q = Quiver(B)
    @test mutate(q, [1, 2, 1]).B == mutate(mutate(mutate(q, 1), 2), 1).B
end

@testset "Quiver accessors" begin
    B = [0 1 0; -1 0 0; 0 0 0]
    q = Quiver(B, 2)

    @test nvertices(q) == 3
    @test labels(q) == ["1", "2", "3"]

    @test is_frozen(q, 1) == false
    @test is_frozen(q, 2) == false
    @test is_frozen(q, 3) == true
    @test_throws InvalidVertex is_frozen(q, 0)
    @test_throws InvalidVertex is_frozen(q, 4)
end

@testset "Label-based Quiver mutation" begin
    B  = [0 1; -1 0]
    q  = Quiver(B, 2, [1, 1], ["u", "v"])
    q1 = mutate(q, "u")
    @test q1.B == mutate(q, 1).B

    @test_throws InvalidArgument mutate(q, "z")
end
