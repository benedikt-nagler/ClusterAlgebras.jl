@testset "Named constructors" begin

    # ─── Symbol constructors ─────────────────────────────────────────────────

    @testset "A_n" begin
        q = Quiver(:A, 1)
        @test q.B == zeros(Int, 1, 1)
        @test q.d == [1]
        @test q.n_mutable == 1

        q2 = Quiver(:A, 2)
        @test q2.B == [0 1; -1 0]
        @test q2.d == [1, 1]

        q3 = Quiver(:A, 3)
        @test q3.B == [0 1 0; -1 0 1; 0 -1 0]
        @test q3.d == [1, 1, 1]
    end

    @testset "B_n" begin
        # B_2: matches the OVERVIEW.md example (with numeric labels instead of ["a","b"])
        q = Quiver(:B, 2)
        @test q.B == [0 1; -2 0]
        @test q.d == [2, 1]

        q3 = Quiver(:B, 3)
        @test q3.B == [0 1 0; -1 0 1; 0 -2 0]
        @test q3.d == [2, 2, 1]

        # Symmetrizability is enforced by constructor; no error means it passed
        @test q3.n_mutable == 3
    end

    @testset "C_n" begin
        q = Quiver(:C, 2)
        @test q.B == [0 2; -1 0]
        @test q.d == [1, 2]

        q3 = Quiver(:C, 3)
        @test q3.B == [0 1 0; -1 0 2; 0 -1 0]
        @test q3.d == [1, 1, 2]
    end

    @testset "D_n" begin
        # D_4: fork at vertex 2
        q = Quiver(:D, 4)
        @test q.B == [0 1 0 0; -1 0 1 1; 0 -1 0 0; 0 -1 0 0]
        @test q.d == [1, 1, 1, 1]
        @test q.n_mutable == 4

        # D_5: fork at vertex 3
        q5 = Quiver(:D, 5)
        @test q5.B[3,4] == 1 && q5.B[3,5] == 1
        @test q5.B[4,3] == -1 && q5.B[5,3] == -1
        @test q5.d == ones(Int, 5)
    end

    @testset "E_n" begin
        q6 = Quiver(:E, 6)
        @test size(q6.B) == (6, 6)
        @test q6.d == ones(Int, 6)
        # Chain edges 1→2→3→4→5
        @test q6.B[1,2] == 1 && q6.B[2,3] == 1 && q6.B[3,4] == 1 && q6.B[4,5] == 1
        # Branch 3→6
        @test q6.B[3,6] == 1 && q6.B[6,3] == -1

        q7 = Quiver(:E, 7)
        @test size(q7.B) == (7, 7)
        @test q7.B[3,7] == 1

        q8 = Quiver(:E, 8)
        @test size(q8.B) == (8, 8)
        @test q8.B[3,8] == 1
    end

    @testset "F_4" begin
        q = Quiver(:F, 4)
        @test q.B == [0 1 0 0; -1 0 1 0; 0 -2 0 1; 0 0 -1 0]
        @test q.d == [2, 2, 1, 1]
        @test q.n_mutable == 4
    end

    @testset "G_2" begin
        q = Quiver(:G, 2)
        @test q.B == [0 1; -3 0]
        @test q.d == [3, 1]
        @test q.n_mutable == 2
    end

    # ─── Symbol constructor errors ────────────────────────────────────────────

    @testset "Symbol constructor errors" begin
        @test_throws InvalidArgument Quiver(:A, 0)
        @test_throws InvalidArgument Quiver(:B, 1)
        @test_throws InvalidArgument Quiver(:C, 1)
        @test_throws InvalidArgument Quiver(:D, 3)
        @test_throws InvalidArgument Quiver(:E, 5)
        @test_throws InvalidArgument Quiver(:E, 9)
        @test_throws InvalidArgument Quiver(:F, 3)
        @test_throws InvalidArgument Quiver(:G, 3)
        @test_throws InvalidArgument Quiver(:Z, 3)
    end

    # ─── String constructor ───────────────────────────────────────────────────

    @testset "String constructor" begin
        @test Quiver("A2").B == Quiver(:A, 2).B
        @test Quiver("B2").B == Quiver(:B, 2).B
        @test Quiver("D4").B == Quiver(:D, 4).B
        @test Quiver("G2").d == Quiver(:G, 2).d

        # Case-insensitive
        @test Quiver("a3").B == Quiver(:A, 3).B
        @test Quiver("g2").B == Quiver(:G, 2).B

        @test_throws InvalidArgument Quiver("")
        @test_throws InvalidArgument Quiver("X3")
        @test_throws InvalidArgument Quiver("A")
        @test_throws InvalidArgument Quiver("Afoo")
    end

    # ─── Structural equality: named vs hand-written ───────────────────────────

    @testset "Equality with matrix constructors" begin
        @test Quiver(:A, 2) == Quiver([0 1; -1 0])
        @test Quiver(:A, 3) == Quiver([0 1 0; -1 0 1; 0 -1 0])
        @test Quiver(:B, 2) == Quiver([0 1; -2 0], 2, [2, 1])
    end

    # ─── Edge-list constructor ─────────────────────────────────────────────────

    @testset "Edge-list constructor" begin
        # Linear A_3 from edge list
        q = Quiver([(1,2), (2,3)])
        @test q.B == Quiver(:A, 3).B
        @test q.n_mutable == 3

        # Single edge
        q2 = Quiver([(1,2)])
        @test q2.B == [0 1; -1 0]
        @test q2.n_mutable == 2

        # Weighted edge (double arrow)
        q3 = Quiver([(1,2,2)])
        @test q3.B == [0 2; -2 0]

        # Reversed arrow
        q4 = Quiver([(2,1)])
        @test q4.B == [0 -1; 1 0]

        # Error: empty list
        @test_throws InvalidArgument Quiver(Tuple{Int,Int}[])
    end

    # ─── Permissive matrix types ──────────────────────────────────────────────

    @testset "Permissive matrix coercion" begin
        B8 = Int8[0 1; -1 0]
        @test Quiver(B8) == Quiver(:A, 2)

        B32 = Int32[0 1 0; -1 0 1; 0 -1 0]
        @test Quiver(B32) == Quiver(:A, 3)
    end

    # ─── Mutation still works after named construction ────────────────────────

    @testset "Mutation after named construction" begin
        q = Quiver(:A, 3)
        q2 = mutate(q, 2)
        @test mutate(q2, 2) == q   # mutation is involution

        s = Seed(Quiver(:G, 2))
        s2 = mutate(s, 1)
        @test s2.cluster[1] != s.cluster[1]
    end

end
