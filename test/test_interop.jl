using Test
using ClusterAlgebras

@testset "dig6" begin

    @testset "golden strings" begin
        # A2: one arc 1 -> 2, so the n^2 adjacency bits are 0100 -> 010000 = 16.
        @test to_dig6(Quiver(:A, 2)) == ("&AO", Pair{Tuple{Int, Int}, Tuple{Int, Int}}[])

        # The Kronecker quiver has the same digraph, with a labelled arc.
        kronecker = Quiver([0 2; -2 0])
        @test to_dig6(kronecker) == ("&AO", [(0, 1) => (2, -2)])

        # B2 is not skew-symmetric: the label records both entries.
        b2 = Quiver([0 1; -2 0], 2, [2, 1])
        @test to_dig6(b2) == ("&AO", [(0, 1) => (1, -2)])

        # Empty and one-vertex quivers.
        @test to_dig6(Quiver(zeros(Int, 1, 1)))[1] == "&@?"
    end

    @testset "round trip is the identity on B" begin
        for q in (Quiver(:A, 3), Quiver(:D, 4), Quiver(:A, 5), Quiver([0 2; -2 0]),
                  Quiver([0 1; -2 0], 2, [2, 1]), Quiver([0 3; -1 0], 2, [1, 3]),
                  Quiver(:E, 6))
            back = from_dig6(to_dig6(q))
            @test back.B == q.B
            @test back.d == q.d
            @test back.n_mutable == nvertices(q)
        end
    end

    @testset "round trip through the canonical form is an isomorphism" begin
        q = Quiver(:D, 4)
        p = permute_vertices(q, [3, 1, 4, 2])
        @test to_dig6(q; canonical = true) == to_dig6(p; canonical = true)
        @test is_isomorphic(from_dig6(to_dig6(q; canonical = true)), q)
    end

    @testset "mutation moves the encoding" begin
        q = Quiver(:A, 3)
        @test to_dig6(mutate(q, 2)) != to_dig6(q)
        @test from_dig6(to_dig6(mutate(q, 2))).B == mutate(q, 2).B
    end

    @testset "frozen vertices are restored via n_mutable" begin
        B = [0 1 0; -1 0 1; 0 -1 0]
        q = Quiver(B, 2)
        back = from_dig6(to_dig6(q); n_mutable = 2)
        @test back.B == B
        @test back.n_frozen == 1
    end

    @testset "accepts Sage-style input" begin
        @test from_dig6(">>digraph6<<&AO").B == Quiver(:A, 2).B
        @test from_dig6("&AO", Dict((0, 1) => (2, -2))).B == [0 2; -2 0]
        @test from_dig6("&AO", [(0, 1) => (2, -2)]).B == [0 2; -2 0]
    end

    @testset "larger vertex counts use the four-character prefix" begin
        n = 70
        B = zeros(Int, n, n)
        for i in 1:(n - 1)
            B[i, i + 1], B[i + 1, i] = 1, -1
        end
        q = Quiver(B)
        s, edges = to_dig6(q)
        @test s[2] == Char(126)
        @test isempty(edges)
        @test from_dig6(s, edges).B == B
    end

    @testset "invalid input throws" begin
        @test_throws InvalidArgument from_dig6("AO")
        @test_throws InvalidArgument from_dig6("&A")
        # "&A~" decodes to a digraph with loops, which no exchange matrix has.
        @test_throws NotSkewSymmetrizable from_dig6("&A~")
    end
end

@testset "qmu" begin

    @testset "section layout matches the applet format" begin
        text = to_qmu(Quiver(:A, 2))
        lines = split(text, '\n')
        @test lines[1:8] == ["//Number of points", "2", "//Vertex radius", "9",
                             "//Labels shown", "1", "//Matrix", "2 2"]
        @test lines[9:10] == ["0 1", "-1 0"]
        @test lines[11] == "//Points"
        @test lines[(end - 4):end] == ["//Historycounter", "-1", "//History", "",
                                       "//Cluster is null"]
    end

    @testset "round trip is the identity on B" begin
        for q in (Quiver(:A, 3), Quiver(:D, 4), Quiver([0 2; -2 0]),
                  Quiver([0 1; -2 0], 2, [2, 1]), Quiver(:E, 6),
                  mutate(Quiver(:A, 4), 2))
            back = from_qmu(to_qmu(q))
            @test back.B == q.B
            @test back.d == q.d
            @test back.n_frozen == 0
        end
    end

    @testset "frozen vertices survive the round trip" begin
        B = [0 1 0 1; -1 0 1 0; 0 -1 0 0; -1 0 0 0]
        q = Quiver(B, 2)
        back = from_qmu(to_qmu(q))
        @test back.n_mutable == 2
        @test back.n_frozen == 2
        @test back.B == B
    end

    @testset "the frozen block is completed the way Sage completes it" begin
        # B[frozen, mutable] is authoritative; the transpose block is rebuilt.
        B = [0 1 7; -1 0 0; -1 0 0]
        q = Quiver(B, 2)
        rows = split(to_qmu(q), '\n')[9:11]
        @test rows == ["0 1 1", "-1 0 0", "-1 0 0"]
        @test from_qmu(to_qmu(q)).B == [0 1 1; -1 0 0; -1 0 0]
    end

    @testset "file round trip" begin
        mktempdir() do dir
            q = Quiver(:D, 4)
            path = write_qmu(joinpath(dir, "d4"), q)
            @test endswith(path, "d4.qmu")
            @test read_qmu(path).B == q.B
            @test write_qmu(joinpath(dir, "d4.qmu"), q) == joinpath(dir, "d4.qmu")
        end
    end

    @testset "malformed input throws" begin
        @test_throws InvalidArgument from_qmu("//Number of points\n2")
        @test_throws InvalidArgument from_qmu("//Matrix\n2 3\n0 1\n-1 0")
        @test_throws InvalidArgument from_qmu("//Matrix\n2 2\n0 1")
        @test_throws InvalidArgument from_qmu("//Matrix\n2 2\n0 1 0\n-1 0 0")
        # Frozen points must come last.
        @test_throws InvalidArgument from_qmu(
            "//Matrix\n2 2\n0 1\n-1 0\n//Points\n9 0 0 1\n9 1 0")
    end
end

@testset "exchange matrix round trip" begin
    for q in (Quiver(:A, 3), Quiver([0 1; -2 0], 2, [2, 1]), Quiver([0 1 0; -1 0 1; 0 -1 0], 2))
        @test Quiver(q.B, q.n_mutable, q.d, labels(q)) == q
    end
end
