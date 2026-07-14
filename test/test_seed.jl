@testset "Seed construction" begin
    B = [0 1; -1 0]
    q = Quiver(B)
    s = Seed(q)

    @test s.quiver === q
    @test length(s.cluster) == 2

    # Custom variable names
    s2 = Seed(q, ["a", "b"])
    @test string(numerator(s2.cluster[1])) == "a"
    @test string(numerator(s2.cluster[2])) == "b"

    # Wrong number of names
    @test_throws ErrorException Seed(q, ["a"])

    # Seed with frozen vertex
    B3 = [0 1 0; -1 0 0; 0 0 0]
    q3 = Quiver(B3, 2)
    s3 = Seed(q3)
    @test length(s3.cluster) == 3
end

@testset "Seed collection interface" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B))

    @test length(s) == 2
    @test s[1] == s.cluster[1]
    @test s[2] == s.cluster[2]

    collected = collect(s)
    @test collected == s.cluster
end

@testset "Seed mutation_path" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B))

    @test s.mutation_path == Int[]

    s1 = mutate(s, 1)
    @test s1.mutation_path == [1]

    s12 = mutate(s1, 2)
    @test s12.mutation_path == [1, 2]

    # Sequence form also accumulates path
    s_seq = mutate(s, [1, 2, 1])
    @test s_seq.mutation_path == [1, 2, 1]
end

@testset "Seed show uses labels" begin
    B = [0 1; -1 0]
    q = Quiver(B, 2, [1, 1], ["alpha", "beta"])
    s = Seed(q)
    buf = IOBuffer()
    show(buf, MIME"text/plain"(), s)
    out = String(take!(buf))
    @test occursin("alpha", out)
    @test occursin("beta", out)
    # must not fall back to hardcoded x1/x2
    @test !occursin("x1 =", out)
    @test !occursin("x2 =", out)
end

@testset "Seed text/latex show" begin
    B = [0 1; -1 0]
    s = Seed(Quiver(B))
    buf = IOBuffer()
    show(buf, MIME"text/latex"(), s)
    out = String(take!(buf))
    @test occursin("\\begin{aligned}", out)
    @test occursin("\\end{aligned}", out)
    @test occursin("x_{1}", out)
    @test occursin("x_{2}", out)
end
