using Test
using ClusterAlgebras

# The cycle on n vertices carrying p arrows one way round and n - p the other:
# the affine type Ã(p, n - p), i.e. the annulus. Degenerate at n = 2, where both
# sides of the cycle join the same pair of vertices, so callers start at n = 3.
function cycle_quiver(n::Int, p::Int)
    @assert n >= 3
    B = zeros(Int, n, n)
    for k in 1:n
        i, j = k, mod1(k + 1, n)
        if k <= p
            B[i, j], B[j, i] = 1, -1
        else
            B[i, j], B[j, i] = -1, 1
        end
    end
    return Quiver(B)
end

# The once-punctured torus: three vertices, double arrows round a cycle.
markov_quiver() = Quiver([0 2 -2; -2 0 2; 2 -2 0])

@testset "block decomposition" begin

    @testset "the block table" begin
        @test length(FST_BLOCKS) == 7
        @test sort([b.label for b in FST_BLOCKS]) ==
              sort([:point, :I, :II, :IIIa, :IIIb, :IV, :V])

        for b in FST_BLOCKS
            @test b.B == -transpose(b.B)            # blocks are skew-symmetric
            @test all(abs.(b.B) .<= 1)              # single arrows inside a block
            @test !isempty(b.outlets)               # a block with no outlet glues to nothing
            @test b.outlets ⊆ 1:size(b.B, 1)
        end

        sizes = Dict(b.label => size(b.B, 1) for b in FST_BLOCKS)
        @test sizes == Dict(:point => 1, :I => 2, :II => 3, :IIIa => 3,
                            :IIIb => 3, :IV => 4, :V => 5)

        outlets = Dict(b.label => length(b.outlets) for b in FST_BLOCKS)
        @test outlets == Dict(:point => 1, :I => 2, :II => 3, :IIIa => 1,
                              :IIIb => 1, :IV => 2, :V => 1)

        # The transcription oracle: the figure's two largest blocks are quivers
        # the literature names, and they come out right.
        by_label = Dict(b.label => b for b in FST_BLOCKS)
        @test mutation_type(Quiver(by_label[:IV].B)) == MutationType(:D, 4)
        @test mutation_type(Quiver(by_label[:V].B)) == MutationType(:D, 4, 1)
        @test mutation_type(Quiver(by_label[:II].B)) == MutationType(:A, 3)

        # IIIa and IIIb are genuinely different blocks, not one drawn twice.
        @test by_label[:IIIa].B != by_label[:IIIb].B
        @test by_label[:IIIa].B == -by_label[:IIIb].B

        # Every block is a quiver of a triangulated piece, so decomposable.
        for b in FST_BLOCKS
            @test is_block_decomposable(Quiver(b.B))
        end
    end

    @testset "roundtrip: reassembling returns the quiver" begin
        for q in [Quiver(:A, 3), Quiver(:A, 5), Quiver(:D, 4), Quiver(:D, 6),
                  cycle_quiver(4, 1), cycle_quiver(5, 2), markov_quiver()]
            d = block_decomposition(q)
            @test d !== nothing
            @test reassemble(d, nvertices(q)) == q.B
        end
    end

    @testset "surface types decompose" begin
        for n in 1:7
            @test is_block_decomposable(Quiver(:A, n))
        end
        for n in 4:7
            @test is_block_decomposable(Quiver(:D, n))       # once-punctured disc
        end
        for (n, p) in [(3, 1), (4, 1), (4, 2), (5, 2), (6, 3)]
            @test is_block_decomposable(cycle_quiver(n, p))  # annulus Ã(p, q)
        end
        @test is_block_decomposable(Quiver([0 2; -2 0]))     # Ã(1, 1), the annulus
        @test is_block_decomposable(markov_quiver())         # once-punctured torus
    end

    @testset "exceptional types do not decompose" begin
        # Felikson-Shapiro-Tumarkin: every minimal non-decomposable quiver is
        # mutation equivalent to X6 or E6, and the eleven exceptional classes
        # are exactly the mutation-finite ones that fail to decompose.
        for n in 6:8
            @test !is_block_decomposable(Quiver(:E, n))
        end
        for n in 6:9
            q = ClusterAlgebras._exceptional_quiver(n)
            q === nothing && continue
            @test !is_block_decomposable(q)
        end
    end

    @testset "decomposable implies mutation finite" begin
        for q in [Quiver(:A, 4), Quiver(:D, 5), cycle_quiver(4, 2), markov_quiver()]
            @test is_block_decomposable(q)
            @test is_mutation_finite(q)
        end
        # The converse fails only on the exceptional classes.
        @test is_mutation_finite(Quiver(:E, 6))
        @test !is_block_decomposable(Quiver(:E, 6))
    end

    @testset "mutation invariance" begin
        # Surface quivers are closed under mutation, so decomposability is a
        # property of the whole class.
        for q in [Quiver(:A, 4), Quiver(:D, 4), markov_quiver()]
            mc = mutation_class(q; up_to_isomorphism = true)
            @test !is_truncated(mc)
            @test all(is_block_decomposable, mc.quivers)
        end
        mc = mutation_class(Quiver(:E, 6); up_to_isomorphism = true)
        @test !is_truncated(mc)
        @test !any(is_block_decomposable, mc.quivers)
    end

    @testset "mutation-infinite quivers do not decompose" begin
        @test !is_block_decomposable(Quiver([0 3; -3 0]))
        @test !is_block_decomposable(Quiver([0 3 -3; -3 0 3; 3 -3 0]))
        # Weight 2 everywhere is within the block bound, but this one is still
        # mutation infinite, so it must not decompose.
        wild = Quiver([0 2 -2 0; -2 0 2 2; 2 -2 0 -2; 0 -2 2 0])
        @test is_block_decomposable(wild) == is_mutation_finite(wild)
    end

    @testset "is_surface_type" begin
        for q in [Quiver(:A, 4), Quiver(:D, 5), markov_quiver(), Quiver(:E, 6)]
            @test is_surface_type(q) == is_block_decomposable(q)
        end
        # The three reasons mutation_type answers `nothing`, told apart.
        m = markov_quiver()
        @test mutation_type(m) === nothing
        @test is_surface_type(m) && is_mutation_finite(m)

        e = ClusterAlgebras._exceptional_quiver(6)     # X6
        @test !is_surface_type(e) && is_mutation_finite(e)

        inf = Quiver([0 3 -3; -3 0 3; 3 -3 0])
        @test !is_surface_type(inf) && !is_mutation_finite(inf)
    end

    @testset "errors and edge cases" begin
        @test is_block_decomposable(Quiver(:A, 1))
        @test_throws InvalidArgument block_decomposition(Quiver(zeros(Int, 0, 0)))
        # Skew-symmetrizable but not skew-symmetric: blocks say nothing there.
        @test_throws InvalidArgument block_decomposition(Quiver(:B, 3))
    end

end
