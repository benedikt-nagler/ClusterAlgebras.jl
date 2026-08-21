@testset "Green sequences" begin

    # ─── Predicates ──────────────────────────────────────────────────────────

    @testset "is_green / is_red / is_all_red" begin
        s = extend(Seed(Quiver(:A, 2)))

        # Initial seed: C = I₂, both vertices green.
        @test is_green(s, 1)
        @test is_green(s, 2)
        @test !is_red(s, 1)
        @test !is_all_red(s)

        # After μ₁: C = [-1 1; 0 1] - c₁ = (−1,0) red, c₂ = (1,1) green.
        s1 = mutate(s, 1)
        @test is_red(s1, 1)
        @test is_green(s1, 2)
        @test !is_all_red(s1)

        # After μ₁,μ₂: C = [0 -1; 1 -1] - c₁ = (0,1) green, c₂ = (−1,−1) red.
        s12 = mutate(s1, 2)
        @test is_green(s12, 1)
        @test is_red(s12, 2)
        @test !is_all_red(s12)

        # After μ₁,μ₂,μ₁: C = [0 -1; -1 0] - both red.
        s121 = mutate(s12, 1)
        @test is_red(s121, 1)
        @test is_red(s121, 2)
        @test is_all_red(s121)
    end

    # ─── A1: single MGS [1] ──────────────────────────────────────────────────

    @testset "A1" begin
        s = Seed(Quiver(:A, 1))
        seqs = maximal_green_sequences(s)
        @test length(seqs) == 1
        @test seqs[1] == [1]
    end

    # ─── A2: exactly two MGS ─────────────────────────────────────────────────

    @testset "A2" begin
        s = Seed(Quiver(:A, 2))
        seqs = maximal_green_sequences(s)

        # Verify count.
        @test length(seqs) == 2

        # Both known MGS must be present.  For B = [0 1; -1 0] (arrow 1 → 2)
        # the sink-first sequence [2,1] has length 2 and the source-first
        # sequence [1,2,1] has length 3.
        @test [2, 1]    ∈ seqs
        @test [1, 2, 1] ∈ seqs
    end

    # ─── Structural checks (A3) ───────────────────────────────────────────────

    @testset "A3 structural" begin
        s = Seed(Quiver(:A, 3))
        seqs = maximal_green_sequences(s)

        # There must be at least one MGS.
        @test !isempty(seqs)

        for seq in seqs
            # Replay the sequence and check every step mutates a green vertex.
            sp = extend(s)
            for k in seq
                @test is_green(sp, k)
                sp = mutate(sp, k)
            end
            # Terminal seed must be all-red.
            @test is_all_red(sp)
        end
    end

    # ─── Auto-extend of TrivialCoefficients seed ─────────────────────────────

    @testset "auto-extend" begin
        s_trivial  = Seed(Quiver(:A, 2))
        s_extended = extend(s_trivial)
        @test maximal_green_sequences(s_trivial)  == maximal_green_sequences(s_extended)
    end

    # ─── D4: finite type, all sequences valid ─────────────────────────────────

    @testset "D4 structural" begin
        s = Seed(Quiver(:D, 4))
        seqs = maximal_green_sequences(s)

        @test !isempty(seqs)
        for seq in seqs
            sp = extend(s)
            for k in seq
                @test is_green(sp, k)
                sp = mutate(sp, k)
            end
            @test is_all_red(sp)
        end
    end

    # ─── Sign tracking (green_sequence_signs) ──────────────────────────────────

    @testset "green_sequence_signs" begin
        s = Seed(Quiver(:A, 2))
        seq = first(maximal_green_sequences(s))    # an actual A₂ MGS

        signs = green_sequence_signs(s, seq)
        @test size(signs) == (length(seq) + 1, 2)  # (length+1) × n_mutable
        @test all(==(:green), signs[1, :])         # start: all green
        @test all(==(:red),   signs[end, :])       # end: all red (MGS theorem)
        # sign-coherent: no seed shows a :mixed c-vector
        @test !any(==(:mixed), signs)

        # matches the step-by-step per-vertex predicates
        sp = extend(s)
        @test signs[1, 1] == (is_green(sp, 1) ? :green : :red)

        # accepts an already-extended seed too, same answer
        @test green_sequence_signs(extend(s), seq) == signs

        # a non-maximal prefix need not end all-red - [1] on A₂ leaves v2 green
        partial = green_sequence_signs(s, [1])
        @test size(partial) == (2, 2)
        @test partial[2, 2] == :green

        # every MGS of A₃ starts all-green, ends all-red
        for seq3 in maximal_green_sequences(Seed(Quiver(:A, 3)))
            m = green_sequence_signs(Seed(Quiver(:A, 3)), seq3)
            @test all(==(:green), m[1, :])
            @test all(==(:red), m[end, :])
        end
    end

    # ─── Certified decision search (mgs_search) ───────────────────────────────

    @testset "mgs_search" begin
        # Acyclic Aₙ: minimal MGS length is exactly n (every vertex must be
        # mutated at least once; the sink-source order achieves n).
        for n in 1:4
            r = mgs_search(Seed(Quiver(:A, n)))
            @test r.status == :found
            @test r.min_length == n
            # The returned sequence is a genuine MGS: replay and check.
            sp = extend(Seed(Quiver(:A, n)))
            for k in r.sequence
                @test is_green(sp, k)
                sp = mutate(sp, k)
            end
            @test is_all_red(sp)
        end

        # A₂ minimal sequence is the sink-first one, [2, 1].
        @test mgs_search(Seed(Quiver(:A, 2))).sequence == [2, 1]

        # D₄: minimal length = rank as well.
        @test mgs_search(Seed(Quiver(:D, 4))).min_length == 4

        # Consistency with the enumerator: min_length equals the shortest
        # enumerated MGS on small finite types.
        for q in (Quiver(:A, 3), Quiver(:D, 4))
            seqs = maximal_green_sequences(Seed(q))
            @test mgs_search(Seed(q)).min_length == minimum(length, seqs)
        end

        # Markov quiver: no MGS (theorem), and by sign coherence every search
        # dead end would be an MGS - so non-existence shows up as an infinite
        # search: a budget verdict at any finite budget, never :found or :none.
        markov = Seed(Quiver([0 2 -2; -2 0 2; 2 -2 0]))
        r = mgs_search(markov; max_length = 50, max_nodes = 2000)
        @test r.status == :unknown                     # node budget binds
        @test r.min_length === nothing && r.sequence === nothing
        # With the length bound binding instead: certified "no MGS of
        # length ≤ 8" (every depth-≤8 state fits in the node budget).
        @test mgs_search(markov; max_length = 8, max_nodes = 10^6).status ==
              :none_within_length

        # Kronecker quiver (double arrow, rank 2): infinite type but an MGS
        # of length 2 exists (sink-first).
        r = mgs_search(Seed(Quiver([0 2; -2 0])))
        @test r.status == :found
        @test r.min_length == 2

        # Length budget below the minimal length → certified "none that short".
        @test mgs_search(Seed(Quiver(:A, 3)); max_length = 2).status ==
              :none_within_length
        # Node budget too small → :unknown, never a wrong certification.
        @test mgs_search(Seed(Quiver(:A, 3)); max_nodes = 2).status == :unknown

        # TrivialCoefficients and extended seeds agree.
        @test mgs_search(Seed(Quiver(:A, 2))) == mgs_search(extend(Seed(Quiver(:A, 2))))

        # Overflow guard: under Int64 this wild quiver silently wrapped its
        # C-entries (they reach ~2^66) and fabricated a length-10 "MGS";
        # BigInt ground truth says no MGS of length ≤ 10 exists. With the
        # Int128 guard the fake :found is gone - pruning may only report
        # :unknown, never certify falsely.
        wild = Seed(Quiver([0 2 0 2 0; -2 0 2 -2 -2; 0 -2 0 2 0;
                            -2 2 -2 0 -2; 0 2 0 2 0]))
        r = mgs_search(wild; max_length = 10, max_nodes = 10^6)
        @test r.status == :unknown
    end

    @testset "verify_mutation_sequence" begin
        a2 = Seed(Quiver(:A, 2))

        # Both A₂ maximal green sequences certify, with the pentagon charges.
        r = verify_mutation_sequence(a2, [2, 1])
        @test r.valid && r.maximal && r.reason === nothing
        @test r.charges == [[0, 1], [1, 0]]
        r = verify_mutation_sequence(a2, [1, 2, 1])
        @test r.valid && r.charges == [[1, 0], [1, 1], [0, 1]]
        # Charges agree with ordered_c_vectors on every enumerated A₃ MGS.
        a3 = Seed(Quiver(:A, 3))
        for seq in maximal_green_sequences(a3)
            @test verify_mutation_sequence(a3, seq).charges == ordered_c_vectors(a3, seq)
        end

        # Red mutation rejected with the offending step named.
        r = verify_mutation_sequence(a2, [2, 1, 1])
        @test !r.valid && r.charges === nothing
        @test r.reason == "vertex 1 not green at step 3"

        # Out-of-range vertex.
        @test occursin("out of range", verify_mutation_sequence(a2, [3]).reason)

        # Green but not maximal: rejected by default, accepted with the flag
        # (maximal = false, prefix charges returned).
        r = verify_mutation_sequence(a2, [1, 2])
        @test !r.valid && occursin("not maximal", r.reason)
        r = verify_mutation_sequence(a2, [1, 2]; require_maximal = false)
        @test r.valid && !r.maximal && r.charges == [[1, 0], [1, 1]]
        @test verify_mutation_sequence(a2, Int[]; require_maximal = false).valid

        # Quiver convenience method and TrivialCoefficients seeds agree.
        @test verify_mutation_sequence(Quiver(:A, 2), [2, 1]).valid
        @test verify_mutation_sequence(extend(a2), [2, 1]).charges ==
              verify_mutation_sequence(a2, [2, 1]).charges

        # Beyond the Int128 guard (the finding-13 quiver): drive the C-entries
        # past 2^62 with an exact BigInt replay that always mutates the green
        # vertex with the largest entries. The verifier's fast path trips its
        # guard and the automatic BigInt fallback must still certify the
        # (genuinely green) word exactly - the verifier has no ceiling.
        Bw = [0 2 0 2 0; -2 0 2 -2 -2; 0 -2 0 2 0; -2 2 -2 0 -2; 0 2 0 2 0]
        B, C = BigInt.(Bw), BigInt.([i == j for i in 1:5, j in 1:5])
        word = Int[]
        expected = Vector{Vector{BigInt}}()
        while maximum(abs, C) <= big(2)^62 && length(word) < 60
            greens = [k for k in 1:5 if all(>=(0), view(C, :, k))]
            @test !isempty(greens)
            k = argmax(Dict(k => maximum(abs, view(C, :, k)) for k in greens))
            push!(word, k)
            push!(expected, C[:, k])
            B, C = ClusterAlgebras._mutate_matrix(B, k),
                   ClusterAlgebras._mutate_C(C, B, k)
        end
        @test maximum(abs, C) > big(2)^62    # the word really exceeds the guard
        r = verify_mutation_sequence(Seed(Quiver(Bw)), word; require_maximal = false)
        @test r.valid && r.charges == expected
        # A red mutation appended after the deep prefix is still refused.
        greens = [k for k in 1:5 if all(>=(0), view(C, :, k))]
        reds = setdiff(1:5, greens)
        @test !isempty(reds)
        r = verify_mutation_sequence(Seed(Quiver(Bw)), vcat(word, reds[1]);
                                     require_maximal = false)
        @test !r.valid && occursin("not green", r.reason)
        # require_green = false accepts the same red step.
        r = verify_mutation_sequence(Seed(Quiver(Bw)), vcat(word, reds[1]);
                                     require_maximal = false, require_green = false)
        @test r.valid
    end

    # Q_{a,b,c}: a arrows 1→2, b arrows 2→3, c arrows 3→1 (Muller arXiv:1503.04675).
    Qabc(a, b, c) = Quiver([0 a -c; -a 0 b; c -b 0])

    @testset "reddening_search" begin
        # Every maximal green sequence is a reddening sequence, so on finite type
        # both exist and the shortest reddening sequence is no longer.
        for q in (Quiver(:A, 2), Quiver(:A, 3), Quiver(:D, 4), Quiver(:B, 3),
                  Quiver(:A, 4))
            m = mgs_search(Seed(q); max_length = 30, max_nodes = 500_000)
            r = reddening_search(Seed(q); max_length = 30, max_nodes = 500_000)
            @test m.status === :found
            @test r.status === :found
            @test r.min_length <= m.min_length
            # The word reaches all-red but need not be green.
            v = verify_mutation_sequence(Seed(q), r.sequence; require_green = false)
            @test v.valid && v.maximal
        end

        # [Mul16] Thm 2.3.1 + Fig. 11: Q_{2,2,3} has NO maximal green sequence and
        # DOES have a reddening sequence. This is the separator between the two
        # notions, and the oracle that catches a search still requiring greenness.
        q223 = Qabc(2, 2, 3)
        @test mgs_search(Seed(q223); max_length = 12,
                         max_nodes = 2_000_000).status === :none_within_length
        r223 = reddening_search(Seed(q223); max_length = 14, max_nodes = 2_000_000)
        @test r223.status === :found
        @test r223.min_length == 6
        # It verifies as a reddening sequence and is refused as a green one.
        @test verify_mutation_sequence(Seed(q223), r223.sequence;
                                       require_green = false).maximal
        vg = verify_mutation_sequence(Seed(q223), r223.sequence; require_green = true)
        @test !vg.valid && occursin("not green", vg.reason)

        # [BDP14]: the Markov quiver Q_{2,2,2} has neither.
        q222 = Qabc(2, 2, 2)
        @test mgs_search(Seed(q222); max_length = 12,
                         max_nodes = 2_000_000).status === :none_within_length
        @test reddening_search(Seed(q222); max_length = 12,
                               max_nodes = 2_000_000).status === :none_within_length

        # [Kel17] Thm 4.7: existence is invariant under mutation. MGS existence is
        # not, so this is a property only the reddening search may have.
        for q0 in (Quiver(:A, 3), Quiver(:D, 4), Qabc(2, 2, 2))
            mc = mutation_class(q0; max_quivers = 200)
            @test !is_truncated(mc)
            sts = unique([reddening_search(Seed(mc[i]); max_length = 14,
                                           max_nodes = 300_000).status
                          for i in 1:length(mc)])
            @test length(sts) == 1
        end

        # Quiver method, and an already-all-red seed as a length-0 sequence. The
        # seed must carry principal coefficients: c-vectors do not exist on a
        # trivial-coefficient seed, where the search restarts from C = I.
        @test reddening_search(Quiver(:A, 2)).status === :found
        red2 = mutate(extend(Seed(Quiver(:A, 2))), [2, 1])
        @test is_all_red(red2)
        @test reddening_search(red2).min_length == 0
    end

end
