using Random

@testset "forks, pre-forks and the forkless part" begin

    A3     = Quiver(:A, 3)
    markov = Quiver([0 2 -2; -2 0 2; 2 -2 0])
    fork3  = Quiver([0 2 -3; -2 0 2; 3 -2 0])
    box(a, b) = Quiver([0 a 0 -b; -a 0 b 0; 0 -b 0 a; b 0 -a 0])

    # ── abundance ─────────────────────────────────────────────────────────────
    @test is_abundant(markov)
    @test is_abundant(fork3)
    @test !is_abundant(A3)
    @test is_abundant(Quiver([0;;]))          # vacuous on one vertex

    # ── the fork definition (Warkentin Def 2.1) ───────────────────────────────
    @test is_fork(fork3)
    @test point_of_return(fork3) == 2
    @test !is_fork(markov)                     # abundant and cyclic, but 2 ≯ 2
    @test !is_fork(A3)                         # not abundant
    @test point_of_return(A3) === nothing
    @test !is_fork(Quiver([0 3; -3 0]))        # rank 2 is always acyclic

    # ── Warkentin Lemma 2.5: mutating a fork away from r gives a fork with the
    #    new vertex as point of return, and strictly more arrows ──────────────
    n_arrows(q) = sum(abs, q.B) ÷ 2
    for F in (fork3, mutate(fork3, 1), mutate(mutate(fork3, 1), 3))
        r = point_of_return(F)
        @test r !== nothing
        for k in 1:nvertices(F)
            k == r && continue
            Fk = mutate(F, k)
            @test point_of_return(Fk) == k
            @test n_arrows(Fk) > n_arrows(F)
        end
    end

    # ── Warkentin Lemma 2.4: a full subquiver of a fork is abundant acyclic, or
    #    a fork with the same point of return ──────────────────────────────────
    let F = mutate(mutate(fork3, 1), 3), r = point_of_return(F), n = nvertices(F)
        for k in 1:n
            k == r && continue
            S = [v for v in 1:n if v != k]
            sub = Quiver(F.B[S, S], length(S))
            rr = point_of_return(sub)
            @test rr === nothing ? (is_abundant(sub) && is_acyclic(sub)) : S[rr] == r
        end
    end

    # ── the point of return is unique ─────────────────────────────────────────
    let rng = MersenneTwister(20260906), n_forks = 0
        for _ in 1:400
            q = random_quiver(rng, 4; max_arrows = 4, density = 1.0)
            r = point_of_return(q)
            r === nothing && continue
            n_forks += 1
            # no second vertex passes the definition
            others = 0
            for s in 1:4
                s == r && continue
                into  = [i for i in 1:4 if i != s && q.B[i, s] > 0]
                outof = [j for j in 1:4 if j != s && q.B[s, j] > 0]
                ok = all(q.B[j, i] > q.B[i, s] && q.B[j, i] > q.B[s, j]
                         for i in into, j in outof)
                ok &= is_acyclic(Quiver(q.B[into, into], length(into)))
                ok &= is_acyclic(Quiver(q.B[outof, outof], length(outof)))
                others += ok
            end
            @test others == 0
        end
        @test n_forks > 0                     # the sampler did produce forks
    end

    # ── Warkentin Thm 3.2: a connected quiver is mutation-infinite iff it is
    #    mutation-equivalent to a fork ───────────────────────────────────────
    has_fork(q; cap = 400) =
        any(is_fork, mutation_class(q; max_quivers = cap).quivers)

    for q in (Quiver(:A, 3), Quiver(:D, 4), Quiver(:A, 2), markov,
              Quiver([0 2; -2 0]))
        @test !is_mutation_finite(q) || !has_fork(q)   # finite ⇒ no fork anywhere
    end
    for q in (Quiver([0 3 -3; -3 0 3; 3 -3 0]), fork3, box(2, 2),
              Quiver([0 2 -1; -2 0 2; 1 -2 0]))
        @test has_fork(q)                              # infinite ⇒ a fork exists
    end

    # ── pre-forks (Warkentin Def 3.7) ─────────────────────────────────────────
    @test !is_pre_fork(A3)
    @test pre_fork_vertices(A3) === nothing
    let P = mutate(box(2, 2), 1)                 # one step off the box quiver
        pf = pre_fork_vertices(P)
        if pf !== nothing
            n = nvertices(P)
            for (k, other) in ((pf.k, pf.k′), (pf.k′, pf.k))
                S = [v for v in 1:n if v != k]
                sub = Quiver(P.B[S, S], length(S))
                @test S[point_of_return(sub)] == pf.r
            end
        end
    end

    # ── mutation-finite ⇒ the forkless part is the whole labeled class ────────
    for q in (A3, markov, Quiver(:D, 4))
        mc = mutation_class(q; max_quivers = 5000)
        fp = forkless_part(q; max_quivers = 5000)
        @test !is_truncated(fp)
        @test length(fp) == length(mc)
        @test has_finite_forkless_part(q) === true
    end

    # ── Ervin Cor 3.6: mutation-invariance.  The part is one set, so it is the
    #    same from every member of it ─────────────────────────────────────────
    let fp = forkless_part(A3)
        keys0 = Set(q.B for q in fp.quivers)
        for i in (2, 5, 9)
            fpi = forkless_part(fp[i])
            @test Set(q.B for q in fpi.quivers) == keys0
        end
    end

    # ── Warkentin Ex 6.11 (via Ervin §5): the (a,b) box quivers with a,b ≥ 2 are
    #    the first family with a finite pre-forkless part and an infinite
    #    forkless part ───────────────────────────────────────────────────────
    for (a, b) in ((2, 2), (2, 3), (3, 3))
        Q = box(a, b)
        @test has_finite_pre_forkless_part(Q; max_quivers = 2000) === true
        @test has_finite_forkless_part(Q; max_quivers = 800) === nothing
        @test is_truncated(forkless_part(Q; max_quivers = 800))
    end
    # min(a,b) = 1 leaves the family: those are mutation-acyclic
    @test has_finite_pre_forkless_part(box(1, 2); max_quivers = 800) === nothing

    # ── Ervin Thm 3.5: hereditary, with the max{M, m} bound ──────────────────
    let Q = box(2, 2), M = length(pre_forkless_part(Q; max_quivers = 2000)), n = 4
        for k in 1:n
            S = [v for v in 1:n if v != k]
            sub = Quiver(Q.B[S, S], length(S))
            fp = pre_forkless_part(sub; max_quivers = 4 * M + 100)
            @test !is_truncated(fp)
            @test length(fp) <= max(M, length(S))
        end
    end

    # ── a finite forkless part implies a finite pre-forkless part ────────────
    @test has_finite_pre_forkless_part(A3) === true

    # ── the ForklessPart object ──────────────────────────────────────────────
    let fp = forkless_part(A3)
        @test fp.kind === :fork
        @test length(fp.adj) == length(fp)
        @test all(length(a) == 3 for a in fp.adj)
        @test all(0 <= j <= length(fp) for a in fp.adj for j in a)
        @test fp[1].B == A3.B
        @test occursin("forkless part", sprint(show, fp))
        @test occursin("pre-forkless part", sprint(show, pre_forkless_part(A3)))
    end

    # ── guards: fork theory is stated for quivers ────────────────────────────
    @test_throws InvalidArgument point_of_return(Quiver([0 1 1; -1 0 0; -1 0 0], 1))
    @test_throws InvalidArgument is_fork(Quiver([0 1; -2 0], 2, [2, 1]))
    @test_throws InvalidArgument forkless_part(A3; max_quivers = 0)
end
