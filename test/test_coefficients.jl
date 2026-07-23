# ─── test_coefficients.jl ─────────────────────────────────────────────────────
# Tests for the principal-coefficient Seed, y-dynamics, separation formula.

# ─── Helper: tropicalize a fraction-field element ─────────────────────────────
#
# For a subtraction-free rational function p/q in y₁,…,yₙ (over ZZ), the
# tropicalization under the min convention is:
#   trop(p/q) = componentwise_min(exponents of p) - componentwise_min(exponents of q)
# where min is taken over all monomials.
function _tropicalize(frac_elem, n::Int)
    p = numerator(frac_elem)
    q = denominator(frac_elem)
    min_exp(poly) = begin
        evs = collect(exponent_vectors(poly))
        isempty(evs) && return zeros(Int, n)
        reduce((a, b) -> min.(a, b), evs)
    end
    return min_exp(p) - min_exp(q)
end

# ─── Backward-compatibility: trivial Seed unchanged ───────────────────────────

@testset "backward compat — Seed{TrivialCoefficients} unchanged" begin
    q  = Quiver(:A, 2)
    s  = Seed(q)
    @test s isa Seed{TrivialCoefficients}

    # Check that mutate still produces the correct A₂ cluster variables
    # Initial cluster: (x1, x2)
    # After μ₁: x1' = (1 + x2)/x1, x2 unchanged
    # After μ₁μ₂: x2' = (x1' + 1)/x2 = (1 + x2 + x1)/(x1*x2)
    # After μ₁μ₂μ₁: x1'' = (x2' + 1)/x1' = (x2 + 1)/x1   (this is x2 originally by 5-period)
    R  = s.ring
    x1, x2 = gens(base_ring(R))
    X1, X2 = R(x1), R(x2)

    s1 = mutate(s, 1)
    @test numerator(s1[1] * s1.ring(x1)) == 1 + X2 ||
          s1[1] == R(1 + x2) // R(x1)

    s12 = mutate(s1, 2)
    s2  = mutate(s, 2)
    # μ₂ on original: x2' = (x1 + 1)/x2
    @test s2[2] == (R(x1) + 1) // R(x2)
end

@testset "backward compat — mutate(Seed, sequence) unchanged" begin
    q = Quiver(:A, 3)
    s = Seed(q)
    # Apply the full 14-seed BFS starting from s; check we get 14 distinct seeds.
    seen  = Set{Any}()
    stack = [s]
    while !isempty(stack)
        si = pop!(stack)
        key = Tuple(sort([denominator_vector(si, k) for k in 1:3]))
        key ∈ seen && continue
        push!(seen, key)
        for k in 1:3; push!(stack, mutate(si, k)); end
    end
    @test length(seen) == 14
end

# ─── Rational y-dynamics: A₂ ─────────────────────────────────────────────────

@testset "rational y-variables — A₂ initial values" begin
    es = extend(Seed(Quiver(:A, 2)))
    ys = y_variables(es; semifield=:rational)
    @test length(ys) == 2
    # Initial y-variables are just y1, y2 in Frac(ZZ[y1,y2])
    yring = es.coeffs.yring
    R     = base_ring(yring)
    y1, y2 = gens(R)
    @test ys[1] == yring(y1)
    @test ys[2] == yring(y2)
end

@testset "rational y-variables — A₂ after μ₁ (FZ oracle)" begin
    # A₂: B = [[0,1],[-1,0]].
    # FZ-IV mutation at k=1 with b_{kj} = B[k,j]: B[1,2] = 1, so
    # y'_1 = y_1^{-1}, y'_2 = y_2 * y_1^{[1]₊} * (1+y_1)^{-1} = y_1*y_2/(1+y_1)
    es1   = mutate(extend(Seed(Quiver(:A, 2))), 1)
    ys    = y_variables(es1; semifield=:rational)
    yring = es1.coeffs.yring
    R     = base_ring(yring)
    y1, y2 = gens(R)
    @test ys[1] == inv(yring(y1))
    @test ys[2] == yring(y1) * yring(y2) // (1 + yring(y1))
end

@testset "rational y-variables — μ_k twice is involutive on y" begin
    # After two mutations at the same vertex, the y-variables should return
    # to their original values (involutivity of mutation).
    es  = extend(Seed(Quiver(:A, 2)))
    es1 = mutate(es, 1)
    es11 = mutate(es1, 1)
    @test y_variables(es11; semifield=:rational) == y_variables(es; semifield=:rational)

    es2  = mutate(es, 2)
    es22 = mutate(es2, 2)
    @test y_variables(es22; semifield=:rational) == y_variables(es; semifield=:rational)
end

# ─── Tropical vs rational agreement ──────────────────────────────────────────

@testset "tropical-vs-rational agreement — A₂ full exchange graph" begin
    es0   = extend(Seed(Quiver(:A, 2)))
    n     = 2
    seen  = Set{Any}()
    stack = [es0]
    while !isempty(stack)
        e = pop!(stack)
        key = Tuple(sort([denominator_vector(e, k) for k in 1:n]))
        key ∈ seen && continue
        push!(seen, key)
        ys_rat  = y_variables(e; semifield=:rational)
        ys_trop = y_variables(e; semifield=:tropical)
        for j in 1:n
            trop_j = _tropicalize(ys_rat[j], n)
            @test trop_j == ys_trop[j]
        end
        for k in 1:n; push!(stack, mutate(e, k)); end
    end
    @test length(seen) == 5
end

@testset "tropical-vs-rational agreement — A₃ full exchange graph" begin
    es0   = extend(Seed(Quiver(:A, 3)))
    n     = 3
    seen  = Set{Any}()
    stack = [es0]
    while !isempty(stack)
        e = pop!(stack)
        key = Tuple(sort([denominator_vector(e, k) for k in 1:n]))
        key ∈ seen && continue
        push!(seen, key)
        ys_rat  = y_variables(e; semifield=:rational)
        ys_trop = y_variables(e; semifield=:tropical)
        for j in 1:n
            trop_j = _tropicalize(ys_rat[j], n)
            @test trop_j == ys_trop[j]
        end
        for k in 1:n; push!(stack, mutate(e, k)); end
    end
    @test length(seen) == 14
end

# ─── Tropical duality: C = (Gᵀ)⁻¹ ────────────────────────────────────────────
#
# For skew-symmetric B, the C- and G-matrices of any seed satisfy Gᵀ·C = I
# (Nakanishi–Zelevinsky tropical duality).  This pins down the convention of
# _mutate_C against the independently-computed G-matrix: a Langlands-dual
# (Bᵀ-pattern) C-matrix fails this identity while still passing the
# tropical-agreement tests above.

@testset "tropical duality GᵀC = I — A₂ and A₃ exchange graphs" begin
    for (type_rank, n_seeds) in ((2, 5), (3, 14))
        n     = type_rank
        es0   = extend(Seed(Quiver(:A, n)))
        seen  = Set{Any}()
        stack = [es0]
        I_n   = [i == j ? 1 : 0 for i in 1:n, j in 1:n]
        while !isempty(stack)
            e = pop!(stack)
            key = Tuple(sort([denominator_vector(e, k) for k in 1:n]))
            key ∈ seen && continue
            push!(seen, key)
            @test transpose(gmatrix(e)) * cmatrix(e) == I_n
            for k in 1:n; push!(stack, mutate(e, k)); end
        end
        @test length(seen) == n_seeds
    end
end

# ─── g-vector = principal grading ────────────────────────────────────────────

@testset "g-vector equals principal grading degree — A₂" begin
    # g-vector of x_k = degree of x_k under deg(xᵢ) = eᵢ, deg(yⱼ) = -B₀[:,j].
    # Oracle values: G[:,k] = g_vector(s,k) = column k of G-matrix (FZ4 convention).
    q  = Quiver(:A, 2)
    s0 = Seed(q)
    es = extend(s0)

    # G = I₂ initially
    @test g_vector(es, 1) == [1, 0]
    @test g_vector(es, 2) == [0, 1]

    es1 = mutate(es, 1)
    # G after [1] == [-1 0; 1 1] → G[:,1]=[-1,1], G[:,2]=[0,1]
    @test g_vector(es1, 1) == [-1, 1]
    @test g_vector(es1, 2) == [0, 1]

    es12 = mutate(es1, 2)
    # G after [1,2] == [-1 -1; 1 0] → G[:,1]=[-1,1], G[:,2]=[-1,0]
    @test g_vector(es12, 1) == [-1, 1]
    @test g_vector(es12, 2) == [-1, 0]

    es121 = mutate(es12, 1)
    # G after [1,2,1] == [0 -1; -1 0] → G[:,1]=[0,-1], G[:,2]=[-1,0]
    @test g_vector(es121, 1) == [0, -1]
    @test g_vector(es121, 2) == [-1, 0]

    es1212 = mutate(es121, 2)
    # G after [1,2,1,2] == [0 1; -1 0] → G[:,1]=[0,-1], G[:,2]=[1,0]
    @test g_vector(es1212, 1) == [0, -1]
    @test g_vector(es1212, 2) == [1, 0]
end

# ─── Separation formula round-trip ───────────────────────────────────────────

@testset "separation formula round-trip — A₂" begin
    q  = Quiver(:A, 2)
    s0 = Seed(q)

    # Check over all 5 seeds of A₂ (BFS)
    seen  = Set{Any}()
    stack_triv = [s0]
    stack_prin = [extend(s0)]

    while !isempty(stack_triv)
        s  = pop!(stack_triv)
        es = pop!(stack_prin)
        key = Tuple(sort([denominator_vector(s, k) for k in 1:2]))
        key ∈ seen && continue
        push!(seen, key)

        for k in 1:2
            sep  = separation_formula_trivial(es, k)
            # The separation formula result lives in Frac(ZZ[x1,x2]) built fresh,
            # while s.cluster[k] lives in the original ring.
            # Compare as strings (since they may be in different ring instances).
            @test string(sep) == string(s[k])
        end

        for kk in 1:2
            push!(stack_triv, mutate(s, kk))
            push!(stack_prin, mutate(es, kk))
        end
    end
    @test length(seen) == 5
end

@testset "separation formula round-trip — A₃" begin
    q  = Quiver(:A, 3)
    s0 = Seed(q)

    seen  = Set{Any}()
    stack_triv = [s0]
    stack_prin = [extend(s0)]

    while !isempty(stack_triv)
        s  = pop!(stack_triv)
        es = pop!(stack_prin)
        key = Tuple(sort([denominator_vector(s, k) for k in 1:3]))
        key ∈ seen && continue
        push!(seen, key)

        for k in 1:3
            sep = separation_formula_trivial(es, k)
            @test string(sep) == string(s[k])
        end

        for kk in 1:3
            push!(stack_triv, mutate(s, kk))
            push!(stack_prin, mutate(es, kk))
        end
    end
    @test length(seen) == 14
end

# ─── Frozen vertices: g-vectors error loudly, c/y-dynamics still work ────────

@testset "frozen vertices — gmatrix/separation error, c- and y-data available" begin
    # Minimal repro from the review: 1 mutable + 1 frozen, B = [0 1; -1 0].
    # x₁' = (1 + x₂)/x₁ contains the frozen variable, so x₁'/F₁(ŷ) is not a
    # Laurent monomial and no g-vector is defined.
    q  = Quiver([0 1; -1 0], 1)
    es = extend(Seed(q))

    @test cvectors(es) == [[1]]
    @test_throws ClusterAlgebraError gmatrix(es)
    @test_throws ClusterAlgebraError separation_formula(es, 1)

    es1 = mutate(es, 1)
    @test cvectors(es1) == [[-1]]
    ys = y_variables(es1; semifield=:rational)
    @test length(ys) == 1
    y1 = gens(base_ring(es1.coeffs.yring))[1]
    @test ys[1] == inv(es1.coeffs.yring(y1))
    @test_throws ClusterAlgebraError gmatrix(es1)
end

# ─── Sign coherence over B₂ and D₄ ──────────────────────────────────────────

@testset "sign coherence — B₂ full exchange graph" begin
    # B₂: exchange matrix for B₂ with d = [2, 1]
    # B = [[0, 1], [-2, 0]]
    B_B2 = [0 1; -2 0]
    q  = Quiver(B_B2, 2, [2, 1])
    es = extend(Seed(q))
    seen  = Set{Any}()
    stack = [es]
    count = 0
    while !isempty(stack)
        e = pop!(stack)
        key = Tuple(sort([denominator_vector(e, k) for k in 1:e.quiver.n_mutable]))
        key ∈ seen && continue
        push!(seen, key)
        count += 1
        @test is_sign_coherent(e)
        for k in 1:e.quiver.n_mutable
            push!(stack, mutate(e, k))
        end
    end
    @test count == 6   # B₂ has 6 distinct seeds
end

@testset "sign coherence — D₄ full exchange graph" begin
    q  = Quiver(:D, 4)
    es = extend(Seed(q))
    seen  = Set{Any}()
    stack = [es]
    count = 0
    while !isempty(stack)
        e = pop!(stack)
        key = Tuple(sort([denominator_vector(e, k) for k in 1:e.quiver.n_mutable]))
        key ∈ seen && continue
        push!(seen, key)
        count += 1
        @test is_sign_coherent(e)
        for k in 1:e.quiver.n_mutable
            push!(stack, mutate(e, k))
        end
    end
    @test count == 50   # D₄ has 50 distinct seeds
end

@testset "ExtendedCoefficients (geometric type)" begin
    # A₂ with principal (identity) frozen block: mutable 1,2; frozen 3,4.
    B = [ 0  1 -1  0;
         -1  0  0 -1;
          1  0  0  0;
          0  1  0  0]
    q  = Quiver(B, 2)                      # n_mutable = 2, n_frozen = 2
    s0 = Seed(q)
    ge = extend_geometric(s0)
    @test ge isa Seed{ExtendedCoefficients}

    R = ge.ring
    x = ge.cluster                          # [x1, x2, x3, x4]

    # Coefficients read off the frozen block: initial C = I₂.
    @test y_variables(ge; semifield = :tropical) == [[1, 0], [0, 1]]
    @test cvectors(ge) == [[1, 0], [0, 1]]
    @test cmatrix(ge) == [1 0; 0 1]
    # Geometric y_j = ∏ frozen x^{B}: y₁ = x₃, y₂ = x₄.
    @test y_variables(ge; semifield = :geometric) == [x[3], x[4]]
    # Full ŷ over all vertices: ŷ₁ = x₃/x₂, ŷ₂ = x₁·x₄.
    @test y_hat(ge) == [x[3] // x[2], x[1] * x[4]]

    # Mutate at the mutable vertex 1: cluster exchange incl. frozen sides.
    g1 = mutate(ge, 1)
    @test g1 isa Seed{ExtendedCoefficients}
    @test g1[1] == (R(gens(base_ring(R))[2]) + R(gens(base_ring(R))[3])) // x[1]  # (x₂+x₃)/x₁

    # Coefficient consistency bridge: the geometric c-vectors (frozen block under
    # matrix mutation) equal the principal C-matrix of the underlying A₂.
    prin1 = mutate(extend(Seed(Quiver(:A, 2))), 1)
    @test y_variables(g1; semifield = :tropical) == cvectors(prin1)
    @test cvectors(g1) == [[-1, 0], [1, 1]]

    # Involution on cluster + quiver (mutation path aside).
    g11 = mutate(g1, 1)
    @test g11.cluster == ge.cluster
    @test g11.quiver  == ge.quiver

    # Guards.
    @test_throws InvalidArgument extend_geometric(Seed(Quiver(:A, 2)))   # no frozen
    @test_throws FrozenVertexMutation mutate(ge, 3)                      # frozen vertex
    @test_throws InvalidArgument y_variables(ge; semifield = :bogus)
end
