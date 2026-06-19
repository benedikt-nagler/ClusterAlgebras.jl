_show_latex(io::IO, x) = try
    show(io, MIME"text/latex"(), x)
catch
    print(io, x)
end

function Base.show(io::IO, q::Quiver)
    n = q.n_mutable + q.n_frozen
    println(io, "Quiver with $n vertices ($(q.n_mutable) mutable, $(q.n_frozen) frozen)")
    println(io, "Exchange matrix B:")
    w = maximum(ndigits(x; pad=1) + (x < 0 ? 1 : 0) for x in q.B; init=1) + 1
    for i in 1:n
        print(io, " ")
        for j in 1:n
            print(io, lpad(q.B[i, j], w))
        end
        tag = i <= q.n_mutable ? "" : "  [frozen]"
        println(io, "   ($(q.labels[i])$tag)")
    end
end

function Base.show(io::IO, ::MIME"text/latex", q::Quiver)
    n = q.n_mutable + q.n_frozen
    print(io, "\\begin{pmatrix}")
    for i in 1:n
        print(io, join(string.(q.B[i, :]), " & "))
        i < n && print(io, " \\\\")
    end
    print(io, "\\end{pmatrix}")
end

# ─── Text display for Seed{TrivialCoefficients} ───────────────────────────────

function Base.show(io::IO, s::Seed{TrivialCoefficients})
    show(io, s.quiver)
    println(io, "Cluster variables:")
    for (i, x) in enumerate(s.cluster)
        tag = i > s.quiver.n_mutable ? "  [frozen]" : ""
        println(io, "  [$i]$tag: $x")
    end
    isempty(s.mutation_path) || println(io, "Mutation path: $(s.mutation_path)")
end

# ─── LaTeX display for Seed{TrivialCoefficients} ──────────────────────────────

function Base.show(io::IO, ::MIME"text/latex", s::Seed{TrivialCoefficients})
    n     = length(s.cluster)
    n_mut = s.quiver.n_mutable
    B     = s.quiver.B

    println(io, "\\begin{aligned}")

    # cluster row vector
    print(io, "  \\mathbf{x} &= \\begin{pmatrix} ")
    for i in 1:n
        _show_latex(io, s.cluster[i])
        i > n_mut && print(io, "^{\\ast}")   # mark frozen slots
        i < n     && print(io, " & ")
    end
    println(io, " \\end{pmatrix} \\\\[6pt]")

    # exchange matrix
    print(io, "  B &= \\begin{pmatrix} ")
    for i in 1:n
        print(io, join(string.(B[i, :]), " & "))
        i < n && print(io, " \\\\ ")
    end
    print(io, " \\end{pmatrix}")

    # optional mutation path
    if !isempty(s.mutation_path)
        path_str = join(string.(s.mutation_path), ",\\,")
        print(io, " \\\\[4pt]\n  \\text{mutations} &= ($(path_str))")
    end

    print(io, "\n\\end{aligned}")
end
