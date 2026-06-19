# Conway–Coxeter SL₂ frieze patterns (Track A, Group D).

"""
    Frieze

An SL₂ Conway–Coxeter frieze pattern associated to a triangulated n-gon.

Fields:
- `n::Int` — polygon size (n ≥ 4); width = n − 3 interior rows.
- `quiddity::Vector{Int}` — quiddity sequence (c₁, …, cₙ), positive integers.
- `entries::Matrix{Int}` — (n+1) × n integer matrix, 1-indexed.
  Row j stores mathematical frieze row j−1:
  row 1 = zeros, row 2 = ones, rows 3..n−1 = interior, row n = ones, row n+1 = zeros.
"""
struct Frieze
    n::Int
    quiddity::Vector{Int}
    entries::Matrix{Int}
end

# ─── Internal recurrence ──────────────────────────────────────────────────────
#
# Diamond rule (0-indexed math): f(r,i)*f(r,i+1) - f(r-1,i)*f(r+1,i+1) = 1
# Rearranged to fill row r+1 from rows r and r-1:
#   f(r+1, i+1) = (f(r,i)*f(r,i+1) - 1) / f(r-1,i)
#
# In 1-indexed Julia (F[j,k] = f(j-1, k-1)):
#   F[r+1, mod1(i+2, n)] = (F[r, i+1]*F[r, mod1(i+2,n)] - 1) / F[r-1, i+1]

function _build_frieze_entries(quiddity::Vector{Int})
    n = length(quiddity)
    F = zeros(Int, n + 1, n)
    F[2, :] .= 1
    F[3, :] .= quiddity
    for r in 3:n
        for i in 0:n-1
            a  = i + 1
            b  = mod1(i + 2, n)
            num = F[r, a] * F[r, b] - 1
            den = F[r - 1, a]
            if den == 0 || num % den != 0
                throw(InvalidArgument(
                    "frieze recurrence: non-integer result $num/$den " *
                    "at matrix row $(r+1), column $b; " *
                    "check that the quiddity sequence is valid"))
            end
            F[r + 1, b] = num ÷ den
        end
    end
    return F
end

# ─── Public constructors ──────────────────────────────────────────────────────

"""
    frieze(n::Int) -> Frieze

Return the Conway–Coxeter SL₂ frieze of the n-gon (n ≥ 4) using the fan
triangulation (all diagonals from vertex 1).

The quiddity sequence is `(n−2, 1, 2, 2, …, 2, 1)`.
The resulting frieze has width `n − 3` interior rows, all positive integers.
"""
function frieze(n::Int)
    n >= 4 || throw(InvalidArgument(
        "frieze requires a polygon with n ≥ 4 vertices, got n = $n"))
    c = fill(2, n)
    c[1] = n - 2
    c[2] = 1
    c[n] = 1
    return Frieze(n, c, _build_frieze_entries(c))
end

"""
    frieze(quiddity::AbstractVector{<:Integer}) -> Frieze

Return the Conway–Coxeter SL₂ frieze with the given quiddity sequence.
The sequence must have length ≥ 4 and consist of positive integers.
"""
function frieze(quiddity::AbstractVector{<:Integer})
    n = length(quiddity)
    n >= 4 || throw(InvalidArgument(
        "quiddity sequence must have length ≥ 4, got $n"))
    all(>(0), quiddity) || throw(InvalidArgument(
        "quiddity sequence must consist of positive integers"))
    c = Vector{Int}(quiddity)
    return Frieze(n, c, _build_frieze_entries(c))
end

# ─── Validation ───────────────────────────────────────────────────────────────

"""
    AbstractAlgebra.is_valid(f::Frieze) -> Bool

Return `true` if every unit diamond in `f` satisfies the unimodular rule
`f(r,i)·f(r,i+1) − f(r−1,i)·f(r+1,i+1) = 1`.

Friezes produced by [`frieze`](@ref) always satisfy this; use `is_valid` to
verify externally supplied patterns.
"""
function AbstractAlgebra.is_valid(f::Frieze)
    n = f.n
    F = f.entries
    # Check diamonds between Julia rows j and j+1 for j = 2..n.
    # In 1-indexed storage, the rule is:
    #   F[j, col] * F[j, nc] - F[j-1, col] * F[j+1, nc] == 1
    # where nc = mod1(col+1, n).
    for j in 2:n
        for col in 1:n
            nc = mod1(col + 1, n)
            F[j, col] * F[j, nc] - F[j - 1, col] * F[j + 1, nc] == 1 || return false
        end
    end
    return true
end

# ─── Display ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, f::Frieze)
    n = f.n
    println(io, "Frieze(n = $n, width = $(n - 3))")
    println(io, "Quiddity: $(f.quiddity)")
    w = maximum(ndigits(x; pad=1) + (x < 0 ? 1 : 0) for x in f.entries; init=1) + 1
    for r in 1:n+1
        println(io, join(lpad(f.entries[r, i], w) for i in 1:n))
    end
end
