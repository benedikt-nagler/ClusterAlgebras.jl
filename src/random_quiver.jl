"""
    random_quiver([rng], rank::Int; n_frozen=0, max_arrows=3, density=0.5, model=:uniform) → Quiver

A random skew-symmetric quiver on `rank` mutable vertices (plus `n_frozen` frozen
ones). Pass an `AbstractRNG` to make a dataset reproducible from its seed.

Models:
- `:uniform` - each entry above the diagonal is nonzero with probability
  `density`, and if nonzero is drawn uniformly from `±(1:max_arrows)`; the
  entries below the diagonal are fixed by skew-symmetry. This is the
  distribution used by the quiver-classification literature.
- `:acyclic` - as `:uniform`, then all arrows are oriented from lower to higher
  vertex index, which makes the quiver acyclic by construction.

Only skew-symmetric quivers are generated (all symmetrizers `1`); there is no
natural uniform measure on skew-symmetrizable ones without first fixing `d`.

# Example
```julia
julia> using Random

julia> q = random_quiver(MersenneTwister(1), 4);

julia> nvertices(q)
4

julia> random_quiver(MersenneTwister(1), 4) == q   # reproducible from the seed
true
```
"""
function random_quiver(rng::AbstractRNG, rank::Int; n_frozen::Int = 0,
                       max_arrows::Int = 3, density::Float64 = 0.5,
                       model::Symbol = :uniform)
    rank >= 0 ||
        throw(InvalidArgument("rank must be non-negative, got $rank"))
    n_frozen >= 0 ||
        throw(InvalidArgument("n_frozen must be non-negative, got $n_frozen"))
    max_arrows >= 1 ||
        throw(InvalidArgument("max_arrows must be at least 1, got $max_arrows"))
    0 <= density <= 1 ||
        throw(InvalidArgument("density must lie in [0, 1], got $density"))
    model in (:uniform, :acyclic) ||
        throw(InvalidArgument("model must be :uniform or :acyclic, got :$model"))

    n = rank + n_frozen
    B = zeros(Int, n, n)
    for i in 1:n, j in (i + 1):n
        rand(rng) < density || continue
        w = rand(rng, 1:max_arrows)
        # :uniform picks an orientation; :acyclic always points i → j, and since
        # i < j every arrow follows the vertex order, so no cycle can close.
        if model === :uniform && rand(rng, Bool)
            w = -w
        end
        B[i, j] =  w
        B[j, i] = -w
    end
    # Arrows between two frozen vertices are not part of the data of a quiver
    # with frozen vertices - mutation never reads that block.
    for i in (rank + 1):n, j in (rank + 1):n
        B[i, j] = 0
    end
    return Quiver(B, rank)
end

random_quiver(rank::Int; kwargs...) = random_quiver(Random.default_rng(), rank; kwargs...)

"""
    random_mutate([rng], q::Quiver, n_steps::Int) → Quiver
    random_mutate([rng], s::Seed, n_steps::Int) → Seed

Apply `n_steps` mutations at uniformly random mutable vertices. Consecutive
repeats are avoided (mutation is an involution, so `μₖμₖ = id` would waste a
step), which makes this a random walk on the exchange graph that never
immediately backtracks.

Every result is mutation-equivalent to the input by construction - useful for
sampling within a known mutation class, and for generating quivers labeled by
the class they were walked from.

# Example
```julia
julia> using Random

julia> q = random_mutate(MersenneTwister(1), Quiver(:A, 3), 10);

julia> is_finite_type(q)   # mutation class is a mutation invariant
true
```
"""
function random_mutate(rng::AbstractRNG, x::Union{Quiver, AbstractSeed}, n_steps::Int)
    n_steps >= 0 ||
        throw(InvalidArgument("n_steps must be non-negative, got $n_steps"))
    n_mutable = x isa Quiver ? x.n_mutable : x.quiver.n_mutable
    n_mutable >= 1 ||
        throw(InvalidArgument("cannot mutate a quiver with no mutable vertices"))

    last_k = 0
    for _ in 1:n_steps
        k = rand(rng, 1:n_mutable)
        while n_mutable > 1 && k == last_k
            k = rand(rng, 1:n_mutable)
        end
        x = mutate(x, k)
        last_k = k
    end
    return x
end

random_mutate(x::Union{Quiver, AbstractSeed}, n_steps::Int) =
    random_mutate(Random.default_rng(), x, n_steps)
