# Denominator vectors (d-vectors) of cluster variables (Track A: B1).
#
# Convention (Fomin–Zelevinsky): the d-vector of cluster variable x is the
# exponent vector of the denominator monomial when x is written as a reduced
# Laurent polynomial in the *initial* cluster.  Initial variables are assigned
# d-vector −e_i (negative simple root) by the FZ convention.
#
# In finite type A_n the multiset of d-vectors over all cluster variables equals
# the set of almost-positive roots (a strong correctness oracle).

"""
    denominator_vector(s::Seed, k::Int) -> Vector{Int}

Return the denominator vector (d-vector) of the `k`-th cluster variable of `s`,
expressed in the initial cluster variables.

Initial cluster variable `x_i` (denominator = 1, numerator = `x_i`) gets the
conventional d-vector `−eᵢ`.  For all other cluster variables the d-vector is
the exponent vector of the (monomial) denominator of the reduced Laurent
expression.
"""
function denominator_vector(s::Seed, k::Int)
    n_total = s.quiver.n_mutable + s.quiver.n_frozen
    1 <= k <= n_total || throw(InvalidVertex(k, n_total))

    xk  = s.cluster[k]
    den = denominator(xk)
    num = numerator(xk)

    # Initial variable: den = 1 and num is a single generator x_i
    if isone(den) && length(num) == 1
        ev = first(exponent_vectors(num))
        sum(ev) == 1 && return -ev   # d-vector of x_i is −eᵢ
    end

    # Non-initial: denominator is a monomial by the Laurent phenomenon
    evs = collect(exponent_vectors(den))
    length(evs) == 1 || throw(InvalidArgument(
        "cluster variable $k has a non-monomial denominator; " *
        "this should not happen for cluster variables (Laurent phenomenon violated?)"))
    return evs[1]
end
