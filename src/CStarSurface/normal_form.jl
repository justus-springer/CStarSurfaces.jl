
@doc raw"""
    beta_plus(X :: CStarSurface)

Return a `β :: DoubleVector` with entries `β[i][j] == slopes(X)[i][j] -
floor(maximum(slopes(X)[i]))`

"""
@attr beta_plus(X :: CStarSurface) = map(ms -> map(m -> m - floor(maximum(ms)), ms), slopes(X))


@doc raw"""
    beta_minus(X :: CStarSurface)

Return a `β :: DoubleVector` with entries `β[i][j] ==
ceil(minimum(slopes(X)[i])) - slopes(X)[i][j]`.

"""
@attr beta_minus(X :: CStarSurface) = map(ms -> map(m -> ceil(minimum(ms)) - m, ms), slopes(X))

############################### Implementation note ###############################
# In a prior implementation of the orientation of P-matrices, beta_plus and
# beta_minus were sorted in a non-standard way. The choice of ordering is
# somewhat arbitrary, but different choices lead to different results for the
# orientation and hence the normal form in some cases. For backward
# compatibilily with existing values in the ldp-database, we use the same
# ordering here that was used before.
###################################################################################

function _custom_lt_vectors(a :: Vector{T}, b :: Vector{T}) where {T} 
    length(a) > length(b) && return true 
    length(a) < length(b) && return false 
    return a > b
end


@doc raw"""
    beta_plus_sorted(X :: CStarSurface)

Return the sorted `beta_plus(X)`. Each vector in `beta_plus(X)` are
individually sorted and the vectors themselves are sorted by first by size and
then lexicographically.

"""
@attr beta_plus_sorted(X :: CStarSurface) = sort(map(xs -> sort(xs, rev=true), beta_plus(X)), lt = _custom_lt_vectors)


@doc raw"""
    beta_minus_sorted(X :: CStarSurface)

Return the sorted `beta_minus(X)`. Each vector in `beta_minus(X)` are
individually sorted and the vectors themselves are sorted by first by size and
then lexicographically.

"""
@attr beta_minus_sorted(X :: CStarSurface) = sort(map(xs -> sort(xs, rev=true), beta_minus(X)), lt = _custom_lt_vectors)


@doc raw"""
    orientation(X :: CStarSurface)

Return the orientation of a C-star surface. This function takes the values `1`,
`0` or `-1`. Note that `orientation` is not an isomorphy invariant, as applying
`InvertLastRow(-1)` inverts the orientation of a C-star surface. All other
types of `AdmissibleOperation`'s leave the orientation invariant.

A C-star surface `X` has orientation `1`, if and only if one of the following
conditions hold:

1. `X.case == :pe`,
2. `X.case ∈ [:ee, :pp]` and `m_plus(X) > m_minus(X)`,
3. `X.case ∈ [:ee, :pp]` and `m_plus(X) == m_minus(X)` and `beta_plus_sorted(X) > beta_minus_sorted(X)`.

Similarly, `X` has orientation `-1` if and only if one of the following
conditions hold:

1. `X.case == :ep`,
2. `X.case ∈ [:ee, :pp]` and `m_plus(X) < m_minus(X)`,
3. `X.case ∈ [:ee, :pp]` and `m_plus(X) == m_minus(X)` and `beta_plus_sorted(X) < beta_minus_sorted(X)`.

The remaining case is that `X.case ∈ [:ee, :pp]` and `m_plus(X) == m_minus(X)`
and `beta_plus_sorted(X) == beta_minus_sorted(X)`, in which case `X` has
orientation `0`.

"""
@attr function orientation(X :: CStarSurface)
    X.case == :pe && return 1
    X.case == :ep && return -1
    m_plus(X) > m_minus(X) && return 1
    m_plus(X) < m_minus(X) && return -1
    beta_plus_sorted(X) > beta_minus_sorted(X) && return 1
    beta_plus_sorted(X) < beta_minus_sorted(X) && return -1
    return 0
end


@doc raw"""
    normal_form(X :: CStarSurface)

Compute the normal form of a C-star surface. Here, a C-star surface `X` is said
to be in normal form, if and only if the following properties hold:

1. `orientation(X) ≠ -1`,
2. `beta_plus(X) == beta_plus_sorted(X)`,
3. `0 ≤ X.d[i][1] < X.l[i][i]` for all `1 ≤ i ≤ r`, where `r+1 == nblocks(X)`.

The third condition can also be phrased as `floor(maximum(slopes(X)[i])) == 0`
for all `1 ≤ i ≤ r`.

The algorithm works by applying an `InvertLastRow` operation to achieve 1.,
then applying `PermutationOfRays` and `PermutationOfBlocks` operations to
achieve 2. and finally, applying an `AdmissibleRowOperation` to achieve 3.
Together, these properties are enough to ensure that `X` and `Y` are isomorphic
if and only if they have the same normal form.

This function then returns a pair `(Y, α)`, where `Y` is a C-star surface in
normal form and `α` is an admissible operation with `α(X) == Y`.

"""
@attr function normal_form(X :: CStarSurface)
    Y = deepcopy(X)
    # Step 1: Bring X into non-negative orientation
    p1 = InvertLastRow(orientation(X) < 0 ? -1 : 1)
    Y = p1(Y)

    # Step 2: Sort the rays in each block by slope
    p2 = PermutationOfRays(map(ms -> inv(perm(sortperm(ms; rev = true))), beta_plus(Y)))
    Y = p2(Y)

    # Step 3: Sort the blocks
    p3 = PermutationOfBlocks(inv(perm(sortperm(Vector(beta_plus(Y)), lt = _custom_lt_vectors))))
    Y = p3(Y)

    # Step 4: Apply row operations to achieve floor(maximum(m_{i1}, ...m_{in_i})) = 0
    # for all i = 1, ..., r
    p4 = AdmissibleRowOperation(map(ms -> -floor(Int, maximum(ms)), slopes(Y)[1:end]))
    Y = p4(Y)

    # Remember that Y is in normal form, so it doesn't have to be
    # computed again
    set_attribute!(Y, :is_normal_form, true)
    set_attribute!(Y, :normal_form, (Y, one(AdmissibleOperation))) 

    return (Y, normalize_admissible_operation(p1 * p2 * p3 * p4))
end


@doc raw"""
    is_normal_form(X :: CStarSurface) = normal_form(X)[1] == X

Check whether a C-star surface is in normal form, see the docstring of
`normal_form`.

"""
@attr is_normal_form(X :: CStarSurface) = normal_form(X)[1] == X


@doc raw"""
    are_isomorphic(X :: CStarSurface, Y :: CStarSurface)

Check whether two C-star surfaces are isomorphic to each other. This function
returns a pair, where the first entry is a boolean and the second entry is
either `nothing` or an admissible operation turning `X` into `Y`.

"""
function are_isomorphic(X :: CStarSurface, Y :: CStarSurface)
    (X_norm, X_op) = normal_form(X)
    (Y_norm, Y_op) = normal_form(Y)
    if X_norm == Y_norm
        return (true, normalize_admissible_operation(X_op * inv(Y_op)))
    else
        return (false, nothing)
    end
end



