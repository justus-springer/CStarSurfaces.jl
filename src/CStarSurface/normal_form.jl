#################################################
# Normal form of C-star surface
#################################################

@attr _beta_plus(X :: CStarSurface) = map(ms -> map(m -> m - floor(maximum(ms)), ms), slopes(X))
@attr _beta_minus(X :: CStarSurface) = map(ms -> map(m -> ceil(minimum(ms)) - m, ms), slopes(X))

###################################### NOTE ##########################################
# In a prior implementation of the orientation of P-matrices, beta_plus and beta_minus
# where sorted in a non-standard way. The choice of ordering is somewhat arbitrary,
# but different choices lead to different results for the orientation and hence the 
# normal form in some cases. For backward compatibilily with existing values in 
# the ldp-database, we use the same ordering here that was used before.
######################################################################################
function _custom_lt_vectors(a :: Vector{T}, b :: Vector{T}) where {T} 
    length(a) > length(b) && return true
    length(a) < length(b) && return false
    return a > b
end

@attr _beta_plus_sorted(X :: CStarSurface) = sort(map(xs -> sort(xs, rev=true), _beta_plus(X)), lt = _custom_lt_vectors)
@attr _beta_minus_sorted(X :: CStarSurface) = sort(map(xs -> sort(xs, rev=true), _beta_minus(X)), lt = _custom_lt_vectors)

@attr orientation(X :: CStarSurface{PE}) = 1
@attr orientation(X :: CStarSurface{EP}) = -1
@attr function orientation(X :: CStarSurface)
    if _m_plus(X) > _m_minus(X)
        return 1
    elseif _m_plus(X) < _m_minus(X)
        return -1
    end
    if _beta_plus_sorted(X) > _beta_minus_sorted(X)
        return 1
    elseif _beta_plus_sorted(X) < _beta_minus_sorted(X)
        return -1
    end
    return 0
end

# Bring a C-star surface `X` into normal form. Returns a tuple (Y,p), where Y is in 
# normal form such that `p(X) = Y`.
function normal_form(X :: CStarSurface)
    Y = deepcopy(X)
    # Step 1: Bring X into non-negative orientation
    p1 = InvertLastRow(orientation(X) < 0 ? -1 : 1)
    Y = p1(Y)

    # Step 2: Sort the rays in each block by slope
    p2 = PermutationOfRays(map(ms -> inv(perm(sortperm(ms; rev = true))), _beta_plus(Y)))
    Y = p2(Y)

    # Step 3: Sort the blocks
    p3 = PermutationOfBlocks(inv(perm(sortperm(Vector(_beta_plus(Y)), lt = _custom_lt_vectors))))
    Y = p3(Y)

    # Step 4: Apply row operations to achieve floor(maximum(m_{i1}, ...m_{in_i})) = 0
    # for all i = 1, ..., r
    p4 = AdmissibleRowOperation(map(ms -> -floor(Int, maximum(ms)), slopes(Y)[1:end]))
    Y = p4(Y)

    return (Y, normalize_admissible_operation(p1 * p2 * p3 * p4))
end

is_normal_form(X :: CStarSurface) = normal_form(X)[1] == X

# Check whether two C-star surfaces are isomorphic to each other
# Returns a pair consisting of a boolean and either an admissible operation turning
# `X` into `Y` or `nothing`.
function are_isomorphic(X :: CStarSurface, Y :: CStarSurface)
    (X_norm, X_op) = normal_form(X)
    (Y_norm, Y_op) = normal_form(Y)
    if X_norm == Y_norm
        return (true, normalize_admissible_operation(X_op * inv(Y_op)))
    else
        return (false, nothing)
    end
end



