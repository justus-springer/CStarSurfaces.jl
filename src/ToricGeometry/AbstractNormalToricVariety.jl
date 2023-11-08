
@doc raw"""
    normal_toric_variety(P :: MatElem; non_redundant :: Bool)

Construct the unique complete normal toric variety of Picard rank one by
providing the rays of the corresponding fan.

"""
function normal_toric_variety(P :: MatElem; non_redundant :: Bool)
    n, r = nrows(P), ncols(P)
    @req n == r + 1 "this constructor only works if the picard rank is one"
    cones = collect(powerset(1 : n, r, r))
    normal_toric_variety(P, cones, non_redundant=non_redundant)
end


@doc raw"""
    affine_toric_chart(X :: NormalToricVarietyType, c :: Vector{Int64})

Return the affine toric chart associated to a cone `c`, given by a list of ray indices.

"""
function affine_toric_chart(X :: NormalToricVarietyType, c :: Vector{Int64})
    P = gen_matrix(X)
    return normal_toric_variety(transpose(P[:, c]), [collect(1:length(c))]; non_redundant=true)
end

@doc raw"""
    remove_torusfactor(X :: NormalToricVarietyType)

Given a normal toric variety with torusfactor $X = X' \times \mathbb{T}^s$,
return $X'$. In terms of fans, this amounts to changing the ambient lattice to
the linear hull of the rays of the fan.

"""
function remove_torusfactor(X :: NormalToricVarietyType)
    P = gen_matrix(X)
    r, A, _ = rcef_rational(P)
    A = A[:, 1 : r]
    # Use `invoke` here to dispatch on the more generic `solve_rational` from 
    # AbstractAlgebra, which doesn't require the first argument to be a square matrix
    new_P, _ = invoke(solve_rational, Tuple{MatElem{ZZRingElem}, MatElem{ZZRingElem}}, A, P)
    normal_toric_variety(transpose(new_P), maximal_cones_indices(X); non_redundant=true)
end


#################################################
# Cones in the rational class group
#################################################


@doc raw"""
    effective_cone(X :: NormalToricVarietyType)

Return the cone of effective divisor classes in the rational vector space
associated to the divisor class group of a normal toric variety.

"""
@attr effective_cone(X :: NormalToricVarietyType) =
positive_hull(transpose(degree_matrix_free_part(X)))


@doc raw"""
    moving_cone(X :: NormalToricVarietyType)

Return the cone of movable divisor classes in the rational vector space
associated to the divisor class group of a normal toric variety.

"""
@attr function moving_cone(X :: NormalToricVarietyType)
    Q0 = degree_matrix_free_part(X)
    n = nrays(X)
    intersect([positive_hull(transpose(Q0[:,c])) for c in powerset(1:n, n-1, n-1)])
end


@doc raw"""
    semiample_cone(X :: NormalToricVarietyType)

Return the cone of semi-ample divisor classes in the rational vector space
associated to the divisor class group of a normal toric variety.

"""
@attr function semiample_cone(X :: NormalToricVarietyType)
    Q0 = degree_matrix_free_part(X)
    intersect([positive_hull(transpose(Q0[:,c])) for c in covering_collection(X)])
end

    

