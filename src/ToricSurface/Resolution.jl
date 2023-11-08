###########################################################################
# This file implements low-level functionality around toric resolutions
# of singularities. Some of this (e.g. Hirzebruch-Jung continued fractions)
# is available in Oscar through polymake. We implement it here in pure
# Julia for increased performance.
###########################################################################

@doc raw"""
    cls_cone_normal_form(A :: ZZMatrix)

Bring a two-dimensional cone into normal form in the sense of Cox, Little,
Schenck Proposition 10.1.1. The result is a triple (d, k, M), where d and k are
the parameters of the cone and M is a 2x2 integral matrix such that M * A = [0
d ; 1 -k]

"""
function cls_cone_normal_form(A :: ZZMatrix)
    u11, u21, u12, u22 = A[1,1], A[2,1], A[1,2], A[2,2]
    # Get integers x.y such that x * u11 + y * u21 == 1 
    _, x, y = gcdx(u11, u21)
    sg, d, l = sign(det(A)), abs(det(A)), x * u12 + y * u22
    s, k = ((l + mod(-l, d)) ÷ d, mod(-l, d))
    M = ZZ[-sg*u21 sg*u11 ; x+sg*s*u21 y-sg*s*u11]
    return (d, k, M)
end


@doc raw"""
    hirzebruch_jung(x :: T, y :: T)

Return the Hirzebruch-Jung continued fraction associated to `x/y`.

"""
function hirzebruch_jung(x :: T, y :: T) where {T <: IntegerUnion}
    res = T[]
    while y > 0
        push!(res, (x + mod(-x, y)) ÷ y)
        x, y = y, mod(-x,y)
    end
    return res;
end    


@doc raw"""
    hilbert_basis_2D(A :: ZZMatrix)

Return the hilbert basis of a two-dimensional cone spanned by the columns of
`A`. 

"""
function hilbert_basis_2D(A :: ZZMatrix)
    d, k, M = cls_cone_normal_form(A)
    hj = hirzebruch_jung(d, k)

    x, y = 0, 1
    a, b = -1, 0
    res = []
    for z in hj
        push!(res, inv(M) * ZZ[y ; -b])
        x, y = y, z * y - x
        a, b = b, z * b - a
    end
    
    return res
end


@doc raw"""
    hilbert_basis_2D(v1 :: Vector{T}, v2 :: Vector{T})

Return the hilbert basis of a two-dimensional cone spanned by given integral
vectors `v1` and `v2`.

"""
function hilbert_basis_2D(v1 :: Vector{T}, v2 :: Vector{T}) where {T <: IntegerUnion}
    hb = hilbert_basis_2D(matrix(ZZ, [v1[1] v2[1] ; v1[2] v2[2]]))
    return [[T(v[1,1]), T(v[2,1])] for v in hb]
end

function primitivize(v :: Vector) 
    v = Int.(lcm(denominator.(v)) .* v)
    return v .÷ gcd(v)
end

@doc raw"""
    hilbert_basis_2D(v1 :: Vector, v2 :: Vector)

Return the hilbert basis of a two-dimensional cone spanned by given vectors
`v1` and `v2`.

"""
hilbert_basis_2D(v1 :: Vector, v2 :: Vector) = hilbert_basis_2D(primitivize(v1), primitivize(v2))


function intersect_lines_2D(M1 :: QQMatrix, M2 :: QQMatrix)
    A = QQ[M1[2,2]-M1[2,1] M1[1,1]-M1[1,2] ; M2[2,2]-M2[2,1] M2[1,1]-M2[1,2]]
    b = QQ[det(M1) ; det(M2)]
    @req det(A) ≠ 0 "No unique intersection"
    (p, z) = can_solve_with_solution(A,b)
    @req p "no solution exists"
    return z
end


@doc raw"""
    intersect_lines_2D(v1 :: Vector, v2 :: Vector, w1 :: Vector, w2 :: Vector)

Given rational vectors `v1`, `v2`, `w1` and `w2` in two-dimensional space,
return the unique intersection point of the line passing through `v1` and `v2`
and the line passing through `w1` and `w2`, if it exists.

"""
function intersect_lines_2D(v1 :: Vector, v2 :: Vector, w1 :: Vector, w2 :: Vector)
    z = intersect_lines_2D(matrix(QQ, [v1 v2]), matrix(QQ, [w1 w2]))
    return [z[1,1], z[2,1]]
end


@doc raw"""
    norm_ratio(v1 :: Vector, v2 :: Vector)   

Return the ratio between the euclidian norms of two rational vectors in
two-dimensional space.

"""
norm_ratio(v1 :: Vector, v2 :: Vector) = v2[2] == 0 ? v1[1] // v2[1] : v1[2] // v2[2]


@doc raw"""
    discrepancy(v1 :: Vector, v2 :: Vector, w :: Vector)

Return the discrepancy of the toric refinement given by introducing the ray
with primitive generator `w` into the cone spanned by the vectors `v1` and
`v2`.

"""
discrepancy(v1 :: Vector, v2 :: Vector, w :: Vector) = 
Rational(norm_ratio(w, intersect_lines_2D(v1, v2, [0, 0], w)) - 1)


@doc raw"""
    toric_affine_surface_resolution(v1 :: Vector, v2 :: Vector)

Return the toric resolution of an affine toric surface given by a cone with
primitive generators `v1` and `v2`. The result is a pair, where the first entry
is the list of rays that need to be inserted and the second is a list of the
discrepancies of the associated exceptional divisors.

"""
function toric_affine_surface_resolution(v1 :: Vector, v2 :: Vector)
    ex_rays = hilbert_basis_2D(v1, v2)
    discrepancies = [discrepancy(v1, v2, w) for w in ex_rays]
    return (ex_rays, discrepancies)
end


