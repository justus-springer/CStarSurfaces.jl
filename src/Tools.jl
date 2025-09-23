function discrepancies(d :: T, k :: T, t :: Rational{T}) where {T <: Integer}
    a0, a1 = 1//t, (k+t) // (d*t)
    res = Rational{T}[]
    push!(res, a0)
    for b ∈ hirzebruch_jung(d, k)
       push!(res, a1)
       a0, a1 = a1, b * a1 - a0
    end
    push!(res, a1)
    return res
end

function discrepancies(A :: Matrix2{T}, t :: Rational{T}) where {T <: Integer}
    d, k, _ = cls_cone_normal_form(A)
    return discrepancies(d, k, t)
end

@doc raw"""
discrepancies(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}, t :: Rational{T}) where {T <: Integer}

A modified version of computing the discrepancies of a two-dimensional cone
which takes an additional parameter `t`. See the `discrepancies` function from
RationalPolygons.jl for the original.

The discrepancy of an element `w` of the hilbert basis of `v1` and `v2` is
defined to be the quotient of the length of `w` by the length of the
intersection point between the ray going through `w` and the line going through
`v1` and `v2`. Here, we add an additional parameter `t`, which modifies this to take the line through `t*v1` and `v2` instead. Hence, for `t = 1`, this behaves
exactly like the unaltered `discrepancies` function.

This extra functionality is needed for computing the discrepancies of elliptic
fixed points of C*-surfaces.

"""
discrepancies(v1 :: LatticePoint{T}, v2 :: LatticePoint{T}, t :: Rational{T}) where {T <: Integer} =
discrepancies(Matrix2{T}(v1[1],v1[2],v2[1],v2[2]), t)
