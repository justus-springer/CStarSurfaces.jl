###########################################################################
# This file implements low-level functionality around toric resolutions
# of singularities. Some of this (e.g. Hirzebruch-Jung continued fractions)
# is available in Oscar through polymake. We implement it here in pure
# Julia for increased performance.
###########################################################################

# Performs integer division in the sense of Cox, Little, Schenk (10.1.1)
# Returns a pair of integers s,k such that l = s*d - k
cls_div(l :: T, d :: T) where {T <: Oscar.IntegerUnion} = ((l + mod(-l, d)) ÷ d, mod(-l, d))

# Brings a two-dimensional cone into normal form in the sense of CLS
# Proposition 10.1.1. Returns a triple (d, k, M), where d and k are 
# the parameters of the cone and M is a 2x2 integral matrix such that
# M * A = [0 d ; 1 -k]
function cls_cone_normal_form(A :: ZZMatrix)
    u11, u21, u12, u22 = A[1,1], A[2,1], A[1,2], A[2,2]
    # Get integers x.y such that x * u11 + y * u21 == 1 
    _, x, y = gcdx(u11, u21)
    sg, d, l = sign(det(A)), abs(det(A)), x * u12 + y * u22
    s, k = cls_div(l, d)
    M = ZZ[-sg*u21 sg*u11 ; x+sg*s*u21 y-sg*s*u11]
    return (d, k, M)
end

function hirzebruch_jung(x :: T, y :: T) where {T <: Oscar.IntegerUnion}
    res = T[]
    while y > 0
        push!(res, (x + (mod(-x, y))) / y)
        x, y = y, mod(-x,y)
    end
    return res;
end    

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

hilbert_basis_2D(v1 :: Vector{T}, v2 :: Vector{T}) where {T <: Oscar.IntegerUnion} =
hilbert_basis_2D(matrix(ZZ, [v1[1] v2[1] ; v1[2] v2[2]]))
