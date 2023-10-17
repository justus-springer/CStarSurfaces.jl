
@doc raw"""
    SurfaceWithTorusActionDivisor = Union{CStarSurfaceDivisor, ToricSurfaceDivisor}

A Weil divisor on a surface with non-trivial torus action.

"""
const SurfaceWithTorusActionDivisor = Union{CStarSurfaceDivisor, ToricSurfaceDivisor}

#################################################
# Intersection numbers
#################################################

function Base.:*(d1 :: SurfaceWithTorusActionDivisor, d2 :: SurfaceWithTorusActionDivisor)
    @req d1.variety === d2.variety "The divisors must be defined on the same variety"
    c1, c2 = coefficients(d1), coefficients(d2)
    M = intersection_matrix(d1.variety)
    n = length(c1)
    return sum([c1[i] * c2[j] * M[i,j] for i = 1 : n, j = 1 : n])
end

