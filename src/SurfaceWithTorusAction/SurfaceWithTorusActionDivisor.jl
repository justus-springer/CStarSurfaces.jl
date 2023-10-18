
@doc raw"""
    SurfaceWithTorusActionDivisor = Union{CStarSurfaceDivisor, ToricSurfaceDivisor}

A Weil divisor on a surface with non-trivial torus action.

"""
const SurfaceWithTorusActionDivisor = Union{CStarSurfaceDivisor, ToricSurfaceDivisor}

#################################################
# Intersection numbers
#################################################

@doc raw"""
    Base.:*(d1 :: SurfaceWithTorusActionDivisor, d2 :: SurfaceWithTorusActionDivisor)

Return the intersection number of two divisors on a surface with torus action.
The divisors have to be defined on the same variety.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> D = cstar_surface_divisor(X, [1, 0, 2, -1])
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), Torus-invariant, non-prime divisor on a normal toric variety)

julia> D * D
4//3
```

"""
function Base.:*(d1 :: SurfaceWithTorusActionDivisor, d2 :: SurfaceWithTorusActionDivisor)
    @req d1.variety === d2.variety "The divisors must be defined on the same variety"
    c1, c2 = coefficients(d1), coefficients(d2)
    M = intersection_matrix(d1.variety)
    n = length(c1)
    return sum([c1[i] * c2[j] * M[i,j] for i = 1 : n, j = 1 : n])
end

