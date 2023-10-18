@doc raw"""
    toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S})

Construct a divisor on a toric surface as a linear combination
of the the torus invariant prime divisors.

# Example

```jldoctest
julia> X = toric_surface([[1,0], [0,1], [-1,-5], [0,-1]])
Normal toric surface

julia> D = toric_surface_divisor(X, [1, 2, -1, 5])
ToricSurfaceDivisor(Normal toric surface, Torus-invariant, non-prime divisor on a normal toric variety)

julia> coefficients(D)
4-element Vector{ZZRingElem}:
 1
 2
 -1
 5
```

"""
toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: Oscar.IntegerUnion} = mori_dream_space_divisor(X, coeffs)
