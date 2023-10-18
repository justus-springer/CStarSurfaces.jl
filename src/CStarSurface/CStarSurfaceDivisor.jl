
@doc raw"""
    cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:Oscar.IntegerUnion})

Construct a divisor on a C-star surface as a linear combination
of the the torus invariant prime divisors. 

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[2, 1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> D = cstar_surface_divisor(X, [0, 1, -1, 3])
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), Torus-invariant, non-prime divisor on a normal toric variety)

julia> coefficients(D)
4-element Vector{ZZRingElem}:
 0
 1
 -1
 3
```

"""
cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:Oscar.IntegerUnion}) = mori_dream_space_divisor(X, coeffs)


@doc raw"""
    cstar_surface_divisor(X :: CStarSurface, coeffs :: DoubleVector{<:Oscar.IntegerUnion})

Construct a divisor on a C-star surface as a linear combination
of the the torus invariant prime divisors, where the coefficients
are given following the double index notation.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[2, 1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> D = cstar_surface_divisor(X, DoubleVector([[0, 1], [-1], [3]]))
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), Torus-invariant, non-prime divisor on a normal toric variety)

julia> coefficients(D)
4-element Vector{ZZRingElem}:
 0
 1
 -1
 3
```

"""
cstar_surface_divisor(X :: CStarSurface, coeffs :: DoubleVector{<:Oscar.IntegerUnion}) = cstar_surface_divisor(X, vcat(coeffs...))
