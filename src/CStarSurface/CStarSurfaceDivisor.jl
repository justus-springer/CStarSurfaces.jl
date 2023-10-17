
@doc raw"""
    cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:Oscar.IntegerUnion})

Construct a divisor on a C-star surface as a linear combination
of the the torus invariant prime divisors. 
"""
cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{S}) where {S <: Oscar.IntegerUnion} = mori_dream_space_divisor(X, coeffs)

cstar_surface_divisor(X :: CStarSurface, coeffs :: DoubleVector{S}) where {S <: Oscar.IntegerUnion} = cstar_surface_divisor(X, vcat(coeffs...))
