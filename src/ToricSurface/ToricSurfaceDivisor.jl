@doc raw"""
    toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S})

Construct a divisor on a toric surface as a linear combination
of the the torus invariant prime divisors.
"""
toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: Oscar.IntegerUnion} = mori_dream_space_divisor(X, coeffs)
