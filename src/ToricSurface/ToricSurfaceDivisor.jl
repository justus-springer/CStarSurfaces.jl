

#################################################
# Julia type for divisors on C-star surfaces
#################################################

const ToricSurfaceDivisor = MoriDreamSpaceDivisor{ToricSurface}

#################################################
# Constructors
#################################################

# Constructor with coefficients given as a single vector
toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: Oscar.IntegerUnion} = mori_dream_space_divisor(X, coeffs)
