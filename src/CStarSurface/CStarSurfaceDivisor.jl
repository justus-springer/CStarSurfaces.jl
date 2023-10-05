
#################################################
# Julia type for divisors on C-star surfaces
#################################################

const CStarSurfaceDivisor{T} = MoriDreamSpaceDivisor{CStarSurface{T}}

#################################################
# Constructors
#################################################

# Constructor with coefficients given as a single vector
cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{S}) where {S <: Oscar.IntegerUnion} = mori_dream_space_divisor(X, coeffs)

# Constructor with coefficients given in double index notation
cstar_surface_divisor(X :: CStarSurface, coeffs :: DoubleVector{S}) where {S <: Oscar.IntegerUnion} = cstar_surface_divisor(X, vcat(coeffs...))
