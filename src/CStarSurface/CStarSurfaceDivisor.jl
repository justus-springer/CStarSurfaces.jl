
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

#################################################
# Intersection numbers
#################################################

function Base.:*(d1 :: CStarSurfaceDivisor, d2 :: CStarSurfaceDivisor)
    @req d1.variety === d2.variety "The divisors must be defined on the same variety"
    c1, c2 = coefficients(d1), coefficients(d2)
    M = intersection_matrix(d1.variety)
    n = length(c1)
    return sum([c1[i] * c2[j] * M[i,j] for i = 1 : n, j = 1 : n])
end

