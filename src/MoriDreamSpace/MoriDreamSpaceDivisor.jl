
mori_dream_space_divisor(X :: T, d :: ToricDivisor) where {T <: MoriDreamSpace} = 
MoriDreamSpaceDivisor(X, d)

function mori_dream_space_divisor(X :: T, coeffs :: Vector{S}) where {T <: MoriDreamSpace, S <: Oscar.IntegerUnion}
    toric_div = toric_divisor(canonical_toric_ambient(X), coeffs)
    MoriDreamSpaceDivisor(X, toric_div)
end

toric_divisor(d :: MoriDreamSpaceDivisor) = d.toric_divisor

coefficients(d :: MoriDreamSpaceDivisor) = coefficients(toric_divisor(d))


#################################################
# Addition, subtraction, scalar multiplication
#################################################

function Base.:+(d1 :: MoriDreamSpaceDivisor{T}, d2 :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    @req d1.variety === d2.variety "The divisors must be defined on the same variety"
    new_coeffiicients = coefficients(d1) + coefficients(d2)
    return mori_dream_space_divisor(d1.variety, new_coeffiicients)
end

function Base.:-(d1 :: MoriDreamSpaceDivisor{T}, d2 :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    @req d1.variety === d2.variety "The divisors must be defined on the same variety"
    new_coeffiicients = coefficients(d1) - coefficients(d2)
    return mori_dream_space_divisor(d1.variety, new_coeffiicients)
end

Base.:-(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} = 
mori_dream_space_divisor(d.variety, -coefficients(d))

Base.:*(c :: S, d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace, S <: Oscar.IntegerUnion} = 
mori_dream_space_divisor(d.variety, [ZZRingElem(c)*x for x in coefficients(d)])



