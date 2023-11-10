
@doc raw"""
    mori_dream_space_divisor(X :: T, coeffs :: Vector{S})

Construct a divisor on a Mori Dream Space as a linear combination of the
(restrictions of) the torus invariant prime divisors of the canonical 
toric ambient variety.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[2, 1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> D = cstar_surface_divisor(X, [0, 1, -1, 3])
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), [0, 1, -1, 3], #undef)

julia> coefficients(D)
4-element Vector{Int64}:
  0
  1
 -1
  3
```

"""
function mori_dream_space_divisor(X :: T, coeffs :: Vector{S}) where {T <: MoriDreamSpace, S <: IntegerUnion}
    MoriDreamSpaceDivisor(X, coeffs)
end


@doc raw"""
    mori_dream_space_divisor(X :: T, d :: ToricDivisor)

Return the Divisor on a Mori Dream Space associated to a toric divisor
on its canonical toric ambient space.

"""
function mori_dream_space_divisor(X :: T, td :: ToricDivisor) where {T <: MoriDreamSpace}
    @req td.toric_variety === canonical_toric_ambient(X) "the given toric divisor must be defined on the canonical toric ambient space of X"
    d = MoriDreamSpaceDivisor(X, coefficients(td))
    set_attribute!(d, :toric_divisor, td)
    return d
end


@doc raw"""
    toric_divisor(d :: MoriDreamSpaceDivisor)

Return the toric divisor on the canonical toric ambient variety associated 
to a divisor on a Mori Dream Space.

"""
@attr toric_divisor(d :: MoriDreamSpaceDivisor) = toric_divisor(canonical_toric_ambient(d.variety), d.coeffs)


@doc raw"""
    coefficients(d :: MoriDreamSpaceDivisor)

Return the coefficients of a divisor on a Mori Dream Space.

"""
coefficients(d :: MoriDreamSpaceDivisor) = d.coeffs


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

Base.:*(c :: S, d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace, S <: IntegerUnion} = 
mori_dream_space_divisor(d.variety, [ZZRingElem(c)*x for x in coefficients(d)])

#################################################
# Basic properties
#################################################

@doc raw"""
    is_prime(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Check whether a given divisor is a prime divisor.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_prime(invariant_divisor(X, 0, 1))
true
```

"""
@attr is_prime(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} = 
sum(coefficients(d)) == 1 && all(c -> c ∈ [0,1], coefficients(d))

function is_prime_with_index(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    !is_prime(d) && return nothing
    cs = coefficients(d)
    for i = 1 : length(cs)
        cs[i] == 1 && return i
    end
    return nothing
end






