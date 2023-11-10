@doc raw"""
    toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: IntegerUnion}

Construct a divisor on a toric surface as a linear combination
of the the torus invariant prime divisors.

# Example

```jldoctest
julia> X = toric_surface([[1,0], [0,1], [-1,-5], [0,-1]])
Normal toric surface

julia> D = toric_surface_divisor(X, [1, 2, -1, 5])
ToricSurfaceDivisor(Normal toric surface, [1, 2, -1, 5], #undef)

julia> coefficients(D)
4-element Vector{Int64}:
  1
  2
 -1
  5
```

"""
toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: IntegerUnion} = mori_dream_space_divisor(X, coeffs)

@doc raw"""
    invariant_divisor(X :: ToricSurface, i :: Int)

Return the $i$-th toric divisor $D^i_X$.

"""
function invariant_divisor(X :: ToricSurface, i :: Int)
    @req 1 ≤ i ≤ nrays(X) "must have 1 ≤ i ≤ nrays(X)"
    coeffs = repeat([0], nrays(X))
    coeffs[i] = 1
    return toric_surface_divisor(X, coeffs)
end

@doc raw"""
    invariant_divisors(X :: ToricSurface)

Return all toric divisors $D^i_X$.

"""
@attr invariant_divisors(X :: ToricSurface) = [invariant_divisor(X, i) for i = 1 : nrays(X)]
