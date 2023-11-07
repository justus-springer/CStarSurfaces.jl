
@doc raw"""
    cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:IntegerUnion})

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
cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:IntegerUnion}) = mori_dream_space_divisor(X, coeffs)

cstar_surface_divisor(X :: CStarSurface{EE}, coeffs :: DoubleVector{<:IntegerUnion}) = 
cstar_surface_divisor(X, vcat(coeffs...))

cstar_surface_divisor(X :: CStarSurface{PE}, coeffs :: DoubleVector{T}, coeff_plus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_plus]))

cstar_surface_divisor(X :: CStarSurface{EP}, coeffs :: DoubleVector{T}, coeff_minus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_minus]))

cstar_surface_divisor(X :: CStarSurface{PP}, coeffs :: DoubleVector{T}, coeff_plus :: T, coeff_minus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_plus, coeff_minus]))

cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{Vector{T}}, coeffs_plus_minus...) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, DoubleVector(coeffs), coeffs_plus_minus...)


@doc raw"""
    double_coefficients(d :: CStarSurfaceDivisor)

    @req is_prime(d) "The given divisor is n
Return the coefficients of a divisor on a $\mathbb{C}^*$-surface, in
double index notation.

"""
@attr function double_coefficients(d :: CStarSurfaceDivisor)
    X = d.variety
    cs = coefficients(d)
    r, ns = nblocks(X) - 1, block_sizes(X)
    res = DoubleVector{Int}(undef, ns)
    N = 0
    for i = 0 : r
        res[i] = cs[N + 1 : N + ns[i]]
        N += ns[i]
    end
    return res
end


@attr function is_prime_with_double_indices(d :: CStarSurfaceDivisor)
    !is_prime(d) && return nothing
    X, cs = d.variety, double_coefficients(d)
    r, ns = nblocks(X) - 1, block_sizes(X)
    ij = [(i,j) for i = 0 : r for j = 1 : ns[i] if cs[i][j] == 1]

    if !isempty(ij)
        i, j = first(ij)
        return (:D_ij, (i, j))
    elseif has_D_plus(X) && d == D_plus(X)
        N = length(coefficients(d))
        k = has_D_minus(X) ? N - 1 : N
        return (:D_plus, k)
    elseif has_D_minus(X) && d == D_minus(X)
        k = length(coefficients(d))
        return (:D_minus, k)
    end

end


@doc raw"""
    contract_prime_divisor(d :: CStarSurfaceDivisor{T}) where {T <: CStarSurfaceCase}

Contract the given prime divisor and return the resulting
$\mathbb{C}^*$-surface. This amounts to deleting the associated ray from the
generator matrix.

"""
function contract_prime_divisor(d :: CStarSurfaceDivisor{T}) where {T <: CStarSurfaceCase}
    @req is_prime(d) "The given divisor is not prime."
    X = d.variety
    d_case, inds = is_prime_with_double_indices(d)

    if d_case == :D_ij
        i, j = inds
        return cstar_surface(deleteat!(deepcopy(X.l), i, j), deleteat!(deepcopy(X.d), i, j), X.case)
    elseif d_case == :D_plus
        new_case = has_D_minus(X) ? :ep : :ee
        return cstar_surface(X.l, X.d, new_case)
    elseif d_case == :D_minus
        new_case = has_D_plus(X) ? :pe : :ee
        return cstar_surface(X.l, X.d, new_case)
    end

end

