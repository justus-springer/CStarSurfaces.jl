
cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{<:IntegerUnion}) = mori_dream_space_divisor(X, coeffs)

@doc raw"""
    cstar_surface_divisor(X :: CStarSurface{EE}, coeffs :: DoubleVector{<:IntegerUnion})

Construct a divisor on a $\mathbb{C}^*$-surface of type (e-e) as a linear
combination of the invariant prime divisors $D^{ij}_X$. The coefficients are
given in double index notation.

"""
cstar_surface_divisor(X :: CStarSurface{EE}, coeffs :: DoubleVector{<:IntegerUnion}) = 
cstar_surface_divisor(X, vcat(coeffs...))


@doc raw"""
    cstar_surface_divisor(X :: CStarSurface{PE}, coeffs :: DoubleVector{T}, coeff_plus :: T) where {T <: IntegerUnion}

Construct a divisor on a $\mathbb{C}^*$-surface of type (p-e) as a linear
combination of the invariant prime divisors $D^{ij}_X$ and $D^+_X$. The
coefficients are given in double index notation.

"""
cstar_surface_divisor(X :: CStarSurface{PE}, coeffs :: DoubleVector{T}, coeff_plus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_plus]))


@doc raw"""
    cstar_surface_divisor(X :: CStarSurface{EP}, coeffs :: DoubleVector{T}, coeff_minus :: T) where {T <: IntegerUnion}

Construct a divisor on a $\mathbb{C}^*$-surface of type (e-p) as a linear
combination of the invariant prime divisors $D^{ij}_X$ and $D^-_X$. The
coefficients are given in double index notation.

"""
cstar_surface_divisor(X :: CStarSurface{EP}, coeffs :: DoubleVector{T}, coeff_minus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_minus]))


@doc raw"""
    cstar_surface_divisor(X :: CStarSurface{PP}, coeffs :: DoubleVector{T}, coeff_plus :: T, coeff_minus :: T) where {T <: IntegerUnion}

Construct a divisor on a $\mathbb{C}^*$-surface of type (p-p) as a linear
combination of the invariant prime divisors $D^{ij}_X$ and $D^+_X, D^-_X$. The
coefficients are given in double index notation.

"""
cstar_surface_divisor(X :: CStarSurface{PP}, coeffs :: DoubleVector{T}, coeff_plus :: T, coeff_minus :: T) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, vcat(coeffs..., [coeff_plus, coeff_minus]))


cstar_surface_divisor(X :: CStarSurface, coeffs :: Vector{Vector{T}}, coeffs_plus_minus...) where {T <: IntegerUnion} = 
cstar_surface_divisor(X, DoubleVector(coeffs), coeffs_plus_minus...)


@doc raw"""
    double_coefficients(d :: CStarSurfaceDivisor)

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
    D_plus(X :: CStarSurface{<:Union{PE,PP}})    

Return the parabolic fixed point curve $D^+$ of a $\mathbb{C}^*$-surface of
type (p-e) or (p-p).

"""
@attr D_plus(X :: CStarSurface{PE}) = 
cstar_surface_divisor(X, [repeat([0], n) for n in _ns(X)], 1)
@attr D_plus(X :: CStarSurface{PP}) = 
cstar_surface_divisor(X, [repeat([0], n) for n in _ns(X)], 1, 0)


@doc raw"""
    D_minus(X :: CStarSurface{<:Union{EP,PP}})    

Return the parabolic fixed point curve $D^-$ of a $\mathbb{C}^*$-surface of
type (e-p) or (p-p).

"""
@attr D_minus(X :: CStarSurface{EP}) = 
cstar_surface_divisor(X, [repeat([0], n) for n in _ns(X)], 1)
@attr D_minus(X :: CStarSurface{PP}) = 
cstar_surface_divisor(X, [repeat([0], n) for n in _ns(X)], 0, 1)


@doc raw"""
    parabolic_fixed_point_curves(X :: CStarSurface)   

Return the parabolic fixed point curves of a $\mathbb{C}^*$-surface.

"""
@attr parabolic_fixed_point_curves(X :: CStarSurface{EE}) = CStarSurfaceDivisor{EE}[]
@attr parabolic_fixed_point_curves(X :: CStarSurface{PE}) = [D_plus(X)]
@attr parabolic_fixed_point_curves(X :: CStarSurface{EP}) = [D_minus(X)]
@attr parabolic_fixed_point_curves(X :: CStarSurface{PP}) = [D_plus(X), D_minus(X)]


@doc raw"""
    invariant_divisor(X :: CStarSurface, i :: Int, j :: Int)

Return the $(i,j)$-th invariant divisor $D^{ij}_X$.

"""
function invariant_divisor(X :: CStarSurface, i :: Int, j :: Int)
    r, ns, m = _r(X), _ns(X), _m(X)
    @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
    @req 1 ≤ j ≤ ns[i] "must have 1 ≤ j ≤ ns[i]"
    coeffs = [repeat([0], n) for n in ns]
    coeffs[i][j] = 1
    m == 0 && return cstar_surface_divisor(X, coeffs)
    m == 1 && return cstar_surface_divisor(X, coeffs, 0)
    return cstar_surface_divisor(X, coeffs, 0, 0)
end

_invariant_divisors_core(X :: CStarSurface) = 
DoubleVector([[invariant_divisor(X, i, j) for j = 1 : _ns(X)[i]] for i = 0 : _r(X)])


@doc raw"""
    invariant_divisors(X :: CStarSurface)

Return all invariant divisors $D^{ij}_X, D^{\pm}_X$. The result is given as a
pair with first entry the divisors of the form $D^{ij}$ and second entry the 
divisors of the form $D^{\pm}$.

"""
@attr invariant_divisors(X :: CStarSurface{EE}) = (_invariant_divisors_core(X), CStarSurfaceDivisor{EE}[])
@attr invariant_divisors(X :: CStarSurface{PE}) = (_invariant_divisors_core(X), [D_plus(X)])
@attr invariant_divisors(X :: CStarSurface{EP}) = (_invariant_divisors_core(X), [D_minus(X)])
@attr invariant_divisors(X :: CStarSurface{PP}) = (_invariant_divisors_core(X), [D_plus(X), D_minus(X)])

@doc raw"""
    contract_prime_divisor(d :: CStarSurfaceDivisor{T}) where {T <: CStarSurfaceCase}

Contract the given prime divisor and return the resulting
$\mathbb{C}^*$-surface. This amounts to deleting the associated ray from the
generator matrix.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ep)
C-star surface of type (e-p)

julia> contract_prime_divisor(D_minus(X))
C-star surface of type (e-e)
```

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

