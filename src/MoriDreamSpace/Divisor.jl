
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


@doc raw"""
    divisor_class(d :: MoriDreamSpaceDivisor)

Return the divisor class of a divisor on a Mori Dream Space `X`, i.e. the
associated element in `class_group(X)`.

"""
function divisor_class(d :: MoriDreamSpaceDivisor)
    f = map_from_torusinvariant_weil_divisor_group_to_class_group(d.variety)
    weil_divisor_group = domain(f)
    return f(weil_divisor_group(coefficients(d)))
end

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

@doc raw"""
    is_movable(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Check whether a given divisor is movable.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_movable(anticanonical_divisor(X))
true
```

"""
@attr function is_movable(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    X = d.variety
    dc = free_part(divisor_class(d))
    cs = [dc.coeff[1,i] for i = 1 : ncols(dc.coeff)]
    return cs ∈ moving_cone(X)
end


@doc raw"""
    is_semiample(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Check whether a given divisor is semiample.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_semiample(anticanonical_divisor(X))
true
```

"""
@attr function is_semiample(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    X = d.variety
    dc = free_part(divisor_class(d))
    cs = [dc.coeff[1,i] for i = 1 : ncols(dc.coeff)]
    return cs ∈ semiample_cone(X)
end


@doc raw"""
    is_ample(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Check whether a given divisor is ample.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_ample(anticanonical_divisor(X))
true
```

"""
@attr function is_ample(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    X = d.variety
    dc = free_part(divisor_class(d))
    cs = [dc.coeff[1,i] for i = 1 : ncols(dc.coeff)]
    Polymake.polytope.contains_in_interior(pm_object(semiample_cone(X)), cs)
end

@doc raw"""
    is_principal(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Checks whether a given divisor is principal.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> D = 3*invariant_divisor(X,1,1) - 2*invariant_divisor(X,2,1)
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), ZZRingElem[0, 0, 3, -2], #undef)

julia> is_principal(D)
true
```

"""
@attr is_principal(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} =
iszero(divisor_class(d))


@doc raw"""
    cartier_index(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}

Return the smallest integer `ι` such that `ι*d` is Cartier.

"""
@attr function cartier_index(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
    X = d.variety
    c = divisor_class(d)
    f = cokernel(map_from_picard_group_to_class_group(X))[2]
    return order(f(c))
end


@doc raw"""
    cartier_index(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace}

Return the smallest integer `ι` such that `ι*d` is Cartier (i.e. principal)
near `x`.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> d = invariant_divisor(X,0,1)
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), [1, 0, 0, 0], #undef)

julia> cartier_index(d, x_plus(X))
18
```

"""
function cartier_index(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace}
    @req d.variety == parent(x) "The divisor and the points must be defined on the same variety"
    c = divisor_class(d)
    f = map_from_class_group_to_local_class_group(x)
    return order(f(c))
end


@doc raw"""
    is_cartier(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} =

Check whether a given divisor is Cartier.

"""
@attr is_cartier(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} =
cartier_index(d) == 1


@doc raw"""
    is_cartier(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace} =

Check whether a given divisor is cartier near a given point.

"""
is_cartier(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace} =
cartier_index(d,x) == 1


@doc raw"""
    is_principal(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace} =

Check whether a given divisor is principal near a given point.

"""
is_principal(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace} =
is_cartier(d,x)


@doc raw"""
    cartier_coefficients(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace}

Given a divisor `d` that is cartier (hence principal) near a point `x`, the
divisor class of `d` can be expressed as an integer combination of invariant
divisor classes w_i, where i ∉ orbit_cone(x). This function returns some
(non-unique) coefficients of this expression.

"""
function cartier_coefficients(d :: MoriDreamSpaceDivisor{T}, x :: MoriDreamSpacePoint) where {T <: MoriDreamSpace}
    @req d.variety == parent(x) "The divisor and the points must be defined on the same variety"
    @req is_cartier(d,x) "The divisor is not Cartier at the point"
    X = parent(x)
    Q = degree_matrix_free_part(X)
    M = Q[:,[i for i = 1 : nrays(X) if i ∉ orbit_cone(x)]]
    dc = free_part(divisor_class(d))
    cs = solve(M, transpose(dc.coeff))
    return [cs[i,1] for i = 1 : nrows(cs)]
end
