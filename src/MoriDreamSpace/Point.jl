@doc raw"""
    parent(x :: MoriDreamSpacePoint)

Return the Mori dream space where `x` lives in.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> x = x_plus(X)
elliptic fixed point x^+

julia> parent(x) === X
true
```

"""
Base.parent(x :: MoriDreamSpacePoint) = x.parent


@doc raw"""
    orbit_cone(x :: MoriDreamSpacePoint)

Given a point $x \in X$ on a Mori dream space, return the index vector of the
cone $\sigma$ of the canonical toric ambient variety such that $x$ is contained
in the toric orbit associated to $\sigma$.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> orbit_cone(x_plus(X))
3-element Vector{Int64}:
 1
 3
 4
```

"""
function orbit_cone end


@doc raw"""
    cox_coordinates(x :: MoriDreamSpacePoint)

Return the Cox coordinates of a point on a Mori dream space.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> cox_coordinates(x_plus(X))
4-element Vector{Int64}:
 0
 1
 0
 0
```

"""
function cox_coordinates end


@doc raw"""
    is_quasismooth(x :: MoriDreamSpacePoint)

Checks whether a point on a Mori dream space is quasismooth.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_quasismooth(x_plus(X))
false
```

"""
function is_quasismooth end


@doc raw"""
    class_group(x :: MoriDreamSpacePoint)

Return the local class group at a given point on a Mori dream space.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> class_group(x_plus(X))
GrpAb: Z/2 x Z/18
```

"""
@attr function class_group(x :: MoriDreamSpacePoint)
    Z = canonical_toric_ambient(parent(x))
    U = affine_toric_chart(Z, orbit_cone(x))
    return class_group(U)
end


@doc raw"""
    class_group_rank(x :: MoriDreamSpacePoint)

Return the rank of the local class group at a point on a Mori Dream Space.

"""
@attr class_group_rank(x :: MoriDreamSpacePoint) = rank(class_group(x))


@doc raw"""
    class_group_torsion(X :: MoriDreamSpacePoint)

Return the list of elementary divisors that make up the torsion part of
the local class group of a point on a Mori Dream Space.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> class_group_torsion(x_plus(X))
2-element Vector{ZZRingElem}:
 2
 18
```

"""
@attr class_group_torsion(x :: MoriDreamSpacePoint) = elementary_divisors(class_group(x))[1 : end - rank(class_group(x))]


@doc raw"""
    class_group_torsion_order(X :: MoriDreamSpacePoint)

Return the order of the torsion part of the class group of a Mori Dream Space.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> class_group_torsion_order(x_plus(X))
36
```

"""
@attr class_group_torsion_order(x :: MoriDreamSpacePoint) = prod(class_group_torsion(x))


@doc raw"""
    map_from_class_group_to_local_class_group(X :: MoriDreamSpacePoint)

Compute the canonical map from the class group of a Mori dream space to 
the local class group at a given point.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> map_from_class_group_to_local_class_group(x_plus(X))
Map with following data
Domain:
=======
Abelian group with structure: Z/2 x Z/6 x Z
Codomain:
=========
Abelian group with structure: Z/2 x Z/18
```

"""
@attr function map_from_class_group_to_local_class_group(x :: MoriDreamSpacePoint)
    Z = canonical_toric_ambient(parent(x))
    K, f = cokernel(map_from_character_lattice_to_torusinvariant_weil_divisor_group(Z))
    c = orbit_cone(x)
    U = affine_toric_chart(Z, c)
    KU, fU = cokernel(map_from_character_lattice_to_torusinvariant_weil_divisor_group(U))
    m = zero_matrix(ZZ, nrays(Z), nrays(U))
    for i in 1:length(c)
        m[c[i],i] = 1
    end
    f = snf(K)[2] * hom(K, KU, m) * inv(snf(KU)[2])
    return hom(class_group(Z), class_group(x), matrix(f))
end


@doc raw"""
    gorenstein_index(X :: MoriDreamSpacePoint)

Return the local gorenstein index of a point on a Mori Dream Space.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> gorenstein_index(x_plus(X))
9
```

"""
@attr gorenstein_index(x :: MoriDreamSpacePoint) = cartier_index(canonical_divisor(parent(x)), x)


@doc raw"""
    gorenstein_coefficients(x :: MoriDreamSpacePoint)

Return the coefficients in the integral linear combination expressing ι*(-K_X)
as a linear combination over the invariant divisor classes w_i with i ∉
orbit_cone(x).

"""
@attr gorenstein_coefficients(x :: MoriDreamSpacePoint) =
cartier_coefficients(gorenstein_index(parent(x)) * anticanonical_divisor(parent(x)), x)


@doc raw"""
    gorenstein_form(x :: MoriDreamSpacePoint)

Return a rational linear form `u` such that `u` evaluates to `a_i` on `v_i` for
all `i ∈ orbit_cone(x)`, where `a_i` is the `i`-th coefficient of the
anticanonical divisor on `parent(x)`.

Note that the gorenstein index at `x` is the smallest integer `ι` such that
`ι*u` is integral.

"""
@attr gorenstein_form(x :: MoriDreamSpacePoint) = cartier_form(anticanonical_divisor(parent(x)), x)


@doc raw"""
    is_factorial(x :: MoriDreamSpacePoint)

Check whether a point on a Mori dream space is factorial, i.e. its local class
group is trivial.

# Example

```jldoctest
julia> X = cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_factorial(x_plus(X))
false
```

"""
@attr is_factorial(x :: MoriDreamSpacePoint) = is_trivial(class_group(x))


@doc raw"""
    is_smooth(x :: MoriDreamSpacePoint)

Check whether a point on a Mori dream space is smooth, i.e. factorial and
quasismooth.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_smooth(x_plus(X))
false
```

"""
@attr is_smooth(x :: MoriDreamSpacePoint) = is_factorial(x) && is_quasismooth(x)
