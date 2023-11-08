
@doc raw"""
    MoriDreamSpacePoint

A point on a Mori dream space. 

Subtypes of `MoriDreamSpacePoint` should at least implement the following
methods: `parent`, [`orbit_cone`](@ref), [`cox_coordinates`](@ref),
[`is_quasismooth`](@ref).

"""
abstract type MoriDreamSpacePoint end


@doc raw"""
    orbit_cone(x :: MoriDreamSpacePoint)

Given a point $x \in X$ on a Mori dream space, return the index vector of the
cone $\sigma$ of the canonical toric ambient variety such that $x$ is contained
in the toric orbit associated to $\sigma$.

This function should be implemented by all subtypes of
[`MoriDreamSpacePoint`](@ref).

"""
function orbit_cone end


@doc raw"""
    cox_coordinates(x :: MoriDreamSpacePoint)

Return the Cox coordinates of a point on a Mori dream space.

This function should be implemented by all subtypes of
[`MoriDreamSpacePoint`](@ref).

"""
function cox_coordinates end


@doc raw"""
    is_quasismooth(x :: MoriDreamSpacePoint)

Checks whether a point on a Mori dream space is quasismooth.

This function should be implemented by all subtypes of
[`MoriDreamSpacePoint`](@ref).

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
@attr function gorenstein_index(x :: MoriDreamSpacePoint)
    K = divisor_class(canonical_divisor_class(parent(x)))
    f = map_from_class_group_to_local_class_group(x)
    return order(f(K))
end

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

"""
@attr is_smooth(x :: MoriDreamSpacePoint) = is_factorial(x) && is_quasismooth(x)
