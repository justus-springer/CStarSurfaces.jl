
@doc raw"""
    MoriDreamSpacePoint

A point on a Mori dream space. 

Subtypes of `MoriDreamSpacePoint` should at least implement the following
methods: [`orbit_cone`](@ref), [`cox_coordinates`](@ref).

"""
abstract type MoriDreamSpacePoint end


@doc raw"""
    orbit_cone(x :: MoriDreamSpacePoint)

Given a point $x \in X$ on a Mori dream space, return the index vector of the
cone $\sigma$ of the canonical toric ambient variety such that $x$ is contained
in the toric orbit associated to $\sigma$.

"""
function orbit_cone end


@doc raw"""
    cox_coordinates(x :: MoriDreamSpacePoint)

Return the Cox coordinates of a point on a Mori dream space.

"""
function cox_coordinates end
