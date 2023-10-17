
@doc raw"""
    SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

Julia type for surfaces with a non-trivial torus action.

"""
const SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

dim(::SurfaceWithTorusAction) = 2


@doc raw"""
    picard_index(X :: SurfaceWithTorusAction)   

Return the index of the Picard group in the class group of a surface with torus
action.

This uses the formula from [Sp23](@cite) and is faster than the more general
implementation.

"""
picard_index(X :: SurfaceWithTorusAction) = prod(values(local_picard_indices(X))) / class_group_torsion_order(X)


@doc raw"""
    anticanonical_self_intersection(X :: SurfaceWithTorusAction)   

Return the self intersection number of the anticanonical divisor on a surface
with torus action.

"""
anticanonical_self_intersection(X :: SurfaceWithTorusAction) = anticanonical_divisor(X) * anticanonical_divisor(X)


@doc raw"""
    exceptional_rays(X :: SurfaceWithTorusAction)   

Return the rays associated to the exceptional divisors in the resolution of
singularities of a surface with torus action.

"""
exceptional_rays(X :: SurfaceWithTorusAction) = canonical_resolution(X)[2]


@doc raw"""
    discrepancies(X :: SurfaceWithTorusAction)

Return the discrepancies associated to the exceptional divisors in the
resolution of singularities of a surface with torus action.

"""
discrepancies(X :: SurfaceWithTorusAction) = canonical_resolution(X)[3]
