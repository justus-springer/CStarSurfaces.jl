
@doc raw"""
    SurfaceWithTorusActionPoint 

The `Union` of [`CStarSurfacePoint`](@ref) and [`ToricSurfacePoint`](@ref).

"""
const SurfaceWithTorusActionPoint = Union{CStarSurfacePoint, ToricSurfacePoint}
