
@doc raw"""
    ToricSurfaceFixedPoint <: ToricSurfacePoint

A toric fixed points on a toric surface.

"""
@attributes mutable struct ToricSurfaceFixedPoint <: ToricSurfacePoint
    parent :: ToricSurface
    i :: Int

    function ToricSurfaceFixedPoint(X :: ToricSurface, i :: Int)
        @req 1 ≤ i ≤ nrays(X) "must have 1 ≤ i ≤ nrays(X)"
        return new(X, i)
    end
end

Base.parent(x :: ToricSurfaceFixedPoint) = x.parent

@attr orbit_cone(x :: ToricSurfaceFixedPoint) = [x.i, _next_ray_index(parent(x), x.i)]

Base.show(io :: IO, x :: ToricSurfaceFixedPoint) = print(io, "toric fixed point x($(x.i), $(_next_ray_index(parent(x),x.i)))")

@attr fixed_points(X :: ToricSurface) = [ToricSurfaceFixedPoint(X, i) for i = 1 : nrays(X)]

@doc raw"""
    toric_surface_fixed_point(X :: ToricSurface, i :: Int)

Return the $i$-th toric fixed point of a toric surface.

"""
function toric_surface_fixed_point(X :: ToricSurface, i :: Int)
    @req 1 ≤ i ≤ nrays(X) "must have 1 ≤ i ≤ nrays(X)"
    return fixed_points(X)[i]
end

@attr is_quasismooth(x :: ToricSurfaceFixedPoint) = true

