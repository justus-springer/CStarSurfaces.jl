
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

@attr function canonical_resolution(x :: ToricSurfaceFixedPoint)
    X, c = parent(x), orbit_cone(x)
    vs = rays(X)
    r = length(vs)

    new_vs = deepcopy(vs)

    v1, v2 = vs[c[1]], vs[c[2]]
    new_rays :: Vector{Vector{Int}}, discrepancies :: Vector{Rational{Int}} = toric_affine_surface_resolution(v1, v2) 
    new_vs = [vs ; new_rays]
    Y = toric_surface(new_vs)

    exceptional_divisors = [invariant_divisor(Y, i) for i = r + 1 : length(new_vs)]

    return (Y, exceptional_divisors, discrepancies)

end

@attr minimal_resolution(x :: ToricSurfaceFixedPoint) = canonical_resolution(x)

@attr is_log_terminal(x :: ToricSurfaceFixedPoint) = true

@attr singularity_type(x :: ToricSurfaceFixedPoint) = SingularityTypeA(length(minimal_resolution(x)[2]))
