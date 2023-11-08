
const SurfaceWithTorusActionFixedPoint = Union{CStarSurfaceFixedPoint, ToricSurfaceFixedPoint}

@attr cox_coordinates(x :: SurfaceWithTorusActionFixedPoint) = 
[i ∈ orbit_cone(x) ? 0 : 1 for i = 1 : nrays(parent(x))]

@attr singularities(X :: SurfaceWithTorusAction) = filter(!is_smooth, fixed_points(X))

@attr number_of_singularities(X :: SurfaceWithTorusAction) = length(singularities(X))

                                            
