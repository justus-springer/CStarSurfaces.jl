
const SurfaceWithTorusActionFixedPoint = Union{CStarSurfaceFixedPoint, ToricSurfaceFixedPoint}

cox_coordinates(x :: SurfaceWithTorusActionFixedPoint) = 
[i ∈ orbit_cone(x) ? 0 : 1 for i = 1 : nrays(parent(x))]
                                            
