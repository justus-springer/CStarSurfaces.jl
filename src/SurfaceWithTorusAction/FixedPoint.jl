
const SurfaceWithTorusActionFixedPoint = Union{CStarSurfaceFixedPoint, ToricSurfaceFixedPoint}

@attr cox_coordinates(x :: SurfaceWithTorusActionFixedPoint) = 
[i ∈ orbit_cone(x) ? 0 : 1 for i = 1 : nrays(parent(x))]

@attr singularities(X :: SurfaceWithTorusAction) = filter(!is_smooth, fixed_points(X))

@attr number_of_singularities(X :: SurfaceWithTorusAction) = length(singularities(X))

@attr log_canonicity(x :: SurfaceWithTorusActionFixedPoint) = 
minimum([[0 // 1] ; canonical_resolution(x)[3]]) + 1

@attr function resolution_graph(x :: SurfaceWithTorusActionFixedPoint)
    (Y, divs, _) = minimal_resolution(x)
    M = Matrix(intersection_matrix(Y))
    inds = map(is_prime_with_index, divs)
    adj_matrix = M[inds, inds]
    nodelabel = [adj_matrix[k, k] for k = 1 : length(divs)]
    return (Graphs.SimpleGraph(adj_matrix), nodelabel)
end


                                            
