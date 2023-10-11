

#############################################################
# Julia type for toric surfaces
#############################################################
@attributes mutable struct ToricSurface <: ToricVarietyMDS
    vs :: Vector{Vector{T}} where {T <: Oscar.IntegerUnion}
    function ToricSurface(vs :: Vector{Vector{T}}) where {T <: Oscar.IntegerUnion}
        @req all(v -> length(v) == 2, vs) "rays must all be two-dimensional"
        new(vs)
    end
end

Base.:(==)(X :: ToricSurface, Y :: ToricSurface) = X.vs == Y.vs 

toric_surface(vs :: Vector{Vector{T}}) where {T <: Oscar.IntegerUnion} = ToricSurface(vs)

function toric_surface(P :: ZZMatrix)
    cols = [[P[j,i] for j = 1 : nrows(P)] for i = 1 : ncols(P)]
    return toric_surface(cols)
end

rays(X :: ToricSurface) = X.vs

#############################################################
# Construction of the canonical toric ambient
# This amounts to determining the unique complete fan 
# structure for the given rays.
#############################################################

# sorts two-dimensional vectors clockwise, where [1,0] is considered
# minimal
function _is_less(v :: Vector{T}, w :: Vector{T}) where {T <: Integer}
    v[2] ≥ 0 && w[2] < 0 && return true
    v[2] < 0 && w[2] ≥ 0 && return false
    v[2] ≥ 0 && w[2] ≥ 0 && return v[1] // v[2] > w[1] // w[2]
    v[2] < 0 && w[2] < 0 && return v[1] // v[2] > w[1] // w[2]
end

# need extra method, since 1 // 0 is not defined for ZZRingElem
_is_less(v :: Vector{ZZRingElem}, w :: Vector{ZZRingElem}) = _is_less(convert(Vector{Int}, v), convert(Vector{Int}, w))

@attr _ordered_ray_indices(X :: ToricSurface) = sortperm(rays(X); lt = _is_less)

# given the index of a ray of X, returns the index of the
# counterclockwise adjacent ray.
function _next_ray_index(X :: ToricSurface, i :: Int) 
    is =_ordered_ray_indices(X)
    r = nrays(X)
    j = indexin(i, is)[1]
    return is[mod(j+1, 1:r)]
end

@attr maximal_cones_indices(X :: ToricSurface) = [[i, _next_ray_index(X, i)] for i = 1 : nrays(X)]

@attr canonical_toric_ambient(X :: ToricSurface) = normal_toric_variety(rays(X), maximal_cones_indices(X))

###############################################################
# Intersection numbers
###############################################################

@attr function intersection_matrix(X :: ToricSurface)
    vs, r = rays(X), nrays(X)
    inds = _ordered_ray_indices(X)

    IM = zeros(Rational{Int}, r, r)

    _det(v, w) = v[1] * w[2] - v[2] * w[1]

    # Intersection numbers of adjacent divisors
    for i = 1 : r
        j = mod(i+1, 1:r)
        v, w = vs[inds[i]], vs[inds[j]]
        IM[inds[i], inds[j]] = 1 // _det(v,w)
        IM[inds[j], inds[i]] = IM[inds[i], inds[j]]
    end

    # self intersection numbers
    for j = 1 : r
        i, k = mod(j-1, 1:r), mod(j+1, 1:r)
        vi, vj, vk = vs[inds[i]], vs[inds[j]], vs[inds[k]]
        IM[inds[j], inds[j]] = - _det(vi, vk) // (_det(vi, vj) * _det(vj, vk))
    end
    
    return matrix(QQ, IM)

end

Base.show(io :: IO, X :: ToricSurface) = Base.print(io, "Normal toric surface")

##############################################################
# Resolution of singularities
##############################################################

@attr function canonical_resolution(X :: ToricSurface)
    vs = rays(X)
    r = length(vs)
    inds = _ordered_ray_indices(X)

    new_vs = deepcopy(vs)
    ex_rays = Vector{Vector{Vector{Int}}}(undef, r)
    discrepancies = Vector{Vector{Rational{Int}}}(undef, r)
    for i = 1 : r
        v1, v2 = vs[inds[i]], vs[inds[mod(i+1, 1:r)]]
        ex_rays[i], discrepancies[i] = toric_affine_surface_resolution(v1, v2) 
        append!(new_vs, ex_rays[i])
    end

    return (toric_surface(new_vs), ex_rays, discrepancies)

end

@attr maximal_log_canonicity(X :: ToricSurface) = minimum(vcat([0], discrepancies(X)...)) + 1
