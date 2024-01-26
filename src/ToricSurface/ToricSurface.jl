
@doc raw"""
    toric_surface(vs :: Vector{Vector{T}}; check_lineality = false) where {T <: IntegerUnion}

Construct a toric surface from a list of integral vectors in two-dimensional
space.

If the keyword argument `check_lineality = true` is given, an error is thrown
when the given vectors do not generate $\mathbb{Q}^2$ as a convex cone. By
default, this is disabled for performance.

# Example

The 5-th Hirzebruch surface.

```jldoctest
julia> toric_surface([[1,0], [0,1], [-1,-5], [0,-1]])
Normal toric surface
```

"""
function toric_surface(vs :: Vector{Vector{T}}; check_lineality = false) where {T <: IntegerUnion} 
    @req all(v -> length(v) == 2, vs) "rays must all be two-dimensional"
    if check_lineality
        @req lineality_dim(positive_hull(vs)) == 2  "The rays must generate QQ^2 as a convex cone"
    end
    return ToricSurface(vs)
end

@doc raw"""
    toric_surface(P :: ZZMatrix; check_lineality = false)
    
Construct a toric surface from an integral matrix, where the columns of
the matrix are the rays of the describing fan.

If the keyword argument `check_lineality = true` is given, an error is thrown
when the given vectors do not generate $\mathbb{Q}^2$ as a convex cone. By
default, this is disabled for performance.

# Example

The 5-th Hirzebruch surface.

```jldoctest
julia> toric_surface(ZZ[1 0 -1 0 ; 0 1 -17 -1])
Normal toric surface
```

"""
function toric_surface(P :: ZZMatrix; check_lineality = false)
    cols = [[P[j,i] for j = 1 : nrows(P)] for i = 1 : ncols(P)]
    return toric_surface(cols; check_lineality)
end

rays(X :: ToricSurface) = X.vs

@attr fano_polytope(X :: ToricSurface) = convex_hull(rays(X), non_redundant = true)


@attr admits_kaehler_einstein_metric(X :: ToricSurface) =
centroid(polarize(fano_polytope(X))) == [0, 0]


@attr is_quasismooth(X :: ToricSurface) = true 

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

@attr canonical_toric_ambient(X :: ToricSurface) = normal_toric_variety(rays(X), maximal_cones_indices(X); non_redundant = true)


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

@doc raw"""
    canonical_resolution(X :: ToricSurface)

Return the canonical resolution of singularities of a toric surface `X`. The
result is a triple `(Y, ex_rays, discrepancies)` where `Y` is the smooth toric
surface in the resolution of singularities of `X`, `ex_rays` contains the rays
of the exceptional divisors in the resolution and `discrepancies` contains
their discrepancies.

# Example

```jldoctest
julia> X = toric_surface(ZZ[1 1 -3 ; 0 4 -7])
Normal toric surface

julia> (Y, ex_div, discr) = canonical_resolution(X);

julia> gen_matrix(Y)
[1   1   -3   1   1   1   0   -1   -2   -1    0]
[0   4   -7   1   2   3   1   -2   -5   -3   -1]

julia> map(E -> E*E, ex_div)
8-element Vector{QQFieldElem}:
 -2
 -2
 -2
 -2
 -3
 -2
 -2
 -3

julia> discr
8-element Vector{Rational{Int64}}:
  0//1
  0//1
  0//1
 -1//5
 -2//5
 -1//7
 -2//7
 -3//7
```

"""
@attr function canonical_resolution(X :: ToricSurface)

    r = nrays(X)

    resolutions = map(canonical_resolution, fixed_points(X))

    new_vs = deepcopy(rays(X))
    discrepancies = Rational{Int}[]
    for (Y, ex_div, discr) in resolutions
        append!(new_vs, rays(Y)[r + 1 : end])
        append!(discrepancies, discr)
    end

    Y = toric_surface(new_vs)

    exceptional_divisors = [invariant_divisor(Y, i) for i = r + 1 : nrays(Y)]

    return (Y, exceptional_divisors, discrepancies)

end


@attr minimal_resolution(X :: ToricSurface) = canonical_resolution(X)


Base.show(io :: IO, X :: ToricSurface) = Base.print(io, "Normal toric surface")
