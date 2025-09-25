
@doc raw"""
    CStarSurface{T<:Integer,C,N,M,R}

The type of a ``\mathbb{C}^*`` surface. It has the following type parameters:

- `T <: Integer`. The integer type to be used, e.g. `Int64` or `BigInt`.
- `C :: CStarSurfaceCase`. One of `EE`, `EP`, `PE` and `PP`. This describes the
  fixed point set of the ``\mathbb{C}^*``-surface, namely the existence of elliptic fixed
  points and parabolic fixed point curves.
- `N :: Int`. The number of rays in the toric ambient variety.
- `M :: Int`. This always equals `2*N` is there purely for technical reasons.
- `R :: Int`. The number of arms of the ``\mathbb{C}^*``-surface.

The type itself has two fields:

- `vertex_matrix :: SMatrix{2,N,T,M}`. This contains the main part of the data:
  A 2xN integral matrix encoding the rays of the toric ambient variety.
- `block_sizes :: SVector{R, T}`. This determines which rays belong to which
  arm in the toric ambient variety. It can be thought of as a partition $n = n_1 + ... + n_r$.

"""
struct CStarSurface{T,C,N,M,R}

    vertex_matrix :: SMatrix{2,N,T,M}

    block_sizes :: SVector{R, T}

    CStarSurface{T,C}(vertex_matrix :: SMatrix{2,N,T,M}, block_sizes :: SVector{R,T}) where {C, T <: Integer, N, M, R} =
    new{T,C,N,M,R}(vertex_matrix, block_sizes)

    function CStarSurface{T,C}(leaf_matrices :: AbstractMatrix{T}...) where {C, T <: Integer}
        R = length(leaf_matrices)
        block_sizes = SVector{R}([size(V,2) for V ∈ leaf_matrices])
        leaf_matrices = map(V ->
            SMatrix{2, size(V,2)}(V[:, sort(1 : size(V,2), by = j -> V[2,j] // V[1,j], rev = true)]),
            leaf_matrices)
        V = hcat(leaf_matrices...)
        return CStarSurface{T,C}(V, block_sizes)
    end
end


@doc raw"""
    vertex_matrix(X :: CStarSurface)

The vertex matrix of a ``\mathbb{C}^*``-surface. This encodes the rays in the toric ambient
variety.

"""
vertex_matrix(X :: CStarSurface) = X.vertex_matrix


@doc raw"""
    number_of_blocks(X :: CStarSurface)

The number of blocks (arms) of a ``\mathbb{C}^*``-surface.

"""
number_of_blocks(:: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R} = R


@doc raw"""
    block_sizes(X :: CStarSurface)

The sizes of the individual blocks of a ``\mathbb{C}^*``-surface.

"""
block_sizes(X :: CStarSurface) = X.block_sizes


@doc raw"""
    case(X :: CStarSurface{T,C}) where {T <: Integer}

The case of the ``\mathbb{C}^*``-surface, as a `CStarSurfaceCase`.

"""
case(:: CStarSurface{T,C}) where {C, T <: Integer} = C


@doc raw"""
    block_sizes(X :: CStarSurface, i :: Int)

The size of the `i`-th block of a ``\mathbb{C}^*``-surface.

"""
block_sizes(X :: CStarSurface, i :: Int) = block_sizes(X)[i+1]


@doc raw"""
    has_elliptic_fixed_point_plus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has an elliptic fixed point x^+ as the source.

"""
has_elliptic_fixed_point_plus(:: CStarSurface{T,C}) where {T<:Integer, C} =
has_elliptic_fixed_point_plus(C)


@doc raw"""
    has_elliptic_fixed_point_minus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has an elliptic fixed point x^- as the sink.

"""
has_elliptic_fixed_point_minus(:: CStarSurface{T,C}) where {T<:Integer, C} =
has_elliptic_fixed_point_minus(C)


@doc raw"""
    has_parabolic_fixed_point_curve_plus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has a curve of parabolic fixed points as the source.

"""
has_parabolic_fixed_point_curve_plus(:: CStarSurface{T,C}) where {T<:Integer,C} =
has_parabolic_fixed_point_curve_plus(C)


@doc raw"""
    has_parabolic_fixed_point_curve_minus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has a curve of parabolic fixed points as the sink.

"""
has_parabolic_fixed_point_curve_minus(:: CStarSurface{T,C}) where {T<:Integer,C} =
has_parabolic_fixed_point_curve_minus(C)


@doc raw"""
    l(X :: CStarSurface, i :: Int, j :: Int)

Return the entry ``l_{ij}`` of the defining triple. By convention, the
indexation of the blocks starts with zero, i.e. ``i`` goes from 0 to
`number_of_blocks(X)-1`. The indexation of the rays in each individual block
starts with one.

"""
function l(X :: CStarSurface, i :: Int, j :: Int)
    V, ns = vertex_matrix(X), block_sizes(X)
    return V[1,sum(ns[1:i])+j]
end


@doc raw"""
    d(X :: CStarSurface, i :: Int, j :: Int)

Return the entry ``l_{ij}`` of the defining triple. By convention, the
indexation of the blocks starts with zero, i.e. ``i`` goes from 0 to
`number_of_blocks(X)-1`. The indexation of the rays in each individual block
starts with one.

"""
function d(X :: CStarSurface, i :: Int, j :: Int)
    V, ns = vertex_matrix(X), block_sizes(X)
    return V[2,sum(ns[1:i])+j]
end


@doc raw"""
    ray(X :: CStarSurface, i :: Int, j :: Int)

Return the `j`-th ray of the `i`-th block of the ``\mathbb{C}^*``-surface. By convention, the
indexation of the blocks starts with zero, i.e. `i` goes from 0 to
`number_of_blocks(X)-1`. The indexation of the rays in each individual block
starts with one.

This returns the ray as a vector with two entries. See also `embedded_ray` for
the embedding into $r+1$-dimensional space, which is the actual ray of the
ambient toric variety.

"""
function ray(X :: CStarSurface, i :: Int, j :: Int)
    V, ns = vertex_matrix(X), block_sizes(X)
    return V[:,sum(ns[1:i])+j]
end



@doc raw"""
    top_ray(X :: CStarSurface, i :: Int)

The topmost ray of the `i`-th block.

"""
function top_ray(X :: CStarSurface, i :: Int)
    V, ns = vertex_matrix(X), block_sizes(X)
    return V[:,sum(ns[1:i])+1]
end


@doc raw"""
    bottom_ray(X :: CStarSurface, i :: Int)

The bottommost ray of the `i`-th block.

"""
function bottom_ray(X :: CStarSurface, i :: Int)
    V, ns = vertex_matrix(X), block_sizes(X)
    return V[:,sum(ns[1:i+1])]
end


@doc raw"""
    embedded_ray(X :: CStarSurface, i :: Int, j :: Int)

The `j`-th ray of the `i`-th block of the toric ambient variety. This is the
embedde ray into r-dimensional space, where `r` is the number of blocks.

"""
function embedded_ray(X :: CStarSurface{T,C,N,M,R}, i :: Int, j :: Int) where {C, T <: Integer, N, M, R}
    l, d = ray(X, i, j)
    if i == 0
        return SVector{R, T}([t == R ? d : -l for t = 1 : R])
    else
        return SVector{R, T}([t == i ? l : t == R ? d : 0 for t = 1 : R])
    end
end


@doc raw"""
    top_embedded_ray(X :: CStarSurface, i :: Int)

The topmost ray of the `i`-th block as an embedded ray.

"""
top_embedded_ray(X :: CStarSurface, i :: Int) = embedded_ray(X, i, 1)


@doc raw"""
    bottom_embedded_ray(X :: CStarSurface, i :: Int)

The bottommost ray of the `i`-th block as an embedded ray.

"""
bottom_embedded_ray(X :: CStarSurface, i :: Int) = embedded_ray(X, i, block_sizes(X, i))


generator_matrix_core(X :: CStarSurface) = 
hcat([hcat([embedded_ray(X, i, j) for j = 1 : block_sizes(X)[i+1]]...) for i = 0 : number_of_blocks(X)-1]...)

@doc raw"""
    generator_matrix(X :: CStarSurface)

The generator matrix of the ambient toric variety of `X`. The columns of this
matrix are the primitive ray generators of the fan of the ambient toric
variety.

"""
generator_matrix(X :: CStarSurface{T,EE}) where {T<:Integer} = generator_matrix_core(X)
generator_matrix(X :: CStarSurface{T,PE,N,M,R}) where {T <: Integer, N, M, R} =
hcat(generator_matrix_core(X),
     SVector{R,Int}([i == R ? 1 : 0 for i = 1 : R]))
generator_matrix(X :: CStarSurface{T,EP,N,M,R}) where {T <: Integer, N, M, R} =
hcat(generator_matrix_core(X),
     SVector{R,Int}([i == R ? -1 : 0 for i = 1 : R]))
generator_matrix(X :: CStarSurface{T,PP,N,M,R}) where {T <: Integer, N, M, R} =
hcat(generator_matrix_core(X),
     SVector{R,Int}([i == R ? 1 : 0 for i = 1 : R]),
     SVector{R,Int}([i == R ? -1 : 0 for i = 1 : R]))


@doc raw"""
    slope(X :: CStarSurface, i :: Int, j :: Int)

The slope of the `j`-th ray of the `i`-th block, i.e. $d_{ij} / l_{ij}$.

"""
function slope(X :: CStarSurface, i :: Int, j :: Int)
    v = ray(X, i, j)
    return  v[2] // v[1]
end


@doc raw"""
    sum_of_maximal_slopes(X :: CStarSurface)

The sum of the maximal slopes over all blocks.

"""
function sum_of_maximal_slopes(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    return sum([slope(X, i, 1) for i = 0 : r])
end


@doc raw"""
    sum_of_minimal_slopes(X :: CStarSurface)

The sum of the minimal slopes over all blocks.

"""
function sum_of_minimal_slopes(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ns = block_sizes(X)
    return sum([slope(X, i, ns[i+1]) for i = 0 : r])
end


@doc raw"""
    l_plus(X :: CStarSurface)

The sum over the reciprocals of $l_{i1}$, minus $(r-1)$.

"""
function l_plus(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ls = [top_ray(X, i)[1] for i = 0 : r]
    return sum(1 .// ls) - r + 1
end


@doc raw"""
    l_minus(X :: CStarSurface)

The sum over the reciprocals of $l_{in_i}$, minus $(r-1)$.

"""
function l_minus(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ls = [bottom_ray(X, i)[1] for i = 0 : r]
    return sum(1 .// ls) - r + 1
end

