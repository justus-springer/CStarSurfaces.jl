
@doc raw"""
    CStarSurface{T<:Integer,C,N,M,R}

The type of a ``\mathbb{C}^*`` surface. It has the following type parameters:

- `T <: Integer`. The integer type to be used, e.g. `Int64` or `BigInt`.
- `C :: CStarSurfaceCase`. One of `EE`, `EP`, `PE` and `PP`. This describes the
  fixed point set of the ``\mathbb{C}^*``-surface, namely the existence of elliptic fixed
  points and parabolic fixed point curves.
- `N :: Int`. The number of rays in the toric ambient variety. This equals ``n_0 + \dots + n_r``
  in the notation of Section ``\ref{sec:background_cstar_surfaces}``.
- `M :: Int`. This always equals `2*N`.
- `R :: Int`. The number of arms of the ``\mathbb{C}^*``-surface. This equals ``r+1`` in the
  notation of Section ``\ref{sec:background_cstar_surfaces}``.

The type itself has two fields:

- `vertex_matrix :: SMatrix{2,N,T,M}`. This contains the main part of the data:
  A 2xN integral matrix encoding the entries ``l_{ij}`` and ``d_{ij}`` of the rays
  of the toric ambient. In the notation of Section ``\ref{sec:background_cstar_surfaces}``,
  the vertex matrix has the form
  ```math
  \begin{bmatrix}
  l_{01} & \dots & l_{0n_0} & \dots & l_{r1} & \dots & l_{rn_r} \\
  d_{01} & \dots & d_{0n_0} & \dots & d_{r1} & \dots & d_{rn_r} \\
  \end{bmatrix}.
  ```
- `block_sizes :: SVector{R, T}`. This determines which rays belong to which
  arm in the toric ambient variety. In the notation of Section ``\ref{sec:background_cstar_surfaces}``,
  it is the tuple ``(n_0, \dots, n_r)``.

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
    block_sizes(X :: CStarSurface, i :: Int)

The size of the `i`-th block of a ``\mathbb{C}^*``-surface.

"""
block_sizes(X :: CStarSurface, i :: Int) = block_sizes(X)[i+1]


@doc raw"""
    case(X :: CStarSurface{T,C}) where {T <: Integer}

The case of the ``\mathbb{C}^*``-surface, as a [`CStarSurfaceCase`](@ref).

"""
case(:: CStarSurface{T,C}) where {C, T <: Integer} = C


function Base.show(io :: IO, X :: CStarSurface{T,C,N,M,R}) where {T<:Integer, C, N, M, R}
    ns = block_sizes(X)
    print(io, "C*-surface of case ", C, " with l = (")
    print(join([ns[i] == 1 ? l(X,i-1,1) : "(" * join([l(X,i-1,j) for j = 1 : ns[i]], ",") * ")"
                for i = 1 : R], ","), ") ")
    print(io, "and d = (")
    print(join([ns[i] == 1 ? d(X,i-1,1) : "(" * join([d(X,i-1,j) for j = 1 : ns[i]], ",") * ")"
                for i = 1 : R], ","), ") ")
end


@doc raw"""
    has_elliptic_fixed_point_plus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has an elliptic fixed point ``x^+`` as the source.
This means the case is either `EE` or `EP`.

"""
has_elliptic_fixed_point_plus(:: CStarSurface{T,C}) where {T<:Integer, C} =
has_elliptic_fixed_point_plus(C)


@doc raw"""
    has_elliptic_fixed_point_minus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has an elliptic fixed point ``x^-`` as the sink.
This means the case is either `EE` or `PE`.

"""
has_elliptic_fixed_point_minus(:: CStarSurface{T,C}) where {T<:Integer, C} =
has_elliptic_fixed_point_minus(C)


@doc raw"""
    has_parabolic_fixed_point_curve_plus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has a curve of parabolic fixed points as the source.
This means the case is either `PE` or `PP`.

"""
has_parabolic_fixed_point_curve_plus(:: CStarSurface{T,C}) where {T<:Integer,C} =
has_parabolic_fixed_point_curve_plus(C)


@doc raw"""
    has_parabolic_fixed_point_curve_minus(X :: CStarSurface)

Whether the ``\mathbb{C}^*``-surface has a curve of parabolic fixed points as the sink.
This means the case is either `EP` or ``PP``.

"""
has_parabolic_fixed_point_curve_minus(:: CStarSurface{T,C}) where {T<:Integer,C} =
has_parabolic_fixed_point_curve_minus(C)


@doc raw"""
    l(X :: CStarSurface, i :: Int, j :: Int)

Return the entry ``l_{ij}`` of the defining triple. By convention, the
indexation of the blocks starts with zero, i.e. ``i`` goes from ``0`` to
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
indexation of the blocks starts with zero, i.e. ``i`` goes from ``0`` to
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
indexation of the blocks starts with zero, i.e. `i` goes from ``0`` to
`number_of_blocks(X)-1`. The indexation of the rays in each individual block
starts with one.

This returns the ray as a vector with two entries. See also [`embedded_ray`](@ref) for
the embedding into `R`-dimensional space, which is the actual ray of the
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
embedded ray into `R`-dimensional space, where `R` is the number of blocks.

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

The generator matrix of the ambient toric variety of ``X``. The columns of this
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

The slope of the `j`-th ray of the `i`-th block, i.e. ``d_{ij} / l_{ij}``.

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

The rational number ``\ell^+ := \frac{1}{l_{01}} + \dots + \frac{1}{l_{r1}} - r + 1``, see
Definition 7.4 of [HaHaSp25](@cite)

"""
function l_plus(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ls = [top_ray(X, i)[1] for i = 0 : r]
    return sum(1 .// ls) - r + 1
end


@doc raw"""
    l_minus(X :: CStarSurface)

The rational number ``\ell^- := \frac{1}{l_{0n_0}} + \dots + \frac{1}{l_{rn_r}} - r + 1``, see
Definition 7.4 of [HaHaSp25](@cite)

"""
function l_minus(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ls = [bottom_ray(X, i)[1] for i = 0 : r]
    return sum(1 .// ls) - r + 1
end

