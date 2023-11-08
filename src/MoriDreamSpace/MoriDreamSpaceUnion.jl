
@doc raw"""
    gen_matrix(X :: MoriDreamSpaceUnion)   

Return the generator matrix of a Mori Dream Space `X`. The columns of this
matrix are the rays of the fan of the canonical toric ambient variety of `X`.

# Example

```jldoctest
julia> gen_matrix(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
[-3   -1   3   0]
[-3   -1   0   2]
[-2   -1   1   1]
```

"""
gen_matrix(X :: MoriDreamSpaceUnion) = transpose(matrix(ZZ, rays(X)))


@doc raw"""
    class_group_rank(X :: MoriDreamSpaceUnion)

Return the rank of the class group of a Mori Dream Space.

# Example

```jldoctest
julia> class_group_rank(cstar_surface([[1, 1], [2], [2]], [[-3, -4], [1], [1]], :pe))
2
```

"""
@attr class_group_rank(X :: MoriDreamSpaceUnion) = rank(class_group(X))


@doc raw"""
    class_group_torsion(X :: MoriDreamSpaceUnion)

Return the list of elementary divisors that make up the torsion part of
the class group of a Mori Dream Space.

# Example

```jldoctest
julia> class_group_torsion(cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee))
2-element Vector{ZZRingElem}:
 2
 6
```

"""
@attr class_group_torsion(X :: MoriDreamSpaceUnion) = elementary_divisors(class_group(X))[1 : end - rank(class_group(X))]


@doc raw"""
    class_group_torsion_order(X :: MoriDreamSpaceUnion)

Return the order of the torsion part of the class group of a Mori Dream Space.

# Example

```jldoctest
julia> class_group_torsion_order(cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee))
12
```

"""
@attr class_group_torsion_order(X :: MoriDreamSpaceUnion) = prod(class_group_torsion(X))


@doc raw"""
    maximal_cones_indices(X :: MoriDreamSpaceUnion)

Return the list of maximal cones of a Mori Dream Space, where each cone
is represented as a list of integers that are the indices of its generating
rays in the fan.

"""
@attr function maximal_cones_indices(X :: MoriDreamSpaceUnion) 
    IM = ray_indices(maximal_cones(X))
    [sort(collect(row(IM,i))) for i = 1 : nrows(IM)]
end


@doc raw"""
    covering_collection(X :: MoriDreamSpaceUnion)

Return the covering collection of a Mori Dream Space. The elements are the
complements of the cones in `maximal_cones_indices(X)`.

"""
@attr covering_collection(X :: MoriDreamSpaceUnion) =
map(c -> [i for i in 1 : nrays(X) if i ∉ c], maximal_cones_indices(X))


@doc raw"""
    cox_ring_weights(X :: MoriDreamSpaceUnion)

Return the weights in the Cox Ring, i.e. the degrees of its generators.

"""
@attr cox_ring_weights(X :: MoriDreamSpaceUnion) = map(degree, gens(cox_ring(X)))


@doc raw"""
    degree_matrix(X :: MoriDreamSpaceUnion)

Return the degree matrix of a Mori Dream Space (often denoted as $Q$).
The columns of this matrix are the degrees of the generator of the Cox
Ring, which are elements of the divisor class group. Note that we write
the torsion parts of these elements *first* (in the upper rows of the 
degree matrix), and the free part after that (in the lower rows of the
degree matrix).

# Example

```jldoctest
julia> degree_matrix(cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee))
[0   1   1   0]
[0   4   1   5]
[1   3   4   2]
```

"""
@attr degree_matrix(X :: MoriDreamSpaceUnion) =
transpose(vcat([w.coeff for w in cox_ring_weights(X)]))


@doc raw"""
    degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)

The torsion part of the degree matrix of a Mori Dream Space. By the
convention in this package, these are the upper rows of `degree_matrix(X)`.

# Example

```jldoctest
julia> degree_matrix_torsion_part(cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee))
[0   1   1   0]
[0   4   1   5]
```

"""
@attr function degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[1 : end - rank(class_group(X)), :]
end


@doc raw"""
    degree_matrix_free_part(X :: MoriDreamSpaceUnion)

The free part of the degree matrix of a Mori Dream Space. By the
convention in this package, these are the lower rows of `degree_matrix(X)`.

# Example

```jldoctest
julia> degree_matrix_free_part(cstar_surface([[2, 2], [2], [4]], [[3, -3], [1], [1]], :ee))
[1   3   4   2]
```

"""
@attr function degree_matrix_free_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[end - rank(class_group(X)) + 1 : end, :]
end

