
@doc raw"""
    gen_matrix(X :: MoriDreamSpaceUnion)   

Return the generator matrix of a Mori Dream Space `X`. The columns of this
matrix are the rays of the fan of the canonical toric ambient variety of `X`.

"""
gen_matrix(X :: MoriDreamSpaceUnion) = transpose(matrix(ZZ, rays(X)))


@doc raw"""
    class_group_rank(X :: MoriDreamSpaceUnion)

Return the rank of the class group of a Mori Dream Space.

"""
@attr class_group_rank(X :: MoriDreamSpaceUnion) = rank(class_group(X))


@doc raw"""
    class_group_torsion(X :: MoriDreamSpaceUnion)

Return the list of elementary divisors that make up the torsion part of
the class group of a Mori Dream Space.

"""
@attr class_group_torsion(X :: MoriDreamSpaceUnion) = elementary_divisors(class_group(X))[1 : end - rank(class_group(X))]


@doc raw"""
    class_group_torsion_order(X :: MoriDreamSpaceUnion)

Return the order of the torsion part of the class group of a Mori Dream Space.

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

"""
@attr degree_matrix(X :: MoriDreamSpaceUnion) =
transpose(vcat([w.coeff for w in cox_ring_weights(X)]))


@doc raw"""
    degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)

The torsion part of the degree matrix of a Mori Dream Space. By the
convention in this package, these are the upper rows of `degree_matrix(X)`.

"""
@attr function degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[1 : end - rank(class_group(X)), :]
end


@doc raw"""
    degree_matrix_free_part(X :: MoriDreamSpaceUnion)

The free part of the degree matrix of a Mori Dream Space. By the
convention in this package, these are the lower rows of `degree_matrix(X)`.

"""
@attr function degree_matrix_free_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[end - rank(class_group(X)) + 1 : end, :]
end


@doc raw"""
    local_gorenstein_indices(X :: MoriDreamSpaceUnion)

Return the local gorenstein indices of a Mori Dream Space. The result is a
dictionary indexed by the maximal cones of `X`.

"""
@attr function local_gorenstein_indices(X :: MoriDreamSpaceUnion)
    gorenstein_indices = Dict{Vector{Int}, ZZRingElem}()
    K = divisor_class(canonical_divisor_class(X))
    for c in maximal_cones_indices(X)
        f = map_from_class_group_to_local_class_group(X, c)
        gorenstein_indices[c] = order(f(K))
    end
    return gorenstein_indices
end


@doc raw"""
    local_gorenstein_index(X :: MoriDreamSpaceUnion, c :: Vector{Int64})

Return the local gorenstein index of a Mori Dream Space at a given cone `c`, 
which must be given as a vector of ray indices.

"""
local_gorenstein_index(X :: MoriDreamSpaceUnion, c :: Vector{Int64}) =
local_gorenstein_indices(X)[c]


@doc raw"""
    local_picard_indices(X :: MoriDreamSpaceUnion)

Return the local picard indices of a Mori Dream Space. The result is a
dictionary indexed by the maximal cones of `X`.

"""
@attr local_picard_indices(X :: MoriDreamSpaceUnion) = 
Dict([(c, order(g)) for (c,g) in local_class_groups(X)])


@doc raw"""
    local_picard_index(X :: MoriDreamSpaceUnion, c :: Vector{Int64})

Return the local picard index of a Mori Dream Space at a given cone `c`, 
which must be given as a vector of ray indices.

"""
local_picard_index(X :: MoriDreamSpaceUnion, c :: Vector{Int64}) =
local_picard_indices(X)[c]


