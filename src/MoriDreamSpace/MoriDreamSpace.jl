
@doc raw"""
    canonical_toric_ambient(X :: MoriDreamSpace)

Construct the canonical toric ambient variety of a Mori Dream Space.

This function should be implemented by all subtypes of `MoriDreamSpace`

"""
function canonical_toric_ambient end


@doc raw"""
    cox_ring_relations(X :: MoriDreamSpace)

Return the list of relations in the Cox Ring of a Mori Dream Space.
Here, a relation is a `RingElem` whose parent is the Cox Ring of
the canonical toric ambient variety.

This function should be implemented by all subtypes of `MoriDreamSpace`.

"""
function cox_ring_relations end


@doc raw"""
    is_toric(X :: MoriDreamSpace)

Check whether a Mori Dream Space is a toric variety.

"""
is_toric(X :: MoriDreamSpace) = isempty(cox_ring_relations(X))


@doc raw"""
    rays(X :: MoriDreamSpace)
    
Return the rays of the canonical toric ambient variety of a Mori Dream Space.

"""
rays(X :: MoriDreamSpace) = rays(canonical_toric_ambient(X))


@doc raw"""
    nrays(X :: MoriDreamSpace)

Return the number of rays of the canonical toric ambient variety of a Mori Dream Space.

"""
nrays(X :: MoriDreamSpace) = length(rays(X))


@doc raw"""
    maximal_cones(X :: MoriDreamSpace)

Return the maximal cones of the fan of the canonical toric ambient variety
of a Mori Dream Space.

"""
maximal_cones(X :: MoriDreamSpace) = maximal_cones(canonical_toric_ambient(X))


@doc raw"""
    class_group(X :: MoriDreamSpace)

Return the class group of a Mori Dream Space.

"""
@attr class_group(X :: MoriDreamSpace) = class_group(canonical_toric_ambient(X))


@doc raw"""
    picard_group(X :: MoriDreamSpace)

Return the Picard group of a Mori Dream Space.

"""
@attr picard_group(X :: MoriDreamSpace) = picard_group(canonical_toric_ambient(X))


@doc raw"""
    map_from_picard_group_to_class_group(X :: MoriDreamSpace)

Return the embedding of the Picard group into the class group of a
Mori Dream Space.

"""
@attr map_from_picard_group_to_class_group(X :: MoriDreamSpace) = map_from_picard_group_to_class_group(canonical_toric_ambient(X))


@doc raw"""
    picard_index(X :: MoriDreamSpace)

Return the index of the Picard group in the class group of a
Mori Dream Space.

"""
@attr picard_index(X :: MoriDreamSpace) = picard_index(canonical_toric_ambient(X))

@doc raw"""
    is_factorial(X :: MoriDreamSpace)

Determine if a Mori Dream Space has at most factorial singularities, i.e.
its canonical toric ambient variety is smooth.

"""
@attr is_factorial(X :: MoriDreamSpace) = is_smooth(canonical_toric_ambient(X))


@doc raw"""
    is_q_factorial(X :: MoriDreamSpace)

Determine if a Mori Dream Space has at most $\mathbb{Q}$-factorial
singularities, i.e. its canonical toric ambient variety is simplicial.

"""
@attr is_q_factorial(X :: MoriDreamSpace) = is_simplicial(canonical_toric_ambient(X))


@doc raw"""
    cox_ring(X :: MoriDreamSpace)

Compute the Cox Ring of a Mori Dream Space.

"""
@attr cox_ring(X :: MoriDreamSpace) = 
isempty(cox_ring_relations(X)) ? 
    cox_ring(canonical_toric_ambient(X)) :
    quo(cox_ring(canonical_toric_ambient(X)), cox_ring_relations(X))[1]


@doc raw"""
    cox_ring_relation_degrees(X :: MoriDreamSpace)

Compute the degrees of the relations in the Cox Ring of a Mori Dream Space.

"""
@attr cox_ring_relation_degrees(X :: MoriDreamSpace) = map(degree, cox_ring_relations(X))


@doc raw"""
    canonical_divisor(X :: MoriDreamSpace)

Return the canonical divisor of a Mori Dream Space.

"""
@attr function canonical_divisor(X :: MoriDreamSpace)
    # disabled for now, since it throws error: "Not a Groebner basis"
    #@req is_complete_intersection(cox_ring(X)) "only avaliable for complete intersection cox rings"
    Z = canonical_toric_ambient(X)
    relation_divisors = [toric_divisor(toric_divisor_class(Z, u)) for u in cox_ring_relation_degrees(X)]
    # we need to treat the treat of no relations in the cox ring seperately, since
    # zero(::Type{ToricDivisor}) can't be defined properly.
    can_divisor = isempty(relation_divisors) ? canonical_divisor(Z) : sum(relation_divisors) + canonical_divisor(Z)
    return mori_dream_space_divisor(X, can_divisor)
end

@doc raw"""
    canonical_divisor(X :: MoriDreamSpace)

Return the anticanonical divisor of a Mori Dream Space.

"""
@attr anticanonical_divisor(X :: MoriDreamSpace) = -canonical_divisor(X)


@doc raw"""
    canonical_divisor_class(X :: MoriDreamSpace)

Return the canonical divisor class of a Mori Dream Space.

"""
@attr canonical_divisor_class(X :: MoriDreamSpace) = toric_divisor_class(toric_divisor(canonical_divisor(X)))


@doc raw"""
    anticanonical_divisor_class(X :: MoriDreamSpace)

Return the anticanonical divisor class of a Mori Dream Space.

"""
@attr anticanonical_divisor_class(X :: MoriDreamSpace) = -canonical_divisor_class(X)


@doc raw"""
    is_fano(X :: MoriDreamSpace)

Check if a Mori Dream Space if fano.

"""
@attr is_fano(X :: MoriDreamSpace) = is_ample(anticanonical_divisor_class(X))

@doc raw"""
    gorenstein_index(X :: MoriDreamSpace)

Return the Gorenstein index of a $\mathbb{Q}$-Gorenstein Mori Dream Space.

"""
@attr function gorenstein_index(X :: MoriDreamSpace)
    c = divisor_class(canonical_divisor_class(X))
    f = cokernel(map_from_picard_group_to_class_group(X))[2]
    order(f(c))
end

@doc raw"""
    local_class_groups(X :: MoriDreamSpace)

Compute the local class groups of a Mori Dream Space `X`. The result is 
a dictionary indexed by the maximal cones of `X`.

"""
@attr local_class_groups(X :: MoriDreamSpace) = local_class_groups(canonical_toric_ambient(X))


@doc raw"""
    local_class_group(X :: MoriDreamSpace, c :: Vector{Int64})

Compute the local class group of a Mori Dream Space `X` at a maximal
cone `c` given by a vector of indices.

"""
local_class_group(X :: MoriDreamSpace, c :: Vector{Int64}) = local_class_groups(X)[c]


@doc raw"""
    maps_from_class_group_to_local_class_groups(X :: MoriDreamSpace)

Compute the canonical maps from the class group of a Mori Dream Space to 
its local class groups. The result is a dictionary indexed by the maximal 
cones of `X`.

"""
@attr maps_from_class_group_to_local_class_groups(X :: MoriDreamSpace) =
maps_from_class_group_to_local_class_groups(canonical_toric_ambient(X))


@doc raw"""
    map_from_class_group_to_local_class_group(X :: MoriDreamSpace, c :: Vector{Int64})

Compute the canonical map from the class group of a Mori Dream Space `X` to
its local class group at a cone `c`, given by a vector of indices.

"""
map_from_class_group_to_local_class_group(X :: MoriDreamSpace, c :: Vector{Int64}) =
maps_from_class_group_to_local_class_groups(X)[c]
