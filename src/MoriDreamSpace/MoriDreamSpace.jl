#################################################
# Julia type for Mori Dream Spaces
#################################################

abstract type MoriDreamSpace end

#################################################
# Subtypes of `MoriDreamSpace` should implement 
# the following methods:
#
# canonical_toric_ambient(X :: MoriDreamSpace)
# cox_ring_relations(X :: MoriDreamSpace)
#################################################


#################################################
# Properties coming from the toric ambient.
#
# The computation of these fields are simply
# passed to the canonical toric ambient. Subtypes
# of `MoriDreamSpace` can and should overwrite
# there, should a more direct approach exist.
#################################################

rays(X :: MoriDreamSpace) = rays(canonical_toric_ambient(X))

nrays(X :: MoriDreamSpace) = length(rays(X))

maximal_cones(X :: MoriDreamSpace) = maximal_cones(canonical_toric_ambient(X))

@attr class_group(X :: MoriDreamSpace) = class_group(canonical_toric_ambient(X))

@attr picard_group(X :: MoriDreamSpace) = picard_group(canonical_toric_ambient(X))

@attr map_from_picard_group_to_class_group(X :: MoriDreamSpace) = map_from_picard_group_to_class_group(canonical_toric_ambient(X))

@attr picard_index(X :: MoriDreamSpace) = picard_index(canonical_toric_ambient(X))

@attr is_factorial(X :: MoriDreamSpace) = is_smooth(canonical_toric_ambient(X))

@attr is_q_factorial(X :: MoriDreamSpace) = is_simplicial(canonical_toric_ambient(X))

@attr cox_ring(X :: MoriDreamSpace) = 
isempty(cox_ring_relations(X)) ? 
    cox_ring(canonical_toric_ambient(X)) :
    quo(cox_ring(canonical_toric_ambient(X)), cox_ring_relations(X))[1]

@attr cox_ring_relation_degrees(X :: MoriDreamSpace) = map(degree, cox_ring_relations(X))

@attr function canonical_divisor(X :: MoriDreamSpace)
    # disabled for now, since it throws error: "Not a Groebner basis"
    #@req is_complete_intersection(cox_ring(X)) "only avaliable for complete intersection cox rings"
    Z = canonical_toric_ambient(X)
    relation_divisors = [toric_divisor(toric_divisor_class(Z, u)) for u in cox_ring_relation_degrees(X)]
    # we need to treat the treat of no relations in the cox ring seperately, since
    # zero(::Type{ToricDivisor}) can't be defined properly.
    can_divisor = isempty(relation_divisors) ? canonical_divisor(Z) : sum(relation_divisors) + canonical_divisor(Z)
    return can_divisor
end

@attr anticanonical_divisor(X :: MoriDreamSpace) = -canonical_divisor(X)

@attr canonical_divisor_class(X :: MoriDreamSpace) = toric_divisor_class(canonical_divisor(X))

@attr anticanonical_divisor_class(X :: MoriDreamSpace) = -canonical_divisor_class(X)

@attr is_fano(X :: MoriDreamSpace) = is_ample(anticanonical_divisor_class(X))

@attr function gorenstein_index(X :: MoriDreamSpace)
    c = divisor_class(canonical_divisor_class(X))
    f = cokernel(map_from_picard_group_to_class_group(X))[2]
    order(f(c))
end


#################################################
# More attributes for Mori dream spaces and toric varieties
#################################################

# Union type of MoriDreamSpace's and Oscar's toric varieties
const MoriDreamSpaceUnion = Union{MoriDreamSpace, Oscar.NormalToricVarietyType}

gen_matrix(X :: MoriDreamSpaceUnion) = transpose(matrix(ZZ, rays(X)))

@attr function maximal_cones_indices(X :: MoriDreamSpaceUnion) 
    IM = ray_indices(maximal_cones(X))
    [sort(collect(row(IM,i))) for i = 1 : nrows(IM)]
end

@attr covering_collection(X :: MoriDreamSpaceUnion) =
map(c -> [i for i in 1 : nrays(X) if i ∉ c], maximal_cones_indices(X))

@attr cox_ring_weights(X :: MoriDreamSpaceUnion) = map(degree, gens(cox_ring(X)))

@attr degree_matrix(X :: MoriDreamSpaceUnion) =
transpose(vcat([w.coeff for w in cox_ring_weights(X)]))

@attr function degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[1 : end - rank(class_group(X)), :]
end

@attr function degree_matrix_free_part(X :: MoriDreamSpaceUnion)
    Q = degree_matrix(X)
    Q[end - rank(class_group(X)) + 1 : end, :]
end

@attr class_group_torsion_order(X :: MoriDreamSpaceUnion) = 
prod(elementary_divisors(class_group(X))[1 : end - rank(class_group(X))])




