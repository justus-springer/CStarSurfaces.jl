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
# Properties
#################################################

@attr xcones_indices(X :: MoriDreamSpace) = cones_indices(canonical_toric_ambient(X))

@attr cox_ring(X :: MoriDreamSpace) = 
quo(cox_ring(canonical_toric_ambient(X)), cox_ring_relations(X))[1]

@attr class_group(X :: MoriDreamSpace) = class_group(canonical_toric_ambient(X))

@attr picard_group(X :: MoriDreamSpace) = picard_group(canonical_toric_ambient(X))

@attr map_from_picard_group_to_class_group(X :: MoriDreamSpace) = map_from_picard_group_to_class_group(canonical_toric_ambient(X))

@attr picard_index(X :: MoriDreamSpace) = picard_index(canonical_toric_ambient(X))

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

@attr function gorenstein_index(X :: MoriDreamSpace)
    c = divisor_class(canonical_divisor_class(X))
    f = cokernel(map_from_picard_group_to_class_group(X))[2]
    order(f(c))
end
