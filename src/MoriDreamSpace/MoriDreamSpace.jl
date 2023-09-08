#################################################
# Julia type for Mori Dream Spaces
#################################################

abstract type MoriDreamSpace <: AbstractNormalToricVariety end

#################################################
# Subtypes of `MoriDreamSpace` should implement 
# the following methods:
#
# canonical_toric_ambient(X :: MoriDreamSpace)
# cox_ring_relations(X :: MoriDreamSpace)
#################################################

pm_object(X :: MoriDreamSpace) = pm_object(canonical_toric_ambient(X))


#################################################
# Properties
#################################################

@attr xcones_indices(X :: MoriDreamSpace) = cones_indices(canonical_toric_ambient(X))

@attr cox_ring(X :: MoriDreamSpace) = 
quo(cox_ring(canonical_toric_ambient(X)), cox_ring_relations(X))[1]

@attr function canonical_divisor(X :: MoriDreamSpace)
    # disabled for now, since it throws error: "Not a Groebner basis"
    #@req is_complete_intersection(cox_ring(X)) "only avaliable for complete intersection cox rings"
    relation_degrees = map(degree, cox_ring_relations(X)) 
    Z = canonical_toric_ambient(X)
    relation_divisors = [toric_divisor(toric_divisor_class(Z, u)) for u in relation_degrees]
    # we need to treat the treat of no relations in the cox ring seperately, since
    # zero(::Type{ToricDivisor}) can't be defined properly.
    mds_canonical_divisor = isempty(relation_divisors) ? canonical_divisor(Z) : sum(relation_divisors) + canonical_divisor(Z)
    # convert to a divisor living on X instead of its ambient toric variety
    mds_canonical_divisor = toric_divisor(X, coefficients(mds_canonical_divisor))
    return mds_canonical_divisor
end

@attr anticanonical_divisor(X :: MoriDreamSpace) = -canonical_divisor(X)



