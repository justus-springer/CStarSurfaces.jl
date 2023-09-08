
@attributes mutable struct SubvarietyOfToricVariety{T<:AbstractNormalToricVariety} <: MoriDreamSpace
    ambient :: T
    relations :: Vector{<:RingElem}

    function SubvarietyOfToricVariety(ambient :: AbstractNormalToricVariety, relations :: Vector{<:RingElem})
        R = cox_ring(ambient)
        @req all(r -> parent(r) == R, relations) "relations must be elements of the cox ring of the ambient toric variety"
        @req all(is_homogeneous, relations) "relations must be homogeneous with respect to the grading of the cox ring of the ambient toric variety"
        new{typeof(ambient)}(ambient, relations)
    end
end

#################################################
# Constructors
#################################################

function subvariety_of_toric_variety(ambient :: AbstractNormalToricVariety, relations :: Vector{<:RingElem}; is_canonical_ambient = false)
    X = SubvarietyOfToricVariety(ambient, relations)
    if is_canonical_ambient
        set_attribute!(X, :canonical_toric_ambient => ambient)
    end
    return X
end

subvariety_of_toric_variety(ambient :: AbstractNormalToricVariety, relation :: RingElem; is_canonical_ambient = false) =
subvariety_of_toric_variety(ambient, [relation]; is_canonical_ambient)

#################################################
# Attributes
#################################################

@attr canonical_toric_ambient(X :: SubvarietyOfToricVariety) = error("not yet implemented")

@attr cox_ring_relations(X :: SubvarietyOfToricVariety) = X.relations

@attr cox_ring(X :: SubvarietyOfToricVariety) = quo(cox_ring(X.ambient), X.relations)[1]

