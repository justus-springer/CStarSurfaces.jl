import Base: 
    +, 
    -, 
    *,
    axes,
    size,
    similar,
    getindex,
    setindex!,
    length,
    firstindex,
    lastindex,
    iterate,
    map,
    resize!,
    one

import Combinatorics:
    powerset

using CustomUnitRanges: 
    filename_for_zerorange


# for now, just import Oscar
# at some point, this should be replaced by a more precice
# import statement, that only imports what is actually needed
using Oscar

import Oscar:
    map_from_torusinvariant_cartier_divisor_group_to_picard_group,
    AbstractNormalToricVariety,
    normal_toric_variety,
    class_group,
    pm_object,
    coefficients,
    cox_ring,
    canonical_divisor,
    anticanonical_divisor,
    is_ample,
    is_fano,
    permuted


