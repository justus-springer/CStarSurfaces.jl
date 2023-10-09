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

using SQLite


# for now, just import Oscar
# at some point, this should be replaced by a more precice
# import statement, that only imports what is actually needed
using Oscar

import Oscar:
    anticanonical_divisor,
    anticanonical_divisor_class,
    canonical_divisor,
    canonical_divisor_class,
    class_group,
    coefficients,
    cox_ring,
    gorenstein_index,
    is_ample,
    is_fano,
    map_from_picard_group_to_class_group,
    normal_toric_variety,
    dim,
    toric_divisor,
    permuted,
    picard_group,
    picard_index,
    pm_object,
    rays,
    nrays,
    maximal_cones


