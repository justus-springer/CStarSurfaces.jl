@doc raw"""
    free_part(dc :: ToricDivisorClass)

Return the coefficients of the free part of a toric divisor class.

"""
@attr function free_part(dc :: ToricDivisorClass)
    r = rank(class_group(dc.toric_variety))
    n = ncols(dc.class.coeff)
    return [dc.class.coeff[1,i] for i = n-r+1 : n]
end


@doc raw"""
    is_movable(dc :: ToricDivisorClass)

Determine whether the toric divisor class `dc` is movable.

"""
@attr is_movable(dc :: ToricDivisorClass) = free_part(dc) ∈ moving_cone(dc.toric_variety)


@doc raw"""
    is_semiample(dc :: ToricDivisorClass)

Determine whether the toric divisor class `dc` is semi-ample.

"""
@attr is_semiample(dc :: ToricDivisorClass) = free_part(dc) ∈ semiample_cone(dc.toric_variety)


@doc raw"""
    is_ample(dc :: ToricDivisorClass)

Determine whether the toric divisor class `dc` is ample.

"""
@attr is_ample(dc :: ToricDivisorClass) =
Polymake.polytope.contains_in_interior(pm_object(semiample_cone(dc.toric_variety)), free_part(dc))



