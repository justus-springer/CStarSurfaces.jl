@attr function free_part(dc :: ToricDivisorClass)
    r = rank(class_group(dc.toric_variety))
    n = ncols(dc.class.coeff)
    return [dc.class.coeff[1,i] for i = n-r+1 : n]
end

@attr is_movable(dc :: ToricDivisorClass) = free_part(dc) ∈ moving_cone(dc.toric_variety)

@attr is_semiample(dc :: ToricDivisorClass) = free_part(dc) ∈ semiample_cone(dc.toric_variety)

@attr is_ample(dc :: ToricDivisorClass) =
Polymake.polytope.contains_in_interior(pm_object(semiample_cone(dc.toric_variety)), free_part(dc))



