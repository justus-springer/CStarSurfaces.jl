#################################################
# Unary minus
#################################################

Base.:-(d :: ToricDivisor) = toric_divisor(d.toric_variety, -coefficients(d))

#################################################
# Properties
#################################################

@attr is_movable(d :: ToricDivisor) = is_movable(toric_divisor_class(d))

@attr is_semiample(d :: ToricDivisor) = is_semiample(toric_divisor_class(d))

# overwritten from Oscar
@attr is_ample(d :: ToricDivisor) = is_ample(toric_divisor_class(d))


