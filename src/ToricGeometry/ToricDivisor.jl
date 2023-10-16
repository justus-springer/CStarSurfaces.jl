#################################################
# Properties
#################################################

@doc raw"""
    is_movable(d :: ToricDivisor)

Determine whether the toric divisor `d` is movable.

"""
@attr is_movable(d :: ToricDivisor) = is_movable(toric_divisor_class(d))


@doc raw"""
    is_semiample(d :: ToricDivisor)

Determine whether the toric divisor `d` is semi-ample.

"""
@attr is_semiample(d :: ToricDivisor) = is_semiample(toric_divisor_class(d))

@doc raw"""
    is_ample(d :: ToricDivisor)

Determine whether the toric divisor `d` is ample.

"""
@attr is_ample(d :: ToricDivisor) = is_ample(toric_divisor_class(d))


