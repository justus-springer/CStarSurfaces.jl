# Basic types

We briefly go over two basic types used in this package:
[`CStarSurfaceCase`](@ref) and [`AbelianGroup`](@ref).

## Cases

```@docs
CStarSurfaceCase
invert_case
has_elliptic_fixed_point_plus(C :: CStarSurfaceCase)
has_elliptic_fixed_point_minus(C :: CStarSurfaceCase)
has_parabolic_fixed_point_curve_plus(C :: CStarSurfaceCase)
has_parabolic_fixed_point_curve_minus(C :: CStarSurfaceCase)
```

## Abelian groups

```@docs
AbelianGroup
rank
elementary_divisors
torsion_order
torsion_part
cokernel
```
