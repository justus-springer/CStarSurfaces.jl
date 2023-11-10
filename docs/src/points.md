# Points

We provide an abstract interface for working with points on Mori dream spaces
in terms of their Cox coordinates and orbit cone. An implementation for fixed
points on $\mathbb{C}^*$-surfaces (elliptic, hyperbolic and parabolic) as well
as toric fixed points on surfaces is provided.

## Types

```@docs
MoriDreamSpacePoint
CStarSurfacePoint
ToricSurfacePoint
SurfaceWithTorusActionPoint
CStarSurfaceFixedPoint
ToricSurfaceFixedPoint
SurfaceWithTorusActionFixedPoint
EllipticFixedPoint
EllipticFixedPointPlus
EllipticFixedPointMinus
HyperbolicFixedPoint
ParabolicFixedPoint
ParabolicFixedPointPlus
ParabolicFixedPointMinus
```

## Constructors

```@docs
x_plus
x_minus
hyperbolic_fixed_point
parabolic_fixed_point_plus
parabolic_fixed_point_minus
```

## Sets of fixed points

```@docs
elliptic_fixed_points
hyperbolic_fixed_points
parabolic_fixed_points_plus
parabolic_fixed_points_minus
parabolic_fixed_points
fixed_points
```

## Attributes

```@docs
parent(x :: MoriDreamSpacePoint)
orbit_cone
cox_coordinates
class_group(x :: MoriDreamSpacePoint)
class_group_rank(x :: MoriDreamSpacePoint)
class_group_torsion(x :: MoriDreamSpacePoint)
class_group_torsion_order(x :: MoriDreamSpacePoint)
map_from_class_group_to_local_class_group(x :: MoriDreamSpacePoint)
gorenstein_index(x :: MoriDreamSpacePoint)
is_quasismooth
is_factorial(x :: MoriDreamSpacePoint)
is_smooth(x :: MoriDreamSpacePoint)
```

## Resolution of singularities

```@docs
canonical_resolution
minimal_resolution
log_canonicity(X :: SurfaceWithTorusActionFixedPoint)
resolution_graph(X :: SurfaceWithTorusActionFixedPoint)
```
