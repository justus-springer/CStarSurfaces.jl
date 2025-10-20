# Fixed points and local properties

Recall that there are three kinds of fixed points of ``\mathbb{C}^*``-surfaces:
elliptic, hyperbolic and parabolic. In `CStarSurfaces.jl`, we model fixed points
by a type hierarchy with an abstract type `FixedPoint` at the top. This allows
us to implement local properties like [`gorenstein_index`](@ref) and
[`log_canonicity`](@ref) separately for each type of fixed point using Julia's
dispatch mechanism. We found that sometimes, hyperbolic and parabolic fixed
point can be treated uniformly, while elliptic fixed points need to be treated
differently. Hence we will refer to hyperbolic and parabolic fixed points
collectively as *bolic* fixed points.

In Subsection ``\ref{doc:Fixed-points}``, we go over the type hierarchy of fixed
points. In Subsection ``\ref{doc:Local-properties}``, we discuss local properties.

## Fixed points

The following graph summarizes the different fixed point types and their subtype
relations. Abstract types have dashed boundaries.

FIXEDPOINTGRAPH

```@docs
FixedPoint
EllipticFixedPoint
EllipticFixedPointPlus
EllipticFixedPointMinus
BolicFixedPoint
HyperbolicFixedPoint
ParabolicFixedPoint
ParabolicFixedPointPlus
ParabolicFixedPointMinus
elliptic_fixed_points
hyperbolic_fixed_points
parabolic_fixed_points
bolic_fixed_points
fixed_points
```

## Local properties

```@docs
toric_chart
class_group
multiplicity
is_factorial
is_quasismooth
is_smooth
gorenstein_index
log_canonicity
is_log_terminal
is_log_canonical
is_terminal
is_canonical
```
