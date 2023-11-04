```@meta
CurrentModule = CStarSurfaces
```

```@setup oscar
using CStarSurfaces, Oscar 
```

# Surfaces with torus action

The two main Julia types in this packages are [`CStarSurface`](@ref) and
[`ToricSurface`](@ref). Some functionality also works for
[`MoriDreamSpace`](@ref)'s, an abstract type of which `CStarSurface` and
`ToricSurface` are subtypes.

## Julia types

```@docs
CStarSurfaceCase
CStarSurface
ToricSurface
SurfaceWithTorusAction
```

## Constructors 

```@docs
cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Symbol)
cstar_surface(P :: ZZMatrix)
toric_surface
gen_matrix
```

## Class group and Picard group

```@docs
class_group(X :: MoriDreamSpace)
class_group_rank
class_group_torsion
class_group_torsion_order
local_class_groups(X :: MoriDreamSpace)
local_class_group(X :: MoriDreamSpace, c :: Vector{Int64})
maps_from_class_group_to_local_class_groups(X :: MoriDreamSpace)
map_from_class_group_to_local_class_group(X :: MoriDreamSpace, c :: Vector{Int64})
picard_group(X :: MoriDreamSpace)
degree_matrix
degree_matrix_free_part
degree_matrix_torsion_part
gorenstein_index(X :: MoriDreamSpace)
local_gorenstein_indices
local_gorenstein_index
picard_index(X :: MoriDreamSpace)
local_picard_indices
local_picard_index
is_factorial
```

## Divisors

```@docs
CStarSurfaceDivisor
cstar_surface_divisor
ToricSurfaceDivisor
toric_surface_divisor
SurfaceWithTorusActionDivisor
canonical_divisor(X :: MoriDreamSpace)
anticanonical_divisor(X :: MoriDreamSpace)
canonical_divisor_class(X :: MoriDreamSpace)
anticanonical_divisor_class(X :: MoriDreamSpace)
Base.:*(d1 :: SurfaceWithTorusActionDivisor, d2 :: SurfaceWithTorusActionDivisor)
```

## Intersection numbers

```@docs
intersection_matrix
anticanonical_self_intersection
```

## Singularities and resolutions

```@docs
is_quasismooth
canonical_resolution
exceptional_rays
discrepancies
maximal_log_canonicity
```

## Kaehler Einstein metrics

```@docs
admits_kaehler_einstein_metric
special_indices
moment_polytopes
```

## Attributes of $\mathbb{C}^*$-surfaces

```@docs
canonical_toric_ambient(X :: CStarSurface)
cox_ring_vars(X :: CStarSurface)
cox_ring_relations(X :: CStarSurface)
has_x_plus(X :: CStarSurface)
has_x_minus(X :: CStarSurface)
has_D_plus(X :: CStarSurface)
has_D_minus(X :: CStarSurface)
x_plus
x_minus
elliptic_fixed_points
hyperbolic_fixed_points
parabolic_fixed_points
D_plus
D_minus
parabolic_fixed_point_curves
number_of_parabolic_fixed_point_curves
nblocks
block_sizes
slopes
is_intrinsic_quadric
```
