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
MoriDreamSpace
CStarSurfaceCase
CStarSurface
ToricSurface
SurfaceWithTorusAction
MoriDreamSpaceUnion
```

## Constructors 

```@docs
cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Symbol)
cstar_surface(P :: ZZMatrix)
toric_surface(vs :: Vector{Vector{T}}) where {T <: IntegerUnion}
toric_surface(P :: ZZMatrix)
```

## Basic attributes

```@docs
gen_matrix(X :: MoriDreamSpaceUnion)
canonical_toric_ambient
cox_ring_relations
cox_ring(X :: MoriDreamSpace)
```

## Class group and Picard group

```@docs
class_group(X :: MoriDreamSpace)
class_group_rank(X :: MoriDreamSpace)
class_group_torsion(X :: MoriDreamSpace)
class_group_torsion_order(X :: MoriDreamSpace)
picard_group(X :: MoriDreamSpace)
degree_matrix(X :: MoriDreamSpaceUnion)
degree_matrix_free_part(X :: MoriDreamSpaceUnion)
degree_matrix_torsion_part(X :: MoriDreamSpaceUnion)
gorenstein_index(X :: MoriDreamSpace)
picard_index(X :: MoriDreamSpace)
```

## Singularities and Resolutions

```@docs
is_quasismooth(X :: MoriDreamSpace)
is_factorial(X :: MoriDreamSpace)
is_smooth(X :: MoriDreamSpace)
is_log_terminal(X :: SurfaceWithTorusAction)
log_canonicity(X :: SurfaceWithTorusAction)
singularities(X :: SurfaceWithTorusAction)
number_of_singularities(X :: SurfaceWithTorusAction)
singularity_types(X :: SurfaceWithTorusAction)
singularity_types_string
canonical_resolution(X :: CStarSurface)
canonical_resolution(X :: ToricSurface)
minimal_resolution(X :: CStarSurface)
```

## Intersection numbers

```@docs
intersection_matrix
anticanonical_self_intersection(X :: SurfaceWithTorusAction)
```

## Kaehler Einstein metrics

```@docs
admits_kaehler_einstein_metric
special_indices(X :: CStarSurface)
moment_polytopes(X :: CStarSurface)
```

## Attributes of $\mathbb{C}^*$-surfaces

```@docs
nblocks
block_sizes
slopes
has_x_plus(X :: CStarSurface)
has_x_minus(X :: CStarSurface)
has_D_plus(X :: CStarSurface)
has_D_minus(X :: CStarSurface)
is_intrinsic_quadric
```
