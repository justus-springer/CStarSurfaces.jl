# Mori Dream Spaces

This package is a about $\mathbb{C}^*$-surfaces. However, some algorithms here
also work for a more general class of varieties admitting a certain embedding
into a toric variety, called Mori dream spaces. To allow for the treatment of
these varieties in the future, this package provides the abstract
[`MoriDreamSpace`](@ref) type, whose only implementations currently are
[`CStarSurface`](@ref) and [`ToricSurface`](@ref).

```@docs
MoriDreamSpace
MoriDreamSpaceUnion
canonical_toric_ambient
cox_ring_relations
cox_ring(X :: MoriDreamSpace)
is_toric(X :: MoriDreamSpace)
rays(X :: MoriDreamSpace)
nrays(X :: MoriDreamSpace)
maximal_cones(X :: MoriDreamSpace)
maximal_cones_indices(X :: MoriDreamSpace)
```
