
# Normal forms

## Normal form of toric surfaces

```@docs
CStarSurfaces.normal_form(X :: ToricSurface)
CStarSurfaces.is_normal_form(X :: ToricSurface)
are_isomorphic(X :: ToricSurface, Y :: ToricSurface)
```

## Types of admissible operations

```@docs
AdmissibleOperation
InvertLastRow
InvertLastRow(factor :: Int)
PermutationOfRays
PermutationOfBlocks
AdmissibleRowOperation
CompositeAdmissibleOperation
normalize_admissible_operation(γ :: CompositeAdmissibleOperation)
```

## Normal form of $\mathbb{C}^*$-surfaces

```@docs
beta_plus
beta_minus
beta_plus_sorted
beta_minus_sorted
orientation
CStarSurfaces.normal_form(X :: CStarSurface)
CStarSurfaces.is_normal_form(X :: CStarSurface)
are_isomorphic(X :: CStarSurface, Y :: CStarSurface)
```
