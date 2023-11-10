
# Divisors

A divisor $D$ on a Mori dream space $X$ is modeled by the integral vector
$a=(a_1, \dots, a_r)$ such that $D = a_1 D^1_X + ... + a_r D^r_X$, where
$D^i_X$ are the restrictions of the torus invariant prime divisors of the
canonical toric ambient variety. For divisors on $\mathbb{C}^*$-surfaces,
double index notation can be used for the coefficients.

## Types

```@docs
MoriDreamSpaceDivisor
CStarSurfaceDivisor
ToricSurfaceDivisor
SurfaceWithTorusActionDivisor
```

## Constuctors

```@docs
mori_dream_space_divisor(X :: T, coeffs :: Vector{S}) where {T <: MoriDreamSpace, S <: IntegerUnion}
mori_dream_space_divisor(X :: T, td :: ToricDivisor) where {T <: MoriDreamSpace}
cstar_surface_divisor(X :: CStarSurface{EE}, coeffs :: DoubleVector{<:IntegerUnion})
cstar_surface_divisor(X :: CStarSurface{PE}, coeffs :: DoubleVector{T}, coeff_plus :: T) where {T <: IntegerUnion}
cstar_surface_divisor(X :: CStarSurface{EP}, coeffs :: DoubleVector{T}, coeff_minus :: T) where {T <: IntegerUnion}
cstar_surface_divisor(X :: CStarSurface{PP}, coeffs :: DoubleVector{T}, coeff_plus :: T, coeff_minus :: T) where {T <: IntegerUnion}
invariant_divisor(X :: CStarSurface, i :: Int, j :: Int)
D_plus
D_minus
toric_surface_divisor(X :: ToricSurface, coeffs :: Vector{S}) where {S <: IntegerUnion}
invariant_divisor(X :: ToricSurface, i :: Int)
```

## (Anti)canonical divisor

```@docs
canonical_divisor(X :: MoriDreamSpace)
anticanonical_divisor(X :: MoriDreamSpace)
canonical_divisor_class(X :: MoriDreamSpace)
anticanonical_divisor_class(X :: MoriDreamSpace)
```

## Attributes

```@docs
coefficients(d :: MoriDreamSpaceDivisor)
toric_divisor(d :: MoriDreamSpaceDivisor)
double_coefficients(d :: CStarSurfaceDivisor)
is_prime(d :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace}
```

## Intersection numbers

```@docs
Base.:*(d1 :: SurfaceWithTorusActionDivisor, d2 :: SurfaceWithTorusActionDivisor)
```

## Contraction

```@docs
contract_prime_divisor(d :: CStarSurfaceDivisor{T}) where {T <: CStarSurfaceCase}
```


