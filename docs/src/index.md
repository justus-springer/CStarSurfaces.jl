```@meta
CurrentModule = CStarSurfaces
```

# CStarSurfaces.jl

A computer algebra package for rational $\mathbb{C}^*$-surfaces in the Julia
programming language. This package makes use of the
[OSCAR](https://www.oscar-system.org) Computer Algebra System.

The approach to $\mathbb{C}^*$-surfaces relies on the general combinatorial
theory of varieties with finitely generated Cox ring developed in
[BeHa07](@cite), [Ha08](@cite) and its specialization to varieties with torus
action initiated in [HaHe13](@cite) and [HaSu10](@cite). As an introductory
reference, we mention [ArDeHaLa15](@cite).

## Features

- Canonical toric ambient variety, divisor class group, local class
  groups, (anti)canonical divisor class,
- Cox Ring, Gorenstein index, Picard index,
- intersection numbers, anticanonical self intersection,
- resolution of singularities, exceptional divisors and discrepancies,
- normal form of defining data, admissible operations, isomorphy test. 

## Installation

[CStarSurfaces.jl](https://github.com/justus-springer/CStarSurfaces.jl) is
available in the [General
Registry](https://github.com/JuliaRegistries/General), hence can be
installed with

```julia
julia> ]add CStarSurfaces
```

## Quick start

We work in the notation of [ArDeHaLa15; Section 5.4](@cite).

Import both [Oscar](https://www.oscar-system.org) and [CStarSurfaces.jl](https://github.com/justus-springer/CStarSurfaces.jl) to get started:

```julia
julia> using Oscar, CStarSurfaces
```

There are essentially two constructors for C-Star surfaces: The first takes the
integral vectors $l_i=(l_{i1}, \dots, l_{in_i})$ and $d_i=(d_{i1}, \dots,
d_{in_i})$ and one of the four symbols `:ee, :pe, :ep, :pp`. The second takes
the generating matrix `P` of the correct shape:

```julia
julia> X = cstar_surface([[1, 1], [4], [4]], [[0, -2], [3], [3]], :ee)
C-star surface of type (e-e)

julia> Y = cstar_surface(ZZ[-1 -1 4 0 ; -1 -1 0 4 ; 0 -2 3 3])
C-star surface of type (e-e)

julia> X == Y
true
```

`gen_matrix` returns the generating matrix (P-Matrix) of a C-star surface:

```julia
julia> gen_matrix(X)
[-1   -1   4   0]
[-1   -1   0   4]
[ 0   -2   3   3]
```

`canonical_toric_ambient` returns the canonical toric ambient variety of a
C-star surface, as an Oscar type:

```julia
julia> Z = canonical_toric_ambient(X)
Normal toric variety
```

We compute some geometric invariants of $X$:

```julia
julia> class_group(X)
GrpAb: Z/8 x Z

julia> cox_ring(X)
Quotient
  of graded multivariate polynomial ring in 4 variables over QQ
  by ideal(T[0][1]*T[0][2] + T[1][1]^4 + T[2][1]^4)

julia> gorenstein_index(X)
3

julia> picard_index(X)
48

julia> K = anticanonical_divisor(X)
CStarSurfaceDivisor{EE}(C-star surface of type (e-e), Torus-invariant, non-prime divisor on a normal toric variety)

julia> K * K # the anticanonical self intersection
2//3

julia> (Xreg, exceptional_divisors, discrepancies) = canonical_resolution(X)
[...]

julia> gen_matrix(Xreg)
[-1   -1   -1   4   1   3   2   1   0   0   0   0   0   0    0]
[-1   -1   -1   0   0   0   0   0   4   1   3   2   1   0    0]
[ 0   -2   -1   3   1   2   1   0   3   1   2   1   0   1   -1]

julia> maximal_log_canonicity(X)
1//3
```

## References

```@bibliography
```

