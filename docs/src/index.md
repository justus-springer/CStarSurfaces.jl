```@meta
CurrentModule = CStarSurfaces
```

```@setup oscar
using CStarSurfaces, Oscar 
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
installed by typing `]add CStarSurfaces` into a Julia REPL.

## Quick start

We work in the notation of [ArDeHaLa15; Section 5.4](@cite).

Import both [Oscar](https://www.oscar-system.org) and [CStarSurfaces.jl](https://github.com/justus-springer/CStarSurfaces.jl) to get started:

```@repl
using Oscar, CStarSurfaces
```

There are essentially two constructors for C-Star surfaces: The first takes the
integral vectors $l_i=(l_{i1}, \dots, l_{in_i})$ and $d_i=(d_{i1}, \dots,
d_{in_i})$ and one of the four symbols `:ee, :pe, :ep, :pp`. The second takes
the generating matrix `P` of the correct shape:

```@repl oscar
X = cstar_surface([[1, 1], [4], [4]], [[0, -2], [3], [3]], :ee)
Y = cstar_surface(ZZ[-1 -1 4 0 ; -1 -1 0 4 ; 0 -2 3 3])
X == Y
```

`gen_matrix` returns the generating matrix (P-Matrix) of a C-star surface:

```@repl oscar
gen_matrix(X)
```

`canonical_toric_ambient` returns the canonical toric ambient variety of a
C-star surface, as an Oscar type:

```@repl oscar
Z = canonical_toric_ambient(X)
```

We compute some geometric invariants of $X$:

```@repl oscar
class_group(X)
cox_ring(X)
gorenstein_index(X)
picard_index(X)
K = anticanonical_divisor(X)
K * K # the anticanonical self intersection
(Xreg, exceptional_divisors, discrepancies) = canonical_resolution(X)
gen_matrix(Xreg)
maximal_log_canonicity(X)
```

## References

```@bibliography
```
