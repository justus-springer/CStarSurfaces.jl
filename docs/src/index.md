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

## Installation and usage

To install the package, type the following into a Julia prompt:

```julia
julia> import Pkg; Pkg.add(url="https://github.com/justus-springer/CStarSurfaces")
```

## Quick start

## References

```@bibliography
```

