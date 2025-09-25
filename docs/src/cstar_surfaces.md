# Surfaces
 
## Cases

We model the four possible cases ``\mathrm{(ee)}, \mathrm{(pe)}, \mathrm{(ep)}``
and ``\mathrm{(pp)}`` by a simple enum type with four values:

```@docs
CStarSurfaceCase
invert_case
has_elliptic_fixed_point_plus(C :: CStarSurfaceCase)
has_elliptic_fixed_point_minus(C :: CStarSurfaceCase)
has_parabolic_fixed_point_curve_plus(C :: CStarSurfaceCase)
has_parabolic_fixed_point_curve_minus(C :: CStarSurfaceCase)
```

## The `CStarSurface` type

Recall from Section ``\ref{sec:background_cstar_surfaces}`` that a
``\mathbb{C}^*``-surface can be described by a defining triple ``(\mathfrak{c},
l, d)`` together with a coefficient matrix ``A \in \mathbb{C}^{2 \times n}``.
In `CStarSurfaces.jl`, we model only the defining triple and ignore the
coefficient matrix. This is because almost all numerical invariants (Gorenstein
index, Picard index, log canonicity etc.) only depend on the defining triple
anyway, so there isn't much use in distinguishing between
``\mathbb{C}^*``-surfaces having equivalent defining triples, but different
coefficient matrices.

Since Julia's indexing of vectors is one-based, the convention for numbering the
blocks differs from the one used throughout Chapter
``\ref{chp:log_del_pezzo_c_star_surfaces}``. In particular, we write ``R`` for
the number of blocks, which equals ``r+1`` in the notation of Construction
``\ref{cns:generator_matrix_cstar_surface}``. For some user-facing functions,
like [`l`](@ref), [`d`](@ref) and [`block_sizes`](@ref), we introduce an offset
to match the convention used in Chapter
``\ref{chp:log_del_pezzo_c_star_surfaces}``.

Like the polygon type in `RationalPolygons.jl`, we use statically sized matrices
for the defining data of a ``\mathbb{C}^*``-surface. This means the type will
depend both on the number of blocks ``r`` and the total number of rays ``n = n_0
+ \dots + n_r`` as a type parameter.

```@docs
CStarSurface
vertex_matrix(X :: CStarSurface)
l
d
number_of_blocks
block_sizes
case
has_elliptic_fixed_point_plus(X :: CStarSurface)
has_elliptic_fixed_point_minus(X :: CStarSurface)
has_parabolic_fixed_point_curve_plus(X :: CStarSurface)
has_parabolic_fixed_point_curve_minus(X :: CStarSurface)
ray
embedded_ray
generator_matrix
slope
sum_of_maximal_slopes
sum_of_minimal_slopes
l_plus
l_minus
```

## Global properties

```@docs
class_group
multiplicity(X :: CStarSurface)
picard_index(X :: CStarSurface)
gorenstein_index(X :: CStarSurface)
is_quasismooth
is_factorial
is_smooth
log_canonicity
is_log_canonical
is_log_terminal
degree(X :: CStarSurface)
```

## Normal form
