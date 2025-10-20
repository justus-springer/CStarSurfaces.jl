# CStar Surfaces
 
This section documents the main functionality around ``\mathbb{C}^*``-surfaces
and their global properties. In Subsection ``\ref{doc:The-CStarSurface-type}``,
we introduce the type `CStarSurface`, which is modeled by defining triples.
Subsection ``\ref{doc:Global-properties}`` goes over functions computing global
invariants. Finally, Subsection
``\ref{doc:Admissible-operations-and-the-normal-form}`` is about the normal
form.

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
blocks used by `CStarSurfaces.jl` differs from the one used throughout Chapter
``\ref{chp:log_del_pezzo_c_star_surfaces}``. In particular, we write ``R`` for the
number of blocks, which equals ``r+1`` in the notation of Construction
``\ref{cns:generator_matrix_cstar_surface}``. For some user-facing functions,
like [`l`](@ref), [`d`](@ref) and [`block_sizes`](@ref), we add an offset to
match the convention used in Chapter
``\ref{chp:log_del_pezzo_c_star_surfaces}``.

Like the polygon type in `RationalPolygons.jl`, we use statically sized matrices
for the defining data of a ``\mathbb{C}^*``-surface. This means the type will
depend both on the number of blocks and the total number of rays ``n = n_0 + \dots + n_r`` as a type parameter.

```@docs
CStarSurface
cstar_surface
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

We describe functions that compute global properties of
``\mathbb{C}^*``-surfaces. For the [`class_group`](@ref), [`multiplicity`](@ref)
and [`picard_index`](@ref), we refer to Section
``\ref{subsec:defining_triple_picard_group}`` and [HaHaSp25](@cite) for more
details. For other properties, like [`gorenstein_index`](@ref),
[`is_quasismooth`](@ref), [`is_log_terminal`](@ref), [`log_canonicity`](@ref),
and [`degree`](@ref), we refer to [HaHaSp25](@cite).

```@docs
class_group(X :: CStarSurface)
multiplicity(X :: CStarSurface)
grading_matrix_free_part(X :: CStarSurface)
grading_matrix_torsion_part(X :: CStarSurface)
grading_matrix(X :: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R}
picard_index(X :: CStarSurface)
gorenstein_index(X :: CStarSurface)
is_quasismooth(X :: CStarSurface)
is_factorial(X :: CStarSurface)
is_smooth(X :: CStarSurface)
log_canonicity(X :: CStarSurface)
is_log_canonical(X :: CStarSurface)
is_log_terminal(X :: CStarSurface)
degree(X :: CStarSurface{T}) where {T <: Integer}
```

## Admissible operations and the normal form

We provide an implementation of the normal form for defining triples from
Section ``\ref{sec:cstar_surface_normal_form}``.

```@docs
mfrak_plus
mfrak_minus
beta_plus
beta_minus
are_equivalent
orientation
is_normal_form
AdmissibleOperation
admissible_operation
inversion
permutation
addition
normal_form_with_operation
normal_form
```
