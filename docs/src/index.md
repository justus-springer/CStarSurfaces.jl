This appendix serves as a reference for `CStarSurfaces.jl`, a Julia package for
computations with rational projective ``\mathbb{C}^*``-surfaces using their
combinatorial description. It is built on top of `RationalPolygons.jl`, which is
described in Appendix ``\ref{apx:julia_rational_polygons}``. 

Let us give an impression of `CStarSurfaces.jl` by means of an example session.
After loading in the package, we can create a `CStarSurface` by providing a
defining triple as in Construction ``\ref{cns:generator_matrix_cstar_surface}``.


```@repl quick_start
using CStarSurfaces
X = cstar_surface(PE, [[1,1],[4,2],[2]], [[-1,-2],[3,-3],[1]])
```

Let us see its generator matrix and grading matrix:

```@repl quick_start
generator_matrix(X)
grading_matrix_free_part(X)
grading_matrix_torsion_part(X)
```

Next, we compute basic global properties, like the class group, Gorenstein index,
Picard index, log canonicity and degree:

```@repl quick_start
class_group(X)
gorenstein_index(X)
picard_index(X)
log_canonicity(X)
degree(X)
```

Let us now consider the fixed points of the ``\mathbb{C}^*``-surface. In this case,
there is the elliptic fixed point ``x^-``, the hyperbolic fixed points
``x_{01}`` and ``x_{11}`` and the parabolic fixed points ``x^+_0, x^+_1`` and
``x^+_2``:

```@repl quick_start
xs = fixed_points(X)
```

We can arrange their local properties into a pretty table as follows:

```@repl quick_start
using PrettyTables
fs = [(X,x) -> x, class_group, gorenstein_index, log_canonicity, is_smooth, is_canonical];
header = ["Point", "Class group", "Gorenstein index", "Log canonicity", "smooth?", "canonical?"];
pretty_table([f(X,x) for x in xs, f in fs]; header)
```

Here, a log canonicity value of `1//0` (infinity by Julia's convention) is used
for smooth points.

This documentation is generated directly from the docstrings of the package's
source code. An online version is available on its GitHub page
[CStarSurfaces_jl](@cite). Every documented item includes a clickable
[`source`](https://github.com/justus-springer/CStarSurfaces.jl) link that directs to
the corresponding line in the source code where it is defined. The content here
refers to version `v1.0.0-thesis` of `CStarSurfaces.jl`. All example sessions have been
tested and verified to work with this version.

```@bibliography
```
