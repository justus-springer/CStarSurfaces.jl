
@doc raw"""
    SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

The `Union` of [`CStarSurface`](@ref) and [`ToricSurface`](@ref).

"""
const SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}


@doc raw"""
    intersection_matrix(X :: SurfaceWithTorusAction)

Return the matrix of intersection numbers of all torus invariant prime divisors
associated of a surface with torus action with each other. The result is a
rational `n` x `n` matrix, where `n = nrays(X)` and the `(i,j)`-th entry is the
intersection number of the prime divisors associated to the `i`-th and `j`-th
ray respectively.

# Examples

```jldoctest
julia> intersection_matrix(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
[1//3   1   2//3   1]
[   1   3      2   3]
[2//3   2   4//3   2]
[   1   3      2   3]
```

```jldoctest
julia> intersection_matrix(toric_surface(ZZ[1 0 -1 0 ; 0 1 -17 -1]))
[0    1   0     1]
[1   17   1     0]
[0    1   0     1]
[1    0   1   -17]
```

"""
function intersection_matrix end


@doc raw"""
    picard_index(X :: SurfaceWithTorusAction)   

Return the index of the Picard group in the class group of a surface with torus
action.

This uses the formula from [Sp23](@cite) and is faster than the more general
implementation.

"""
picard_index(X :: SurfaceWithTorusAction) = 
prod([order(class_group(x)) for x in fixed_points(X)]) / class_group_torsion_order(X)


@doc raw"""
    anticanonical_self_intersection(X :: SurfaceWithTorusAction)   

Return the self intersection number of the anticanonical divisor on a surface
with torus action.

# Example

```jldoctest
julia> anticanonical_self_intersection(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
3
```

"""
anticanonical_self_intersection(X :: SurfaceWithTorusAction) = anticanonical_divisor(X) * anticanonical_divisor(X)


@doc raw"""
    singularities(X :: SurfaceWithTorusAction)

Return the list of singular points of a given surface with torus action.

# Example

The $E_6$ singular cubic surface.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> singularities(X)
1-element Vector{CStarSurfaceFixedPoint{EE}}:
 elliptic fixed point x^+
```

"""
@attr singularities(X :: SurfaceWithTorusAction) = filter(!is_smooth, fixed_points(X))


@doc raw"""
    number_of_singularities(X :: SurfaceWithTorusAction)

Return the number of singularities of a given surface with torus action.

# Example

The $E_6$ singular cubic surface.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> number_of_singularities(X)
1
```

"""
@attr number_of_singularities(X :: SurfaceWithTorusAction) = length(singularities(X))


@doc raw"""
    singularity_types(X :: SurfaceWithTorusAction)

Return the list of singularity types of all singularities of a surface with
torus action.

# Examples

A $\mathbb{C}^*$-surface with two Gorenstein singularities.

```jldoctest
julia> X = cstar_surface([[1, 2], [3], [3]], [[-1, -3], [2], [2]], :ee)
C-star surface of type (e-e)

julia> singularity_types(X)
2-element Vector{SingularityTypeADE}:
 A2
 E6
```

A non-log terminal $\mathbb{C}*$-surface.

```jldoctest
julia> X = cstar_surface([[5, 7],[3],[2]], [[-1, 2], [1], [-1]], :ee)
C-star surface of type (e-e)

julia> singularity_types(X)
3-element Vector{SingularityType}:
 Non log terminal singularity
 E8
 A2
```

"""
@attr singularity_types(X :: SurfaceWithTorusAction) = map(singularity_type, singularities(X))


@doc raw"""
    is_log_terminal(X :: SurfaceWithTorusAction)

Check whether a surface with torus action has at most log terminal
singularities.

# Example

The $E_6$ singular cubic.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_log_terminal(X)
true
```

A non-log terminal $\mathbb{C}^*$-surface.

```jldoctest
julia> X = cstar_surface([[5, 7], [3], [2]], [[-1, 2], [1], [-1]], :ee)
C-star surface of type (e-e)

julia> is_log_terminal(X)
false
```

"""
@attr is_log_terminal(X :: SurfaceWithTorusAction) = all(is_log_terminal, fixed_points(X))


@doc raw"""
    log_canonicity(X :: SurfaceWithTorusAction)

Given a surface with torus action $X$, return the maximal rational number
$\varepsilon$ such that $X$ is $\varepsilon$-log canonical. By definition,
this is the minimal discrepancy in the resolution of singularities plus one.

# Example

```jldoctest
julia> log_canonicity(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
1//1
```

"""
@attr log_canonicity(X :: SurfaceWithTorusAction) = minimum(map(log_canonicity, fixed_points(X)))


@doc raw"""
    admits_kaehler_einstein_metric(X :: SurfaceWithTorusAction)

Checks whether a surface with torus action admits a Kaehler-Einstein metric.

# Examples

```jldoctest
julia> admits_kaehler_einstein_metric(cstar_surface([[1,1], [4], [4]], [[-1,-2], [3], [3]], :ee))
true
```

```jldoctest
julia> admits_kaehler_einstein_metric(toric_surface([[1,0], [1,5], [-2,-5]]))
true
```

"""
function admits_kaehler_einstein_metric end
