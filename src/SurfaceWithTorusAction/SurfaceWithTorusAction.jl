
@doc raw"""
    SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

Julia type for surfaces with a non-trivial torus action.

"""
const SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

dim(::SurfaceWithTorusAction) = 2


@doc raw"""
    picard_index(X :: SurfaceWithTorusAction)   

Return the index of the Picard group in the class group of a surface with torus
action.

This uses the formula from [Sp23](@cite) and is faster than the more general
implementation.

"""
picard_index(X :: SurfaceWithTorusAction) = prod(values(local_picard_indices(X))) / class_group_torsion_order(X)


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
    exceptional_divisors(X :: SurfaceWithTorusAction)   

Return the exceptional divisors in the resolution of singularities of a surface
with torus action.

"""
exceptional_divisors(X :: SurfaceWithTorusAction) = canonical_resolution(X)[2]


@doc raw"""
    discrepancies(X :: SurfaceWithTorusAction)

Return the discrepancies associated to the exceptional divisors in the
resolution of singularities of a surface with torus action.

"""
discrepancies(X :: SurfaceWithTorusAction) = canonical_resolution(X)[3]


@doc raw"""
    maximal_log_canonicity(X :: SurfaceWithTorusAction)

Given a surface with torus action $X$, return the maximal rational number
$\varepsilon$ such that $X$ is $\varepsilon$-log canonical. By definition,
this is the minimal discrepancy in the resolution of singularities plus one.

# Example

```jldoctest
julia> maximal_log_canonicity(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
1//1
```

"""
@attr function maximal_log_canonicity(X :: SurfaceWithTorusAction) 
    # we add a superficial zero into the list of discrepancies to ensure a
    # well-defined (and correct) result in case there are no exceptional rays
    # (i.e. the surface is already smooth).
    ds = vcat([0], [d for (_,d) in discrepancies(X)]...)
    # the maximal log canonicity equals the minimal discrepancy plus one
    return minimum(ds) + 1
end
