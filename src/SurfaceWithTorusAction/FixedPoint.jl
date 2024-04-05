
@doc raw"""
    SurfaceWithTorusActionFixedPoint 

The `Union` of [`CStarSurfaceFixedPoint`](@ref) and
[`ToricSurfaceFixedPoint`](@ref).

"""
const SurfaceWithTorusActionFixedPoint = Union{CStarSurfaceFixedPoint, ToricSurfaceFixedPoint}


@doc raw"""
    canonical_resolution(x :: SurfaceWithTorusActionFixedPoint)

Return the canonical resolution of singularities of a given fixed point on a
surface with torus action. The result is a triple `(Y, ex_div, discr)` where
`Y` is the resulting surface after the resolution step, `ex_div` contains the
exceptional divisors over `x` and `discrepancies` contains their discrepancies.

# Example

Resolving the elliptic fixed point $x^+$ of the $E_6$ singular cubic.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> (Y, ex_div, discr) = canonical_resolution(x_plus(X));

julia> gen_matrix(Y)
[-3   -1   -2   -1   3   2   1   0   0   0]
[-3   -1   -2   -1   0   0   0   2   1   0]
[-2   -1   -1    0   1   1   1   1   1   1]

julia> map(E -> E*E, ex_div)
6-element Vector{QQFieldElem}:
 -2
 -2
 -2
 -2
 -2
 -2

julia> discr
6-element Vector{Rational{Int64}}:
 0//1
 0//1
 0//1
 0//1
 0//1
 0//1

```

"""
function canonical_resolution end


@doc raw"""
    singularity_type(x :: SurfaceWithTorusActionFixedPoint)

Return the singularity type of a fixed point on a surface with torus action.

# Examples

The $E_6$ singular cubic.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> singularity_type(x_plus(X))
E6
```

A non-log terminal singularity.

```jldoctest
julia> X = cstar_surface([[5, 7],[3],[2]], [[-1, 2], [1], [-1]], :ee)
C-star surface of type (e-e)

julia> singularity_type(x_plus(X))
Non log terminal singularity
```


"""
function singularity_type end


@doc raw"""
    number_of_exceptional_prime_divisors(x :: SurfaceWithTorusActionFixedPoint)

Return the number of exceptional prime divisors over a given point in the
minimal resolution of singularities.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> number_of_exceptional_prime_divisors(x_plus(X))
6
```

"""
@attr number_of_exceptional_prime_divisors(x :: SurfaceWithTorusActionFixedPoint) =
number_of_exceptional_prime_divisors(singularity_type(x))

@doc raw"""
    singularity_kind(x :: SurfaceWithTorusActionFixedPoint)

Return one of the symbols :A, :D or :E, depending on the singularity type
of the point `x`.

"""
function singularity_kind end
 

@doc raw"""
    is_log_terminal(x :: SurfaceWithTorusActionFixedPoint)

Check whether a point on a surface with torus action is at most a log terminal
singularity.

# Example

The $E_6$ singular cubic.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> is_log_terminal(x_plus(X))
true
```

A non-log terminal singularity.

```jldoctest
julia> X = cstar_surface([[5, 7],[3],[2]], [[-1, 2], [1], [-1]], :ee)
C-star surface of type (e-e)

julia> is_log_terminal(x_plus(X))
false
```

"""
function is_log_terminal end


@attr cox_coordinates(x :: SurfaceWithTorusActionFixedPoint) = 
[i ∈ orbit_cone(x) ? 0 : 1 for i = 1 : nrays(parent(x))]


@doc raw"""
    log_canonicity(X :: SurfaceWithTorusActionFixedPoint)

Return the maximal rational number $\varepsilon$ such that a given point on a
surface with torus action is $\varepsilon$-log canonical. By definition, this
is the minimal discrepancy in the resolution of singularities plus one.

# Example

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> log_canonicity(x_plus(X))
1//1
```

"""
@attr log_canonicity(x :: SurfaceWithTorusActionFixedPoint) = 
minimum([[0 // 1] ; canonical_resolution(x)[3]]) + 1


@doc raw"""
    resolution_graph(x :: SurfaceWithTorusActionFixedPoint)

Return the resolution graph of the minimal resolution at a given fixed point of
a surface with torus action. The result is a pair with first entry a
`Graphs.SimpleGraph` and second entry the list of self intersection numbers of
the exceptional divisors, which serve as node labels of the graph.

# Example

The $E_6$ singular cubic surface.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> (graph, nodelabel) = resolution_graph(x_plus(X))
(Graphs.SimpleGraphs.SimpleGraph{Int64}(11, [[1, 2], [1, 2, 6], [3, 4], [3, 4, 6], [5, 6], [2, 4, 5, 6]]), Nemo.QQFieldElem[-2, -2, -2, -2, -2, -2])
```

The resolution graph can be visualized with `GraphPlot.jl`:

```julia
julia> using GraphPlot

julia> gplothtml(graph, nodelabel = nodelabel)

```

"""
@attr function resolution_graph(x :: SurfaceWithTorusActionFixedPoint)
    (Y, divs, _) = minimal_resolution(x)
    M = Matrix(intersection_matrix(Y))
    inds = map(is_prime_with_index, divs)
    adj_matrix = M[inds, inds]
    nodelabel = [adj_matrix[k, k] for k = 1 : length(divs)]
    return (Graphs.SimpleGraph(adj_matrix), nodelabel)
end


                                            
