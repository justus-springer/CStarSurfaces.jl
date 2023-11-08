
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
    minimal_resolution_exceptional_divisors(X :: SurfaceWithTorusAction)   

Return the exceptional divisors in the minimal resolution of singularities of a surface
with torus action.

"""
minimal_resolution_exceptional_divisors(X :: SurfaceWithTorusAction) = minimal_resolution(X)[2]


@doc raw"""
    minimal_resolution_discrepancies(X :: SurfaceWithTorusAction)

Return the discrepancies associated to the exceptional divisors in the
minimal resolution of singularities of a surface with torus action.

"""
minimal_resolution_discrepancies(X :: SurfaceWithTorusAction) = minimal_resolution(X)[3]


@doc raw"""
    canonical_resolution_exceptional_divisors(X :: SurfaceWithTorusAction)   

Return the exceptional divisors in the canonical resolution of singularities of a surface
with torus action.

"""
canonical_resolution_exceptional_divisors(X :: SurfaceWithTorusAction) = canonical_resolution(X)[2]


@doc raw"""
    canonical_resolution_discrepancies(X :: SurfaceWithTorusAction)

Return the discrepancies associated to the exceptional divisors in the
canonical resolution of singularities of a surface with torus action.

"""
canonical_resolution_discrepancies(X :: SurfaceWithTorusAction) = canonical_resolution(X)[3]


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
    resolution_graphs(X :: SurfaceWithTorusAction)

Return the resolution graphs of the minimal resolution of a surface with
torus action. The result is a dictionary indexed by the fixed points of $X$,
where each value is a pair with first entry a `Graphs.SimpleGraph` and second
entry the list of self intersection numbers of the exceptional divisors,
which serve as node labels of the graph.

# Example

Drawing the resolution graph of the $E_6$ singular cubic surface using 
`GraphPlot.jl`.

```jldoctest
julia> import Graphs, GraphPlot

julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> res_graphs = resolution_graphs(X);

julia> (graph, nodelabel) = res_graphs[x_plus(X)];

julia> Graphs.nv(graph)
6

julia> GraphPlot.gplothtml(graph; nodelabel = nodelabel); # Opens a browser window displaying the graph

```

"""
@attr function resolution_graphs(X :: SurfaceWithTorusAction)
    (Y, ex_div, _) = minimal_resolution(X)
    M = Matrix(intersection_matrix(Y))
    res_graphs = Dict{Vector{Int}, Tuple{Graphs.SimpleGraph, Vector}}()
    for (x, divs) in ex_div
        inds = map(is_prime_with_index, divs)
        adj_matrix = M[inds, inds]
        nodelabel = [adj_matrix[k, k] for k = 1 : length(divs)]
        res_graphs[x] = (Graphs.SimpleGraph(adj_matrix), nodelabel)
    end
    return res_graphs
end


@doc raw"""
    resolution_graph(X :: SurfaceWithTorusAction, x :: Vector{Int})   

Return the resolution graph of the minimal resolution of singularities at a
given point of a surface with torus action. The point must be given as an index
vector of the corresponding maximal cone. The result is a pair with first entry
a `Graphs.SimpleGraph` and second entry the list of self intersection numbers
of the exceptional divisors, which serve as node labels of the graph.


# Example

Drawing the resolution graph of the $E_6$ singular cubic surface using 
`GraphPlot.jl`.

```jldoctest
julia> import Graphs, GraphPlot

julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> (graph, nodelabel) = resolution_graph(X, x_plus(X))

julia> Graphs.nv(graph)
6

julia> GraphPlot.gplothtml(graph; nodelabel = nodelabel); # Opens a browser window displaying the graph

```

"""
resolution_graph(X :: SurfaceWithTorusAction, x :: Vector{Int}) = resolution_graphs(X)[x]
