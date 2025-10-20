@doc raw"""
    FixedPoint

Abstract supertype of a fixed point on a ``\mathbb{C}^*``-surface. Here,
a fixed point is understood to be a formal symbol associated to the
defining triple, see Definition ``\ref{def:defining_triple_fixed_points}``.

This type has the two subtypes [`EllipticFixedPoint`](@ref) and [`BolicFixedPoint`](@ref).

"""
abstract type FixedPoint end


@doc raw"""
    fixed_points(X :: CStarSurface)

Return all fixed points of ``X``, see Definition ``\ref{def:defining_triple_fixed_points}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> fixed_points(X)
6-element Vector{FixedPoint}:
 x⁻
 x⁺₀
 x⁺₁
 x⁺₂
 x₀₁
 x₁₁
```

"""
function fixed_points(X :: CStarSurface)
    xs = FixedPoint[]
    append!(xs, elliptic_fixed_points(X))
    append!(xs, bolic_fixed_points(X))
    return xs
end


@doc raw"""
    toric_chart(X :: CStarSurface, x :: FixedPoint)

Return the generator matrix of the toric chart around a fixed point.
This is the local generator matrix, as in Definition ``\ref{def:defining_triple_local_generator_matrix}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> toric_chart(X, EllipticFixedPointMinus())
3×3 SMatrix{3, 3, Int64, 9} with indices SOneTo(3)×SOneTo(3):
 -1   1  0
 -1   0  2
 -4  -1  1
```

"""
function toric_chart end


@doc raw"""
    gorenstein_index(X :: CStarSurface, x :: FixedPoint)

Return the local Gorenstein index at the point ``x``. See Propositions 8.8
and 8.9 of [HaHaSp25](@cite) for their formulas in terms of defining data.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> gorenstein_index(X, EllipticFixedPointMinus())
3
```

"""
function gorenstein_index end


@doc raw"""
    multiplicity(X :: CStarSurface, x :: FixedPoint)

Return the order of the local class group at the point ``x``, see Definition
``\ref{def:defining_triple_local_multiplicity}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> multiplicity(X, EllipticFixedPointMinus())
9
```

"""
function multiplicity end

multiplicity(X :: CStarSurface, x :: FixedPoint) =
abs(det_bareiss(toric_chart(X, x)))


@doc raw"""
    is_factorial(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is factorial, i.e. its multiplicity is one.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> filter(x -> is_factorial(X,x), fixed_points(X))
4-element Vector{FixedPoint}:
 x⁺₀
 x⁺₁
 x₀₁
 x₁₁
```

"""
function is_factorial end

is_factorial(X :: CStarSurface, x :: FixedPoint) =
multiplicity(X, x) == 1


@doc raw"""
    is_quasismooth(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is quasismooth. Bolic fixed points are always
quasismooth. See Summary 8.1 of [HaHaSp25](@cite).

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> all(x -> is_quasismooth(X,x), fixed_points(X))
true
```

"""
function is_quasismooth end


@doc raw"""
    is_smooth(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is smooth. This is equivalent to being
factorial and quasismooth.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> filter(x -> is_smooth(X,x), fixed_points(X))
4-element Vector{FixedPoint}:
 x⁺₀
 x⁺₁
 x₀₁
 x₁₁
```

"""
function is_smooth end

is_smooth(X :: CStarSurface, x :: FixedPoint) =
is_factorial(X, x) && is_quasismooth(X, x)


@doc raw"""
    class_group(X :: CStarSurface, x :: FixedPoint)

Return the local class group at the point ``x``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> class_group(X, EllipticFixedPointMinus())
Z/9
```

"""
function class_group end

class_group(X :: CStarSurface, x :: FixedPoint) =
cokernel(toric_chart(X, x))


@doc raw"""
    log_canonicity(X :: CStarSurface, x :: FixedPoint)

Return the maximal ``\varepsilon > 0`` such that the singularity at
``x`` is ``\varepsilon``-log canonical. By definition, this is set to
be `1//0` (infinity) for smooth points.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> [log_canonicity(X,x) for x in fixed_points(X)]
6-element Vector{Rational{Int64}}:
 1//3
 1//0
 1//0
  1
 1//0
 1//0
```

"""
function log_canonicity end

function log_canonicity(X :: CStarSurface, x :: FixedPoint)
    ds = log_canonicities(X, x)
    if isempty(ds)
        return 1//0
    else
        d = minimum(ds)
        # For any value above 1, fall back to infinity
        return d > 1 ? 1//0 : d
    end
end


@doc raw"""
    is_log_terminal(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0)

Check whether a point on a ``\mathbb{C}^*``-surface is ``\varepsilon``-log terminal.

"""
function is_log_terminal end

is_log_terminal(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0) =
log_canonicity(X, x) > ε


@doc raw"""
    is_log_canonical(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0)

Check whether a point on a ``\mathbb{C}^*``-surface is ``\varepsilon``-log canonical.

"""
function is_log_canonical end

is_log_canonical(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0) =
log_canonicity(X, x) ≥ ε


@doc raw"""
    is_terminal(X :: CStarSurface, x :: FixedPoint)

Check whether a point on a ``\mathbb{C}^*``-surface is terminal.
This is equivalent to being smooth.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> filter(x -> is_terminal(X,x), fixed_points(X))
4-element Vector{FixedPoint}:
 x⁺₀
 x⁺₁
 x₀₁
 x₁₁
```

"""
function is_terminal end

is_terminal(X :: CStarSurface, x :: FixedPoint) =
is_log_terminal(X, x, 1)


@doc raw"""
    is_canonical(X :: CStarSurface, x :: FixedPoint)

Check whether a point on a ``\mathbb{C}^*``-surface is canonical.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[1,1],[2]], [[-3,-4],[0,-1],[1]])
C*-surface of case PE with l = ((1,1),(1,1),2) and d = ((-3,-4),(0,-1),1) 

julia> filter(x -> is_canonical(X,x), fixed_points(X))
5-element Vector{FixedPoint}:
 x⁺₀
 x⁺₁
 x⁺₂
 x₀₁
 x₁₁
```

"""
function is_canonical end

is_canonical(X :: CStarSurface, x :: FixedPoint) =
is_log_canonical(X, x, 1)
