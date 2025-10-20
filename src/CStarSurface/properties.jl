@doc raw"""
    class_group(X :: CStarSurface)

Return the divisor class group of ``X``.
See Definition ``\ref{def:defining_triple_class_group}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> class_group(X)
Z^2 x Z/3
```

"""
class_group(X :: CStarSurface) = cokernel(transpose(generator_matrix(X)))


@doc raw"""
    multiplicity(X :: CStarSurface)

Return the order of the torsion part of the class group of ``X``.

"""
multiplicity(X :: CStarSurface) = torsion_order(class_group(X))


@doc raw"""
    picard_index(X :: CStarSurface)

Return the index of the Picard group inside the divisor class group of ``X``.
This equals the product of all local multiplicities divided by the global multiplicity, see
Theorem ``\ref{thm:picard_index_formula_defining_triples}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> picard_index(X)
90
```

"""
picard_index(X :: CStarSurface) = prod([multiplicity(X,x) for x in fixed_points(X)]) ÷ multiplicity(X)


@doc raw"""
    gorenstein_index(X :: CStarSurface)

Return the Gorenstein index of ``X``. This is the least common multiple of the
local Gorenstein indices, see also [`gorenstein_index`](@ref).

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> gorenstein_index(X)
15
```

"""
gorenstein_index(X :: CStarSurface) = lcm([gorenstein_index(X,x) for x in fixed_points(X)])


@doc raw"""
    is_quasismooth(X :: CStarSurface)   

Check whether ``X`` is quasismooth, i.e. its characteristic space is
smooth. See Summary 8.1 of [HaHaSp25](@cite) for a description in terms of the defining data.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> is_quasismooth(X)
true
```

"""
is_quasismooth(X :: CStarSurface) = all(x -> is_quasismooth(X, x), elliptic_fixed_points(X))


@doc raw"""
    is_factorial(X :: CStarSurface)   

Check whether ``X`` is locally factorial, i.e. all local class groups are trivial.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> is_factorial(X)
false
```

"""
is_factorial(X :: CStarSurface) = all(x -> is_factorial(X, x), fixed_points(X))


@doc raw"""
    is_smooth(X :: CStarSurface)

Check whether ``X`` is smooth. This is equivalent to being
factorial and quasismooth.

"""
is_smooth(X :: CStarSurface) = is_factorial(X) && is_quasismooth(X)


@doc raw"""
    log_canonicity(X :: CStarSurface)

Return the maximal ``\varepsilon > 0`` such that the surface has ``\varepsilon``-log canonical
singularities. For smooth surfaces, this returns `1//0`, which is infinity by Julia's convention.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> log_canonicity(X)
2//5
```

"""
log_canonicity(X :: CStarSurface) = minimum([log_canonicity(X, x) for x in fixed_points(X)])


@doc raw"""
    is_log_terminal(X :: CStarSurface, ε :: Real = 0)   

Check whether ``X`` is ``\varepsilon``-log terminal.

"""
is_log_terminal(X :: CStarSurface, ε :: Real = 0) = log_canonicity(X) > ε


@doc raw"""
    is_log_canonical(X :: CStarSurface, ε :: Real = 0)

Check whether ``X`` is ``\varepsilon``-log canonical.

"""
is_log_canonical(X :: CStarSurface, ε :: Real = 0) = log_canonicity(X) ≥ ε


@doc raw"""
    is_terminal(X :: CStarSurface)

Check whether ``X`` is terminal. This is equivalent to being smooth.

"""
is_terminal(X :: CStarSurface) = is_log_terminal(X, 1)


@doc raw"""
    is_canonical(X :: CStarSurface, x :: FixedPoint)

Check whether a ``\mathbb{C}^*``-surface is canonical.

"""
is_canonical(X :: CStarSurface) = is_log_canonical(X, 1)


@doc raw"""
    degree(X :: CStarSurface)

The self-intersection number of an anticanonical divisor of ``X``.
See Proposition 7.9 of [HaHaSp25](@cite) for a formula in terms of the defining data.

"""
function degree(X :: CStarSurface{T}) where {T <: Integer}
    # We use the formula given in Proposition 7.9 of [HHS23].
    r = number_of_blocks(X) - 1
    l = [[ray(X, i, j)[1] for j = 1 : block_sizes(X, i)] for i = 0 : r]
    λ = [Rational{T}[2 - l[i+1][j+1] // l[i+1][j] - l[i+1][j] // l[i+1][j+1] for j = 1 : block_sizes(X, i) - 1] for i = 0 : r]
    res = sum([λ[i+1][j] // multiplicity(X, HyperbolicFixedPoint(i, j))
               for i = 0 : r for j = 1 : block_sizes(X,i)-1])

    if has_elliptic_fixed_point_plus(X)
        res += l_plus(X)^2 // sum_of_maximal_slopes(X)
    else
        res += 2 * l_plus(X) - sum_of_maximal_slopes(X)
    end

    if has_elliptic_fixed_point_minus(X)
        res -= l_minus(X)^2 // sum_of_minimal_slopes(X)
    else
        res += 2 * l_minus(X) + sum_of_minimal_slopes(X)
    end

    return res

end


@doc raw"""
    grading_matrix(X :: CStarSurface)

Return a tuple ``(Q_0,Q_1)`` where ``Q_0`` is the free part and ``Q_1`` is the
torsion part of the grading matrix associated to ``X``.
See Definition ``\ref{def:defining_triple_class_group}``.

"""
function grading_matrix(X :: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R}
    P = generator_matrix(X)
    S, U, _ = snf_with_transform(transpose(P))
    i0 = findfirst(i -> S[i,i] > 1, 1 : R) |> x -> x !== nothing ? x : R + 1
    
    n = size(P,2)
    Q_free = SMatrix{n - R, n, T}(U[R+1:end,:])
    Q_torsion = SMatrix{R-i0+1, n, T}(vcat([mod.(U[[i],:], S[i,i]) for i = i0 : R]...))
    return Q_free, Q_torsion
end


@doc raw"""
    grading_matrix_free_part(X :: CStarSurface)

Return the free part of the grading matrix associated to ``X``.
See Definition ``\ref{def:defining_triple_class_group}``.


# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> grading_matrix_free_part(X)
2×5 SMatrix{2, 5, Int64, 10} with indices SOneTo(2)×SOneTo(5):
  5  1   2   2  0
 -3  0  -1  -1  1
```

"""
grading_matrix_free_part(X :: CStarSurface) = grading_matrix(X)[1]


@doc raw"""
    grading_matrix_torsion_part(X :: CStarSurface)

Return the torsion part of the grading matrix associated to ``X``.
See Definition ``\ref{def:defining_triple_class_group}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[3],[3]], [[-1,-3],[2],[2]])
C*-surface of case PE with l = ((1,1),3,3) and d = ((-1,-3),2,2) 

julia> grading_matrix_torsion_part(X)
1×5 SMatrix{1, 5, Int64, 5} with indices SOneTo(1)×SOneTo(5):
 0  0  1  2  0
```

"""
grading_matrix_torsion_part(X :: CStarSurface) = grading_matrix(X)[2]
