@doc raw"""
    class_group(X :: CStarSurface)

The divisor class group of a C*-surface.

"""
class_group(X :: CStarSurface) = cokernel(transpose(generator_matrix(X)))


@doc raw"""
    multiplicity(X :: CStarSurface)

The order of the torsion part of the class group of `X`.

"""
multiplicity(X :: CStarSurface) = torsion_order(class_group(X))


@doc raw"""
    picard_index(X :: CStarSurface)

The index of the Picard group inside the divisor class group of a C*-surface.

"""
picard_index(X :: CStarSurface) = prod([multiplicity(X,x) for x in fixed_points(X)]) ÷ multiplicity(X)


@doc raw"""
    gorenstein_index(X :: CStarSurface)

The Gorenstein index of a C*-surface.

"""
gorenstein_index(X :: CStarSurface) = lcm([gorenstein_index(X,x) for x in fixed_points(X)])


@doc raw"""
    is_quasismooth(X :: CStarSurface)   

Check whether the C*-surface is quasismooth, i.e. its characteristic space is
smooth.

"""
is_quasismooth(X :: CStarSurface) = all(x -> is_quasismooth(X, x), elliptic_fixed_points(X))


@doc raw"""
    is_factorial(X :: CStarSurface)   

Check whether the C*-surface is locally factorial, i.e. all local class groups
are trivial.

"""
is_factorial(X :: CStarSurface) = all(x -> is_factorial(X, x), fixed_points(X))


@doc raw"""
    is_smooth(X :: CStarSurface)

Check whether the C*-surface is smooth.

"""
is_smooth(X :: CStarSurface) = is_factorial(X) && is_quasismooth(X)


@doc raw"""
    log_canonicity(X :: CStarSurface)

Return the maximal `ε` such that the surface has at most ε-log canonical
singularities. For smooth surfaces, this returns `1//0` which is infinity.

"""
log_canonicity(X :: CStarSurface) = minimum([log_canonicity(X, x) for x in fixed_points(X)])


@doc raw"""
    is_log_terminal(X :: CStarSurface, ε :: Real = 0)   

Check whether a C*-surface is ε-log terminal. The default value of ε is zero,
which is the usual notion of log terminality.

"""
is_log_terminal(X :: CStarSurface, ε :: Real = 0) = log_canonicity(X) > ε


@doc raw"""
    is_log_canonical(X :: CStarSurface, ε :: Real = 0)

Check whether a C*-surface is ε-log canonical. The default value of ε is zero,
which is the usual notion of log canonical..

"""
is_log_canonical(X :: CStarSurface, ε :: Real = 0) = log_canonicity(X) ≥ ε


@doc raw"""
    is_terminal(X :: CStarSurface)

Check whether a C*-surface is terminal. This is equivalent to being smooth.

"""
is_terminal(X :: CStarSurface) = is_log_terminal(X, 1)


@doc raw"""
    is_canonical(X :: CStarSurface, x :: FixedPoint)

Check whether a C*-surface is canonical.

"""
is_canonical(X :: CStarSurface) = is_log_canonical(X, 1)


@doc raw"""
    degree(X :: CStarSurface)

The self-intersection number of an anticanonical divisor on the C*-surface.

"""
function degree(X :: CStarSurface{T,C}) where {C, T <: Integer}
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

function degree_matrix(X :: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R}
    P = generator_matrix(X)
    S, U, _ = snf_with_transform(transpose(P))
    i0 = findfirst(i -> S[i,i] > 1, 1 : R) |> x -> x !== nothing ? x : R + 1
    
    n = size(P,2)
    Q_free = SMatrix{n - R, n, T}(U[R+1:end,:])
    Q_torsion = SMatrix{R-i0+1, n, T}(vcat([mod.(U[[i],:], S[i,i]) for i = i0 : R]...))
    return Q_free, Q_torsion
end


degree_matrix_free_part(P :: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R} =
degree_matrix(P)[1]

degree_matrix_torsion_part(P :: CStarSurface{T,C,N,M,R}) where {C, T <: Integer, N, M, R} =
degree_matrix(P)[2]
