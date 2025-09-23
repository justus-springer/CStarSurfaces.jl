@doc raw"""
    FixedPoint

Abstract supertype of a fixed point on a C*-surface. It has two subtypes
`EllipticFixedPoint` and `BolicFixedPoint`, which itself has subtypes
`HyperbolicFixedPoint` and `ParabolicFixedPoint`.

"""
abstract type FixedPoint end


@doc raw"""
    fixed_points(X :: CStarSurface)

Return all fixed point of a C*-surface.

"""
function fixed_points(X :: CStarSurface)
    xs = FixedPoint[]
    append!(xs, elliptic_fixed_points(X))
    append!(xs, parabolic_fixed_points(X))
    append!(xs, hyperbolic_fixed_points(X))
    return xs
end


@doc raw"""
    toric_chart(X :: CStarSurface, x :: FixedPoint)

Return the generator matrix of the toric chart around a fixed point. This will
be the submatrix of the generator matrix of `X`, where we only take those rays
that are contained in the orbit cone of `x`.

"""
function toric_chart end


@doc raw"""
    gorenstein_index(X :: CStarSurface, x :: FixedPoint)

Return the local Gorenstein index at the point `x`.

"""
function gorenstein_index end


@doc raw"""
    multiplicity(X :: CStarSurface, x :: FixedPoint)

Return the order of the local class group at the point `x`.

"""
multiplicity(X :: CStarSurface, x :: FixedPoint) =
abs(det_bareiss(toric_chart(X, x)))


@doc raw"""
    is_factorial(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is factorial, i.e. has trivial local class group.

"""
is_factorial(X :: CStarSurface, x :: FixedPoint) =
multiplicity(X, x) == 1


@doc raw"""
    is_quasismooth(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is quasismooth.

"""
function is_quasismooth end


@doc raw"""
    is_smooth(X :: CStarSurface, x :: FixedPoint)

Check whether a fixed point is smooth.

"""
is_smooth(X :: CStarSurface, x :: FixedPoint) =
is_factorial(X, x) && is_quasismooth(X, x)


@doc raw"""
    class_group(X :: CStarSurface, x :: FixedPoint)

Return the local class group at the point `x`.

"""
class_group(X :: CStarSurface, x :: FixedPoint) =
cokernel(toric_chart(X, x))


@doc raw"""
    log_canonicity(X :: CStarSurface, x :: FixedPoint)

Return the maximal `ε` such that the singularity at `x` is ε-log canonical. By
definition, this is set to be `1//0` (infinity) for smooth points.

"""
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

Check whether a point on a C*-surface is ε-log terminal. The default value
of ε is zero, which is the usual notion of log terminality.

"""
is_log_terminal(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0) =
log_canonicity(X, x) > ε


@doc raw"""
    is_log_canonical(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0)

Check whether a point on a C*-surface is ε-log canonical. The default value
of ε is zero, which is the usual notion of log canonical..

"""
is_log_canonical(X :: CStarSurface, x :: FixedPoint, ε :: Real = 0) =
log_canonicity(X, x) ≥ ε


@doc raw"""
    is_terminal(X :: CStarSurface, x :: FixedPoint)

Check whether a point on a C*-surface is terminal. This is equivalent to being
smooth.

"""
is_terminal(X :: CStarSurface, x :: FixedPoint) =
is_log_terminal(X, x, 1)


@doc raw"""
    is_canonical(X :: CStarSurface, x :: FixedPoint)

Check whether a point on a C*-surface is canonical.

"""
is_canonical(X :: CStarSurface, x :: FixedPoint) =
is_log_canonical(X, x, 1)
