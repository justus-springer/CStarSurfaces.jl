
@doc raw"""
    BolicFixedPoint <: FixedPoint

Abstract supertype of [`ParabolicFixedPoint`](@ref) and [`HyperbolicFixedPoint`](@ref).

"""
abstract type BolicFixedPoint <: FixedPoint end

@doc raw"""
    bolic_fixed_points(X :: CStarSurface)

Return all bolic (hyperbolic and parabolic) fixed point of ``X``,
see Definition ``\ref{def:defining_triple_fixed_points}``.

"""
function bolic_fixed_points(X :: CStarSurface)
    xs = BolicFixedPoint[]
    append!(xs, parabolic_fixed_points(X))
    append!(xs, hyperbolic_fixed_points(X))
    return xs
end

log_canonicities(X :: CStarSurface, x :: BolicFixedPoint) =
discrepancies(toric_chart(X, x))[2:end-1]

function gorenstein_index(X :: CStarSurface, x :: BolicFixedPoint)
    U = toric_chart(X, x)
    return abs(det(U)) ÷ gcd(U[1,1] - U[1,2], U[2,1] - U[2,2])
end

is_quasismooth(X :: CStarSurface, x :: BolicFixedPoint) = true

