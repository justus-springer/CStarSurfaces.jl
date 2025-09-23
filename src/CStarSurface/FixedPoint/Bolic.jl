
@doc raw"""
    BolicFixedPoint <: FixedPoint

Abstract supertype of `ParabolicFixedPoint` and `HyperbolicFixedPoint`.

"""
abstract type BolicFixedPoint <: FixedPoint end

log_canonicities(X :: CStarSurface, x :: BolicFixedPoint) =
discrepancies(toric_chart(X, x))[2:end-1]

function gorenstein_index(X :: CStarSurface, x :: BolicFixedPoint)
    U = toric_chart(X, x)
    return abs(det(U)) ÷ gcd(U[1,1] - U[1,2], U[2,1] - U[2,2])
end

is_quasismooth(X :: CStarSurface, x :: BolicFixedPoint) = true

