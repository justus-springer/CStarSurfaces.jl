@doc raw"""
    EllipticFixedPoint <: FixedPoint

Abstract supertype of elliptic fixed points of a ``\mathbb{C}^*``-surface.
It has two subtypes [`EllipticFixedPointPlus`](@ref) and [`EllipticFixedPointMinus`](@ref).

"""
abstract type EllipticFixedPoint <: FixedPoint end


@doc raw"""
    EllipticFixedPointPlus <: EllipticFixedPoint

The elliptic fixed point ``x^+``. This struct has no fields.

"""
struct EllipticFixedPointPlus <: EllipticFixedPoint end


@doc raw"""
    EllipticFixedPointMinus <: EllipticFixedPoint

The elliptic fixed point ``x^-``. This struct has no fields.

"""
struct EllipticFixedPointMinus <: EllipticFixedPoint end


@doc raw"""
    elliptic_fixed_points(X :: CStarSurface)

Return all elliptic fixed points of ``X``, see Definition ``\ref{def:defining_triple_fixed_points}``.

"""
function elliptic_fixed_points(X :: CStarSurface)
    xs = EllipticFixedPoint[]
    has_elliptic_fixed_point_plus(X) && push!(xs, EllipticFixedPointPlus())
    has_elliptic_fixed_point_minus(X) && push!(xs, EllipticFixedPointMinus())
    return xs
end

Base.show(io :: IO, ::MIME"text/plain", :: EllipticFixedPointPlus) =
print(io, "x⁺")

Base.show(io :: IO, :: EllipticFixedPointMinus) =
print(io, "x⁻")

toric_chart(X :: CStarSurface, :: EllipticFixedPointPlus) =
hcat([top_embedded_ray(X,i) for i = 0 : number_of_blocks(X)-1]...)

toric_chart(X :: CStarSurface, :: EllipticFixedPointMinus) =
hcat([bottom_embedded_ray(X,i) for i = 0 : number_of_blocks(X)-1]...)

function gorenstein_index(X :: CStarSurface{T,C,N,M,R}, :: EllipticFixedPointPlus) where {C, T <: Integer, N, M, R}
    ls = SVector{R}([top_ray(X, i-1)[1] for i = 1 : R])
    ds = SVector{R}([top_ray(X, i-1)[2] for i = 1 : R])

    ls_prods = SVector{R}([prod([ls[k] for k = 1 : R if k ≠ i]) for i = 1 : R])
    us = SVector{R-1}([(R-2) * ds[i] * ls_prods[i] +
        sum([(ds[j] - ds[i]) * 
            prod([ls[k] for k = 1 : R if k ≠ i && k ≠ j]) 
        for j = 1 : R if j ≠ i])
    for i = 2 : R])
    u_last = sum(ls_prods) - (R-2) * prod(ls)

    num = sum([ds[i] * ls_prods[i] for i = 1 : R])

    return num ÷ gcd(u_last, gcd(us))
end

function gorenstein_index(X :: CStarSurface{T,C,N,M,R}, :: EllipticFixedPointMinus) where {C, T <: Integer, N, M, R}
    ls = SVector{R}([bottom_ray(X, i-1)[1] for i = 1 : R])
    ds = SVector{R}([bottom_ray(X, i-1)[2] for i = 1 : R])

    ls_prods = SVector{R}([prod([ls[k] for k = 1 : R if k ≠ i]) for i = 1 : R])
    us = SVector{R-1}([(R-2) * ds[i] * ls_prods[i] +
        sum([(ds[j] - ds[i]) * 
            prod([ls[k] for k = 1 : R if k ≠ i && k ≠ j]) 
        for j = 1 : R if j ≠ i])
    for i = 2 : R])
    u_last = sum(ls_prods) - (R-2) * prod(ls)

    num = sum([ds[i] * ls_prods[i] for i = 1 : R])

    return -num ÷ gcd(u_last, gcd(us))
end

function log_canonicities(X :: CStarSurface{T,C} , :: EllipticFixedPointPlus) where {C, T <: Integer}
    ds = Rational{T}[]
    d_plus = sum_of_maximal_slopes(X) // l_plus(X)
    for i = 0 : number_of_blocks(X) - 1
        append!(ds, discrepancies(LatticePoint{T}(0,1), top_ray(X, i), d_plus)[1:end-1])
    end
    return ds
end

function log_canonicities(X :: CStarSurface{T,C} , :: EllipticFixedPointMinus) where {C, T <: Integer}
    ds = Rational{T}[]
    d_minus = -sum_of_minimal_slopes(X) // l_minus(X)
    for i = 0 : number_of_blocks(X) - 1
        append!(ds, discrepancies(LatticePoint{T}(0,-1), bottom_ray(X, i), d_minus)[1:end-1])
    end
    return ds
end

function is_quasismooth(X :: CStarSurface{T,C,N,M,R}, :: EllipticFixedPointPlus) where {C,T<:Integer,N,M,R}
    ls = SVector{R}([top_ray(X, i-1)[1] for i = 1 : R])
    return length(filter(l -> l > 1, ls)) <= 2
end

function is_quasismooth(X :: CStarSurface{T,C,N,M,R}, :: EllipticFixedPointMinus) where {C,T<:Integer,N,M,R}
    ls = SVector{R}([bottom_ray(X, i-1)[1] for i = 1 : R])
    return length(filter(l -> l > 1, ls)) <= 2
end
