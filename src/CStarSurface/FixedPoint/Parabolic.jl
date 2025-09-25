@doc raw"""
    ParabolicFixedPoint <: BolicFixedPoint

Abstract supertype of parabolic fixed points of a ``\mathbb{C}^*``-surface.
It has two subtypes [`ParabolicFixedPointPlus`](@ref) and [`ParabolicFixedPointMinus`](@ref).

"""
abstract type ParabolicFixedPoint <: BolicFixedPoint end


@doc raw"""
    ParabolicFixedPointPlus <: ParabolicFixedPoint

A parabolic fixed point ``x^+_i``. It has one field `i :: Int`, which
is the index of the arm of the fixed point.

"""
struct ParabolicFixedPointPlus <: ParabolicFixedPoint
    i :: Int
end


@doc raw"""
    ParabolicFixedPointMinus <: ParabolicFixedPoint

A parabolic fixed point ``x^-_i``. It has one field `i :: Int`, which
is the index of the arm of the fixed point.

"""
struct ParabolicFixedPointMinus <: ParabolicFixedPoint
    i :: Int
end


@doc raw"""
    parabolic_fixed_points(X :: CStarSurface)

Return all parabolic fixed points of a ``\mathbb{C}^*``-surface,
see Definition ``\ref{def:defining_triple_fixed_points}``.

"""
function parabolic_fixed_points(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    xs = ParabolicFixedPoint[]

    if has_parabolic_fixed_point_curve_plus(X)
        append!(xs, [ParabolicFixedPointPlus(i) for i = 0 : r])
    end

    if has_parabolic_fixed_point_curve_minus(X)
        append!(xs, [ParabolicFixedPointMinus(i) for i = 0 : r])
    end

    return xs
end

Base.show(io :: IO, x :: ParabolicFixedPointPlus) =
print(io, "Parabolic fixed point x^+($(x.i))")

Base.show(io :: IO, x :: ParabolicFixedPointMinus) =
print(io, "Parabolic fixed point x^-($(x.i))")

toric_chart(X :: CStarSurface{T,C}, x :: ParabolicFixedPointPlus) where {C, T <: Integer} =
hcat(top_ray(X, x.i), LatticePoint{T}(0,1))

toric_chart(X :: CStarSurface{T,C}, x :: ParabolicFixedPointMinus) where {C, T <: Integer} =
hcat(LatticePoint{T}(0,-1), bottom_ray(X, x.i))
