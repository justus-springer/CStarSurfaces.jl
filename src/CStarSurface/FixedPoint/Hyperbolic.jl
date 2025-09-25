@doc raw"""
    HyperbolicFixedPoint <: BolicFixedPoint

A struct describing a hyperbolic fixed point ``x_{ij}`` of a ``\mathbb{C}^*``-surface.
It has two fields `i :: Int` and `j :: Int`, which are the indices of the block and
the ray inside the block that this fixed point is associated to.

"""
struct HyperbolicFixedPoint <: BolicFixedPoint
    i :: Int
    j :: Int
end

Base.show(io :: IO, x :: HyperbolicFixedPoint) =
print(io, "Hyperbolix fixed point x($(x.i), $(x.j))")

toric_chart(X :: CStarSurface, x :: HyperbolicFixedPoint) =
hcat(ray(X, x.i, x.j), ray(X, x.i, x.j+1))

@doc raw"""
    hyperbolic_fixed_points(X :: CStarSurface)

Return all hyperbolic fixed points of a ``\mathbb{C}^*``-surface,
see Definition ``\ref{def:defining_triple_fixed_points}``.

"""
function hyperbolic_fixed_points(X :: CStarSurface)
    r = number_of_blocks(X) - 1
    ns = block_sizes(X)

    xs = HyperbolicFixedPoint[]
    for i = 0 : number_of_blocks(X) - 1
        append!(xs, [HyperbolicFixedPoint(i,j) for j = 1 : ns[i+1] - 1])
    end

    return xs
end
