
@doc raw"""
    CStarSurfaceFixedPoint{T <: CStarSurfaceCase} <: CStarSurfacePoint{T}

A fixed point on a $\mathbb{C}^*$-surface.

"""
abstract type CStarSurfaceFixedPoint{T} <: CStarSurfacePoint{T} end


#################################################
# Elliptic fixed points
#################################################


@doc raw"""
    EllipticFixedPoint{T <: CStarSurfaceCase} <: CStarSurfaceFixedPoint{T}

An elliptic fixed point on a $\mathbb{C}^*$-surface.

"""
abstract type EllipticFixedPoint{T} <: CStarSurfaceFixedPoint{T} end


@doc raw"""
    EllipticFixedPointPlus{T <: Union{EE,EP}} <: EllipticFixedPoint{T}

An elliptic fixed point $x^+$ on a $\mathbb{C}$-surface of type (e-e) or (e-p).
There should only ever be one instance of this type for any given
`CStarSurface`, which is accessible via [`x_plus`](@ref).

"""
@attributes mutable struct EllipticFixedPointPlus{T <: Union{EE,EP}} <: EllipticFixedPoint{T}
    parent :: CStarSurface{T}
    
    EllipticFixedPointPlus(X :: CStarSurface{T}) where {T <: Union{EE,EP}} = new{T}(X)
end

Base.parent(x :: EllipticFixedPointPlus) = x.parent

Base.show(io :: IO, x :: EllipticFixedPointPlus) = print(io, "elliptic fixed point x^+")

@attr orbit_cone(x :: EllipticFixedPointPlus) = _sigma_plus(parent(x))


@doc raw"""
    x_plus(X :: CStarSurface{<:Union{EE,EP}})

Return the elliptic fixed point $x^+$ of a $\mathbb{C}^*$-surface of type
(e-e) or (e-p).

"""
@attr x_plus(X :: CStarSurface{<:Union{EE,EP}}) = EllipticFixedPointPlus(X)

@attr is_quasismooth(x :: EllipticFixedPointPlus) =
_almost_all_one(first.(_slope_ordered_l(parent(x))))


@attr function canonical_resolution(x :: EllipticFixedPointPlus)
    X = parent(x)
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_l(X), _slope_ordered_d(X)

    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    discrepancies = Rational{Int}[]
    for i = 0 : r
        v = [first(l[i]), first(d[i])]
        new_rays, new_discrepancies = toric_affine_surface_resolution(v, [0, _d_plus(X)])
        append!(new_l[i], map(first, new_rays))
        append!(new_d[i], map(last, new_rays))
        append!(discrepancies, new_discrepancies)
    end

    new_case = X.case == :ee ? :pe : :pp
    Y = cstar_surface(new_l, new_d, new_case)

    exceptional_divisors = CStarSurfaceDivisor[]
    for i = 0 : r
        append!(exceptional_divisors, [invariant_divisor(Y, i, j) for j = ns[i] + 1 : length(new_l[i])])
    end

    push!(exceptional_divisors, D_plus(Y))
    push!(discrepancies, _l_plus(X) // m_plus(X) - 1)

    return (Y, exceptional_divisors, discrepancies)

end


@doc raw"""
    EllipticFixedPointMinus{T <: Union{EE,PE}} <: EllipticFixedMinus{T}

An elliptic fixed point $x^-$ on a $\mathbb{C}$-surface of type (e-e) or (p-e).
There should only ever be one instance of this type for any given
`CStarSurface`, which is accessible via [`x_minus`](@ref).

"""
@attributes mutable struct EllipticFixedPointMinus{T <: Union{EE,PE}} <: EllipticFixedPoint{T}
    parent :: CStarSurface{T}

    EllipticFixedPointMinus(X :: CStarSurface{T}) where {T <: Union{EE,PE}} = new{T}(X)
end

Base.parent(x :: EllipticFixedPointMinus) = x.parent

Base.show(io :: IO, x :: EllipticFixedPointMinus) = print(io, "elliptic fixed point x^-")

@attr orbit_cone(x :: EllipticFixedPointMinus) = _sigma_minus(parent(x))


@doc raw"""
    x_minus(X :: CStarSurface{<:Union{EE,PE}})

Return the elliptic fixed point $x^-$ of a $\mathbb{C}^*$-surface of type
(e-e) or (p-e).

"""
@attr x_minus(X :: CStarSurface{<:Union{EE,PE}}) = EllipticFixedPointMinus(X)


@attr is_quasismooth(x :: EllipticFixedPointMinus) =
_almost_all_one(last.(_slope_ordered_l(parent(x))))


@attr function canonical_resolution(x :: EllipticFixedPointMinus)
    X = parent(x)
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_l(X), _slope_ordered_d(X)

    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    discrepancies = Rational{Int}[]
    for i = 0 : r
        v = [last(l[i]), last(d[i])]
        new_rays, new_discrepancies = toric_affine_surface_resolution(v, [0, -_d_minus(X)])
        append!(new_l[i], map(first, new_rays))
        append!(new_d[i], map(last, new_rays))
        append!(discrepancies, new_discrepancies)
    end

    new_case = X.case == :ee ? :ep : :pp
    Y = cstar_surface(new_l, new_d, new_case)

    exceptional_divisors = CStarSurfaceDivisor[]
    for i = 0 : r
        append!(exceptional_divisors, [invariant_divisor(Y, i, j) for j = ns[i] + 1 : length(new_l[i])])
    end

    push!(exceptional_divisors, D_minus(Y))
    push!(discrepancies, _l_minus(X) // m_minus(X) - 1)

    return (Y, exceptional_divisors, discrepancies)

end


@attr function minimal_resolution(x :: EllipticFixedPoint)
    (Y, divs, discrepancies) = deepcopy(canonical_resolution(x))

    contractible_curves = [(k, divs[k]) 
        for k = 1 : length(divs)
        if divs[k] * divs[k] == -1]

    while !isempty(contractible_curves)
        k, d = first(contractible_curves)
        Y = contract_prime_divisor(d)
        deleteat!(divs, k)
        deleteat!(discrepancies, k)

        # adjust the coefficients of the remaining exceptional divisors
        d_case, inds = is_prime_with_double_indices(d)
        m = number_of_parabolic_fixed_point_curves(Y)
        for l = 1 : length(divs)
            if d_case == :D_ij
                i, j = inds
                new_coeffs = deleteat!(double_coefficients(divs[l]), i, j)
                last_coeffs = coefficients(divs[l])[end-m+1 : end]
                new_d = cstar_surface_divisor(Y, new_coeffs, last_coeffs...)
            else
                k = inds
                new_coeffs = deleteat!(coefficients(divs[l]), k)
                new_d = cstar_surface_divisor(Y, new_coeffs)
            end
            divs[l] = new_d
        end

        contractible_curves = [(k, divs[k]) 
            for k = 1 : length(divs)
            if divs[k] * divs[k] == -1]

    end

    return (Y, divs, discrepancies)

end


@doc raw"""
    elliptic_fixed_points(X :: CStarSurface)

Return the elliptic fixed points of a $\mathbb{C}^*$-surface.

"""
@attr elliptic_fixed_points(X :: CStarSurface{EE}) = [x_plus(X), x_minus(X)]
@attr elliptic_fixed_points(X :: CStarSurface{PE}) = [x_minus(X)]
@attr elliptic_fixed_points(X :: CStarSurface{EP}) = [x_plus(X)]
@attr elliptic_fixed_points(X :: CStarSurface{PP}) = EllipticFixedPoint{PP}[]


#################################################
# Hyperbolic fixed points
#################################################


@doc raw"""
    HyperbolicFixedPoint{T <: CStarSurfaceCase} <: CStarSurfaceFixedPoint{T}

A hyperbolic fixed point on a $\mathbb{C}$^*-surface.

"""
@attributes mutable struct HyperbolicFixedPoint{T} <: CStarSurfaceFixedPoint{T}
    parent :: CStarSurface{T}
    i :: Int
    j :: Int
    
    function HyperbolicFixedPoint(X :: CStarSurface{T}, i :: Int, j :: Int) where {T <: CStarSurfaceCase}
        r, ns = nblocks(X) - 1, block_sizes(X)
        @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
        @req 1 ≤ j ≤ ns[i] - 1 "must have 1 ≤ j ≤ ns[i] - 1"
        return new{T}(X, i, j)
    end

end

Base.parent(x :: HyperbolicFixedPoint) = x.parent

Base.show(io :: IO, x :: HyperbolicFixedPoint) = print(io, "hyperbolic fixed point x($(x.i), $(x.j))")

@attr orbit_cone(x :: HyperbolicFixedPoint) = _tau(parent(x), x.i, x.j)


@doc raw"""
    hyperbolic_fixed_points(X :: CStarSurface)   

Return the hyperbolic fixed points of a $\mathbb{C}^*$-surface as a
`DoubleVector`.

"""
@attr hyperbolic_fixed_points(X :: CStarSurface) = 
DoubleVector([[HyperbolicFixedPoint(X, i, j) for j = 1 : _ns(X)[i] - 1] for i = 0 : _r(X)])


@doc raw"""
    hyperbolic_fixed_point(X :: CStarSurface, i :: Int, j :: Int)

Return the hyperbolic fixed point $x_{ij}$ of a $\mathbb{C}^*$-surface, where
$0 ≤ i ≤ r$ and $1 ≤ j ≤ n_i - 1$.

"""
function hyperbolic_fixed_point(X :: CStarSurface, i :: Int, j :: Int)
    r, ns = nblocks(X) - 1, block_sizes(X)
    @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
    @req 1 ≤ j ≤ ns[i] - 1 "must have 1 ≤ j ≤ ns[i] - 1"
    return hyperbolic_fixed_points(X)[i][j]
end


@attr is_quasismooth(x :: HyperbolicFixedPoint) = true


@attr function canonical_resolution(x :: HyperbolicFixedPoint)
    X, i, j = parent(x), x.i, x.j
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_l(X), _slope_ordered_d(X)
    v1, v2 = [l[i][j], d[i][j]], [l[i][j+1], d[i][j+1]]
    
    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    new_rays, discrepancies = toric_affine_surface_resolution(v1, v2)
    append!(new_l[i], map(first, new_rays))
    append!(new_d[i], map(last, new_rays))

    Y = cstar_surface(new_l, new_d, X.case)

    exceptional_divisors = [invariant_divisor(Y, i, j) for j = ns[i] + 1 : length(new_l[i])]

    return (Y, exceptional_divisors, discrepancies)

end


@attr minimal_resolution(x :: HyperbolicFixedPoint) = canonical_resolution(x)


#################################################
# Parabolic fixed points
#################################################


@doc raw"""
    ParabolicFixedPoint{T <: CStarSurfaceCase} <: CStarSurfaceFixedPoint{T}

A parabolic fixed point on a $\mathbb{C}^*$-surface.

"""
abstract type ParabolicFixedPoint{T} <: CStarSurfaceFixedPoint{T} end

@attr is_quasismooth(x :: ParabolicFixedPoint) = true

@doc raw"""
    ParabolicFixedPointPlus{T<:Union{PE,PP}} <: ParabolicFixedPoint{T}

A parabolic fixed point $x^+_i$ on a $\mathbb{C}^*$-surface of type (p-e) or
(p-p)

"""
@attributes mutable struct ParabolicFixedPointPlus{T<:Union{PE,PP}} <: ParabolicFixedPoint{T}
    parent :: CStarSurface{T}
    i :: Int

    function ParabolicFixedPointPlus(X :: CStarSurface{T}, i :: Int) where {T <: Union{PE,PP}}
        r = nblocks(X) - 1
        @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
        return new{T}(X, i)
    end
end

Base.parent(x :: ParabolicFixedPointPlus) = x.parent

Base.show(io :: IO, x :: ParabolicFixedPointPlus) = print(io, "parabolic fixed point x^+($(x.i))")

@attr orbit_cone(x :: ParabolicFixedPointPlus) = _tau_plus(parent(x), x.i)


@doc raw"""
    parabolic_fixed_points_plus(X :: CStarSurface{T}) where {T <: Union{PE,PP}}

Return the parabolic fixed points $x_i^+$ of a $\mathbb{C}^*$-surface.

"""
@attr parabolic_fixed_points_plus(X :: CStarSurface{T}) where {T <: Union{PE,PP}} =
ZeroVector([ParabolicFixedPointPlus(X, i) for i = 0 : nblocks(X) - 1])


@doc raw"""
    parabolic_fixed_point_plus(X :: CStarSurface{<:Union{PE,PP}}, i :: Int)

Return the parabolic fixed point $x_i^+$ of a $\mathbb{C}^*$-surface, where $0
≤ i ≤ r$.


"""
function parabolic_fixed_point_plus(X :: CStarSurface{<:Union{PE,PP}}, i :: Int)
    r = nblocks(X) - 1
    @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
    return parabolic_fixed_points_plus(X)[i]
end


@attr function canonical_resolution(x :: ParabolicFixedPointPlus)
    X, i = parent(x), x.i
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_l(X), _slope_ordered_d(X)

    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    v = [first(l[i]), first(d[i])]
    new_rays, discrepancies = toric_affine_surface_resolution(v, [0, 1])
    append!(new_l[i], map(first, new_rays))
    append!(new_d[i], map(last, new_rays))

    Y = cstar_surface(new_l, new_d, X.case)

    exceptional_divisors = [invariant_divisor(Y, i, j) for j = ns[i] + 1 : length(new_l[i])]

    return (Y, exceptional_divisors, discrepancies)

end

@doc raw"""
    ParabolicFixedPointMinus{T<:Union{EP,PP}} <: ParabolicFixedPoint{T}

A parabolic fixed point $x^-_i$ on a $\mathbb{C}^*$-surface of type (e-p) or
(p-p).

"""
@attributes mutable struct ParabolicFixedPointMinus{T<:Union{EP,PP}} <: ParabolicFixedPoint{T}
    parent :: CStarSurface{T}
    i :: Int

    function ParabolicFixedPointMinus(X :: CStarSurface{T}, i :: Int) where {T <: Union{EP,PP}}
        r = nblocks(X) - 1
        @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
        return new{T}(X, i)
    end
end

Base.parent(x :: ParabolicFixedPointMinus) = x.parent

Base.show(io :: IO, x :: ParabolicFixedPointMinus) = print(io, "parabolic fixed point x^-($(x.i))")

@attr orbit_cone(x :: ParabolicFixedPointMinus) = _tau_minus(parent(x), x.i)


@doc raw"""
    parabolic_fixed_points_minus(X :: CStarSurface{T}) where {T <: Union{EP,PP}}

Return the parabolic fixed points $x_i^-$ of a $\mathbb{C}^*$-surface.

"""
@attr parabolic_fixed_points_minus(X :: CStarSurface{T}) where {T <: Union{EP,PP}} =
ZeroVector([ParabolicFixedPointMinus(X, i) for i = 0 : nblocks(X) - 1])


@doc raw"""
    parabolic_fixed_point_minus(X :: CStarSurface{<:Union{EP,PP}}, i :: Int)

Return the parabolic fixed point $x_i^-$ of a $\mathbb{C}^*$-surface, where $0
≤ i ≤ r$.

"""
function parabolic_fixed_point_minus(X :: CStarSurface{<:Union{EP,PP}}, i :: Int)
    r = nblocks(X) - 1
    @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
    return parabolic_fixed_points_minus(X)[i]
end


@attr function canonical_resolution(x :: ParabolicFixedPointMinus)
    X, i = parent(x), x.i
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_l(X), _slope_ordered_d(X)

    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    v = [last(l[i]), last(d[i])]
    new_rays, discrepancies = toric_affine_surface_resolution(v, [0, -1])
    append!(new_l[i], map(first, new_rays))
    append!(new_d[i], map(last, new_rays))

    Y = cstar_surface(new_l, new_d, X.case)

    exceptional_divisors = [invariant_divisor(Y, i, j) for j = ns[i] + 1 : length(new_l[i])]

    return (Y, exceptional_divisors, discrepancies)

end

@attr minimal_resolution(x :: ParabolicFixedPoint) = canonical_resolution(x)

@doc raw"""
    parabolic_fixed_points(X :: CStarSurface)

Return the parabolic fixed points of a $\mathbb{C}^*$-surface.

"""
@attr parabolic_fixed_points(X :: CStarSurface{EE}) = ParabolicFixedPoint{EE}[] 
@attr parabolic_fixed_points(X :: CStarSurface{PE}) = Vector(parabolic_fixed_points_plus(X))
@attr parabolic_fixed_points(X :: CStarSurface{EP}) = Vector(parabolic_fixed_points_minus(X))
@attr parabolic_fixed_points(X :: CStarSurface{PP}) = [Vector(parabolic_fixed_points_plus(X)) ; Vector(parabolic_fixed_points_minus(X))]


@doc raw"""
    fixed_points(X :: CStarSurface)

Return all fixed points of a $\mathbb{C}^*$-action. This is the union of
[`elliptic_fixed_points`](@ref), [`hyperbolic_fixed_points`](@ref) and
[`parabolic_fixed_points`](@ref).

"""
@attr fixed_points(X :: CStarSurface{T}) where {T <: CStarSurfaceCase} = 
CStarSurfaceFixedPoint{T}[
    elliptic_fixed_points(X) ; 
    vcat(hyperbolic_fixed_points(X)...) ; 
    parabolic_fixed_points(X)
]
