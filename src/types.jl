
######################################################################
# Types of varieties
######################################################################

@doc raw"""
    MoriDreamSpace

Julia type for Mori Dream Spaces.

Subtypes of `MoriDreamSpace` should at least implement the following
methods: [`canonical_toric_ambient`](@ref), [`cox_ring_relations`](@ref),
[`is_quasismooth`](@ref).

"""
abstract type MoriDreamSpace end


@doc raw"""
    MoriDreamSpaceUnion

The Union of `MoriDreamSpace` and Oscar's `NormalToricVarietyType`.

"""
const MoriDreamSpaceUnion = Union{MoriDreamSpace, NormalToricVarietyType}


@doc raw"""
    ToricVarietyMDS

Julia type for toric Mori Dream Spaces.
"""
abstract type ToricVarietyMDS <: MoriDreamSpace end


@doc raw"""
    CStarSurfaceCase

The abstract supertype of the possible possible configurations of source and
sink of a $\mathbb{C}^*$-surface.

This type has the four subtypes `EE`, `EP`, `PE` and `PP`, named after the
existence of elliptic fixed points and parabolic fixed point curves in the
source and sink of a $\mathbb{C}^*$-surface respectively.

"""
abstract type CStarSurfaceCase end

struct EE <: CStarSurfaceCase end
struct PE <: CStarSurfaceCase end
struct EP <: CStarSurfaceCase end
struct PP <: CStarSurfaceCase end


@doc raw"""
    ZeroVector{T}

A zero-indexed vector. 
"""
const ZeroVector{T, S} = OffsetVector{T, S}


@doc raw"""
    DoubleVector{T}

A zero-indexed vector of one-indexed vectors.

This type of indexing is very common when working with $\mathbb{C}^*$-
surfaces.

"""
const DoubleVector{T, S} = ZeroVector{Vector{T}, S}


@doc raw"""
    CStarSurface{T<:CStarSurfaceCase} <: MoriDreamSpace

A $\mathbb{C}^*$-surface of case `T <: CStarSurfaceCase`. As a Julia type, it
gets modeled by a stuct with fields `l`, `d` and `case`, where `l` and `d` are
zero-indexed vectors of one-indexed vectors of integers and `case` is one of
the four symbols `:ee`, `:pe`, `:ep` and `:pp`.

"""
@attributes mutable struct CStarSurface{T<:CStarSurfaceCase} <: MoriDreamSpace
    l :: DoubleVector{Int64}
    d :: DoubleVector{Int64}
    case :: Symbol
    CStarSurface{T}(l, d, case) where {T <: CStarSurfaceCase} = new{T}(l, d, case)
end

Base.:(==)(X :: CStarSurface, Y :: CStarSurface) = X.l == Y.l && X.d == Y.d && X.case == Y.case


@doc raw"""
    ToricSurface <: MoriDreamSpace

A toric surface. As a Julia type, it gets modeled by a struct with a single
field `vs` that stores the primitive generator of the rays of the
two-dimensional complete fan describing the toric surface.

"""
@attributes mutable struct ToricSurface <: ToricVarietyMDS
    vs :: Vector{Vector{T}} where {T <: IntegerUnion}
    ToricSurface(vs) = new(vs)
end

Base.:(==)(X :: ToricSurface, Y :: ToricSurface) = X.vs == Y.vs 

######################################################################
# Types of divisors
######################################################################

@doc raw"""
    MoriDreamSpaceDivisor{T <: MoriDreamSpace}

A Weil divisor on a Mori Dream Space of type `T`.

"""
@attributes mutable struct MoriDreamSpaceDivisor{T <: MoriDreamSpace}
    variety :: T
    coeffs :: Vector{<:IntegerUnion}
    MoriDreamSpaceDivisor(variety :: T, coeffs :: Vector{<:IntegerUnion}) where {T <: MoriDreamSpace} = new{T}(variety, coeffs)
end

Base.:(==)(d1 :: MoriDreamSpaceDivisor{T}, d2 :: MoriDreamSpaceDivisor{T}) where {T <: MoriDreamSpace} = 
d1.variety === d2.variety && d1.coeffs == d2.coeffs


@doc raw"""
    CStarSurfaceDivisor{T} = MoriDreamSpaceDivisor{CStarSurface{T}}

A Weil divisor on a $\mathbb{C}^*$-surface of type `T <: CStarSurfaceCase`.
"""
const CStarSurfaceDivisor{T} = MoriDreamSpaceDivisor{CStarSurface{T}}


@doc raw"""
    ToricSurfaceDivisor = MoriDreamSpaceDivisor{ToricSurface} 

A Weil divisor on a toric surface.
"""
const ToricSurfaceDivisor = MoriDreamSpaceDivisor{ToricSurface}


