
######################################################################
# Types of varieties
######################################################################

@doc raw"""
    MoriDreamSpace

Julia type for Mori Dream Spaces.

Subtypes of `MoriDreamSpace` should at least implement the following
methods:

`canonical_toric_ambient(X :: MoriDreamSpace)`

`cox_ring_relations(X :: MoriDreamSpace)`
"""
abstract type MoriDreamSpace end


@doc raw"""
    MoriDreamSpaceUnion 

The union of `MoriDreamSpace` and Oscar's `NormalToricVarietyType`.

"""
const MoriDreamSpaceUnion = Union{MoriDreamSpace, Oscar.NormalToricVarietyType}


@doc raw"""
    ToricVarietyMDS

Julia type for toric Mori Dream Spaces.
"""
abstract type ToricVarietyMDS <: MoriDreamSpace end


@doc raw"""
    CStarSurfaceCase

The abstract supertype of the four possible cases `EE`, `PE`, `EP` and
`PP` of $\mathbb{C}^*$-surfaces.
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
struct ZeroVector{T} <: AbstractVector{T}
    parent :: Vector{T}
end


@doc raw"""
    DoubleVector{T}

A zero-indexed vector of one-indexed vectors.

This type of indexing is very common when working with $\mathbb{C}^*$-
surfaces.
"""
const DoubleVector{T} = ZeroVector{Vector{T}}


@doc raw"""
    CStarSurface{T<:CStarSurfaceCase}

A $\mathbb{C}^*$-surface of case `T <: CStarSurfaceCase`.
"""
@attributes mutable struct CStarSurface{T<:CStarSurfaceCase} <: MoriDreamSpace
    l :: DoubleVector{Int64}
    d :: DoubleVector{Int64}
    case :: Symbol
    CStarSurface{T}(l, d, case) where {T <: CStarSurfaceCase} = new{T}(l, d, case)
end


@doc raw"""
    ToricSurface

A toric surface.
"""
@attributes mutable struct ToricSurface <: ToricVarietyMDS
    vs :: Vector{Vector{T}} where {T <: Oscar.IntegerUnion}
    ToricSurface(vs) = new(vs)
end


######################################################################
# Types of divisors
######################################################################

@doc raw"""
    MoriDreamSpaceDivisor{T <: MoriDreamSpace}

A Weil divisor on a Mori Dream Space of type `T`.

Since there is a 1-to-1 correspondence between (classes of) divisors
on a Mori Dream Space and (classes of) divisors on its canonical 
toric ambient variety, this type is modeled as a wrapper around
the `ToricDivisor` from Oscar.
"""
mutable struct MoriDreamSpaceDivisor{T <: MoriDreamSpace}
    variety :: T
    toric_divisor :: ToricDivisor
    function MoriDreamSpaceDivisor(variety :: T, toric_divisor :: ToricDivisor) where {T <: MoriDreamSpace}
        @req canonical_toric_ambient(variety) == toric_divisor.toric_variety "the toric divisor must be defined on the canonical toric ambient variety of X"
        new{T}(variety, toric_divisor)
    end
end


@doc raw"""
    CStarSurfaceDivisor{T} 

A Weil divisor on a $\mathbb{C}^*$-surface of type `T <: CStarSurfaceCase`.
"""
const CStarSurfaceDivisor{T} = MoriDreamSpaceDivisor{CStarSurface{T}}


@doc raw"""
    ToricSurfaceDivisor 

A Weil divisor on a toric surface.
"""
const ToricSurfaceDivisor = MoriDreamSpaceDivisor{ToricSurface}
