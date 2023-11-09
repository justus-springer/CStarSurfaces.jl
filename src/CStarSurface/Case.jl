
@doc raw"""
    invert_case(c :: Union{Symbol, Type{<:CStarSurfaceCase}} [, invert :: Bool])   

Flips source and sink for a given C-star surface case. It can either take one
of the four symbols `:ee, :pe, :ep, :pp` or one of the four types `EE, PE, EP, PP`.

If the optional parameter `invert` is set to `false`, this function does nothing.

"""
invert_case(::Type{EE}) = EE
invert_case(::Type{PE}) = EP
invert_case(::Type{EP}) = PE
invert_case(::Type{PP}) = PP
function invert_case(c :: Symbol)
    c == :ee && return :ee
    c == :pe && return :ep
    c == :ep && return :pe
    c == :pp && return :pp
    throw(DomainError(c, "symbol must be one of :ee, :pe, :ep and :pp"))
end
invert_case(c::Union{Symbol, Type{<:CStarSurfaceCase}}, invert :: Bool) = invert ? invert_case(c) : c


@doc raw"""
    has_x_plus(::Type{<:CStarSurfaceCase})

Checks whether a given `CStarSurfaceCase` has an elliptic fixed point in the
source, commonly denoted $x^+$.

"""
has_x_plus(::Type{EE}) = true
has_x_plus(::Type{PE}) = false
has_x_plus(::Type{EP}) = true
has_x_plus(::Type{PP}) = false


@doc raw"""
    has_x_minus(::Type{<:CStarSurfaceCase})

Checks whether a given `CStarSurfaceCase` has an elliptic fixed point in the
sink, commonly denoted $x^-$.

"""
has_x_minus(c::Type{T}) where {T <: CStarSurfaceCase} = has_x_plus(invert_case(c))


@doc raw"""
    has_D_plus(::Type{<:CStarSurfaaceCase})

Checks whether a given `CStarSurfaceCase` has a parabolic fixed point curve
source, commonly denoted $D^+$.

"""
has_D_plus(c::Type{T}) where {T <: CStarSurfaceCase} = !has_x_plus(c)


@doc raw"""
    has_D_minus(::Type{<:CStarSurfaaceCase})

Checks whether a given `CStarSurfaceCase` has a parabolic fixed point curve
sink, commonly denoted $D^-$.

"""
has_D_minus(c::Type{T}) where {T <: CStarSurfaceCase} = !has_x_minus(c)

function _case_sym_to_type(c :: Symbol)
    c == :ee && return EE
    c == :pe && return PE 
    c == :ep && return EP
    c == :pp && return PP
    throw(DomainError(c, "symbol must be one of :ee, :pe, :ep and :pp"))
end

_case_type_to_sym(::Type{EE}) = :ee
_case_type_to_sym(::Type{PE}) = :pe
_case_type_to_sym(::Type{EP}) = :ep
_case_type_to_sym(::Type{PP}) = :pp


Base.:+(::Type{EE}, ::Type{EE}) = EE
Base.:+(::Type{EE}, ::Type{PE}) = PE
Base.:+(::Type{EE}, ::Type{EP}) = EP
Base.:+(::Type{EE}, ::Type{PP}) = PP
Base.:+(::Type{PE}, ::Type{EE}) = PE
Base.:+(::Type{PE}, ::Type{PE}) = PE
Base.:+(::Type{PE}, ::Type{EP}) = PP
Base.:+(::Type{PE}, ::Type{PP}) = PP
Base.:+(::Type{EP}, ::Type{EE}) = EP
Base.:+(::Type{EP}, ::Type{PE}) = PP
Base.:+(::Type{EP}, ::Type{EP}) = EP
Base.:+(::Type{EP}, ::Type{PP}) = PP
Base.:+(::Type{PP}, ::Type{EE}) = PP
Base.:+(::Type{PP}, ::Type{PE}) = PP
Base.:+(::Type{PP}, ::Type{EP}) = PP
Base.:+(::Type{PP}, ::Type{PP}) = PP
