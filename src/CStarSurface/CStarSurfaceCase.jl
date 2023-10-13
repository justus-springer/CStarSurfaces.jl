
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

has_x_plus(::Type{EE}) = true
has_x_plus(::Type{PE}) = false
has_x_plus(::Type{EP}) = true
has_x_plus(::Type{PP}) = false

has_x_minus(c::Type{T}) where {T <: CStarSurfaceCase} = has_x_plus(invert_case(c))

has_D_plus(c::Type{T}) where {T <: CStarSurfaceCase} = !has_x_plus(c)

has_D_minus(c::Type{T}) where {T <: CStarSurfaceCase} = !has_x_minus(c)

function _case_sym_to_type(c :: Symbol)
    c == :ee && return EE
    c == :pe && return PE
    c == :ep && return EP
    c == :pp && return PP
    throw(DomainError(c, "symbol must be one of :ee, :pe, :ep and :pp"))
end

