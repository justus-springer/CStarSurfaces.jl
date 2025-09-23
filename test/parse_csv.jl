# some auxiliary methods to parse the data files.
# These methods tell CSV.read how to parse the columns of the tsv file correctly
Base.tryparse(::Type{Matrix{T}}, s::String) where {T} = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{Vector{T}}, s::String) where {T} = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{Rational{Int}}, s::String) = Meta.eval(Meta.parse(s))
function Base.tryparse(::Type{CStarSurfaceCase}, s::String)
    if s == "EE"
        return EE
    elseif s == "PE"
        return PE
    elseif s == "EP"
        return EP
    elseif s == "PP"
        return PP
    else
        error("must be one of EE, PE, EP or PP")
    end
end

const TSV_COLUMN_TYPES = Dict([
    :V => Matrix{Int},
    :ns => Vector{Int},
    :case => CStarSurfaceCase,
    :class_group_rank => Int,
    :class_group_torsion => Vector{Int},
    :orientation => Int,
    :gorenstein_index => Int,
    :picard_index => Int,
    :anticanonical_self_intersection => Rational{Int},
    :maximal_log_canonicity => Rational{Int}
])

