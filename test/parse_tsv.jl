# some auxiliary methods to parse the data files 
# "ToricSurface/examples.tsv" and "CStarSurface/examples.tsv"

# These methods tell CSV.read how to parse the columns of the tsv file correctly
Base.tryparse(::Type{Vector{T}}, s::String) where {T} = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{Matrix}, s::String) = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{ZZMatrix}, s::String) = matrix(ZZ, Meta.eval(Meta.parse(s)))
Base.tryparse(::Type{QQMatrix}, s::String) = matrix(QQ, Meta.eval(Meta.parse(s)))
Base.tryparse(::Type{Rational{Int}}, s::String) = Meta.eval(Meta.parse(s))

const TSV_COLUMN_TYPES = Dict([
    :P => ZZMatrix, 
    :case => Symbol,
    :classGroupRank => Int,
    :classGroupTorsion => Vector{Int},
    :orientation => Int,
    :gorensteinIndex => Int,
    :picardIndex => Int,
    :intersectionMatrix => QQMatrix,
    :anticanonicalSelfIntersection => Rational{Int},
    :maximalLogCanonicity => Rational{Int}
])


