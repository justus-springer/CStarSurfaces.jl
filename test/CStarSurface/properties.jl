# using CSV, DataFrames, Test

# Loads in a list of 200 randomly chosen examples from the ldp-database
EXAMPLES_FILE = joinpath(dirname(@__FILE__), "examples.tsv")
# EXAMPLES_FILE = "test/CStarSurface/examples.tsv"

# These methods tell CSV.read how to parse the columns of the tsv file correctly
Base.tryparse(::Type{Vector{T}}, s::String) where {T} = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{Matrix}, s::String) = Meta.eval(Meta.parse(s))
Base.tryparse(::Type{ZZMatrix}, s::String) = matrix(ZZ, Meta.eval(Meta.parse(s)))
Base.tryparse(::Type{QQMatrix}, s::String) = matrix(QQ, Meta.eval(Meta.parse(s)))
Base.tryparse(::Type{Rational{Int}}, s::String) = Meta.eval(Meta.parse(s))

col_types = Dict([
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

df = CSV.read(EXAMPLES_FILE, DataFrame, delim='\t', types=col_types)

Xs = map(df_row -> cstar_surface(df_row.P), eachrow(df))

@testset "examples - class group" begin
    for i = 1 : length(Xs)
        X = Xs[i]
        @test rank(class_group(X)) == df[i, "classGroupRank"]
        @test class_group(X).snf[1 : end - rank(class_group(X))] == df[i, "classGroupTorsion"]
    end
end

#@testset "examples - is fano" begin
    #for i = 1 : length(Xs)
        #@test is_fano(Xs[i])
    #end
#end

@testset "examples - gorenstein index" begin
    for i = 1 : length(Xs)
        @test gorenstein_index(Xs[i]) == df[i, "gorensteinIndex"]
    end
end

@testset "examples - picard index" begin
    for i = 1 : length(Xs)
        @test picard_index(Xs[i]) == df[i, "picardIndex"]
    end
end

@testset "examples - orientation" begin
    for i = 1 : length(Xs)
        @test orientation(Xs[i]) == df[i, "orientation"]
    end
end

@testset "examples - intersection matrix" begin
    for i = 1 : length(Xs)
        @test intersection_matrix(Xs[i]) == df[i, "intersectionMatrix"]
    end
end

@testset "examples - anticanonical self intersection" begin
    for i = 1 : length(Xs)
        @test anticanonical_self_intersection(Xs[i]) == df[i, "anticanonicalSelfIntersection"]
    end
end

@testset "examples - maximal log canonicity" begin
    for i = 1 : length(Xs)
        @test maximal_log_canonicity(Xs[i]) == df[i, "maximalLogCanonicity"]
    end
end

@testset "examples - is normal form" begin
    for i = 1 : length(Xs)
        # the non-toric entries should all be in normal form
        if nblocks(Xs[i]) > 2
            @test is_normal_form(Xs[i])
        end
    end
end

