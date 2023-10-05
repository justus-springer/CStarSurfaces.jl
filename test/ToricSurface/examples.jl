# Loads in a list of 200 randomly chosen examples from the ldp-database
filepath = joinpath(dirname(@__FILE__), "toricsurfaces.tsv")
# EXAMPLES_FILE = "test/CStarSurface/cstarsurfaces.tsv"

df = CSV.read(filepath, DataFrame, delim='\t', types=TSV_COLUMN_TYPES)

Xs = map(df_row -> toric_surface(df_row.P), eachrow(df))

@testset "examples (toric) - class group" begin
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

@testset "examples (toric) - gorenstein index" begin
    for i = 1 : length(Xs)
        @test gorenstein_index(Xs[i]) == df[i, "gorensteinIndex"]
    end
end

@testset "examples (toric) - picard index" begin
    for i = 1 : length(Xs)
        @test picard_index(Xs[i]) == df[i, "picardIndex"]
    end
end

@testset "examples (toric) - intersection matrix" begin
    for i = 1 : length(Xs)
        @test intersection_matrix(Xs[i]) == df[i, "intersectionMatrix"]
    end
end

@testset "examples (toric) - anticanonical self intersection" begin
    for i = 1 : length(Xs)
        @test anticanonical_self_intersection(Xs[i]) == df[i, "anticanonicalSelfIntersection"]
    end
end

@testset "examples (toric) - maximal log canonicity" begin
    for i = 1 : length(Xs)
        @test maximal_log_canonicity(Xs[i]) == df[i, "maximalLogCanonicity"]
    end
end

