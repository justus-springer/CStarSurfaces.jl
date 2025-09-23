using CStarSurfaces, CSV, Test, StaticArrays, DataFrames

include("parse_csv.jl")

filepath = joinpath(dirname(@__FILE__), "examples.tsv")

df = CSV.read(filepath, DataFrame, delim='\t', types=TSV_COLUMN_TYPES)

function to_CStarSurface(row :: DataFrameRow)
    V, ns, c = row.V, row.ns, row.case
    return CStarSurface{Int64,c}(SMatrix{size(V,1),size(V,2)}(V), SVector{length(ns)}(ns))
end

Xs = map(to_CStarSurface, eachrow(df))

@testset "Class group" begin
    for i = 1 : length(Xs)
        G = class_group(Xs[i])
        @test rank(G) == df[i, "class_group_rank"]
        @test elementary_divisors(G) == df[i, "class_group_torsion"]
    end
end

@testset "Picard index" begin
    for i = 1 : length(Xs)
        @test picard_index(Xs[i]) == df[i, "picard_index"]
    end
end

@testset "Gorenstein index" begin
    for i = 1 : length(Xs)
        @test gorenstein_index(Xs[i]) == df[i, "gorenstein_index"]
    end
end

@testset "Log canonicity" begin
    for i = 1 : length(Xs)
        @test log_canonicity(Xs[i]) == df[i, "maximal_log_canonicity"]
    end
end

@testset "Degree" begin
    for i = 1 : length(Xs)
        @test degree(Xs[i]) == df[i, "anticanonical_self_intersection"]
    end
end


