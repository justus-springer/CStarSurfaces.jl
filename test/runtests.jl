using CStarSurfaces
using Oscar
using CSV
using DataFrames
using Test


include("CStarSurface/constructors.jl")
include("CStarSurface/admissible_operations.jl")
include("CStarSurface/isomorphy_test.jl")

include("parse_tsv.jl")
include("CStarSurface/examples.jl")
include("ToricSurface/examples.jl")
