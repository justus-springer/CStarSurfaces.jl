using CStarSurfaces
using Oscar
using CSV
using DataFrames
using Test
using SQLite


include("CStarSurface/constructors.jl")
include("CStarSurface/admissible_operations.jl")

include("parse_tsv.jl")
include("CStarSurface/examples.jl")
include("ToricSurface/examples.jl")

include("CStarSurface/isomorphy_test.jl")
include("ToricSurface/isomorphy_test.jl")

include("Database/sqlite_database.jl")
