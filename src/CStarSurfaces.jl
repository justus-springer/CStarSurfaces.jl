module CStarSurfaces

include("imports.jl")
include("exports.jl")
include("Tools.jl")

include("types.jl")

include("SingularityType.jl")

include("MoriDreamSpace/MoriDreamSpace.jl")
include("MoriDreamSpace/MoriDreamSpaceUnion.jl")
include("MoriDreamSpace/Divisor.jl")
#include("MoriDreamSpace/SubvarietyOfToricVariety.jl")
include("MoriDreamSpace/ToricVarietyMDS.jl")
include("MoriDreamSpace/Point.jl")

include("ToricGeometry/AbstractNormalToricVariety.jl")
include("ToricGeometry/ToricDivisor.jl")
include("ToricGeometry/ToricDivisorClass.jl")

include("ToricSurface/Resolution.jl")
include("ToricSurface/ToricSurface.jl")
include("ToricSurface/normal_form.jl")
include("ToricSurface/ToricSurfaceDivisor.jl")
include("ToricSurface/Point.jl")
include("ToricSurface/FixedPoint.jl")

include("CStarSurface/ZeroVector.jl")
include("CStarSurface/DoubleVector.jl")
include("CStarSurface/Case.jl")
include("CStarSurface/CStarSurface.jl")
include("CStarSurface/Divisor.jl")
include("CStarSurface/AdmissibleOperation.jl")
include("CStarSurface/normal_form.jl")
include("CStarSurface/Point.jl")
include("CStarSurface/FixedPoint.jl")

include("SurfaceWithTorusAction/SurfaceWithTorusAction.jl")
include("SurfaceWithTorusAction/Divisor.jl")
include("SurfaceWithTorusAction/Point.jl")
include("SurfaceWithTorusAction/FixedPoint.jl")

include("Database/DatabaseAdapter.jl")
include("Database/SQLiteAdapter.jl")
include("Database/SQLiteAdapterSurfaces.jl")

end 
