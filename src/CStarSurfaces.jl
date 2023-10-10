module CStarSurfaces

include("imports.jl")
include("exports.jl")
include("Tools.jl")


include("MoriDreamSpace/MoriDreamSpace.jl")
include("MoriDreamSpace/MoriDreamSpaceDivisor.jl")
#include("MoriDreamSpace/SubvarietyOfToricVariety.jl")
include("MoriDreamSpace/ToricVarietyMDS.jl")

include("ToricGeometry/AbstractNormalToricVariety.jl")
include("ToricGeometry/ToricDivisor.jl")
include("ToricGeometry/ToricDivisorClass.jl")

include("ToricSurface/ToricSurfaceResolution.jl")
include("ToricSurface/ToricSurface.jl")
include("ToricSurface/normal_form.jl")
include("ToricSurface/ToricSurfaceDivisor.jl")

include("CStarSurface/ZeroVector.jl")
include("CStarSurface/DoubleVector.jl")
include("CStarSurface/CStarSurface.jl")
include("CStarSurface/CStarSurfaceDivisor.jl")
include("CStarSurface/AdmissibleOperation.jl")
include("CStarSurface/normal_form.jl")

include("SurfaceWithTorusAction/SurfaceWithTorusAction.jl")
include("SurfaceWithTorusAction/SurfaceWithTorusActionDivisor.jl")

include("Database/DatabaseAdapter.jl")
include("Database/SQLiteAdapter.jl")
include("Database/SQLiteAdapterSurfaces.jl")

end 
