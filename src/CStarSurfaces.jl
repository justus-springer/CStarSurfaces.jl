module CStarSurfaces

include("imports.jl")
include("exports.jl")
include("Tools.jl")

include("ToricGeometry/AbstractNormalToricVariety.jl")
include("ToricGeometry/ToricDivisor.jl")
include("ToricGeometry/ToricDivisorClass.jl")

include("MoriDreamSpace/MoriDreamSpace.jl")
#include("MoriDreamSpace/SubvarietyOfToricVariety.jl")

include("CStarSurface/ZeroVector.jl")
include("CStarSurface/DoubleVector.jl")
include("CStarSurface/CStarSurface.jl")
include("CStarSurface/AdmissibleOperation.jl")
include("CStarSurface/normal_form.jl")

end 
