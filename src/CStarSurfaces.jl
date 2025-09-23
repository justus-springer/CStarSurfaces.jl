module CStarSurfaces

include("imports.jl")
include("exports.jl")

include("Tools.jl")
include("snf.jl")
include("AbelianGroup.jl")

include("CStarSurface/Case.jl")
include("CStarSurface/CStarSurface.jl")
include("CStarSurface/properties.jl")
include("CStarSurface/normal_form.jl")
include("CStarSurface/io.jl")

include("CStarSurface/FixedPoint/FixedPoint.jl")
include("CStarSurface/FixedPoint/Elliptic.jl")
include("CStarSurface/FixedPoint/Bolic.jl")
include("CStarSurface/FixedPoint/Parabolic.jl")
include("CStarSurface/FixedPoint/Hyperbolic.jl")

include("Classifications/picard_index.jl")

end
