
ZeroVector{T}(v :: AbstractVector{T}) where {T} = OffsetVector{T}(v, 0:length(v)-1)

ZeroVector(v :: AbstractVector) = OffsetVector(v, 0:length(v)-1)

ZeroVector{T}(::UndefInitializer, m :: Integer) where {T} = ZeroVector(Vector{T}(undef, m))

ZeroVector(::UndefInitializer, m :: Integer) = ZeroVector(Vector(undef, m))

Vector(v :: ZeroVector{T, Vector{T}}) where {T} = v.parent

permuted(v :: ZeroVector, x :: PermGroupElem) = ZeroVector(permuted(v.parent, x))
