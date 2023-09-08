
const DoubleVector{T} = ZeroVector{Vector{T}}

DoubleVector(v :: Vector{Vector{T}}) where {T} = ZeroVector{Vector{T}}(v)

DoubleVector{T}(::UndefInitializer, ns :: AbstractVector{Int}) where {T} = ZeroVector([Vector{T}(undef, n) for n in ns])
DoubleVector(::UndefInitializer, ns :: AbstractVector{Int}) = ZeroVector([Vector(undef, n) for n in ns])

map2(f, vs :: DoubleVector...) = map(function(xs...) map(f, xs...) end, vs...)

all2(v :: DoubleVector) = all(all, v)
all2(f, vs :: DoubleVector...) = all2(map2(f, vs...))
