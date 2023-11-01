
include(filename_for_zerorange)

ZeroVector(v :: OffsetArrays.OffsetVector) = ZeroVector(v.parent)

Vector{T}(v :: ZeroVector{T}) where {T} = v.parent

ZeroVector{T}(::UndefInitializer, m::Integer) where {T} = ZeroVector(Vector{T}(undef, m))
ZeroVector(::UndefInitializer, m::Integer) = ZeroVector(Vector(undef, m))

Base.axes(v :: ZeroVector) = (ZeroRange(length(v.parent)),)

Base.size(v :: ZeroVector) = size(v.parent)

function Base.similar(A::AbstractArray, T::Type, shape::Tuple{ZeroRange,Vararg{ZeroRange}})
    similar(A, T, map(r -> Base.OneTo(length(r)), shape))
end

getindex(a :: ZeroVector, i :: Int) = a.parent[i+1]

function setindex!(a :: ZeroVector{T}, v :: T, i :: Int) where T
    a.parent[i+1] = v
end

firstindex(a :: ZeroVector) = 0

lastindex(a :: ZeroVector) = lastindex(a.parent) - 1

iterate(a :: ZeroVector) = (a[0], 1)
iterate(a :: ZeroVector, s :: Int) = s < size(a,1) ? (a[s], s+1) : nothing

length(a :: ZeroVector) = length(a.parent)

map(f, vs :: ZeroVector...) = ZeroVector(map(f, map(v -> v.parent, vs)...))

function resize!(v :: ZeroVector, n :: Integer)
    resize!(v.parent, n)
    return v
end

permuted(v :: ZeroVector, x :: PermGroupElem) = ZeroVector(permuted(v.parent, x))

