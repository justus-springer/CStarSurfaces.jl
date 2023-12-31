

DoubleVector{T}(v :: AbstractVector{Vector{T}}) where {T} = ZeroVector{T}(v)

DoubleVector(v :: AbstractVector{Vector{T}}) where {T} = ZeroVector(v)

DoubleVector{T}(::UndefInitializer, ns :: AbstractVector{Int}) where {T} = ZeroVector([Vector{T}(undef, n) for n in ns])

DoubleVector(::UndefInitializer, ns :: AbstractVector{Int}) = ZeroVector([Vector(undef, n) for n in ns])

function deleteat!(d :: DoubleVector, i :: Integer, j :: Integer)
    deleteat!(d[i], j)
    return d
end


@doc raw"""
    map2(f, vs :: DoubleVector...)   

Apply a function `f` to all inner elements of a `DoubleVector`.

"""
map2(f, vs :: DoubleVector...) = map(function(xs...) map(f, xs...) end, vs...)

@doc raw"""
    all2(v :: DoubleVector)   

Test whether all inner elements of a `DoubleVector` are `true`.

"""
all2(v :: DoubleVector) = all(all, v)
all2(f, vs :: DoubleVector...) = all2(map2(f, vs...))
