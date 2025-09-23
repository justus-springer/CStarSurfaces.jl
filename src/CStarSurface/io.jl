
function parse_cstar_surfaces(files :: AbstractVector{String}, _ :: Type{T} = Int64) where {T <: Integer}
    Xs = CStarSurface{T}[]
    for fp ∈ files
        f = open(fp, "r")
        for str ∈ readlines(f)
            substrings = strip.(split(str, ";"))
            C = parse_case(substrings[1])
            block_sizes = parse.(T, split(substrings[2], ","))
            R = length(block_sizes)
            N = sum(block_sizes)
            block_sizes = SVector{R}(block_sizes)
            vertex_matrix = transpose(SMatrix{N,2,T}(parse.(T, split(substrings[3], ","))))
            push!(Xs, CStarSurface{T,C}(vertex_matrix, block_sizes))
        end
        close(f)
    end
    return Xs
end

parse_cstar_surfaces(file :: String, _ :: Type{T} = Int64) where {T <: Integer} =
parse_cstar_surfaces([file], T)

function write_cstar_surfaces(Xs :: Vector{<:CStarSurface{T}}, filepath :: String, mode :: String = "w") where {T <: Integer}
    f = open(filepath, mode) 
    for X ∈ Xs
        println(f, case(X), ";", join(block_sizes(X), ","), ";", join(transpose(vertex_matrix(X)), ","))
    end
    close(f)
end
    
