is_almost_free(xs :: AbstractVector) =
all([gcd(deleteat!(copy(xs), k)) == 1 for k = 1 : length(xs)])

@doc raw"""
    decompositions(n :: T, len :: Int) where {T <: Integer}

Return all decompositions of `n` into a product of `len` positive integers.

"""
function decompositions(n :: T, len :: Int) where {T <: Integer}
    len == 1 && return [[n]]
    res = Vector{T}[]
    for x = 1 : n
        m, m_rem = divrem(n, x)
        m_rem == 0 || continue
        for xs in decompositions(m, len - 1)
            push!(res, push!(xs, x))
        end
    end

    return res
end


@doc raw"""
    VarBound{T <: Integer}

A simple struct with two fields `val :: T` and `eq :: Bool`. If `eq` is true, it is
understood to be the condition that a variable equals `val`. If `eq` is false, it is
the condition that a variable is at least `val`.
"""
struct VarBound{T <: Integer}
    val :: T
    eq :: Bool
    VarBound(x :: T, eq :: Bool) where {T <: Integer} = new{T}(x, eq)
end

eq(x :: Integer) = VarBound(x, true)
geq(x :: Integer) = VarBound(x, false)

function Base.show(io::IO, vb::VarBound{T}) where T <: Integer
    if vb.eq
        print(io, "=", vb.val)
    else
        print(io, "≥", vb.val)
    end
end


@enum SingularityType eAeA eAeD eDeD eAeE eDeE eEeE eDp eEp

function classify_by_picard_index(p :: T, ::Val{eAeA}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    append!(res, classify_ee_from_l0(p, 1, 1, geq(2), geq(2))) # ((1, 1), ≥2, ≥2)
    return res
end

function classify_by_picard_index(p :: T, ::Val{eAeD}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    append!(res, classify_ee_from_ls(p, eq(1), geq(2), 2, 2)) # ((1, ≥2), 2, 2)
    append!(res, classify_ee_from_l0(p, 1, 2, geq(3), eq(2))) # ((1, 2), ≥3, 2)
    return res
end

function classify_by_picard_index(p :: T, ::Val{eDeD}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    append!(res, classify_ee_from_ls(p, geq(2), geq(2), 2, 2)) # ((≥2, ≥2), 2, 2)
    append!(res, classify_ee_from_l0(p, 2, 2, geq(3), eq(2))) # ((2, 2), ≥3, 2)
    append!(res, classify_ee_from_l0(p, 1, 1, geq(2), eq(2), eq(2))) # ((1, 1), ≥2, 2, 2)
    return res
end

function classify_by_picard_index(p :: T, ::Val{eAeE}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    for (l01, l02, l1, l2) in [(1,3,3,2), (1,4,3,2), (1,5,3,2), (1,3,4,2),
                               (1,3,5,2), (1,2,3,3), (1,2,4,3), (1,2,5,3)]
        append!(res, classify_ee_from_l0(p, l01, l02, eq(l1), eq(l2)))
    end
    return res
end

function classify_by_picard_index(p :: T, ::Val{eDeE}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    for (l01, l02, l1, l2) in [(2,3,3,2), (2,3,4,2), (2,3,5,2), (2,4,3,2), (2,5,3,2)]
        append!(res, classify_ee_from_l0(p, l01, l02, eq(l1), eq(l2)))
    end
    return res
end

function classify_by_picard_index(p :: T, ::Val{eEeE}) where {T <: Integer}
    res = CStarSurface{T, EE}[]
    for (l01, l02, l1, l2) in [(2,2,3,3), (2,2,4,3), (2,2,5,3), (3,3,3,2), (3,4,3,2), (3,5,3,2), 
                               (4,4,3,2), (4,5,3,2), (5,5,3,2), (3,3,4,2), (3,3,5,2)]
        append!(res, classify_ee_from_l0(p, l01, l02, eq(l1), eq(l2)))
    end
    for (l01, l02, l1, l2, l3) in [(1,1,3,3,2), (1,1,4,3,2), (1,1,5,3,2)]
        append!(res, classify_ee_from_l0(p, l01, l02, eq(l1), eq(l2), eq(l3)))
    end

    return res
end

function classify_by_picard_index(p :: T, ::Val{eDp}) where {T <: Integer}
    res = CStarSurface{T, PE}[]
    append!(res, classify_pe(p, geq(2), eq(2), eq(2)))
    return res
end


function classify_by_picard_index(p :: T, ::Val{eEp}) where {T <: Integer}
    res = CStarSurface{T, PE}[]
    for (l0, l1, l2) in [(5,3,2), (4,3,2), (3,3,2)]
        append!(res, classify_pe(p, eq(l0), eq(l1), eq(l2)))
    end
    return res
end


function classify_by_picard_index(p :: T) where {T <: Integer}
    res = CStarSurface{T}[]
    for c in instances(SingularityType)
        append!(res, classify_by_picard_index(p, Val(c)))
    end
    return res
end


function classify_ee_from_l0(p :: T, l01 :: T, l02 :: T, ls_bounds :: VarBound{T}...) where {T <: Integer}

    R = length(ls_bounds)
    N = R+2
    res = CStarSurface{T, EE, N, 2N, R+1}[]

    for μ = 1 : p
        p1, p_rem = divrem(p, μ)
        p_rem == 0 || continue
        for μ01 = 1 : p1
            p2, p_rem = divrem(p1, μ01)
            p_rem == 0 || continue
            for w01 = 1 : p2
                w02, w02_rem = divrem(p2, w01)
                w02_rem == 0 || continue

                relation_degree = (l01 * w01 + l02 * w02)

                all([relation_degree % ls_bounds[i].val == 0 for i = 1 : R if ls_bounds[i].eq]) || continue

                ls_ranges = [ls_bounds[i].eq ? [ls_bounds[i].val] :
                             filter(l -> relation_degree % l == 0, ls_bounds[i].val : relation_degree) for i = 1 : R]

                for ls in Iterators.product(ls_ranges...)
                    ws = relation_degree .÷ ls
                    is_almost_free(vcat([w01, w02], ws...)) || continue
                    ds_ranges = [filter(d -> gcd(ls[i], d) == 1, 1 : ls[i]-1) for i = 1 : R]
                    for ds in Iterators.product(ds_ranges...)
                        d01, d01_rem = divrem(μ * w02 - l01 * sum(ds[i] * prod(ls[j] for j = 1 : R if j ≠ i) for i = 1 : R), prod(ls))
                        d01_rem == 0 || continue
                        gcd(l01, d01) == 1 || continue
                        d02, d02_rem = divrem(l02 * d01 - μ01, l01)
                        d02_rem == 0 || continue
                        gcd(l02, d02) == 1 || continue
                        d01 * w01 + d02 * w02 + sum(ds[i] * ws[i] for i = 1 : R) == 0 || continue
                        X = CStarSurface{T,EE}([l01 l02 ; d01 d02], [[ls[i] ds[i]]' for i = 1 : R]...)
                        all(Y -> !are_equivalent(X, Y), res) || continue
                        push!(res, X)
                    end
                end
            end
        end
    end

    return res

end


function classify_ee_from_ls(p :: T, l01_bound :: VarBound{T}, l02_bound :: VarBound{T}, ls :: T...) where {T <: Integer}

    R = length(ls)
    N = R+2
    res = CStarSurface{T, EE, N, 2N, R+1}[]

    for μ = 1 : p
        p1, p_rem = divrem(p, μ)
        p_rem == 0 || continue
        for μ01 = 1 : p1
            p2, p_rem = divrem(p1, μ01)
            p_rem == 0 || continue
            for w01 = 1 : p2
                w02, w02_rem = divrem(p2, w01)
                w02_rem == 0 || continue

                l01_range = l01_bound.eq ? [l01_bound.val] : (l01_bound.val : fld(prod(ls) * μ01, μ * w01))
                for l01 in l01_range

                    l02, l02_rem = divrem(prod(ls) * μ01 - l01 * μ * w01, μ * w02)
                    l02_rem == 0 || continue
                    if l02_bound.eq
                        l02 == l02_bound.val || continue
                    else
                        l02 ≥ l02_bound.val || continue
                    end

                    relation_degree = (l01 * w01 + l02 * w02)
                    all([relation_degree % ls[i] == 0 for i = 1 : R]) || continue
                    ws = relation_degree .÷ ls
                    is_almost_free(vcat([w01, w02], ws...)) || continue

                    ds_ranges = [filter(d -> gcd(ls[i], d) == 1, 1 : ls[i]-1) for i = 1 : R]

                    for ds in Iterators.product(ds_ranges...)
                        d01, d01_rem = divrem(μ * w02 - l01 * sum(ds[i] * prod(ls[j] for j = 1 : R if j ≠ i) for i = 1 : R), prod(ls))
                        d01_rem == 0 || continue
                        gcd(l01, d01) == 1 || continue
                        d02, d02_rem = divrem(l02 * d01 - μ01, l01)
                        d02_rem == 0 || continue
                        gcd(l02, d02) == 1 || continue
                        d01 * w01 + d02 * w02 + sum(ds[i] * ws[i] for i = 1 : R) == 0 || continue
                        X = CStarSurface{T, EE}([l01 l02 ; d01 d02], [[ls[i] ds[i]]' for i = 1 : R]...)
                        all(Y -> !are_equivalent(X, Y), res) || continue
                        push!(res, X)
                    end
                end
            end
        end
    end

    return res

end

function classify_pe(p :: T, ls_bounds :: VarBound{T}...) where {T <: Integer}
    R = length(ls_bounds)-1
    res = CStarSurface{T, PE, R+1, 2*(R+1), R+1}[]

    nonfixed_indices = filter(i -> !ls_bounds[i].eq, 1 : R+1)

    # We have that the picard index p equals l0 * ... * l1 * wplus
    # Hence the product of the already fixed ls must divide p
    p0, p0_rem = divrem(p, prod([l.val for l in ls_bounds if l.eq]))
    p0_rem == 0 || return res
    
    # We now decompose into the remaining factors and wplus
    for decomp in decompositions(p0, length(nonfixed_indices) + 1)

        all(j -> decomp[j] ≥ ls_bounds[nonfixed_indices[j]].val, 1 : length(nonfixed_indices)) || continue
        wplus = last(decomp)

        ls = [l.val for l in ls_bounds]
        for j = 1 : length(nonfixed_indices)
            ls[nonfixed_indices[j]] = decomp[j]
        end

        μ_max = gcd([prod(ls) ÷ l for l in ls])
        ds_ranges = [filter(d -> gcd(ls[i], d) == 1, 1 : ls[i]-1) for i = 2 : R+1]

        for μ = 1 : μ_max
            μ_max % μ == 0 || continue
            ws = [prod(ls) ÷ (ls[i] * μ) for i = 1 : R+1]
            is_almost_free(vcat([wplus], ws)) || continue

            for ds in Iterators.product(ds_ranges...)
                d0, d0_rem = divrem(-μ*wplus - ls[1] * sum(ds[i] * prod(ls[j] for j = 2 : R+1 if j ≠ i+1) for i = 1 : R), prod(ls[2:end]))
                d0_rem == 0 || continue

                gcd(d0, ls[1]) == 1 || continue

                d0 * ws[1] + sum(ds[i] * ws[i+1] for i = 1 : R) + wplus == 0 || continue

                X = CStarSurface{T, PE}([ls[1] d0]', [[ls[i+1] ds[i]]' for i = 1 : R]...)
                all(Y -> !are_equivalent(X, Y), res) || continue
                push!(res, X)

            end
        end
    end

    return res

end
