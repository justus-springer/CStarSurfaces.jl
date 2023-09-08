abstract type CStarSurfaceCase end
struct EE <: CStarSurfaceCase end
struct PE <: CStarSurfaceCase end
struct EP <: CStarSurfaceCase end
struct PP <: CStarSurfaceCase end

invert_case(c :: EE) = EE()
invert_case(c :: PE) = EP()
invert_case(c :: EP) = PE()
invert_case(c :: PP) = PP()
invert_case(c :: CStarSurfaceCase, invert :: Bool) = invert ? invert_case(c) : c

@attributes mutable struct CStarSurface{T<:CStarSurfaceCase} <: MoriDreamSpace
    l :: DoubleVector{Int64}
    d :: DoubleVector{Int64}
    case :: T

    CStarSurface{T}(l, d) where {T<:CStarSurfaceCase} = new{T}(l,d,T())
end

Base.:(==)(X :: CStarSurface, Y :: CStarSurface) = X.l == Y.l && X.d == Y.d && X.case == Y.case

#################################################
# Constructors
#################################################

_symbol_to_case_dict = Dict([(:ee, EE()), (:pe, PE()), (:ep, EP()), (:pp, PP())])
_to_case(s :: Symbol) = _symbol_to_case_dict[s]
_to_case(c :: CStarSurfaceCase) = c

function cstar_surface(ls :: DoubleVector{Int64}, ds :: DoubleVector{Int64}, case :: Union{Symbol, CStarSurfaceCase})
    @req length(ls) == length(ds) "ls and ds must be of the same length"
    r = length(ls)
    @req all(i -> length(ls[i]) == length(ds[i]), axes(ls,1)) "ls[i] and ds[i] must be of the same length for all i"
    @req all2((l,d) -> gcd(l,d) == 1, ls, ds) "ls[i][j] and ds[i][j] must be coprime for all i and j"

    c = _to_case(case)
    X = CStarSurface{typeof(c)}(ls, ds)

    _set_coordinate_names_cstar(X)
    return X
end

cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Union{Symbol, CStarSurfaceCase}) = cstar_surface(DoubleVector(ls), DoubleVector(ds), case)

function _bas_vec(r :: Int, i :: Int)
    @req i >= 0 "index cannot be negative"
    @req i <= r+1 "index cannot exceed r+1"
    v = fill(0, r+1)
    if i == 0
        for j = 1 : r
            v[j] = -1
        end
    else
        v[i] = 1
    end
    return v
end

function _is_cstar_column(v :: Vector{T}) where {T <: Oscar.IntegerUnion}
    r = length(v) - 1
    v0 = v[begin : end-1]
    l = maximum(abs.(v0))
    if l == 0 return nothing end
    d = last(v)
    for i = 0 : r
        if v0 == l * basis_vector(T, r, i)
            return (i, Int(l), Int(d))
        end
    end
    return nothing
end

function cstar_surface(P :: ZZMatrix)
    ls = DoubleVector{Int}(undef, 0)
    ds = DoubleVector{Int}(undef, 0)
    cols = [[P[j,i] for j = 1 : nrows(P)] for i = 1 : ncols(P)]
    r = nrows(P) - 1 
    last_i = -1
    for col in cols
        ild = _is_cstar_column(col)
        
        if isnothing(ild)
            break
        end

        (i, l, d) = ild
        if i == last_i
            push!(ls[i], l)
            push!(ds[i], d)
        elseif i == last_i + 1
            push!(ls, [l])
            push!(ds, [d])
        else
            throw(ArgumentError("given matrix is not in P-Matrix shape"))
        end

        last_i = i
    end
    @req last_i == r "given matrix is not in P-Matrix shape"

    m = length(cols) - sum(map(length, ls))
    if m == 0
        case = :ee
    elseif m == 1 && last(cols) == _bas_vec(r,r+1)
        case = :pe
    elseif m == 1 && last(cols) == -_bas_vec(r,r+1)
        case = :ep
    elseif m == 2 && cols[end-1] == _bas_vec(r,r+1) && last(cols) == -_bas_vec(r,r+1)
        case = :pp
    else
        throw(ArgumentError("given matrix is not in P-Matrix shape"))
    end

    return cstar_surface(ls, ds, case)
end


#################################################
# Basic attributes
#################################################


@attr nblocks(X :: CStarSurface) = length(X.l)

@attr slopes(X :: CStarSurface) = map2(//, X.d, X.l)

_r(X :: CStarSurface) = nblocks(X) - 1

_n(X :: CStarSurface, i :: Int) = length(X.l[i])
@attr _ns(X :: CStarSurface) = map(length, X.l)
@attr _n(X :: CStarSurface) = sum(_ns(X))

_m(X :: CStarSurface{EE}) = 0
_m(X :: CStarSurface{PE}) = 1
_m(X :: CStarSurface{EP}) = 1
_m(X :: CStarSurface{PP}) = 2


#################################################
# Construction of canonical toric ambient
#################################################

function _slope_ordered_ray_indices(X :: CStarSurface) 
    is = map(v -> sortperm(v, rev=true), slopes(X))
    N = 0
    for i in axes(is,1)
        is[i] .+= N
        N += length(is[i])
    end
    return is
end

_sigma_plus(X :: CStarSurface) = Vector(map(first, _slope_ordered_ray_indices(X)))

_sigma_minus(X :: CStarSurface) = Vector(map(last, _slope_ordered_ray_indices(X)))

_taus(X :: CStarSurface) = vcat(map(vs -> [[vs[j], vs[j+1]] for j = 1 : length(vs) - 1], _slope_ordered_ray_indices(X))...)

_vplus_index(X :: CStarSurface{PE}) = sum(map(length, X.l)) + 1
_vplus_index(X :: CStarSurface{PP}) = sum(map(length, X.l)) + 1

_taus_plus(X :: CStarSurface) = Vector(map(vs -> [first(vs), _vplus_index(X)], _slope_ordered_ray_indices(X)))

_vminus_index(X :: CStarSurface{EP}) = sum(map(length, X.l)) + 1
_vminus_index(X :: CStarSurface{PP}) = sum(map(length, X.l)) + 2

_taus_minus(X :: CStarSurface) = Vector(map(vs -> [last(vs), _vminus_index(X)], _slope_ordered_ray_indices(X)))

_max_cones_indices(X :: CStarSurface{EE}) = append!(_taus(X), [_sigma_plus(X), _sigma_minus(X)])
_max_cones_indices(X :: CStarSurface{PE}) = append!(_taus(X), [_sigma_minus(X)], _taus_plus(X))
_max_cones_indices(X :: CStarSurface{EP}) = append!(_taus(X), [_sigma_plus(X)], _taus_minus(X))
_max_cones_indices(X :: CStarSurface{PP}) = append!(_taus(X), _taus_plus(X), _taus_minus(X))

_ray(X :: CStarSurface, i :: Int, j :: Int) = 
X.l[i][j] * _bas_vec(_r(X),i) + X.d[i][j] * _bas_vec(_r(X),nblocks(X))

_rays_core(X :: CStarSurface) = [_ray(X, i, j) for i in axes(X.l, 1) for j in axes(X.l[i], 1)]

_vplus(X :: CStarSurface) = _bas_vec(_r(X), nblocks(X))
_vminus(X :: CStarSurface) = -_bas_vec(_r(X), nblocks(X))

_rays(X :: CStarSurface{EE}) = _rays_core(X)
_rays(X :: CStarSurface{PE}) = push!(_rays_core(X), _vplus(X))
_rays(X :: CStarSurface{EP}) = push!(_rays_core(X), _vminus(X))
_rays(X :: CStarSurface{PP}) = push!(_rays_core(X), _vplus(X), _vminus(X))

@attr canonical_toric_ambient(X :: CStarSurface) = normal_toric_variety(_rays(X), _max_cones_indices(X))

# Unfortunately, polymake reorders the rays in unpredictable ways, destroying the
# ordering presrcibed by the double index notation of c-star surfaces.
# This function finds the indices matching to a given ray vector.
function _find_ray_index(X :: CStarSurface, ray :: RayVector{QQFieldElem})
    int_ray = map(Int, lcm(denominator.(ray)) * ray)
    if int_ray == _vplus(X)
        return _vplus_index(X) - _n(X)
    elseif int_ray == _vminus(X)
        return _vminus_index(X) - _n(X)
    end

    (i,l,d) = _is_cstar_column(int_ray)
    j = findfirst(ld -> ld == (l,d), collect(zip(X.l[i], X.d[i])))
    return (i,j)

end

@attr _sorted_index_pairs(X :: CStarSurface) = map(ray -> _find_ray_index(X, ray), rays(X))

#################################################
# Cox Ring
#################################################

_coord_name(i :: Int, j :: Int) = "T[$(i)][$(j)]"
_coord_name(k :: Int) = "S[$(k)]"

# Sets the coordinate names in the cox ring according to the double index notation
function _set_coordinate_names_cstar(X :: CStarSurface)
    coord_names = map(ij -> _coord_name(ij...), _sorted_index_pairs(X))
    set_coordinate_names(X, coord_names)
    set_coordinate_names(canonical_toric_ambient(X), coord_names)
end

# Returns the variables in the Cox Ring of the canonical toric ambient, according
# to the double index notation of C-Star surfaces.
# The result is a tuple, whose first entry is a DoubleVector consisting of the 
# variables T[i][j] and whose second entry is the Vector of variables (S[1], ... S[m])
@attr function cox_ring_vars(X :: CStarSurface)
    cox_gens = gens(cox_ring(canonical_toric_ambient(X)))
    Ts = DoubleVector{MPolyDecRingElem}(undef, _ns(X))
    Ss = Vector{MPolyDecRingElem}(undef, _m(X))
    indices = _sorted_index_pairs(X)
    for k in axes(indices, 1)
        if length(indices[k]) == 1
            Ss[indices[k]] = cox_gens[k]
        else
            (i,j) = indices[k]
            Ts[i][j] = cox_gens[k]
        end
    end
    return (Ts,Ss)
end

function _monomial(X :: CStarSurface, i :: Int)
    T = cox_ring_vars(X)[1]
    return prod([T[i][j]^X.l[i][j] for j = 1 : _n(X,i)])
end

_trinomial(X :: CStarSurface, i :: Int) = _monomial(X, i) + _monomial(X, i+1) + _monomial(X, i+2)

@attr cox_ring_relations(X :: CStarSurface) = [_trinomial(X,i) for i = 0 : _r(X) - 2]


#################################################
# Printing
#################################################


Base.show(io :: IO, X :: CStarSurface{EE}) = print(io, "C-star surface of type (e-e)")
Base.show(io :: IO, X :: CStarSurface{PE}) = print(io, "C-star surface of type (p-e)")
Base.show(io :: IO, X :: CStarSurface{EP}) = print(io, "C-star surface of type (e-p)")
Base.show(io :: IO, X :: CStarSurface{PP}) = print(io, "C-star surface of type (p-p)")

