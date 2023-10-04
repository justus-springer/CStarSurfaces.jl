abstract type CStarSurfaceCase end
struct EE <: CStarSurfaceCase end
struct PE <: CStarSurfaceCase end
struct EP <: CStarSurfaceCase end
struct PP <: CStarSurfaceCase end

invert_case(:: EE) = EE()
invert_case(:: PE) = EP()
invert_case(:: EP) = PE()
invert_case(:: PP) = PP()
invert_case(c :: CStarSurfaceCase, invert :: Bool) = invert ? invert_case(c) : c

has_x_plus(:: EE) = true
has_x_plus(:: PE) = false
has_x_plus(:: EP) = true
has_x_plus(:: PP) = false

has_x_minus(c :: CStarSurfaceCase) = has_x_plus(invert_case(c))

has_D_plus(c :: CStarSurfaceCase) = !has_x_plus(c)
has_D_minus(c :: CStarSurfaceCase) = !has_x_minus(c)

@attributes mutable struct CStarSurface{T<:CStarSurfaceCase} <: MoriDreamSpace
    l :: DoubleVector{Int64}
    d :: DoubleVector{Int64}
    case :: T

    CStarSurface{T}(l, d) where {T<:CStarSurfaceCase} = new{T}(l,d,T())
end

Base.:(==)(X :: CStarSurface, Y :: CStarSurface) = X.l == Y.l && X.d == Y.d && X.case == Y.case

has_x_plus(X :: CStarSurface) = has_x_plus(X.case)
has_x_minus(X :: CStarSurface) = has_x_minus(X.case)
has_D_plus(X :: CStarSurface) = has_D_plus(X.case)
has_D_minus(X :: CStarSurface) = has_D_minus(X.case)

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
    return CStarSurface{typeof(c)}(ls, ds)

end

"""
    cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Symbol)

Construct a C-Star surface from the given data
"""
cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Union{Symbol, CStarSurfaceCase}) = cstar_surface(DoubleVector(ls), DoubleVector(ds), case)

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
    elseif m == 1 && last(cols) == basis_vector(r+1,r+1)
        case = :pe
    elseif m == 1 && last(cols) == -basis_vector(r+1,r+1)
        case = :ep
    elseif m == 2 && cols[end-1] == basis_vector(r+1,r+1) && last(cols) == -basis_vector(r+1,r+1)
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

_r(X :: CStarSurface) = nblocks(X) - 1

_n(X :: CStarSurface, i :: Int) = length(X.l[i])
@attr _ns(X :: CStarSurface) = map(length, X.l)
@attr _n(X :: CStarSurface) = sum(_ns(X))

_m(X :: CStarSurface{EE}) = 0
_m(X :: CStarSurface{PE}) = 1
_m(X :: CStarSurface{EP}) = 1
_m(X :: CStarSurface{PP}) = 2

@attr slopes(X :: CStarSurface) = map2(//, X.d, X.l)

@attr _m_plus(X :: CStarSurface) = sum(map(maximum, slopes(X)))
@attr _m_minus(X :: CStarSurface) = -sum(map(minimum, slopes(X)))

@attr _slope_ordering_permutations(X :: CStarSurface) = map(v -> sortperm(v, rev=true), slopes(X))

@attr _slope_ordered_l(X :: CStarSurface) = map(getindex, X.l, _slope_ordering_permutations(X))

@attr _slope_ordered_d(X :: CStarSurface) = map(getindex, X.d, _slope_ordering_permutations(X))

@attr _l_plus(X :: CStarSurface) = sum([1 // first(_slope_ordered_l(X)[i]) for i = 0 : _r(X)]) - _r(X) + 1
@attr _l_minus(X :: CStarSurface) = sum([1 // last(_slope_ordered_l(X)[i]) for i = 0 : _r(X)]) - _r(X) + 1



#################################################
# Construction of canonical toric ambient
#################################################


@attr function _slope_ordered_ray_indices(X :: CStarSurface) 
    ns, r = _ns(X), _r(X)
    is = ZeroVector([sum(ns[0 : i-1]) .+ collect(1 : ns[i]) for i = 0 : r])
    return map(getindex, is, _slope_ordering_permutations(X))
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

@attr maximal_cones_indices(X :: CStarSurface{EE}) = append!(_taus(X), [_sigma_plus(X), _sigma_minus(X)])
@attr maximal_cones_indices(X :: CStarSurface{PE}) = append!(_taus(X), [_sigma_minus(X)], _taus_plus(X))
@attr maximal_cones_indices(X :: CStarSurface{EP}) = append!(_taus(X), [_sigma_plus(X)], _taus_minus(X))
@attr maximal_cones_indices(X :: CStarSurface{PP}) = append!(_taus(X), _taus_plus(X), _taus_minus(X))

_ray(X :: CStarSurface, i :: Int, j :: Int) = 
X.l[i][j] * [basis_vector(_r(X),i) ; 0] + X.d[i][j] * basis_vector(_r(X)+1, _r(X)+1)

_rays_core(X :: CStarSurface) = [_ray(X, i, j) for i in axes(X.l, 1) for j in axes(X.l[i], 1)]

_vplus(X :: CStarSurface) = basis_vector(_r(X)+1, _r(X)+1)
_vminus(X :: CStarSurface) = -basis_vector(_r(X)+1, _r(X)+1)

rays(X :: CStarSurface{EE}) = _rays_core(X)
rays(X :: CStarSurface{PE}) = push!(_rays_core(X), _vplus(X))
rays(X :: CStarSurface{EP}) = push!(_rays_core(X), _vminus(X))
rays(X :: CStarSurface{PP}) = push!(_rays_core(X), _vplus(X), _vminus(X))

_coordinate_names_core(X :: CStarSurface) = ["T[$(i)][$(j)]" for i in axes(X.l, 1) for j in axes(X.l[i], 1)]

_coordinate_names(X :: CStarSurface{EE}) = _coordinate_names_core(X)
_coordinate_names(X :: CStarSurface{PE}) = push!(_coordinate_names_core(X), "S[1]")
_coordinate_names(X :: CStarSurface{EP}) = push!(_coordinate_names_core(X), "S[1]")
_coordinate_names(X :: CStarSurface{PP}) = push!(_coordinate_names_core(X), "S[1]", "S[2]")

@attr function canonical_toric_ambient(X :: CStarSurface) 
    # passing non_redundant=true does two important things here:
    # 1. It skips time-consuming checks on the input rays,
    # 2. It prevents polymake from reordering the rays, making the interface to
    # toric geometry in Oscar much more reliable.
    Z = normal_toric_variety(rays(X), maximal_cones_indices(X); non_redundant=true)
    set_coordinate_names(Z, _coordinate_names(X))
    return Z
end

#################################################
# Cox Ring
#################################################

# Returns the variables in the Cox Ring of the canonical toric ambient, according
# to the double index notation of C-Star surfaces.
# The result is a tuple, whose first entry is a DoubleVector consisting of the 
# variables T[i][j] and whose second entry is the Vector of variables (S[1], ... S[m])
@attr function cox_ring_vars(X :: CStarSurface)
    r, ns, m = _r(X), _ns(X), _m(X)
    Z = canonical_toric_ambient(X)
    cox_gens = gens(cox_ring(Z))

    Ts = DoubleVector{MPolyDecRingElem}(undef, ns)
    k = 1
    for i = 0 : r, j = 1 : ns[i]
        Ts[i][j] = cox_gens[k]
        k += 1
    end

    Ss = cox_gens[end - (m-1) : end]

    return (Ts,Ss)
end

function _monomial(X :: CStarSurface, i :: Int)
    T = cox_ring_vars(X)[1]
    return prod([T[i][j]^X.l[i][j] for j = 1 : _n(X,i)])
end

_trinomial(X :: CStarSurface, i :: Int) = _monomial(X, i) + _monomial(X, i+1) + _monomial(X, i+2)

@attr cox_ring_relations(X :: CStarSurface) :: 
    Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}} = 
[_trinomial(X,i) for i = 0 : _r(X) - 2]

#################################################
# Intersection numbers
#
# We compute intersection numbers following
# Chapter 7 of "Classifying log del Pezzo surfaces
# with torus action" (arXiv:2302.03095)
#################################################

function _mcal(X :: CStarSurface, i :: Int, j :: Int)
    @req 0 ≤ j ≤ _ns(X)[i] "j must be between 0 and n_i"

    j == 0 && return has_x_plus(X) ? -1/_m_plus(X) : 0
    j == _ns(X)[i] && return has_x_minus(X) ? -1/_m_minus(X) : 0

    ms = sort(slopes(X)[i], rev=true)
    return 1 / (ms[j] - ms[j+1])
end

@attr _mcals(X :: CStarSurface) = ZeroVector([[_mcal(X, i, j) for j = 0 : _ns(X)[i]] for i = 0 : _r(X)])

@attr function intersection_matrix(X :: CStarSurface)
    r, n, m, ns = _r(X), _n(X), _m(X), _ns(X)

    # The resulting intersection matrix
    IM = zeros(Rational{Int}, n + m, n + m)

    # The entries ls sorted by slope
    ls = ZeroVector(map(i -> X.l[i][sortperm(slopes(X)[i], rev=true)], 0 : r))

    inds = _slope_ordered_ray_indices(X)

    # 1. Intersection numbers of adjacent rays in the leaf cones
    for i = 0 : r
        if has_D_plus(X) 
            IM[_vplus_index(X), inds[i][1]] = 1 // ls[i][1] 
        end

        for j = 1 : ns[i]-1
            IM[inds[i][j], inds[i][j+1]] = _mcal(X, i, j) // (ls[i][j] * ls[i][j+1])
        end

        if has_D_minus(X) 
            IM[inds[i][ns[i]], _vminus_index(X)] = 1 // ls[i][1] 
        end
    end

    # 2. Intersection numbers of top and bottom most rays with each other
    for i = 0 : r, k = i+1 : r
        if has_x_plus(X)
            IM[inds[i][1], inds[k][1]] = ns[i] * ns[k] == 1 ? 
                -(_mcal(X, i, 0) + _mcal(X, i, ns[i])) // (ls[i][1] * ls[k][1]) :
                -_mcal(X, i, 0) // (ls[i][1] * ls[k][1])
        end
        if has_x_minus(X)
            IM[inds[i][ns[i]], inds[k][ns[k]]] = ns[i] * ns[k] == 1 ? 
                -(_mcal(X, i, 0) + _mcal(X, i, ns[i])) // (ls[i][ns[i]] * ls[k][ns[k]]) :
                -_mcal(X, i, ns[i]) // (ls[i][ns[i]] * ls[k][ns[k]])
        end
    end

    # 3. Self intersection numbers
    if has_D_plus(X)
        IM[_vplus_index(X), _vplus_index(X)] = - _m_plus(X)
    end
    if has_D_minus(X)
        IM[_vminus_index(X), _vminus_index(X)] = - _m_minus(X)
    end
    for i = 0 : r, j = 1 : ns[i]
        IM[inds[i][j], inds[i][j]] = - 1 // ls[i][j]^2 * (_mcal(X, i, j-1) + _mcal(X, i, j))
    end

    # fill in numbers to make the matrix symmetric
    for i = 1 : n+m, j = 1 : n+m
        if IM[i,j] == 0 
            IM[i,j] = IM[j,i]
        end
    end

    return matrix(QQ, IM)

end

anticanonical_self_intersection(X :: CStarSurface) = anticanonical_divisor(X) * anticanonical_divisor(X)


#################################################
# Resolution of singularities
#
# The resolution of singularities for C-star surfaces
# is obtained in two steps: In the tropical step,
# we add the rays [0,...,0,1] and [0,...,0,-1] to 
# the fan, if not already present. The second step
# is toric: We perform regular subdivision all cones
# in the fan, which are now all two-dimensional
#################################################

# See Section 9 of "Classifying log del Pezzo surfaces with 
# torus action" (arxiv:2302.03095) for the definition of 
# these numbers
_d_plus(X :: CStarSurface{EE}) = _m_plus(X) // _l_plus(X)
_d_plus(X :: CStarSurface{EP}) = _m_plus(X) // _l_plus(X)
_d_plus(X :: CStarSurface{PE}) = 1
_d_plus(X :: CStarSurface{PP}) = 1

_d_minus(X :: CStarSurface{EE}) = _m_minus(X) // _l_minus(X)
_d_minus(X :: CStarSurface{EP}) = 1
_d_minus(X :: CStarSurface{PE}) = _m_minus(X) // _l_minus(X)
_d_minus(X :: CStarSurface{PP}) = 1


# We extend the ls and ds on both sides to make it easier to loop through all the 
# cones of the fan. Together, these numbers form the outer vertices of the 
# simplicial complex \mathcal{A}_P from arxiv:2302.03095.
_slope_ordered_extended_l(X :: CStarSurface) = map(ls -> [0 ; ls ; 0], _slope_ordered_l(X))
_slope_ordered_extended_d(X :: CStarSurface) = map(ds -> [_d_plus(X) ; ds ; -_d_minus(X)], _slope_ordered_d(X))

# Computes the canonical resolution of singularities of a C-star surface
# together with the exceptional rays and discrepancies. The output format is:
# (Y, ex_rays, discrepancies), where `Y` is the smooth C-star surface that 
# is the resolution of X. `ex_rays` and `discrepancies` are both given as pairs, 
# where the first entry comes from the toric step and the second entry from the 
# tropical step. The toric ones are given as a `DoubleVector`, adhering the double
# index notation.
@attr function canonical_resolution(X :: CStarSurface)
    ns, r = _ns(X), _r(X) 
    l, d = _slope_ordered_extended_l(X), _slope_ordered_extended_d(X)

    # The toric step
    new_l, new_d = deepcopy(X.l), deepcopy(X.d)
    ex_rays_toric = DoubleVector{Vector{Vector{Int}}}(undef, Vector(ns) .+ 1)
    discrepancies_toric = DoubleVector{Vector{Rational{Int}}}(undef, Vector(ns) .+ 1)
    for i = 0 : r, j = 1 : ns[i] + 1
        v1, v2 = [l[i][j], d[i][j]], [l[i][j+1], d[i][j+1]]
        ex_rays_toric[i][j], discrepancies_toric[i][j] = toric_affine_surface_resolution(v1, v2)
        append!(new_l[i], map(first, ex_rays_toric[i][j]))
        append!(new_d[i], map(last, ex_rays_toric[i][j]))
    end

    # The tropical step
    ex_rays_tropic = Vector{Int}[]
    discrepancies_tropic = Rational{Int}[]
    if has_x_plus(X)
        push!(ex_rays_tropic, [0,1])
        push!(discrepancies_tropic, _l_plus(X) // _m_plus(X) - 1)
    end
    if has_x_minus(X)
        push!(ex_rays_tropic, [0,-1])
        push!(discrepancies_tropic, _l_minus(X) // _m_minus(X) - 1)
    end

    ex_rays = (ex_rays_toric, ex_rays_tropic)
    discrepancies = (discrepancies_toric, discrepancies_tropic)
    
    return (cstar_surface(new_l, new_d, :pp), ex_rays, discrepancies)
end

exceptional_rays(X :: CStarSurface) = canonical_resolution(X)[2]

discrepancies(X :: CStarSurface) = canonical_resolution(X)[3]

@attr function maximal_log_canonicity(X :: CStarSurface) 
    (ds_tor, ds_tro) = discrepancies(X) 
    # we add a superficial zero into the list of discrepancies to ensure a
    # well-defined (and correct) result in case there are no exceptional rays
    # (i.e. the surface is already smooth).
    ds = vcat([0], ds_tro, map(d -> vcat(d...), ds_tor)...)
    # the maximal log canonicity equals the minimal discrepancy plus one
    return minimum(ds) + 1
end


#################################################
# Printing
#################################################


Base.show(io :: IO, X :: CStarSurface{EE}) = print(io, "C-star surface of type (e-e)")
Base.show(io :: IO, X :: CStarSurface{PE}) = print(io, "C-star surface of type (p-e)")
Base.show(io :: IO, X :: CStarSurface{EP}) = print(io, "C-star surface of type (e-p)")
Base.show(io :: IO, X :: CStarSurface{PP}) = print(io, "C-star surface of type (p-p)")

