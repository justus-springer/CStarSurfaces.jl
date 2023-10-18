
#################################################
# Constructors
#################################################

@doc raw"""
    cstar_surface(ls :: DoubleVector{Int64}, ds :: DoubleVector{Int64}, case :: Symbol)

Construct a C-star surface from the integral vectors $l_i=(l_{i1}, ...,
l_{in_i})$ and $d_i=(d_{i1}, ..., d_{in_i})$ and a given C-star surface case.
The parameters `ls` and `ds` are given both given as a `DoubleVector`, which is
a zero-indexed vector of one-indexed vectors. They must be of the same length
and satisfy `gcd(ls[i][j], ds[i][j]) == 1` for all `i` and `j`. The parameter
`case` can be one of the four symbols `:ee, :pe, :ep, :pp`.

# Example

The ``E_6`` singular cubic.

```jldoctest
julia> X = cstar_surface(DoubleVector([[3, 1], [3], [2]]), DoubleVector([[-2, -1], [1], [1]]), :ee)
C-star surface of type (e-e)
julia> gen_matrix(X)
[-3   -1   3   0]
[-3   -1   0   2]
[-2   -1   1   1]

```

"""
function cstar_surface(ls :: DoubleVector{Int64}, ds :: DoubleVector{Int64}, case :: Symbol)
    @req length(ls) == length(ds) "ls and ds must be of the same length"
    r = length(ls)
    @req all(i -> length(ls[i]) == length(ds[i]), axes(ls,1)) "ls[i] and ds[i] must be of the same length for all i"
    @req all2((l,d) -> gcd(l,d) == 1, ls, ds) "ls[i][j] and ds[i][j] must be coprime for all i and j"

    return CStarSurface{_case_sym_to_type(case)}(ls, ds, case)

end


@doc raw"""
    cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Symbol)

Construct a C-star surface from the integral vectors $l_i=(l_{i1}, ...,
l_{in_i})$ and $d_i=(d_{i1}, ..., d_{in_i})$ and a given C-star surface case.
The parameters `ls` and `ds` are given both given as a vector of vectors. They
must be of the same length and satisfy `gcd(ls[i][j], ds[i][j]) == 1` for all
`i` and `j`. The parameter `case` can be one of the four symbols `:ee, :pe,
:ep, :pp`.

# Example

The ``E_6`` singular cubic.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)
julia> gen_matrix(X)
[-3   -1   3   0]
[-3   -1   0   2]
[-2   -1   1   1]

```

"""
cstar_surface(ls :: Vector{Vector{Int64}}, ds :: Vector{Vector{Int64}}, case :: Symbol) = cstar_surface(DoubleVector(ls), DoubleVector(ds), case)

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


@doc raw"""
    cstar_surface(P :: ZZMatrix)

Construct a C-star surface from a generator matrix of the correct format. That
is, $P$ must be of one of the following forms:

```math
\begin{array}{lcclc}
\text{(e-e)} & 
\begin{bmatrix}
L \\
d
\end{bmatrix} &
\qquad &
\text{(p-e)} &
\begin{bmatrix}
L & 0 \\
d & 1
\end{bmatrix} \\

\text{(e-p)} &
\begin{bmatrix}
L & 0 \\
d & -1
\end{bmatrix} &
\qquad &
\text{(p-p)} &
\begin{bmatrix}
L & 0 & 0 \\
d & 1 & -1
\end{bmatrix}
\end{array},
```

where for some integral vectors ``l_i=(l_{i1}, ..., l_{in_i}) \in
\mathbb{Z}^{n_i}_{\geq 0}`` and ``d_i=(d_{i1}, ..., d_{in_i}) \in
\mathbb{Z}^{n_i}`` with ``\gcd(l_{ij}, d_{ij}) = 1``, we have

```math
L = \begin{bmatrix}
-l_0 & l_1 & \dots & 0 \\
\vdots & & \ddots & 0 \\
-l_0 & 0 & \dots & l_r
\end{bmatrix}, 
\qquad
d = \begin{bmatrix}
d_0 & \dots & d_r
\end{bmatrix}.
```

# Example

The ``E_6`` singular cubic.

```jldoctest
julia> cstar_surface(ZZ[-3 -1 3 0 ; -3 -1 0 2 ; -2 -1 1 1])
C-star surface of type (e-e)

```

"""
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

@doc raw"""
    has_x_plus(X :: CStarSurface{<:CStarSurfaceCase})

Checks whether a given C-star surface has an elliptic fixed point in the
source, commonly denoted $x^+$.

"""
has_x_plus(X :: CStarSurface{T}) where {T <: CStarSurfaceCase} = has_x_plus(T)


@doc raw"""
    has_x_minus(X :: CStarSurface{<:CStarSurfaceCase})

Checks whether a given C-star surface has an elliptic fixed point in the
sink, commonly denoted $x^-$.

"""
has_x_minus(X :: CStarSurface{T}) where {T <: CStarSurfaceCase} = has_x_minus(T)


@doc raw"""
    has_D_plus(X :: CStarSurface{<:CStarSurfaceCase})

Checks whether a given C-star surface has a parabolic fixed point curve in the
source, commonly denoted $D^+$.

"""
has_D_plus(X :: CStarSurface{T}) where {T <: CStarSurfaceCase} = has_D_plus(T)


@doc raw"""
    has_D_minus(X :: CStarSurface{<:CStarSurfaceCase})

Checks whether a given C-star surface has a parabolic fixed point curve in the
sink, commonly denoted $D^-$.

"""
has_D_minus(X :: CStarSurface{T}) where {T <: CStarSurfaceCase} = has_D_minus(T)

@doc raw"""
    nblocks(X :: CStarSurface)

Returns the number of blocks in the generator matrix of a C-star surface.

"""
@attr nblocks(X :: CStarSurface) = length(X.l)

_r(X :: CStarSurface) = nblocks(X) - 1

_n(X :: CStarSurface, i :: Int) = length(X.l[i])
@attr _ns(X :: CStarSurface) = map(length, X.l)
@attr _n(X :: CStarSurface) = sum(_ns(X))


@doc raw"""
    block_sizes(X :: CStarSurface)   

Returns the sizes of the blocks in the generator matrix of a C-star surface.
The result is a zero-indexed vector of type `ZeroVector`.

"""
block_sizes(X :: CStarSurface) = _ns(X)

_m(X :: CStarSurface{EE}) = 0
_m(X :: CStarSurface{PE}) = 1
_m(X :: CStarSurface{EP}) = 1
_m(X :: CStarSurface{PP}) = 2


@doc raw"""
    number_of_parabolic_fixed_point_curves(X :: CStarSurface)

Returns the number of parabolic fixed point curves of a C-star surface.

"""
number_of_parabolic_fixed_point_curves(X :: CStarSurface) = _m(X)


@doc raw"""
    slopes(X :: CStarSurface)

Returns the `DoubleVector` of slopes of a C-star surface, i.e a `DoubleVector`
with `slopes(X)[i][j] = X.d[i][j] // X.l[i][j]`.

"""
@attr slopes(X :: CStarSurface) = map2(//, X.d, X.l)

@attr m_plus(X :: CStarSurface) = sum(map(maximum, slopes(X)))
@attr m_minus(X :: CStarSurface) = -sum(map(minimum, slopes(X)))

@attr _slope_ordering_permutations(X :: CStarSurface) = map(v -> sortperm(v, rev=true), slopes(X))

@attr _slope_ordered_l(X :: CStarSurface) = map(getindex, X.l, _slope_ordering_permutations(X))

@attr _slope_ordered_d(X :: CStarSurface) = map(getindex, X.d, _slope_ordering_permutations(X))

@attr _l_plus(X :: CStarSurface) = sum([1 // first(_slope_ordered_l(X)[i]) for i = 0 : _r(X)]) - _r(X) + 1
@attr _l_minus(X :: CStarSurface) = sum([1 // last(_slope_ordered_l(X)[i]) for i = 0 : _r(X)]) - _r(X) + 1


@doc raw"""
    is_intrinsic_quadric(X :: CStarSurface)

Checks whether a C-star surface is an intrinsic quadric, i.e its Cox Ring
has a single quadratic relation.

"""
@attr is_intrinsic_quadric(X :: CStarSurface) = 
nblocks(X) == 3 && all(ls -> sum(ls) == 2, X.l)


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

@doc raw"""
    cox_ring_vars(X :: CStarSurface)

Return the variables in the Cox Ring of the canonical toric ambient of a C-star
surface, following the double index notation. The result
is a tuple, whose first entry is a DoubleVector consisting of the variables
T[i][j] and whose second entry is the Vector of variables [S[1], ... S[m]], 
where `m` is the number of parabolic fixed point curves.

"""
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

    j == 0 && return has_x_plus(X) ? -1/m_plus(X) : 0
    j == _ns(X)[i] && return has_x_minus(X) ? -1/m_minus(X) : 0

    ms = sort(slopes(X)[i], rev=true)
    return 1 / (ms[j] - ms[j+1])
end

@attr _mcals(X :: CStarSurface) = ZeroVector([[_mcal(X, i, j) for j = 0 : _ns(X)[i]] for i = 0 : _r(X)])


@doc raw"""
    intersection_matrix(X :: CStarSurface)

Return the matrix of intersection numbers of all restrictions of toric prime
divisors of a C-star surface `X` with each other. The result is a rational `n` x
`n` matrix, where `n = nrays(X)` and the `(i,j)`-th entry is the intersection
number of the toric prime divisors associated to the `i`-th and `j`-th ray
respectively.

# Example

```jldoctest
julia> intersection_matrix(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
[1//3   1   2//3   1]
[   1   3      2   3]
[2//3   2   4//3   2]
[   1   3      2   3]
```

"""
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
        IM[_vplus_index(X), _vplus_index(X)] = - m_plus(X)
    end
    if has_D_minus(X)
        IM[_vminus_index(X), _vminus_index(X)] = - m_minus(X)
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
_d_plus(X :: CStarSurface{EE}) = m_plus(X) // _l_plus(X)
_d_plus(X :: CStarSurface{EP}) = m_plus(X) // _l_plus(X)
_d_plus(X :: CStarSurface{PE}) = 1
_d_plus(X :: CStarSurface{PP}) = 1

_d_minus(X :: CStarSurface{EE}) = m_minus(X) // _l_minus(X)
_d_minus(X :: CStarSurface{EP}) = 1
_d_minus(X :: CStarSurface{PE}) = m_minus(X) // _l_minus(X)
_d_minus(X :: CStarSurface{PP}) = 1


# We extend the ls and ds on both sides to make it easier to loop through all the 
# cones of the fan. Together, these numbers form the outer vertices of the 
# simplicial complex \mathcal{A}_P from arxiv:2302.03095.
_slope_ordered_extended_l(X :: CStarSurface) = map(ls -> [0 ; ls ; 0], _slope_ordered_l(X))
_slope_ordered_extended_d(X :: CStarSurface) = map(ds -> [_d_plus(X) ; ds ; -_d_minus(X)], _slope_ordered_d(X))


@doc raw"""
    canonical_resolution(X :: CStarSurface)

Return the canonical resolution of singularities of a C-star surface surface
`X`. The result is a triple `(Y, ex_rays, discrepancies)` where `Y` is the
smooth C-star surface in the resolution of singularities of `X`, `ex_rays`
contains the rays of the exceptional divisors in the resolution and
`discrepancies` contains their discrepancies. `ex_rays` and `discrepancies`
itself are given as pairs, where the first entry comes from the toric step and
the second entry from the tropical step. The first entry is given as a
`DoubleVector`, adhering to the double index notation of C-star surfaces.

# Example

Resolution of singularities of the $E_6$ singular cubic surface.

```jldoctest
julia> X = cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee)
C-star surface of type (e-e)

julia> (Y, (ex_toric, ex_tropic), (discr_toric, discr_tropic)) = canonical_resolution(X);

julia> gen_matrix(Y)
[-3   -1   -1   -2   3   1   2   1   0   0   0   0    0]
[-3   -1   -1   -2   0   0   0   0   2   1   1   0    0]
[-2   -1    0   -1   1   1   1   0   1   1   0   1   -1]

julia> ex_toric
3-element DoubleVector{Vector{Vector{Int64}}} with indices ZeroRange(3):
 [[[1, 0], [2, -1]], [], []]
 [[[1, 1], [2, 1]], [[1, 0]]]
 [[[1, 1]], [[1, 0]]]

julia> ex_tropic
2-element Vector{Vector{Int64}}:
 [0, 1]
 [0, -1]

julia> discr_toric
3-element DoubleVector{Vector{Rational{Int64}}} with indices ZeroRange(3):
 [[0//1, 0//1], [], []]
 [[0//1, 0//1], [1//1]]
 [[0//1], [2//1]]

julia> discr_tropic
2-element Vector{Rational{Int64}}:
 0//1
 4//1
```

"""
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
        push!(discrepancies_tropic, _l_plus(X) // m_plus(X) - 1)
    end
    if has_x_minus(X)
        push!(ex_rays_tropic, [0,-1])
        push!(discrepancies_tropic, _l_minus(X) // m_minus(X) - 1)
    end

    ex_rays = (ex_rays_toric, ex_rays_tropic)
    discrepancies = (discrepancies_toric, discrepancies_tropic)
    
    return (cstar_surface(new_l, new_d, :pp), ex_rays, discrepancies)
end


@doc raw"""
    maximal_log_canonicity(X :: CStarSurface)

Given a C-star surface $X$, return the maximal rational number $\varepsilon$
such that $X$ is $\varepsilon$-log canonical. By definition, this is the
minimal discrepancy in the resolution of singularities plus one.

# Example

```jldoctest
julia> maximal_log_canonicity(cstar_surface([[3, 1], [3], [2]], [[-2, -1], [1], [1]], :ee))
1//1
```

"""
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

