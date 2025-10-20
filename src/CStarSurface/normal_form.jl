@doc raw"""
    mfrak_plus(X :: CStarSurface)

The sum of the rounded down slopes ``\lfloor m_{i1} \rfloor``, see
Definition ``\ref{def:defining_triple_magic_numbers}``.

"""
function mfrak_plus(X :: CStarSurface{T,C}) where {C, T <: Integer}
    r = number_of_blocks(X) - 1
    return sum([floor(T, slope(X, i, 1)) for i = 0 : r])
end


@doc raw"""
    mfrak_minus(X :: CStarSurface)

Minus the sum of the rounded up slopes ``\lceil m_{in_i} \rceil``, see
Definition ``\ref{def:defining_triple_magic_numbers}``.

"""
function mfrak_minus(X :: CStarSurface{T,C}) where {C, T <: Integer}
    r = number_of_blocks(X) - 1
    ns = block_sizes(X)
    return -sum([ceil(T, slope(X, i, ns[i+1])) for i = 0 : r])
end


@doc raw"""
    beta_plus(X :: CStarSurface)

The nested vector with entries ``m_{ij} - \lfloor m_{i1} \rfloor``, see
Definition ``\ref{def:defining_triple_magic_numbers}``.
It is not yet sorted.

"""
function beta_plus(X :: CStarSurface{T,C}) where {C, T <: Integer}
    r = number_of_blocks(X) - 1
    ns = block_sizes(X)
    return [[slope(X, i, j) - floor(T, slope(X, i, 1))
             for j = 1 : ns[i+1]]
             for i = 0 : r]
end

@doc raw"""
    beta_minus(X :: CStarSurface)

The nested vector with entries ``\lceil m_{in_i} \rceil - m_{ij}``, see
Definition ``\ref{def:defining_triple_magic_numbers}``.
It is not yet sorted.

"""
function beta_minus(X :: CStarSurface{T,C}) where {C, T <: Integer}
    r = number_of_blocks(X) - 1
    ns = block_sizes(X)
    return [[ceil(T, slope(X, i, ns[i+1])) - slope(X, i, ns[i+1]-j+1)
             for j = 1 : ns[i+1]]
             for i = 0 : r]
end


@doc raw"""
    are_equivalent(X :: CStarSurface{T}, Y :: CStarSurface{T}) where {T <: Integer}

Check whether two ``\mathbb{C}^*``-surfaces are isomorphic, i.e. whether the
defining triples are equivalent. This is done by comparing ``(\mathfrak{m}^+, \beta^+)`` and
``(\mathfrak{m}^-, \beta^-)``, see
Lemma ``\ref{lem:admissible_operations_magic_numbers}``.

"""
are_equivalent(X :: CStarSurface{T,C1}, Y :: CStarSurface{T,C2}) where {C1, C2, T <: Integer} =
(C1 == C2 && mfrak_plus(X) == mfrak_plus(Y) && sort(beta_plus(X)) == sort(beta_plus(Y))) ||
(C1 == invert_case(C2) && mfrak_plus(X) == mfrak_minus(Y) && sort(beta_plus(X)) == sort(beta_minus(Y)))


@doc raw"""
    orientation(X :: CStarSurface)

The orientation of a ``\mathbb{C}^*``-surface, see Definition
``\ref{def:defining_triple_orientation}``. This is either ``-1``, ``0``, or ``1``.

"""
function orientation(X :: CStarSurface{T,C}) where {C, T <: Integer}
    C == PE && return 1
    C == EP && return -1
    mfrak_plus(X) > mfrak_minus(X) && return 1
    mfrak_plus(X) < mfrak_minus(X) && return -1
    sort(beta_plus(X)) > sort(beta_minus(X)) && return 1
    sort(beta_plus(X)) < sort(beta_minus(X)) && return -1
    return 0
end


@doc raw"""
    is_normal_form(X :: CStarSurface)

Check whether a (defining triple of a) ``\mathbb{C}^*``-surface is in normal form.
This holds if the following three conditions are satisfied, see Definition ``\ref{def:defining_triple_normal_form}``.

- The [`orientation`](@ref) is non-negative,
- [`beta_plus`](@ref) is sorted lexicographically,
- We have ``0 \leq d_{i1} < l_{i1}`` for all ``i = 1, \dots, r``.

"""
is_normal_form(X :: CStarSurface{T,C,N,M,R}) where {T <: Integer, C, N, M, R} =
orientation(X) ≥ 0 && beta_plus(X) == sort(beta_plus(X)) && all(i -> 0 ≤ d(X, i, 1) < l(X, i, 1), 1 : R-1)


@doc raw"""
    AdmissibleOperation{T<:Integer,R}

An admissible operation of a ``\mathbb{C}^*``-surface. It consists of three fields:

- `invert :: Bool`: Whether the operation contains an inversion,
- `perm :: SVector{R,Int}`: The permutation to apply to the blocks,
- `c :: SVector{R,T}`: The coefficients to apply in the addition.

See also Definition ``\ref{def:admissible_operations}`` and Lemma
``\ref{lem:admissible_operations_normal_form}``. Note that both `perm` and `c` are static vectors of length equal to the number of blocks.
In particular, we require that ``c_1 = -(c_2 + \dots + c_r)``
(in contrast to Section ``\ref{subsec:defining_triples_normal_form}``, the indexation here
is one-based).

Admissible operations can be applied to ``\mathbb{C}^*``-surfaces using standard Julia call syntax.
See for instance the example at [`inversion`](@ref).

"""
struct AdmissibleOperation{T<:Integer,R}
    invert :: Bool
    perm :: SVector{R,Int}
    c :: SVector{R,T}
end


@doc raw"""
    admissible_operation(invert :: Bool, perm :: SVector{R,Int}, c :: SVector{S,T}) where {R,S,T<:Integer}

Construct an admissible operation from an inversion, a permutation and an addition.
Here, we must have ``S = R-1``.

"""
function admissible_operation(invert :: Bool, perm :: SVector{R,Int}, c :: SVector{S,T}) where {R,S,T<:Integer}
    R == S+1 || error("length of perm must equal length of c plus one")
    return AdmissibleOperation{T,R}(invert, perm, SVector{R,T}(-sum(c), c...))
end


@doc raw"""
    inversion(R :: Int, T :: Type{<:Integer} = Int)

Return an inversion as an admissible operation, see Definition ``\ref{def:admissible_operations}``.
It takes the number of blocks of the
``\mathbb{C}^*``-surface and optionally an integer type as input.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[2],[4]], [[-2,-3],[1],[3]])
C*-surface of case PE with l = ((1,1),2,4) and d = ((-2,-3),1,3) 

julia> α = inversion(3)
AdmissibleOperation{Int64, 3}(true, [1, 2, 3], [0, 0, 0])

julia> generator_matrix(α(X))
3×5 SMatrix{3, 5, Int64, 15} with indices SOneTo(3)×SOneTo(5):
 -1  -1   2   0   0
 -1  -1   0   4   0
  3   2  -1  -3  -1
```

"""
inversion(R :: Int, T :: Type{<:Integer} = Int) =
AdmissibleOperation{T,R}(true, SVector{R,Int}(1:R), SVector{R,T}(zeros(T,R)))


@doc raw"""
    permutation(perm :: SVector{R,Int}, T :: Type{<:Integer} = Int) where {R}

Return a permutation as an admissible operation, see Definition ``\ref{def:admissible_operations}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[2],[4]], [[-2,-3],[1],[3]])
C*-surface of case PE with l = ((1,1),2,4) and d = ((-2,-3),1,3) 

julia> α = permutation(@SVector [2,3,1])
AdmissibleOperation{Int64, 3}(false, [2, 3, 1], [0, 0, 0])

julia> generator_matrix(α(X))
3×5 SMatrix{3, 5, Int64, 15} with indices SOneTo(3)×SOneTo(5):
 -2  4   0   0  0
 -2  0   1   1  0
  1  3  -2  -3  1
```

"""
permutation(perm :: SVector{R,Int}, T :: Type{<:Integer} = Int) where {R} =
AdmissibleOperation{T,R}(false, perm, SVector{R,T}(zeros(T,R)))


@doc raw"""
    addition(c :: SVector{R,T}) where {R, T <: Integer}

Return an addition as an admissible operation, see Definition ``\ref{def:admissible_operations}``.

# Example:

```jldoctest
julia> X = cstar_surface(PE, [[1,1],[2],[4]], [[-2,-3],[1],[3]])
C*-surface of case PE with l = ((1,1),2,4) and d = ((-2,-3),1,3) 

julia> α = addition(@SVector [2,-1])
AdmissibleOperation{Int64, 3}(false, [1, 2, 3], [-1, 2, -1])

julia> generator_matrix(α(X))
3×5 SMatrix{3, 5, Int64, 15} with indices SOneTo(3)×SOneTo(5):
 -1  -1  2   0  0
 -1  -1  0   4  0
 -3  -4  5  -1  1
```

"""
addition(c :: SVector{R,T}) where {R, T <: Integer} =
AdmissibleOperation{T,R+1}(false, SVector{R+1,Int}(1:R+1), SVector{R+1}(-sum(c), c...))


@doc raw"""
    (ψ :: AdmissibleOperation{T,R})(X :: CStarSurface{T,C,N,M,R}) where {T <: Integer, C, N, M, R}

Apply an admissible operation to a ``\mathbb{C}^*``-surface.

"""
function (ψ :: AdmissibleOperation{T,R})(X :: CStarSurface{T,C,N,M,R}) where {T <: Integer, C, N, M, R}
    new_C = invert_case(C, ψ.invert)
    new_ns = SVector{R,T}([X.block_sizes[ψ.perm[i]] for i = 1 : R])
    new_V = hcat([SVector{2,T}(l(X, ψ.perm[i]-1, ψ.invert ? new_ns[i]-j+1 : j),
                               (ψ.invert ? -1 : 1) * d(X, ψ.perm[i]-1, ψ.invert ? new_ns[i]-j+1 : j) +
                               ψ.c[i] * l(X, ψ.perm[i]-1, ψ.invert ? new_ns[i]-j+1 : j))
                  for i = 1 : R for j = 1 : new_ns[i]]...)
    return CStarSurface{T,new_C}(new_V, new_ns)
end


@doc raw"""
    normal_form_with_operation(X :: CStarSurface)

Return a pair ``(Y, \psi)``, where ``\psi`` is an admissible operation turning ``X`` into normal form
and ``Y = \psi(X)``, see Proposition ``\ref{prp:defining_triple_normal_form_unique}`` 

# Example:

```jldoctest
julia> X = cstar_surface(EP, [[2],[1,1],[4]], [[-5],[2,1],[9]])
C*-surface of case EP with l = (2,(1,1),4) and d = (-5,(2,1),9) 

julia> Y, α = normal_form_with_operation(X)
(C*-surface of case PE with l = ((1,1),2,4) and d = ((-2,-3),1,3) , AdmissibleOperation{Int64, 3}(true, [2, 1, 3], [-1, -2, 3]))

julia> α(X) == Y
true
```

"""
function normal_form_with_operation(X :: CStarSurface{T,C,N,M,R}) where {T <: Integer, C, N, M, R}
    invert = orientation(X) < 0
    perm = SVector{R,T}(sortperm(invert ? beta_minus(X) : beta_plus(X)))
    c = SVector{R-1,T}([invert ?
                        ceil(slope(X, perm[i]-1, block_sizes(X)[perm[i]])) :
                        -floor(slope(X, perm[i]-1, 1)) for i = 2 : R])
    ψ = admissible_operation(invert, perm, c)
    return (ψ(X), ψ)
end


@doc raw"""
    normal_form(X :: CStarSurface)

Return the normal form of a defining triple of a ``\mathbb{C}^*``-surface, see Proposition ``\ref{prp:defining_triple_normal_form_unique}``.

"""
normal_form(X :: CStarSurface) = normal_form_with_operation(X)[1]
