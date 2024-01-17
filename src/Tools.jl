#################################################
# Reset all attributes
#################################################

function reset_attributes!(x)
    try
        empty!(x.__attrs)
    finally
        return x
    end
end

reset_attributes!(xs :: AbstractVector) = map(reset_attributes!, xs)

# i-th canonical basis vector, with the convention that the zeroth 
# one has minus ones everywhere
function basis_vector(::Type{T}, n, i) where {T <: IntegerUnion}
    @req i >= 0 "index cannot be negative"
    @req i <= n "index cannot exceed n"
    i == 0 && return fill(T(-1), n)
    v = fill(T(0), n)
    v[i] = T(1)
    return v
end
basis_vector(n, i) = basis_vector(Int, n, i)


#################################################
# Reduced row echelon form
#################################################


function rcef(M :: MatrixElem{T}) where {T <: RingElement}
    r, A = rref(transpose(M))
    (r, transpose(A))
end

function rcef_rational(M :: MatrixElem{T}) where {T <: RingElement}
    r, A, d = rref_rational(transpose(M))
    (r, transpose(A), d)
end

is_rcef(M :: MatrixElem{T}) where {T <: RingElement} = is_rref(transpose(M))

is_rcef(M :: MatrixElem{T}) where {T <: FieldElement} = is_rref(transpose(M))


#################################################
# Check for complete intersection
#################################################

is_complete_intersection(I :: MPolyIdeal) = 
length(minimal_generating_set(I)) == codim(I)

is_complete_intersection(R :: MPolyQuoRing) = is_complete_intersection(R.I)

is_complete_intersection(R :: MPolyRing) = true


#################################################
# missing polymake interfaace functions
#################################################

centroid(P :: Polyhedron{QQFieldElem}) = P.pm_polytope.CENTROID[begin+1 : end]


# Checks whether at most two entries in a list are greater than one
_almost_all_one(v :: AbstractVector) = length(filter(x -> x > 1, v)) <= 2

# Checks if a list of positive integers is a platonic tuple
_is_platonic_tuple(v :: AbstractVector) = sum([1 // v[i] for i = 1 : length(v)]) > length(v) - 2

function _platonicity_type(v :: AbstractVector)
    @req _is_platonic_tuple(v) "Not a platonic tuple"
    !_is_platonic_tuple(v) && th
    q = sort(v, rev=true)
    q0, q1, q2 = sort(v, rev=true)[1 : 3]
    q2 == 1 && return :A
    (q1, q2) == (2,2) && return :D
    (q0, q1, q2) == (3,3,2) && return :E6
    (q0, q1, q2) == (4,3,2) && return :E7
    (q0, q1, q2) == (5,3,2) && return :E8
end

#################################################
# Methods for Hecke abelian groups
#################################################

free_part(A :: GrpAbFinGen) = free_abelian_group(rank(A))

free_part(x :: GrpAbFinGenElem, A :: GrpAbFinGen) = A(x.coeff[1,end - rank(A) + 1:end])

free_part(x :: GrpAbFinGenElem) = free_part(x, free_part(x.parent))





