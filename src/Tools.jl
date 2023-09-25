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
function basis_vector(::Type{T}, n, i) where {T <: Oscar.IntegerUnion}
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

