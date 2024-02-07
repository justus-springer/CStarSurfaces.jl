@doc raw"""
    AdmissibleOperation

Abstract Julia type for admissible operations on C-star surfaces. All subtypes
`T <: AdmissibleOperation` should at least implement the following functions:

```julia
(α :: T)(X :: CStarSurface)
Base.:*(α :: T, β :: T)
Base.one(α :: T)
Base.one(::Type{T})
Base.inv(α :: T)
```

"""
abstract type AdmissibleOperation end


##################################################################################
# Julia type for admissible operations inverting the last row.
#
# It is encoded by an integer `factor` that can be +1 or -1.
##################################################################################

@doc raw"""
    InvertLastRow <: AdmissibleOperation 

Julia type for an admissible operation that may or may not invert the last row
in the generator matrix of a C-star surface. Geometrically, this means swapping
source and sink of the C-star action.

"""
struct InvertLastRow <: AdmissibleOperation 
    factor :: Int # only 1 or -1

    @doc raw"""
        InvertLastRow(factor :: Int)

    Construct an admissible operation that multiplies the last row of the
    generator matrix of a C-star surface by `factor`, which can take the values
    `1` and `-1`.
    
    """
    function InvertLastRow(factor :: Int)
        @req factor == 1 || factor == -1 "factor must be 1 or -1"
        new(factor)
    end
end


@doc raw"""
    InvertLastRow()   

Construct an admissible operation that multiplies the last row of the 
generator matrix of a C-star surface by `-1`.

"""
InvertLastRow() = InvertLastRow(-1)

Base.:(==)(α :: InvertLastRow, β :: InvertLastRow) = α.factor == β.factor

Base.one(::InvertLastRow) = InvertLastRow(1)

Base.one(::Type{InvertLastRow}) = InvertLastRow(1)

Base.:*(α :: InvertLastRow, β :: InvertLastRow) = InvertLastRow(α.factor * β.factor)

Base.inv(α :: InvertLastRow) = InvertLastRow(1 ÷ α.factor)

(α :: InvertLastRow)(X :: CStarSurface) = cstar_surface(X.l, map(d -> α.factor * d, X.d), invert_case(X.case, α.factor == -1))

function Base.show(io :: IO, α :: InvertLastRow) 
    print(io, "Multiply the last row by $(α.factor).")
end


@doc raw"""
    PermutationOfRays <: AdmissibleOperation

Julia type for an admissible operation that permutes rays within the blocks of
the generator matrix of C-star surfaces. Only permutations that are local to
each block are permitted.

"""
struct PermutationOfRays <: AdmissibleOperation
    ray_perms :: ZeroVector{PermGroupElem}

    @doc raw"""
        PermutationOfRays(ray_perms :: ZeroVector{PermGroupElem})   

    Construct an admissible operation that permutes the rays in the `i`-th block
    by the permutation `ray_perms[i]` for all `i`.
    
    """
    PermutationOfRays(ray_perms :: ZeroVector{PermGroupElem}) = new(ray_perms)
end

PermutationOfRays(ray_perms :: Vector{PermGroupElem}) = PermutationOfRays(ZeroVector(ray_perms))

@doc raw"""
    PermutationOfRays(i :: Int, r :: Int, ray_perm :: PermGroupElem)

Construct an admissible operation that permutes the rays in the `i`-th block of
a C-star surface with `r+1` blocks by the permutation `ray_perm`. Here, we 
must have `0 ≤ i ≤ r`.

"""
function PermutationOfRays(i :: Int, r :: Int, ray_perm :: PermGroupElem) 
    @req 0 ≤ i ≤ r "must have 0 ≤ i ≤ r"
    PermutationOfRays(append!(repeat([perm([1])], i), [ray_perm], repeat([perm([1])], r-i)))
end

Base.:(==)(α :: PermutationOfRays, β :: PermutationOfRays) = α.ray_perms == β.ray_perms

function (α :: PermutationOfRays)(X :: CStarSurface)
    @req length(α.ray_perms) == nblocks(X) "number of permutations must equal number of blocks"
    cstar_surface(map(permuted, X.l, α.ray_perms), map(permuted, X.d, α.ray_perms), X.case)
end

Base.one(α :: PermutationOfRays) = PermutationOfRays(map(one, α.ray_perms))

Base.one(::Type{PermutationOfRays}) = PermutationOfRays([perm([1])])

Base.:*(α :: PermutationOfRays, β :: PermutationOfRays) = PermutationOfRays(map(*, α.ray_perms, β.ray_perms))

Base.inv(α :: PermutationOfRays) = PermutationOfRays(map(inv, α.ray_perms))

function Base.show(io :: IO, α :: PermutationOfRays)
    print(io, "Permute the rays in the blocks by: ")
    for i in 0 : length(α.ray_perms) - 2
        print(io, α.ray_perms[i])
        print(", ")
    end
    print(last(α.ray_perms))
    print(".")
end


@doc raw"""
    PermutationOfBlocks <: AdmissibleOperation

Julia type for an admissible operation that permutes the blocks of the
generator matrix of a C-star surface.

"""
struct PermutationOfBlocks <: AdmissibleOperation
    block_perm :: PermGroupElem

    @doc raw"""
        PermutationOfBlocks(block_perm :: PermGroupElem)   

    Construct an admissible operation that permutes the blocks of the 
    generator matrix of a C-star surface by `block_perm`.
    
    """
    PermutationOfBlocks(block_perm :: PermGroupElem) = new(block_perm)
end

Base.:(==)(α :: PermutationOfBlocks, β :: PermutationOfBlocks) = α.block_perm == β.block_perm

Base.one(α :: PermutationOfBlocks) = PermutationOfBlocks(one(α.block_perm))
Base.one(::Type{PermutationOfBlocks}) = PermutationOfBlocks(perm([1]))

(p :: PermutationOfBlocks)(X :: CStarSurface) =
cstar_surface(permuted(X.l, p.block_perm), permuted(X.d, p.block_perm), X.case)

Base.:*(p1 :: PermutationOfBlocks, p2 :: PermutationOfBlocks) = PermutationOfBlocks(p1.block_perm * p2.block_perm)

Base.inv(p :: PermutationOfBlocks) = PermutationOfBlocks(inv(p.block_perm))

function Base.show(io :: IO, α :: PermutationOfBlocks)
    print(io, "Permute the blocks by $(α.block_perm).")
end


@doc raw"""
    AdmissibleRowOperation <: AdmissibleOperation

Julia type for an admissible operation that adds integral multiples of the 
upper rows of the generator matrix of a C-star surface to its last row.

"""
struct AdmissibleRowOperation <: AdmissibleOperation
    factors :: ZeroVector{Int}

    function AdmissibleRowOperation(factors :: ZeroVector{Int})
        @req factors[0] == -sum(factors[1:end]) "The zeroth factor must equal the negative sum of all other factors"
        new(factors)
    end

end

@doc raw"""
    AdmissibleRowOperation(factors :: Vector{Int})

Construct an admissible operation that adds `factors[i]` times the `i`-th
row to the last row of the generator matrix of a C-star surface for all `i`.

"""
AdmissibleRowOperation(factors :: Vector{Int}) = 
AdmissibleRowOperation(ZeroVector(append!([-sum(factors)], factors)))

Base.:(==)(α :: AdmissibleRowOperation, β :: AdmissibleRowOperation) = α.factors == β.factors

Base.one(α :: AdmissibleRowOperation) = AdmissibleRowOperation(ZeroVector(repeat([0], length(α.factors))))

Base.one(α :: Type{AdmissibleRowOperation}) = AdmissibleRowOperation([0])

function (p :: AdmissibleRowOperation)(X :: CStarSurface)
    @req length(p.factors) == nblocks(X) "length of factors does not match the format of the C-star surface"
    new_d = DoubleVector(map(i -> X.d[i] + p.factors[i] * X.l[i], 0 : nblocks(X) - 1))
    cstar_surface(X.l, new_d, X.case)
end

Base.:*(p1 :: AdmissibleRowOperation, p2 :: AdmissibleRowOperation) = AdmissibleRowOperation(map(+, p1.factors, p2.factors))

Base.inv(p :: AdmissibleRowOperation) = AdmissibleRowOperation(map(-, p.factors))

function Base.show(io :: IO, α :: AdmissibleRowOperation)
    print(io, "Add multiples of the upper rows to the last row with factors $(α.factors[1:end]).")
end


@doc raw"""
    CompositeAdmissibleOperation <: AdmissibleOperation

Julia type for an admissible operation that arises as the composition of
multiple admissible operations.

"""
struct CompositeAdmissibleOperation <: AdmissibleOperation
    ops :: Vector{AdmissibleOperation}

    @doc raw"""
        CompositeAdmissibleOperation(ops :: Vector{AdmissibleOperation})    

    Construct an admissible operation that applies the admissible operations `ops`
    sequentially (first `ops[1]`, then `ops[2]` etc).
    
    """
    CompositeAdmissibleOperation(ops :: Vector) = new(ops)
end

Base.:(==)(α :: CompositeAdmissibleOperation, β :: CompositeAdmissibleOperation) = α.ops == β.ops

Base.getindex(α :: CompositeAdmissibleOperation, i :: Int) = α.ops[i]

Base.one(α :: CompositeAdmissibleOperation) = CompositeAdmissibleOperation(map(one, α.ops))

Base.one(::Type{CompositeAdmissibleOperation}) = CompositeAdmissibleOperation([])

function (p :: CompositeAdmissibleOperation)(X)
    Y = deepcopy(X)
    for op in p.ops
        Y = op(Y)
    end
    return Y
end

Base.:*(p1 :: CompositeAdmissibleOperation, p2 :: CompositeAdmissibleOperation) = CompositeAdmissibleOperation([p1.ops ; p2.ops])

Base.inv(p :: CompositeAdmissibleOperation) = CompositeAdmissibleOperation(map(inv, reverse(p.ops)))

function Base.show(io :: IO, α :: CompositeAdmissibleOperation)
    if isempty(α.ops)
        println("Trivial admissible operation")
    else
        println("Do the following steps:")
        for i in axes(α.ops, 1)
            print("$i. ")
            println(α.ops[i])
        end
    end
end 

# fallback methods for admissible operations
Base.one(::Type{AdmissibleOperation}) = CompositeAdmissibleOperation([])
Base.:*(p1 :: CompositeAdmissibleOperation, p2 :: AdmissibleOperation) = CompositeAdmissibleOperation([p1.ops ; [p2]])
Base.:*(p1 :: AdmissibleOperation, p2 :: CompositeAdmissibleOperation) = CompositeAdmissibleOperation([[p1] ; p2.ops])
Base.:*(p1 :: AdmissibleOperation, p2 :: AdmissibleOperation) = CompositeAdmissibleOperation([p1, p2])


#################################################################################
# Normalization of admissible operations
#
# We provide an algorithm for turning an admissible operation into normal form.
# It consists of two phases: Swapping and merging. In the swapping phase, all
# pairs of admissible operations violating the desired order (see docstring of
# `normalize_admissible_operation`) are swapped and replaced by an equivalent
# pair in the correct order. Essentially, this is the bubble sort algorithm,
# where during swapping, the admissible operations might be replaced by
# different ones. The replacement is done by various implementations of the
# `_swap_ops` function, satisfying the following condition:
#
# If `(β',α') = _swap_ops(α, β)`, we have `typeof(α) == typeof(α') && typeof(β)
# = typeof(β')` and `α * β` has the same effect on a C-star surface as `β' * α'`.
#
# In the merging phase, admissible operations of the same type are composed
# into a single operation for each type. Additionally, admissible operations
# that are trivial are removed.
###################################################################################

_swap_ops(α :: PermutationOfBlocks, β :: InvertLastRow) = (β, α)

_swap_ops(α :: PermutationOfRays, β :: InvertLastRow) = (β, α)

_swap_ops(α :: AdmissibleRowOperation, β :: InvertLastRow) = (β, inv(α))

_swap_ops(α :: PermutationOfBlocks, β :: PermutationOfRays) = (PermutationOfRays(permuted(β.ray_perms, inv(α.block_perm))), α)

_swap_ops(α :: AdmissibleRowOperation, β :: PermutationOfRays) = (β, α)

_swap_ops(α :: AdmissibleRowOperation, β :: PermutationOfBlocks) = (β, AdmissibleRowOperation(permuted(α.factors, β.block_perm)))

_admissible_operations_ordering = Dict([(InvertLastRow, 1), (PermutationOfRays, 2), (PermutationOfBlocks, 3), (AdmissibleRowOperation, 4)])

@doc raw"""
    normalize_admissible_operation(γ :: CompositeAdmissibleOperation)

Normalize a composite admissible operation by reordering the steps and merging
multiple operations of the same kind into one. The resulting
`CompositeAdmissibleOperation` is guaranteed to consist of at most four
admissible operations in the following order:

`InvertLastRow`, `PermutationOfRays`, `PermutationOfBlocks`,
`AdmissibleRowOperation`.

If only one type of operation remains, that operation is returned instead.

"""
function normalize_admissible_operation(γ :: CompositeAdmissibleOperation)
    val(α :: AdmissibleOperation) = _admissible_operations_ordering[typeof(α)]

    # Phase 1: Swapping
    ops = deepcopy(γ.ops)
    for i = 1 : length(ops) - 1
        for j = 2 : length(ops)
            if val(ops[j-1]) > val(ops[j])
                (β, α) = _swap_ops(ops[j-1], ops[j])
                ops[j-1] = β
                ops[j] = α
            end
        end
    end

    # Phase 2: Merging
    inv_op = prod(filter(α -> α isa InvertLastRow, ops))
    permute_rays_op = prod(filter(α -> α isa PermutationOfRays, ops))
    permute_blocks = prod(filter(α -> α isa PermutationOfBlocks, ops))
    row_op = prod(filter(α -> α isa AdmissibleRowOperation, ops))

    # only keep the non-trivial admissible operations
    new_ops = filter(!isone, [inv_op, permute_rays_op, permute_blocks, row_op])

    if length(new_ops) == 1
        return first(new_ops)
    else
        return CompositeAdmissibleOperation(new_ops)
    end

end

# fallback method
normalize_admissible_operation(γ :: AdmissibleOperation) = normalize_admissible_operation(CompositeAdmissibleOperation([γ]))

 
