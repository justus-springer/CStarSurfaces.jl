@doc raw"""
    AbelianGroup{T <: Integer}

A struct reprensting a finitely generated abelian group. It consists of two
fields:

- `rank :: T`: The rank of the group.
- `elementary_divisors :: T`: The elementary divisors, i.e. the orders of the
  torsion part.

"""
struct AbelianGroup{T <: Integer}
    rank :: T
    elementary_divisors :: Vector{T}
end


function Base.show(io :: IO, G :: AbelianGroup{T}) where {T <: Integer}
    r, ds = G.rank, G.elementary_divisors
    free_string = r == 1 ? "Z" : "Z^$r"
    torsion_string = join(["Z/$d" for d in G.elementary_divisors], " x ")
    if r == 0
        if isempty(ds)
            print(io, "Z/1")
        else
            print(io, torsion_string)
        end
    else
        if isempty(ds)
            print(io, free_string)
        else
            print(io, free_string * " x " * torsion_string)
        end
    end
end


@doc raw"""
    rank(G :: AbelianGroup)

The rank of an abelian group.

"""
rank(G :: AbelianGroup) = G.rank


@doc raw"""
    elementary_divisors(G :: AbelianGroup)

The elementary divisors of the abelian group, i.e. the orders of the torsion
part.

"""
elementary_divisors(G :: AbelianGroup) = G.elementary_divisors


@doc raw"""
    torsion_order(G :: AbelianGroup)

The order of the torsion part of the abelian group. This equals the product
of its elementary divisors.

"""
torsion_order(G :: AbelianGroup) = prod(elementary_divisors(G))


@doc raw"""
    torsion_part(G :: AbelianGroup)

The torsion part of an abelian group.

"""
torsion_part(G :: AbelianGroup) = AbelianGroup(0, elementary_divisors(G))


@doc raw"""
    cokernel(A :: AbstractMatrix)

The cokernel of a matrix as an abelian group.

"""
function cokernel(A :: AbstractMatrix)
    r = size(A,1) - size(A,2)
    ds = elementary_divisors(A)
    filter!(d -> d > 1, ds)
    return AbelianGroup(r, ds)
end
