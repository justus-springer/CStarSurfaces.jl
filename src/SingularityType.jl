@doc raw"""
    SingularityType

Abstract supertype of various "types" of singularities.

"""
abstract type SingularityType end

@doc raw"""
    SingularityTypeADE <: SingularityType

Abstract supertype of ADE singularities, i.e. log terminal ones.

"""
abstract type SingularityTypeADE <: SingularityType end


@doc raw"""
    SingularityTypeA <: SingularityTypeADE

The singularity type $A_n$. It has a single field `n` that holds the number
of nodes in the resolution graph.

"""
struct SingularityTypeA <: SingularityTypeADE
    n :: Integer
end

Base.show(io :: IO, st :: SingularityTypeA) = print(io, "A$(st.n)")


@doc raw"""
    SingularityTypeD <: SingularityTypeADE

The singularity type $D_n$. It has a single field `n` that holds the number
of nodes in the resolution graph.

"""
struct SingularityTypeD <: SingularityTypeADE
    n :: Integer
end

Base.show(io :: IO, st :: SingularityTypeD) = print(io, "D$(st.n)")


@doc raw"""
    SingularityTypeE6 <: SingularityTypeADE

The singularity type $E_6$.

"""
struct SingularityTypeE6 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE6) = print(io, "E6")


@doc raw"""
    SingularityTypeE7 <: SingularityTypeADE

The singularity type $E_7$.

"""
struct SingularityTypeE7 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE7) = print(io, "E7")


@doc raw"""
    SingularityTypeE8 <: SingularityTypeADE

The singularity type $E_8$.

"""
struct SingularityTypeE8 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE8) = print(io, "E8")


@doc raw"""
    SingularityTypeNonLogTerminal <: SingularityType

The singularity type of non log terminal singularities.

"""
struct SingularityTypeNonLogTerminal <: SingularityType end

Base.show(io :: IO, st :: SingularityTypeNonLogTerminal) = print(io, "Non log terminal singularity")
