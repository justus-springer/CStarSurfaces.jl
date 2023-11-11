abstract type SingularityType end

abstract type SingularityTypeADE <: SingularityType end

struct SingularityTypeA <: SingularityTypeADE
    n :: Integer
end

Base.show(io :: IO, st :: SingularityTypeA) = print(io, "A$(st.n)")

struct SingularityTypeD <: SingularityTypeADE
    n :: Integer
end

Base.show(io :: IO, st :: SingularityTypeD) = print(io, "D$(st.n)")

struct SingularityTypeE6 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE6) = print(io, "E6")

struct SingularityTypeE7 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE7) = print(io, "E7")

struct SingularityTypeE8 <: SingularityTypeADE end

Base.show(io :: IO, st :: SingularityTypeE8) = print(io, "E8")

struct SingularityTypeNonLogTerminal <: SingularityType end

Base.show(io :: IO, st :: SingularityTypeNonLogTerminal) = print(io, "Non log terminal singularity")
