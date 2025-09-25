@doc raw"""
    CStarSurfaceCase

An enum type with four values: `EE`, `PE`, `EP`, `PP`.

"""
Base.Enums.@enum CStarSurfaceCase EE PE EP PP


@doc raw"""
    has_elliptic_fixed_point_plus(C :: CStarSurfaceCase)

Return true if `C` is `EE` or `EP`, false otherwise.

"""
has_elliptic_fixed_point_plus(C :: CStarSurfaceCase) =
C == EE || C == EP


@doc raw"""
    has_elliptic_fixed_point_minus(C :: CStarSurfaceCase)

Return true if `C` is `EE` or `PE`, false otherwise.

"""
has_elliptic_fixed_point_minus(C :: CStarSurfaceCase) =
C == EE || C == PE


@doc raw"""
    has_parabolic_fixed_point_curve_plus(C :: CStarSurfaceCase)

Return true if `C` is `PP` or `PE`, false otherwise.

"""
has_parabolic_fixed_point_curve_plus(C :: CStarSurfaceCase) =
C == PP || C == PE


@doc raw"""
    has_parabolic_fixed_point_curve_minus(C :: CStarSurfaceCase)

Return true if `C` is `PP` or `EP`, false otherwise.

"""
has_parabolic_fixed_point_curve_minus(C :: CStarSurfaceCase) =
C == PP || C == EP


@doc raw"""
    invert_case(C :: CStarSurfaceCase, invert :: Bool = true)

Invert a case, see Definition ``\ref{def:admissible_operations}``.
It sends `EE` to `EE`, `PE` to `EP`, `EP` to `PE` and `PP` to `PP`.
If `invert` is set to false, it just returns the given case.

"""
function invert_case(C :: CStarSurfaceCase, invert :: Bool = true)
    !invert && return C
    C == EE && return EE
    C == PE && return EP
    C == EP && return PE
    C == PP && return PP
end

function parse_case(str :: AbstractString)
    str == "EE" && return EE
    str == "EP" && return EP
    str == "PE" && return PE
    str == "PP" && return PP
end
