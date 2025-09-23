Base.Enums.@enum CStarSurfaceCase EE PE EP PP

has_elliptic_fixed_point_plus(C :: CStarSurfaceCase) =
C == EE || C == EP

has_elliptic_fixed_point_minus(C :: CStarSurfaceCase) =
C == EE || C == PE

has_parabolic_fixed_point_curve_plus(C :: CStarSurfaceCase) =
C == PP || C == PE

has_parabolic_fixed_point_curve_minus(C :: CStarSurfaceCase) =
C == PP || C == EP

function swap_case(C :: CStarSurfaceCase, swap :: Bool = true)
    !swap && return C
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
