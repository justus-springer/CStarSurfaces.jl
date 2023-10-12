
#####################################################
# Julia type for an SQLite connection to a database
# of surfaces with torus action
#####################################################

const SQLiteAdapterSurfaces = SQLiteAdapter{SurfaceWithTorusAction}

##################################################
# Default functions to compute the columns in the 
# SQLite database
##################################################

_db_is_toric(X :: SurfaceWithTorusAction) = is_toric(X) ? 1 : 0

function _db_gen_matrix(X :: SurfaceWithTorusAction)
    P = gen_matrix(X)
    return string([[Int(P[i,j]) for j = 1 : ncols(P)] for i = 1 : nrows(P)])
end

_db_rays(X :: SurfaceWithTorusAction) = "[" * join(["[" * join(r, ", ") * "]" for r in rays(X)], ", ") * "]"

_db_nrays(X :: SurfaceWithTorusAction) = nrays(X)

_db_lss(X :: ToricSurface) = missing
_db_lss(X :: CStarSurface) = string(Vector(X.l))

_db_dss(X :: ToricSurface) = missing
_db_dss(X :: CStarSurface) = string(Vector(X.d))

_db_case_(X :: ToricSurface) = missing
_db_case_(X :: CStarSurface) = string(X.case)

_db_block_sizes(X :: ToricSurface) = missing
_db_block_sizes(X :: CStarSurface) = string(Vector(block_sizes(X)))

_db_nblocks(X :: ToricSurface) = missing
_db_nblocks(X :: CStarSurface) = nblocks(X)

_db_number_of_parabolic_fixed_point_curves(X :: ToricSurface) = missing
_db_number_of_parabolic_fixed_point_curves(X :: CStarSurface) = number_of_parabolic_fixed_point_curves(X)

_db_orientation(X :: ToricSurface) = missing
_db_orientation(X :: CStarSurface) = orientation(X)

_db_is_intrinsic_quadric(X :: ToricSurface) = missing
_db_is_intrinsic_quadric(X :: CStarSurface) = is_intrinsic_quadric(X) ? 1 : 0

_db_class_group_rank(X :: SurfaceWithTorusAction) = class_group_rank(X)

_db_class_group_torsion(X :: SurfaceWithTorusAction) = "[" * join(class_group_torsion(X), ", ") * "]"

_db_class_group_torsion_order(X :: SurfaceWithTorusAction) = Int(class_group_torsion_order(X))

function _db_degree_matrix(X :: SurfaceWithTorusAction)
  Q = degree_matrix(X)
  return string([[Int(Q[i,j]) for j = 1 : ncols(Q)] for i = 1 : nrows(Q)])
end

function _db_canonical_divisor_class(X :: SurfaceWithTorusAction)
  c = divisor_class(canonical_divisor_class(X)).coeff
  return string([Int(c[1,i]) for i = 1 : ncols(c)])
end

_db_gorenstein_index(X :: SurfaceWithTorusAction) = Int(gorenstein_index(X))

_db_picard_index(X :: SurfaceWithTorusAction) = Int(picard_index(X))

_db_maximal_log_canonicity_numerator(X :: SurfaceWithTorusAction) = Int(numerator(maximal_log_canonicity(X)))

_db_maximal_log_canonicity_denominator(X :: SurfaceWithTorusAction) = Int(denominator(maximal_log_canonicity(X)))

_db_anticanonical_self_intersection_numerator(X :: SurfaceWithTorusAction) = Int(numerator(anticanonical_self_intersection(X)))

_db_anticanonical_self_intersection_denominator(X :: SurfaceWithTorusAction) = Int(denominator(anticanonical_self_intersection(X)))

default_column_functions(::Type{<:SurfaceWithTorusAction}) = Dict([
:is_toric => _db_is_toric,
:gen_matrix => _db_gen_matrix,
:rays => _db_rays,
:nrays => _db_nrays,
:lss => _db_lss,
:dss => _db_dss,
:case_ => _db_case_,
:block_sizes => _db_block_sizes,
:nblocks => _db_nblocks,
:number_of_parabolic_fixed_point_curves => _db_number_of_parabolic_fixed_point_curves,
:orientation => _db_orientation,
:is_intrinsic_quadric => _db_is_intrinsic_quadric,
:class_group_rank => _db_class_group_rank,
:class_group_torsion => _db_class_group_torsion,
:class_group_torsion_order => _db_class_group_torsion_order,
:degree_matrix => _db_degree_matrix,
:canonical_divisor_class => _db_canonical_divisor_class,
:gorenstein_index => _db_gorenstein_index,
:picard_index => _db_picard_index,
:maximal_log_canonicity_numerator => _db_maximal_log_canonicity_numerator,
:maximal_log_canonicity_denominator => _db_maximal_log_canonicity_denominator,
:anticanonical_self_intersection_numerator => _db_anticanonical_self_intersection_numerator,
:anticanonical_self_intersection_denominator => _db_anticanonical_self_intersection_denominator,
])

######################################################################
# Find a surface in a database.
# Returns nothing when the surface is not present.
######################################################################

function find_in_database(db :: SQLiteAdapterSurfaces, X :: SurfaceWithTorusAction)
    X = is_toric(X) ? normal_form(X) : normal_form(X)[1]
    is_toric_str, gen_matrix_str = _db_is_toric(X), _db_gen_matrix(X)
    sql = "SELECT $(db.primary_key) FROM $(db.table_name) WHERE \
             is_toric == $is_toric_str AND \
             gen_matrix == \"$gen_matrix_str\""
    res = DBInterface.execute(db.db, sql) |> rowtable
    isempty(res) && return nothing
    return res[1][Symbol(db.primary_key)]
end

######################################################################
# Import a single surface from an SQLite row
######################################################################

function sqlite_import_row(::Type{SurfaceWithTorusAction}, row :: Union{SQLite.Row, NamedTuple})
    P = matrix(ZZ, eval(Meta.parse(row[:gen_matrix])))
    return row[:is_toric] == 1 ? toric_surface(P) : cstar_surface(P)
end








