
@doc raw"""
    SQLiteAdapterSurfaces = SQLiteAdapter{SurfaceWithTorusAction}

An adapter to an SQLite database holding objects of type
`SurfaceWithTorusAction`.

"""
const SQLiteAdapterSurfaces = SQLiteAdapter{SurfaceWithTorusAction}



# The list of column names together with their SQLite definition. The third
# entry in each tuple says if the column is explicitly exported. It is false
# for columns that are automatically computed from other columns, such as
# `log_canonicity`. The fourth entry says whether the column is used to
# cache attributes when imported.
#
# If the third entry is `true`, there should be a function named 
# `_db_<column_name>` specifying how the column is exported (see below).
#
# If the fourth entry is `true`, there should be a function named
# `_import_db_<column_name>` specifying how to import the column (see below).
const _db_column_defs = [
    (:is_toric, "INTEGER NOT NULL", true, false),
    (:gen_matrix, "TEXT NOT NULL", true, false),
    (:rays, "TEXT", true, false),
    (:nrays, "INTEGER", true, false),
    (:lss, "TEXT", true, false),
    (:dss, "TEXT", true, false),
    (:case_, "TEXT", true, false),
    (:block_sizes, "TEXT", true, false),
    (:nblocks, "INTEGER", true, false),
    (:number_of_parabolic_fixed_point_curves, "INTEGER", true, false),
    (:orientation, "INTEGER", true, false),
    (:is_intrinsic_quadric, "INTEGER", true, false),
    (:class_group_rank, "INTEGER", true, true),
    (:class_group_torsion, "TEXT", true, true),
    (:class_group_torsion_order, "INTEGER", true, true),
    (:degree_matrix, "TEXT", true, false),
    (:canonical_divisor_class, "TEXT", true, false),
    (:gorenstein_index, "INTEGER", true, true),
    (:picard_index, "INTEGER", true, true),
    (:log_canonicity_numerator, "INTEGER", true, false),
    (:log_canonicity_denominator, "INTEGER", true, false),
    (:log_canonicity, "REAL AS (CAST(log_canonicity_numerator AS FLOAT) / CAST(log_canonicity_denominator AS FLOAT))", false, true),
    (:anticanonical_self_intersection_numerator, "INTEGER", true, false),
    (:anticanonical_self_intersection_denominator, "INTEGER", true, false),
    (:anticanonical_self_intersection, "REAL AS (CAST(anticanonical_self_intersection_numerator AS FLOAT) / CAST(anticanonical_self_intersection_denominator AS FLOAT))", false, true),
    (:admits_kaehler_ricci_soliton, "INTEGER", false, false),
    (:admits_kaehler_einstein_metric, "INTEGER", true, true),
    (:admits_sasaki_einstein_metric, "INTEGER", false, false),
    (:is_quasismooth, "INTEGER", true, true),
    (:is_factorial, "INTEGER", true, true),
    (:is_smooth, "INTEGER", true, true),
    (:number_of_singularities, "INTEGER", true, true),
    (:singularity_types_string, "TEXT", true, true),
    (:number_of_exceptional_prime_divisors, "INTEGER", true, true),
    (:singularity_kind_x_plus, "TEXT", true, false),
    (:number_of_exceptional_prime_divisors_x_plus, "INTEGER", true, false),
    (:singularity_kind_x_minus, "TEXT", true, false),
    (:number_of_exceptional_prime_divisors_x_minus, "INTEGER", true, false),
    (:is_combinatorially_minimal, "INTEGER", true, true)
]

@doc raw"""
    create_table(db :: SQLiteAdapterSurfaces; temp = false, ifnotexists = true)

Create a new SQLite table with default column definitions holding
`SurfaceWithTorusAction`s.

# Example

Create a new SQLite database for holding `SurfaceWithTorusAction`s. The
resulting table can be inspected with
[`SQLite.jl`](https://juliadatabases.org/SQLite.jl/stable/):

```jldoctest
julia> db = SQLiteAdapterSurfaces("my_database.db", "surfaces", "surface_id")
SQLiteAdapterSurfaces(SQLite.DB("my_database.db"), "my_database.db", "surfaces", "surface_id")

julia> create_table(db);

julia> using SQLite

julia> SQLite.tables(db.db)
1-element Vector{SQLite.DBTable}:
 SQLite.DBTable("surfaces", Tables.Schema:
 :surface_id                                   Union{Missing, Int64}
 :is_toric                                     Union{Missing, Int64}
 :gen_matrix                                   Union{Missing, String}
 :rays                                         Union{Missing, String}
 :nrays                                        Union{Missing, Int64}
 :lss                                          Union{Missing, String}
 :dss                                          Union{Missing, String}
 :case_                                        Union{Missing, String}
 :block_sizes                                  Union{Missing, String}
 :nblocks                                      Union{Missing, Int64}
 :number_of_parabolic_fixed_point_curves       Union{Missing, Int64}
 :orientation                                  Union{Missing, Int64}
 :is_intrinsic_quadric                         Union{Missing, Int64}
 :class_group_rank                             Union{Missing, Int64}
 :class_group_torsion                          Union{Missing, String}
 :class_group_torsion_order                    Union{Missing, Int64}
 :degree_matrix                                Union{Missing, String}
 :canonical_divisor_class                      Union{Missing, String}
 :gorenstein_index                             Union{Missing, Int64}
 :picard_index                                 Union{Missing, Int64}
 :log_canonicity_numerator                     Union{Missing, Int64}
 :log_canonicity_denominator                   Union{Missing, Int64}
 :log_canonicity                               Union{Missing, Float64}
 :anticanonical_self_intersection_numerator    Union{Missing, Int64}
 :anticanonical_self_intersection_denominator  Union{Missing, Int64}
 :anticanonical_self_intersection              Union{Missing, Float64}
 :admits_kaehler_ricci_soliton                 Union{Missing, Int64}
 :admits_kaehler_einstein_metric               Union{Missing, Int64}
 :admits_sasaki_einstein_metric                Union{Missing, Int64}
 :is_quasismooth                               Union{Missing, Int64}
 :is_factorial                                 Union{Missing, Int64}
 :is_smooth                                    Union{Missing, Int64}
 :number_of_singularities                      Union{Missing, Int64}
 :singularity_types_string                      Union{Missing, String}
 :number_of_exceptional_prime_divisors          Union{Missing, Int64}
 :singularity_kind_x_plus                       Union{Missing, String}
 :number_of_exceptional_prime_divisors_x_plus   Union{Missing, Int64}
 :singularity_kind_x_minus                      Union{Missing, String}
 :number_of_exceptional_prime_divisors_x_minus  Union{Missing, Int64}
 :is_combinatorially_minimal                    Union{Missing, Int64})

```

"""
function create_table(db :: SQLiteAdapterSurfaces; temp = false, ifnotexists = true)
    temp = temp ? "TEMP" : ""
    ifnotexists = ifnotexists ? "IF NOT EXISTS" : ""
    columns = [string(db.primary_key, ' ', "INTEGER PRIMARY KEY ASC")]
    for (column_name, column_def, _, _) in _db_column_defs
        push!(columns, string(column_name, ' ', column_def))
    end
    sql = "CREATE $temp TABLE $ifnotexists $(db.table_name) ($(join(columns, ',')))"
    DBInterface.execute(db.db, sql)
end



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
  c = canonical_divisor_class(X).coeff
  return string([Int(c[1,i]) for i = 1 : ncols(c)])
end

_db_gorenstein_index(X :: SurfaceWithTorusAction) = Int(gorenstein_index(X))

_db_picard_index(X :: SurfaceWithTorusAction) = Int(picard_index(X))

_db_log_canonicity_numerator(X :: SurfaceWithTorusAction) = Int(numerator(log_canonicity(X)))

_db_log_canonicity_denominator(X :: SurfaceWithTorusAction) = Int(denominator(log_canonicity(X)))

_db_anticanonical_self_intersection_numerator(X :: SurfaceWithTorusAction) = Int(numerator(anticanonical_self_intersection(X)))

_db_anticanonical_self_intersection_denominator(X :: SurfaceWithTorusAction) = Int(denominator(anticanonical_self_intersection(X)))

_db_admits_kaehler_einstein_metric(X :: SurfaceWithTorusAction) = 
admits_kaehler_einstein_metric(X) ? 1 : 0 

_db_is_quasismooth(X :: SurfaceWithTorusAction) = 
is_quasismooth(X) ? 1 : 0

_db_is_factorial(X :: SurfaceWithTorusAction) = 
is_factorial(X) ? 1 : 0

_db_is_smooth(X :: SurfaceWithTorusAction) = 
is_smooth(X) ? 1 : 0

_db_number_of_singularities(X :: SurfaceWithTorusAction) = 
number_of_singularities(X)

_db_singularity_types_string(X :: SurfaceWithTorusAction) =
singularity_types_string(X)

_db_number_of_exceptional_prime_divisors(X :: SurfaceWithTorusAction) =
number_of_exceptional_prime_divisors(X)

_db_singularity_kind_x_plus(X :: CStarSurface{<:Union{EE,EP}}) =
"$(singularity_kind(x_plus(X)))"
_db_singularity_kind_x_plus(X :: CStarSurface{<:Union{PE,PP}}) = missing
_db_singularity_kind_x_plus(X :: ToricSurface) = missing

_db_number_of_exceptional_prime_divisors_x_plus(X :: CStarSurface{<:Union{EE,EP}}) =
number_of_exceptional_prime_divisors(x_plus(X))
_db_number_of_exceptional_prime_divisors_x_plus(X :: CStarSurface{<:Union{PE,PP}}) = missing
_db_number_of_exceptional_prime_divisors_x_plus(X :: ToricSurface) = missing

_db_singularity_kind_x_minus(X :: CStarSurface{<:Union{EE,PE}}) =
"$(singularity_kind(x_minus(X)))"
_db_singularity_kind_x_minus(X :: CStarSurface{<:Union{EP,PP}}) = missing
_db_singularity_kind_x_minus(X :: ToricSurface) = missing

_db_number_of_exceptional_prime_divisors_x_minus(X :: CStarSurface{<:Union{EE,PE}}) =
number_of_exceptional_prime_divisors(x_minus(X))
_db_number_of_exceptional_prime_divisors_x_minus(X :: CStarSurface{<:Union{EP,PP}}) = missing
_db_number_of_exceptional_prime_divisors_x_minus(X :: ToricSurface) = missing

_db_is_combinatorially_minimal(X :: SurfaceWithTorusAction) = 
is_combinatorially_minimal(X) ? 1 : 0



@doc raw"""
    default_column_functions(::Type{<:SurfaceWithTorusAction})   

The default columns names and how to compute them when exporting objects
of type `SurfaceWithTorusAction` to an SQLite database. 

The function names all have the prefix `_db_` followed by the name of the
column, for instance `CStarSurfaces._db_gen_matrix` (they are not exported by
default). They basically wrap the corresponding attribute function, giving the
result as a language-agnostic string for database storage instead of a Julia
type.

"""
default_column_functions(::Type{<:SurfaceWithTorusAction}) = 
Dict([column_name => eval(Symbol(:_db_, column_name)) 
      for (column_name, _, exported, _) in _db_column_defs if exported])


@doc raw"""
    find_in_database(db :: SQLiteAdapterSurfaces, X :: SurfaceWithTorusAction)

Tries to find a `SurfaceWithTorusAction` in an SQLite database. If the surface
is in the database, this function returns its `id`. Otherwise, it returns
`nothing`.

"""
function find_in_database(db :: SQLiteAdapterSurfaces, X :: SurfaceWithTorusAction)
    X = is_toric(X) ? normal_form(X) : normal_form(X)[1]
    is_toric_str, gen_matrix_str = _db_is_toric(X), _db_gen_matrix(X)
    sql = "SELECT $(db.primary_key) FROM $(db.table_name) WHERE is_toric == $is_toric_str AND gen_matrix == \"$(gen_matrix_str)\""
    res = DBInterface.execute(db.db, sql) |> rowtable
    isempty(res) && return nothing
    return res[1][Symbol(db.primary_key)]
end


##################################################
# Importer functions for attribute caching
##################################################

_import_db_class_group_rank(row :: Union{SQLite.Row, NamedTuple}) =
row[:class_group_rank]

_import_db_class_group_torsion(row :: Union{SQLite.Row, NamedTuple}) =
eval(Meta.parse(row[:class_group_torsion]))

_import_db_class_group_torsion_order(row :: Union{SQLite.Row, NamedTuple}) =
row[:class_group_torsion_order]

_import_db_gorenstein_index(row :: Union{SQLite.Row, NamedTuple}) =
row[:gorenstein_index]

_import_db_picard_index(row :: Union{SQLite.Row, NamedTuple}) =
row[:picard_index]

_import_db_log_canonicity(row :: Union{SQLite.Row, NamedTuple}) =
row[:log_canonicity_numerator] // row[:log_canonicity_denominator]

_import_db_anticanonical_self_intersection(row :: Union{SQLite.Row, NamedTuple}) =
row[:anticanonical_self_intersection_numerator] // row[:anticanonical_self_intersection_denominator]

_import_db_admits_kaehler_einstein_metric(row :: Union{SQLite.Row, NamedTuple}) =
row[:admits_kaehler_einstein_metric] == 1

_import_db_is_quasismooth(row :: Union{SQLite.Row, NamedTuple}) =
row[:is_quasismooth] == 1

_import_db_is_factorial(row :: Union{SQLite.Row, NamedTuple}) =
row[:is_factorial] == 1

_import_db_is_smooth(row :: Union{SQLite.Row, NamedTuple}) =
row[:is_smooth] == 1

_import_db_number_of_singularities(row :: Union{SQLite.Row, NamedTuple}) =
row[:number_of_singularities]

_import_db_singularity_types_string(row :: Union{SQLite.Row, NamedTuple}) =
row[:singularity_types_string]

_import_db_number_of_exceptional_prime_divisors(row :: Union{SQLite.Row, NamedTuple}) =
row[:number_of_exceptional_prime_divisors]

_import_db_is_combinatorially_minimal(row :: Union{SQLite.Row, NamedTuple}) =
row[:is_combinatorially_minimal] == 1


@doc raw"""
    sqlite_import_row(::Type{SurfaceWithTorusAction}, row :: Union{SQLite.Row, NamedTuple})

Imports a single row from an `SQLiteAdapterSurfaces` into a
`SurfaceWithTorusAction`.

"""
function sqlite_import_row(::Type{SurfaceWithTorusAction}, row :: Union{SQLite.Row, NamedTuple}; import_attributes = false)
    P = matrix(ZZ, eval(Meta.parse(row[:gen_matrix])))
    X = row[:is_toric] == 1 ? toric_surface(P) : cstar_surface(P)

    if import_attributes
        attribute_names = [attr for (attr, _, _, imported) in _db_column_defs if imported]
        for attr_name in attribute_names
            f = eval(Symbol(:_import_db_, attr_name))
            set_attribute!(X, attr_name, f(row))
        end
    end

    return X
end
