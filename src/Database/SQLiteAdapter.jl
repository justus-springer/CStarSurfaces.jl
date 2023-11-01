@doc raw"""
    SQLiteAdapter{T} <: DatabaseAdapter{T}

An adapter to an SQLite database holding objects of type `T` where `T <:
MoriDreamSpace`. The type `T` should at least implement the following methods
(see their docstrings for more information):

`default_column_functions`, `find_in_database`, `sqlite_import_row`.

"""
struct SQLiteAdapter{T} <: DatabaseAdapter{T}
    db :: SQLite.DB
    file_path :: AbstractString
    table_name :: AbstractString
    primary_key :: AbstractString

    SQLiteAdapter{T}(f :: AbstractString, table_name :: AbstractString, primary_key :: AbstractString) where {T <: MoriDreamSpace} = new{T}(SQLite.DB(f), f, table_name, primary_key)
end


@doc raw"""
    default_column_functions(::Type{T}) where {T <: MoriDreamSpace}

Returns a `Dict{Symbol, Function}` that serves as a default for the 
names of the columns to export and how to export them for a given 
subtype of `MoriDreamSpace`. Should be implemented by all subtypes
of `MoriDreamSpace` where database functionality is desired.

The fallback definition for a general `T` returns an empty dictionary.

"""
default_column_functions(::Type{T}) where {T <: MoriDreamSpace} = 
Dict{Symbol, Function}([])


@doc raw"""
    default_insert_predicate(::Type{T}) where {T <: MoriDreamSpace}

Returns a function of type `(db :: SQLiteAdapter{T}, X :: T) -> Bool` that
serves as the default insert predicate when exporting objects of type `T` to an
SQLite database. The default implementation returns `true` if and only if
`find_in_database(db, X)` returns nothing, hence avoiding duplicate entries in
the database. Note that there is no default implementation for
`find_in_database`, it has to be implemented for each subtype of
`MoriDreamSpace` where database functionality is desired.

"""
default_insert_predicate(::Type{T}) where {T <: MoriDreamSpace} =
function(db :: SQLiteAdapter{T}, X :: T)
    return isnothing(find_in_database(db, X))
end


@doc raw"""
    export_to_database(db_adapter :: SQLiteAdapter{T}, Xs :: AbstractVector; kwargs...) where {T <: MoriDreamSpace}

Export a list `Xs` of varieties of type `T` to an SQLite database. The
following keyword arguments are supported:

- `column_functions`: Defaults to `default_column_functions(T)`. A dictionary of
  type `Dict{Symbol, Function}` containing for each exported column name a function
  that tells `export_to_database` how to compute that column. Each function takes
  an object of type `T` and returns an SQLite compatible data type (`Int`, `Float`, `String` or `Nothing`).
- `insert_predicate`: Defaults to `default_insert_predicate(T)`. A function of 
  type `(db :: SQLiteAdapter{T}, X :: T) -> Bool`. Only objects where this
  predicate evaluates to `true` are exported.

"""
function export_to_database(
        db_adapter :: SQLiteAdapter{T}, 
        Xs :: AbstractVector;
        column_functions = default_column_functions(T),
        insert_predicate = default_insert_predicate(T)) where {T <: MoriDreamSpace}

    @info "Starting export to SQLite database file $(db_adapter.file_path)."
    @info "Number of objects: $(length(Xs))."

    db = db_adapter.db
    table_name = db_adapter.table_name
    if isempty(column_functions)
        @info "No column functions provided. Inserting empty rows."
        stmt = SQLite.Stmt(db, "INSERT INTO $table_name DEFAULT VALUES")
    else
        @info "Provided column functions: " keys(column_functions)
        columns = join(keys(column_functions), ",")
        values = join(":" .* string.(keys(column_functions)), ",")
        stmt = SQLite.Stmt(db, "INSERT INTO $table_name ($columns) VALUES ($values)")
    end

    @info """
    Starting to insert objects into database. 
    Use a `TerminalLogger` from `TerminalLoggers.jl` to see a nice progress bar.
    """

    skip_count = 0
    @progress for i in axes(Xs, 1)
        X = Xs[i]

        if !insert_predicate(db_adapter, X)
            skip_count += 1
            continue
        end

        val_dict = Dict([k => f(X) for (k,f) in column_functions])
        DBInterface.execute(stmt, val_dict)

    end

    inserted_count = length(Xs) - skip_count
    @info """
    Inserted $inserted_count objects.
    Skipped $skip_count objects.
    """
    
end


@doc raw"""
    import_from_database(db :: SQLiteAdapter{T}, sql :: String = "TRUE") where {T <: MoriDreamSpace}

Return a list of objects of type `T` that match a given sql query. The string
`sql` is used after the WHERE in a SELECT statement, hence can contain
restrictions on the columns as well as an ORDER BY and a LIMIT clause.

"""
function import_from_database(db :: SQLiteAdapter{T}, sql :: String = "TRUE") where {T <: MoriDreamSpace}

    @info "Importing from SQLite database file $(db.file_path)..."

    select_stmt = SQLite.Stmt(db.db, "SELECT * FROM $(db.table_name) WHERE $sql")
    rows = DBInterface.execute(select_stmt) |> rowtable
    count = length(rows)

    @info "Number of imported objects: $count."

    @info "Converting to Julia type `$T`."

    Xs = T[]
    @progress for i = 1 : count
        push!(Xs, sqlite_import_row(T, rows[i]))
    end
    return Xs
end


@doc raw"""
    import_from_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}) where {T <: MoriDreamSpace}   

Return the list of objects of type `T` from an SQLite database with the given
`ids`.

"""
import_from_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}) where {T <: MoriDreamSpace} = 
import_from_database(db, "$(db.primary_key) IN ($(join(ids, ", ")))")


@doc raw"""
    import_from_database(db :: SQLiteAdapter{T}, id :: Int) where {T <: MoriDreamSpace}

Return the object of type `T` from an SQLite database with the given `id`.

"""
import_from_database(db :: SQLiteAdapter{T}, id :: Int) where {T <: MoriDreamSpace} = first(import_from_database(db, [id]))


_argsym_to_arg(T :: Type{<:MoriDreamSpace}, row :: Union{SQLite.Row, NamedTuple}, argsym :: Symbol) = 
argsym == :variety ? sqlite_import_row(T, row) : row[argsym]

@doc raw"""
    update_in_database(db :: SQLiteAdapter{T}, column_functions :: Dict{Symbol, <:Function}; sql :: String, column_function_args :: AbstractVector{Symbol}) where {T <: MoriDreamSpace}

Update all rows in an SQLite database matching a given SQL query by recomputing
the columns given by `column_functions`. See also `export_to_database` and
`import_from_database`.

Keyword arguments:

- `sql :: String`, defaults to `"TRUE"`. The SQL expression used for filtering the 
  rows in the database that should be updated. It is inserted into a SELECT statement
  after the WHERE, hence can contain restrictions on the columns as well as 
  an ORDER BY and a LIMIT clause.
- `column_function_args :: AbstractVector{String}`, defaults to `[:variety]`. The 
  list of arguments that each column function receives. These can either be names of
  columns in the database, or the special symbol `:variety`, in which case the 
  variety corresponding to the row is imported using `sqlite_import_row` and then 
  passed as an argument.

"""
function update_in_database(
        db :: SQLiteAdapter{T}, 
        column_functions :: Dict{Symbol, <:Function};
        sql :: String = "TRUE",
        column_function_args :: AbstractVector{Symbol} = [:variety]) where {T <: MoriDreamSpace}

    table_name, primary_key = db.table_name, db.primary_key

    @info "Updating columns in SQLite database file $(db.file_path)."
    @info "Provided column functions: " keys(column_functions)
    @info "Importing..."

    select_stmt = SQLite.Stmt(db.db, "SELECT * FROM $(db.table_name) WHERE $sql")
    rows = DBInterface.execute(select_stmt) |> rowtable
    count = length(rows)

    @info "Number of affected rows: $count."

    @info "Updating columns..."

    set_equations = join(["$k = :$k" for k in keys(column_functions)], ", ")
    update_stmt = SQLite.Stmt(db.db, "UPDATE $table_name SET $set_equations WHERE $primary_key = :$primary_key")

    @progress for i = 1 : count
        row = rows[i]
        arglist = [_argsym_to_arg(T, row, argsym) for argsym in column_function_args]
        val_dict = Dict{Symbol,Any}([k => f(arglist...) for (k,f) in column_functions])
        push!(val_dict, Symbol(primary_key) => row[Symbol(primary_key)])
        DBInterface.execute(update_stmt, val_dict)
    end

end

@doc raw"""
    update_in_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}, column_functions :: Dict{Symbol, <:Function}; column_function_args :: AbstractVector{Symbol}) where {T <: MoriDreamSpace})

Update all rows in an SQLite database with the given `ids` by recomputing the
columns given by `column_functions`.

"""
update_in_database(
        db :: SQLiteAdapter{T}, 
        column_functions :: Dict{Symbol, <:Function},
        ids :: AbstractVector{Int}; 
        column_function_args :: AbstractVector{Symbol} = [:variety]) where {T <: MoriDreamSpace} =
update_in_database(db, column_functions; sql = "$(db.primary_key) IN ($(join(ids, ", ")))", column_function_args = column_function_args)


@doc raw"""
    update_in_database(db :: SQLiteAdapter{T}, id :: Int, column_functions :: Dict{Symbol, <:Function}, column_function_args :: AbstractVector{Symbol}) where {T <: MoriDreamSpace}))

Update all rows in an SQLite database with the given `id` by recomputing the
columns given by `column_functions`.

"""
update_in_database(
        db :: SQLiteAdapter{T}, 
        id :: Int, 
        column_functions :: Dict{Symbol, <:Function};
        column_function_args :: AbstractVector{Symbol} = [:variety]) where {T <: MoriDreamSpace} =
update_in_database(db, column_functions, [id]; column_function_args = column_function_args)
