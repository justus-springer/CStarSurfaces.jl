struct SQLiteAdapter{T} <: DatabaseAdapter{T}
    db :: SQLite.DB
    file_path :: AbstractString
    table_name :: AbstractString
    primary_key :: AbstractString

    SQLiteAdapter{T}(f :: AbstractString, table_name :: AbstractString, primary_key :: AbstractString) where {T <: MoriDreamSpace} = new{T}(SQLite.DB(f), f, table_name, primary_key)
end

######################################################################
# Exporting to a SQLite database
######################################################################

# Fallback definition. Subtypes of `MoriDreamSpace` should implement this,
# describing the default column names the functions to compute them
# when inserting into a database..
default_column_functions(::Type{T}) where {T <: MoriDreamSpace} = 
Dict{Symbol, Function}([])

default_insert_predicate(::Type{T}) where {T <: MoriDreamSpace} =
function(X...) true end

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

        if !insert_predicate(X)
            @info "Skipping object no. $i"
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

######################################################################
# Importing from an SQLite database
######################################################################

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

import_from_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}) where {T <: MoriDreamSpace} = 
import_from_database(db, "$(db.primary_key) IN ($(join(ids, ", ")))")

import_from_database(db :: SQLiteAdapter{T}, id :: Int) where {T <: MoriDreamSpace} = first(import_from_database(db, [id]))

######################################################################
# Updating columns in an SQLite database
######################################################################

function update_in_database(
        db :: SQLiteAdapter{T}, 
        column_functions :: Dict{Symbol, <:Function};
        sql :: String = "TRUE") where {T <: MoriDreamSpace}

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
        X = sqlite_import_row(T, row)
        val_dict = Dict{Symbol,Any}([k => f(X) for (k,f) in column_functions])
        push!(val_dict, Symbol(primary_key) => row[Symbol(primary_key)])
        DBInterface.execute(update_stmt, val_dict)
    end

end

update_in_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}, column_functions :: Dict{Symbol, <:Function}) where {T <: MoriDreamSpace} = 
update_in_database(db, "$(db.primary_key) IN ($(join(ids, ", ")))", column_functions)

update_in_database(db :: SQLiteAdapter{T}, id :: Int, column_functions :: Dict{Symbol, <:Function}) where {T <: MoriDreamSpace} = update_in_database(db, [id], column_functions)
