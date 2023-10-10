struct SQLiteAdapter{T} <: DatabaseAdapter{T}
    db :: SQLite.DB
    table_name :: AbstractString
    primary_key :: AbstractString

    SQLiteAdapter{T}(f :: AbstractString, table_name :: AbstractString, primary_key :: AbstractString) where {T <: MoriDreamSpace} = new{T}(SQLite.DB(f), table_name, primary_key)
end

# Fallback definition. Subtypes of `MoriDreamSpace` should implement this,
# describing the default column names the functions to compute them
# when inserting into a database..
default_column_functions(::Type{T}) where {T <: MoriDreamSpace} = 
Dict{Symbol, Function}([])

default_insert_predicate(::Type{T}) where {T <: MoriDreamSpace} =
function(X...) true end

function export_to_database(
        db_adapter :: SQLiteAdapter{T}, 
        table_name :: AbstractString, 
        Xs :: AbstractVector;
        column_functions = default_column_functions(T),
        insert_predicate = default_insert_predicate(T)) where {T <: MoriDreamSpace}

    db = db_adapter.db
    if isempty(column_functions)
        stmt = SQLite.Stmt(db, "INSERT INTO $table_name DEFAULT VALUES")
    else
        columns = join(keys(column_functions), ",")
        values = join(":" .* string.(keys(column_functions)), ",")
        stmt = SQLite.Stmt(db, "INSERT INTO $table_name ($columns) VALUES ($values)")
    end

    for X in Xs
        !insert_predicate(X) && continue
        val_dict = Dict([k => f(X) for (k,f) in column_functions])
        DBInterface.execute(stmt, val_dict)
    end
    
end

function import_from_database(db :: SQLiteAdapter{T}, str :: AbstractString) where {T <: MoriDreamSpace}
    stmt = SQLite.Stmt(db.db, str)
    Xs = T[]
    for row in DBInterface.execute(stmt)
        push!(Xs, sqlite_import_row(T, row))
    end
    return Xs
end

import_from_database(db :: SQLiteAdapter{T}, ids :: AbstractVector{Int}) where {T <: MoriDreamSpace} = 
import_from_database(db, "SELECT * FROM $(db.table_name) WHERE $(db.primary_key) IN ($(join(ids, ", ")))")

import_from_database(db :: SQLiteAdapter{T}, id :: Int) where {T <: MoriDreamSpace} = import_from_database(db, [id])



