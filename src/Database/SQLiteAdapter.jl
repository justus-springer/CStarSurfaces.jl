struct SQLiteAdapter{T} <: DatabaseAdapter{T}
    db :: SQLite.DB
    SQLiteAdapter{T}(f :: AbstractString) where {T <: MoriDreamSpace} = new{T}(SQLite.DB(f))
    SQLiteAdapter{T}() where {T <: MoriDreamSpace} = new{T}(SQLite.DB())
end

function export_into_database(
        db :: SQLiteAdapter{T}, 
        table_name :: AbstractString, 
        Xs :: AbstractVector{T}, 
        column_functions :: Dict{Symbol, Function}) where {T <: MoriDreamSpace}

    str = "INSERT INTO $table_name ("
    for col_name in keys(column_functions)
        str *= "$col_name,"
    end
    str = str[begin : end-1] # remove trailing comma
    str *= ") VALUES ("
    for col_name in keys(column_functions)
        str *= ":$col_name,"
    end
    str = str[begin : end-1] # remove trailing comma
    str *= ")"

    stmt = SQLite.Stmt(db, str)

    for X in Xs
        val_dict = Dict([k => f(X) for (k,f) in column_functions])
        DBInterface.execute(stmt, val_dict)
    end
    
end

