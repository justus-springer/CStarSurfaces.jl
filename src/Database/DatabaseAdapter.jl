@doc raw"""
    DatabaseAdapter{T <: MoriDreamSpace}

Abstract type of a database adapter holding objects of a common subtype of
`MoriDreamSpace`. Subtypes of `DatabaseAdapter{T}` should at least implement
the following functions:

`import_from_database`, `find_in_database`.

If the database is writeable, the following functions should be implemented as
well:

`export_to_database`, `update_in_database`.

"""
abstract type DatabaseAdapter{T <: MoriDreamSpace} end
