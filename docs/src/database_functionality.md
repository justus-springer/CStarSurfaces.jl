# Database functionality

CStarSurfaces.jl has the ability to interact with database connections for
saving and retrieving Mori dream spaces and their properties. The interface here
is generic, but currently only an SQLite database holding
[`SurfaceWithTorusAction`](@ref)'s is supported.

A connection to the online
[ldp-database](https://www.math.uni-tuebingen.de/forschung/algebra/ldp-database/)
is planned.

```@docs
DatabaseAdapter
SQLiteAdapter
SQLiteAdapterSurfaces
import_from_database
find_in_database
default_column_functions
default_insert_predicate
export_to_database
update_in_database
```
