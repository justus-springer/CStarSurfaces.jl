
DB_FILE_PATH = joinpath(dirname(@__FILE__), "test.db")
# DB_FILE_PATH = "test/Database/test.db"
CREATE_TABLE_FILE_PATH = joinpath(dirname(@__FILE__), "create.sql")
# CREATE_TABLE_FILE_PATH = "test/Database/create.sql"
TABLE_NAME = "surfaces"
PRIMARY_KEY = "surface_id"

@testset "Insert surfaces into test database" begin 

    rm(DB_FILE_PATH, force=true)
    db = SQLiteAdapter{SurfaceWithTorusAction}(DB_FILE_PATH, TABLE_NAME, PRIMARY_KEY)
    DBInterface.execute(db.db, read(CREATE_TABLE_FILE_PATH, String))

    @testset "Insert C-star surface" begin

        X = cstar_surface([[1,3], [2], [3]], [[1,-4], [1], [1]], :ee)
        export_to_database(db, [X]) 
        row = first(DBInterface.execute(db.db, "SELECT * FROM $TABLE_NAME WHERE $PRIMARY_KEY = 1"))

        @test row[:surface_id] == 1
        @test row[:is_toric] == 0
        @test row[:gen_matrix] == "[[-1, -3, 2, 0], [-1, -3, 0, 3], [1, -4, 1, 1]]"
        @test row[:rays] ==  "[[-1, -1, 1], [-3, -3, -4], [2, 0, 1], [0, 3, 1]]"
        @test row[:nrays] == 4
        @test row[:lss] == "[[1, 3], [2], [3]]"
        @test row[:dss] == "[[1, -4], [1], [1]]" 
        @test row[:case_] == "ee"
        @test row[:block_sizes] == "[2, 1, 1]"
        @test row[:nblocks] == 3
        @test row[:number_of_parabolic_fixed_point_curves] == 0
        @test row[:orientation] == 1
        @test row[:class_group_rank] == 1 
        @test row[:class_group_torsion] == "[]"
        @test row[:class_group_torsion_order] == 1
        @test row[:degree_matrix] == "[[9, 11, 21, 14]]"
        @test row[:canonical_divisor_class] == "[-13]"
        @test row[:gorenstein_index] == 693
        @test row[:picard_index] == 693
        @test row[:log_canonicity_numerator] == 1
        @test row[:log_canonicity_denominator] == 3
        @test row[:anticanonical_self_intersection_numerator] == 169
        @test row[:anticanonical_self_intersection_denominator] == 693
        @test row[:admits_kaehler_einstein_metric] == 0
        @test row[:is_quasismooth] == 0
        @test row[:is_factorial] == 0
        @test row[:is_smooth] == 0

        @test import_from_database(db, row[:surface_id]) == X

        # Exporting again should do nothing
        export_to_database(db, [X]) 
        @test first(DBInterface.execute(db.db, "SELECT COUNT(*) FROM $TABLE_NAME"))[1] == 1

    end

    @testset "Insert toric surface" begin

        X = toric_surface(ZZ[1 1 -5 ; 0 3 -12])
        export_to_database(db, [X]) 
        row = first(DBInterface.execute(db.db, "SELECT * FROM SURFACES WHERE surface_id = 2"))

        @test row[:surface_id] == 2
        @test row[:is_toric] == 1
        @test row[:gen_matrix] == "[[1, 1, -5], [0, 3, -12]]"
        @test row[:rays] == "[[1, 0], [1, 3], [-5, -12]]"
        @test row[:nrays] == 3
        @test ismissing(row[:lss])
        @test ismissing(row[:dss])
        @test ismissing(row[:case_])
        @test ismissing(row[:block_sizes])
        @test ismissing(row[:nblocks])
        @test ismissing(row[:number_of_parabolic_fixed_point_curves])
        @test ismissing(row[:orientation])
        @test row[:class_group_rank] == 1 
        @test row[:class_group_torsion] == "[3]"
        @test row[:class_group_torsion_order] == 3
        @test row[:degree_matrix] == "[[0, 2, 1], [1, 4, 1]]"
        @test row[:canonical_divisor_class] == "[0, -6]"
        @test row[:gorenstein_index] == 2
        @test row[:picard_index] == 36
        @test row[:log_canonicity_numerator] == 1
        @test row[:log_canonicity_denominator] == 2
        @test row[:anticanonical_self_intersection_numerator] == 3
        @test row[:anticanonical_self_intersection_denominator] == 1
        @test row[:admits_kaehler_einstein_metric] == 0
        @test row[:is_quasismooth] == 1
        @test row[:is_factorial] == 0
        @test row[:is_smooth] == 0

        @test import_from_database(db, row[:surface_id]) == X

        # Exporting again should do nothing
        export_to_database(db, [X]) 
        @test first(DBInterface.execute(db.db, "SELECT COUNT(*) FROM $TABLE_NAME"))[1] == 2

    end

    rm(DB_FILE_PATH)

end
