
DB_FILE_PATH = joinpath(dirname(@__FILE__), "test.db")
# DB_FILE_PATH = "test/Database/test.db"
CREATE_TABLE_FILE_PATH = joinpath(dirname(@__FILE__), "create.sql")
# CREATE_TABLE_FILE_PATH = "test/Database/create.sql"
TABLE_NAME = "surfaces"

@testset "Insert surfaces into test database" begin 

    rm(DB_FILE_PATH, force=true)
    db = SQLiteAdapter{SurfaceWithTorusAction}(DB_FILE_PATH)
    DBInterface.execute(db.db, read(CREATE_TABLE_FILE_PATH, String))

    @testset "Insert C-star surface" begin

        X = cstar_surface([[1,3], [2], [3]], [[1,-4], [1], [1]], :ee)
        export_to_database(db, TABLE_NAME, [X]) 
        row = first(DBInterface.execute(db.db, "SELECT * FROM SURFACES WHERE surface_id = 1"))

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
        @test row[:maximal_log_canonicity_numerator] == 1
        @test row[:maximal_log_canonicity_denominator] == 3
        @test row[:anticanonical_self_intersection_numerator] == 169
        @test row[:anticanonical_self_intersection_denominator] == 693

    end

    @testset "Insert toric surface" begin

        X = toric_surface(ZZ[-5 1 1 ; -3 0 3])
        export_to_database(db, TABLE_NAME, [X]) 
        row = first(DBInterface.execute(db.db, "SELECT * FROM SURFACES WHERE surface_id = 2"))

        @test row[:surface_id] == 2
        @test row[:is_toric] == 1
        @test row[:gen_matrix] == "[[-5, 1, 1], [-3, 0, 3]]"
        @test row[:rays] == "[[-5, -3], [1, 0], [1, 3]]"
        @test row[:nrays] == 3
        @test row[:lss] == Missing
        @test row[:dss] == Missing
        @test row[:case_] == Missing
        @test row[:block_sizes] == Missing
        @test row[:nblocks] == Missing
        @test row[:number_of_parabolic_fixed_point_curves] == Missing
        @test row[:orientation] == Missing
        @test row[:class_group_rank] == 1 
        @test row[:class_group_torsion] == "[3]"
        @test row[:class_group_torsion_order] == 3
        @test row[:degree_matrix] == "[[0, 2, 1], [1, 4, 1]]"
        @test row[:canonical_divisor_class] == "[0, -6]"
        @test row[:gorenstein_index] == 2
        @test row[:picard_index] == 36
        @test row[:maximal_log_canonicity_numerator] == 1
        @test row[:maximal_log_canonicity_denominator] == 2
        @test row[:anticanonical_self_intersection_numerator] == 3
        @test row[:anticanonical_self_intersection_denominator] == 1

    end

    rm(DB_FILE_PATH)

end
