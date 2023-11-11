CREATE TABLE IF NOT EXISTS surfaces (
    -- the unique id of a surface in the database.
    surface_id INTEGER PRIMARY KEY ASC,

    -- 1 if the surface toric, 0 otherwise.
    is_toric INTEGER NOT NULL,

    -- the generator matrix (P-Matrix) of the surface. The columns of 
    -- this matrix are the rays of the canonical toric ambient variety.
    -- This field is stored as a list of rows with each row a list of 
    -- integers. For example: 
    -- "[[-1, -1, 2, 0], [-1, -1, 0, 2], [0, -2, 1, 1]]"
    gen_matrix TEXT NOT NULL,

    -- The rays of the canonical toric ambient, stored as a list of lists
    -- of integers. For example:
    -- "[[-1, -1, 0], [-1, -1, -2], [2, 0, 1], [0, 2, 1]]"
    rays TEXT,

    -- the number of rays of the canonical toric ambient.
    nrays INTEGER,

    -- The vectors l_i=(l_{i1}, ..., l_{in_i}) that make up the upper
    -- rows of the generator matrix of a C-star surface. For example:
    -- [[1, 1], [2], [2]]
    lss TEXT, -- NULL for toric surfaces

    -- The vectors d_i=(d_{i1}, ..., d_{in_i}) that make up the lower
    -- rows of the generator matrix of a C-star surface. For example:
    -- [[0, -2], [1], [1]]
    dss TEXT, -- NULL for toric surfaces

    -- The case of the C-star surface, i.e. one of the strings
    -- "ee", "pe", "ep" and "pp".
    case_ TEXT, -- NULL for toric surfaces

    -- The sizes of the blocks in the generator matrix of a C-star
    -- surface. For example:
    -- [2, 1, 1].
    block_sizes TEXT,-- NULL for toric surfaces

    -- The number of blocks in the generator matrix of a C-star surface
    nblocks INTEGER, -- NULL for toric surfaces

    -- The number of parabolic fixed point curves of a C-star surface.
    -- This is always either 0, 1 or 2.
    number_of_parabolic_fixed_point_curves INTEGER,-- NULL for toric surfaces

    -- The orientation of the generator matrix of a C-star surface.
    -- This is always either -1, 0 or 1.
    orientation INTEGER, -- NULL for toric surfaces

    -- Whether the C-star surface is an intrinsic quadric, i.e. its Cox
    -- Ring has a single quadratic equation.
    is_intrinsic_quadric INTEGER, --NULL for toric surfaces

    -- The rank of the divisor class group.
    class_group_rank INTEGER,

    -- The elementary divisors of the torsion part of the class group,
    -- as a list of integers.
    class_group_torsion TEXT,
    
    -- The order of the torsion part of the class group. This equals the
    -- product of the entries of `class_group_torsion`.
    class_group_torsion_order INTEGER,

    -- A (non-unique) choice of the degree matrix of a surface with torus
    -- action. This is a gale dual of the generator matrix.
    degree_matrix TEXT,

    -- A representation of a canoniacal divisor class as a list of
    -- integer coefficients. By convention, the free part comes first, 
    -- then the torsion part.
    canonical_divisor_class TEXT,

    -- The gorenstein index of a surface with torus action.
    gorenstein_index INTEGER,

    -- The Picard index of a surface with torus action.
    picard_index INTEGER,

    -- the numerator of the macimal log canonicity.
    log_canonicity_numerator INTEGER,

    -- the denominator of the macimal log canonicity.
    log_canonicity_denominator INTEGER,

    -- The maximal rational number $\epsilon$ such that the surface is
    -- $\epsilon$-log terminal, represented by a floating point number.
    -- This column is automatically calculated from numerator and denominator.
    log_canonicity REAL AS (CAST(log_canonicity_numerator AS FLOAT) / CAST(log_canonicity_denominator AS FLOAT)),

    -- the numerator of the anticanonical self intersection.
    anticanonical_self_intersection_numerator INTEGER,

    -- the denominator of the anticanonical self intersection.
    anticanonical_self_intersection_denominator INTEGER,

    -- the anticanonical self intersection as a floating point number.
    anticanonical_self_intersection REAL AS (CAST(anticanonical_self_intersection_numerator AS FLOAT) / CAST(anticanonical_self_intersection_denominator AS FLOAT)),

    -- whether the surface admits a Kaehler-Ricci soliton 
    admits_kaehler_ricci_soliton INTEGER,

    -- whether the surface admits a Kaehler Einstein metric
    admits_kaehler_einstein_metric INTEGER,

    -- whether the surface admits a Sasaki Einstein metric
    admits_sasaki_einstein_metric INTEGER,

    -- whether the surface is quasi-smooth
    is_quasismooth INTEGER,

    -- whether the surface is factorial
    is_factorial INTEGER,

    -- whether the surface is smooth
    is_smooth INTEGER

)













