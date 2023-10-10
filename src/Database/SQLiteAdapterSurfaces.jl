
#####################################################
# Julia type for an SQLite connection to a database
# of surfaces with torus action
#####################################################

const SQLiteAdapterSurfaces = SQLiteAdapter{SurfaceWithTorusAction}

##################################################
# Default functions to compute the columns in the 
# SQLite database
##################################################

default_column_functions(::Type{<:SurfaceWithTorusAction}) = Dict([

  :is_toric => is_toric,

  :gen_matrix => function(X)
      P = gen_matrix(X)
      return string([[Int(P[i,j]) for j = 1 : ncols(P)] for i = 1 : nrows(P)])
  end,

  :rays => X -> string(rays(X)),

  :nrays => nrays,

  :lss => X -> is_toric(X) ? Missing : string(Vector(X.l)),

  :dss => X -> is_toric(X) ? Missing : string(Vector(X.d)),

  :case_ => X -> string(X.case),

  :block_sizes => X -> string(Vector(block_sizes(X))),

  :nblocks => nblocks,

  :number_of_parabolic_fixed_point_curves => number_of_parabolic_fixed_point_curves,

  :orientation => orientation,

  :class_group_rank => class_group_rank,

  :class_group_torsion => X -> "[" * join(class_group_torsion(X), ", ") * "]",

  :class_group_torsion_order => X -> Int(class_group_torsion_order(X)),

  :degree_matrix => function(X)
      Q = gen_matrix(X)
      return string([[Int(Q[i,j]) for j = 1 : ncols(Q)] for i = 1 : nrows(Q)])
  end,

  :canonical_divisor_class => function(X)
      c = divisor_class(canonical_divisor_class(X)).coeff
      return string([Int(c[1,i]) for i = 1 : ncols(c)])
  end,

  :gorenstein_index => X -> Int(gorenstein_index(X)),

  :picard_index => X -> Int(picard_index(X)),

  :maximal_log_canonicity_numerator => X -> Int(numerator(maximal_log_canonicity(X))),

  :maximal_log_canonicity_denominator => X -> Int(denominator(maximal_log_canonicity(X))),

  :anticanonical_self_intersection_numerator => X -> Int(numerator(anticanonical_self_intersection(X))),

  :anticanonical_self_intersection_denominator => X -> Int(denominator(anticanonical_self_intersection(X)))

])





