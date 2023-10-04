
const SurfaceWithTorusAction = Union{CStarSurface, ToricSurface}

dim(::SurfaceWithTorusAction) = 2

#################################################
# Picard Index
#
# For surfaces with torus action, we can use the formula for
# the picard index in terms of the local class group,
# which is more efficient than the general algorithm
#################################################

picard_index(X :: SurfaceWithTorusAction) = prod(values(local_picard_indices(X))) / class_group_torsion_order(X)


#################################################
# Anticanonical self intersection
#################################################

anticanonical_self_intersection(X :: SurfaceWithTorusAction) = anticanonical_divisor(X) * anticanonical_divisor(X)
