_cycles(n :: Int) = [map(x -> mod(x+i, 1:n), 1:n) for i = 1:n]

# TODO: keep track of and return the unimodular transformation and 
# column permutation turning X into its normal form
@doc raw"""
    normal_form(X :: ToricSurface)

Bring a toric surface into normal form. Here, we call a toric surface with
generator matrix $P$ in normal form, if and only if:

1.  $P$ in Hermite normal form,
2.  neighbouring columns in $P$ are adjacent rays (i.e. the columns are 
    sorted either clockwise or counterclockwise),
3.  $P$ is lexicographically minimal among all generator matrices satisfying 1. and 2.

These properties are enough to ensure that that two toric surfaces `X` and `Y`
are isomorphic if and only if `gen_matrix(normal_form(X)) ==
gen_matrix(normal_form(Y))`.

"""
@attr function normal_form(X :: ToricSurface)
    r = nrays(X)
    # sort the rays counterclockwise
    ordered_rays = rays(X)[_ordered_ray_indices(X)]
    # take the hermite normal form of all possible permutations of the
    # rays keeping the counterclockwise ordering
    cycls = [_cycles(r) ; map(reverse, _cycles(r))]
    As = map(is -> hnf(matrix(ZZ, hcat(ordered_rays[is]...))), cycls)
    # take the lexicographical minumum of all the hermite normal forms
    _lt(A, B) = vcat(A...) < vcat(B...)
    A = sort(As; lt = _lt)[1]

    Y = toric_surface(A)
    # Remember that Y is in normal form, so it doesn't have to be
    # computed again
    set_attribute!(Y, :is_normal_form, true)
    set_attribute!(Y, :normal_form, Y)

    return Y
end


@doc raw"""
    is_normal_form(X :: ToricSurface)

Check whether a toric surface is in normal form.

"""
@attr is_normal_form(X :: ToricSurface) = normal_form(X) == X


@doc raw"""
    are_isomorphic(X :: ToricSurface, Y :: ToricSurface)

Check whether two toric surfaces are isomorphic to each other.

"""
are_isomorphic(X :: ToricSurface, Y :: ToricSurface) = normal_form(X) == normal_form(Y)
