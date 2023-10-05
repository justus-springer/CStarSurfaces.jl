
@testset "normal form and isomorphy test (C-star)" begin

    ##########################################################
    # A list of c-star surfaces, all isomorphic to each other
    ##########################################################

    Xs = [
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; -1 -2 2 2 1]),
    cstar_surface(ZZ[-3 3 0 0 0 ; -3 0 1 1 0 ; -3 0 0 0 2 ; 2 2 -1 -2 1]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; -2 -1 2 2 1]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; -3 -4 2 2 5]),
    cstar_surface(ZZ[-2 3 0 0 0 ; -2 0 3 0 0 ; -2 0 0 1 1 ; 1 2 2 -2 -1]),
    cstar_surface(ZZ[-2 3 0 0 0 ; -2 0 3 0 0 ; -2 0 0 1 1 ; -1 -2 -2 2 1]),
    cstar_surface(ZZ[-2 3 0 0 0 ; -2 0 3 0 0 ; -2 0 0 1 1 ; -7 1 -2 4 3]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; 2 1 -1 -1 -1])]

    for i = 1 : length(Xs), j = i+1 : length(Xs)
        X, Y = Xs[i], Xs[j]
        @test are_isomorphic(X,Y)[1]
        α = are_isomorphic(X,Y)[2]
        @test α(X) == Y
    end


    ##############################################################
    # A list of c-star surfaces, all non-isomorphic to each other
    ##############################################################
    Xs = [
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; -1 -2 2 2 1]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; 0 -2 2 2 1]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; 3 -2 2 2 1]),
    cstar_surface(ZZ[-1 -1 3 0 0 ; -1 -1 0 2 0 ; -1 -1 0 0 3 ; -1 -2 2 1 1]),
    cstar_surface(ZZ[-1 -1 4 0 0 ; -1 -1 0 3 0 ; -1 -1 0 0 2 ; -1 -2 3 2 1])]

    for i = 1 : length(Xs), j = i+1 : length(Xs)
        X, Y = Xs[i], Xs[j]
        @test !are_isomorphic(X,Y)[1]
    end

end
    
    
    


