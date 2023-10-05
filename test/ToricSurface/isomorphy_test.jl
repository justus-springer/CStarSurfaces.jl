
@testset "normal form and isomorphy test (toric)" begin

    # A list of toric surfaces, all isomorphic to each other
    Xs = [
          toric_surface(ZZ[-1 -2 -1 7 0 1 ; 0 -3 -3 3 1 0]),
          toric_surface(ZZ[-1 7 -1 -2 1 0 ; 0 3 -3 -3 0 1]),
          toric_surface(ZZ[-1 7 -1 -2 1 0 ; -1 10 -4 -5 1 1]),
          toric_surface(ZZ[-3 0 -3 3 1 0 ; -2 -1 -1 7 0 1]),
          toric_surface(ZZ[-3 0 -3 3 1 0 ; 1 -1 2 4 -1 1]),
          toric_surface(ZZ[1 0 -3 -3 0 3 ; 0 1 2 1 -1 -7])
    ]
        
    for i = 1 : length(Xs), j = i+1 : length(Xs)
        X, Y = Xs[i], Xs[j]
        @test are_isomorphic(X,Y)
    end

    # A list of toric surfaces, all non-isomorphic to each other
    Xs = [
          toric_surface(ZZ[-1 -2 -1 7 0 1 ; 0 -3 -3 3 1 0]),
          toric_surface(ZZ[-2 -2 -1 7 0 1 ; 0 -3 -3 3 1 0]),
          toric_surface(ZZ[-1 -2 -1 7 0 ; 0 -3 -3 3 1])
    ]

    for i = 1 : length(Xs), j = i+1 : length(Xs)
        X, Y = Xs[i], Xs[j]
        @test !are_isomorphic(X,Y)
    end

end
