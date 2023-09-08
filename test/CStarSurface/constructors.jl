@testset "Constructors for CStarSurface" begin

    X1 = cstar_surface([[1,2,3,1],[5],[1,1]], [[3,3,-1,-2],[2],[0,-1]], :ee)
    X2 = cstar_surface(DoubleVector([[1,2,3,1],[5],[1,1]]), DoubleVector([[3,3,-1,-2],[2],[0,-1]]), :ee)
    X3 = cstar_surface(ZZ[-1 -2 -3 -1 5 0 0; -1 -2 -3 -1 0 1 1; 3 3 -1 -2 2 0 -1])
    
    @test X1.l.parent == [[1,2,3,1],[5],[1,1]]
    @test X1.d.parent == [[3,3,-1,-2],[2],[0,-1]]
    @test X1.case == EE()
    @test gen_matrix(X1) == ZZ[-1 -2 -3 -1 5 0 0; -1 -2 -3 -1 0 1 1; 3 3 -1 -2 2 0 -1] broken=true

    @test X1 == X2
    @test X2 == X3
    @test X1 == X3

    @test cstar_surface(ZZ[-1 -3 2 0 ; -1 -3 0 3 ; 1 -4 1 1]).case == EE()
    @test cstar_surface(ZZ[-1 -3 2 0 0 ; -1 -3 0 3 0 ; 1 -4 1 1 1]).case == PE()
    @test cstar_surface(ZZ[-1 -3 2 0 0 ; -1 -3 0 3 0 ; 1 -4 1 1 -1]).case == EP()
    @test cstar_surface(ZZ[-1 -3 2 0 0 0 ; -1 -3 0 3 0 0 ; 1 -4 1 1 1 -1]).case == PP()

    @test_throws "must be of the same length" cstar_surface([[1,3],[1,3],[2]], [[0,5],[0,-7],[1],[1]], :ee)
    @test_throws "must be of the same length" cstar_surface([[1,3],[1,3],[2]], [[0,5],[-7],[1]], :ee)
    @test_throws "must be coprime" cstar_surface([[1,3],[1,3],[2]], [[0,6],[0,-7],[1]], :ee)

    @test_throws "given matrix is not in P-Matrix shape" cstar_surface(ZZ[-1 -3 2 0 ; -1 -2 0 3 ; 1 -4 1 1])
    @test_throws "given matrix is not in P-Matrix shape" cstar_surface(ZZ[2 -1 -3 0 ; 0 1 -3 3 ; 1 1 -4 1])
    @test_throws "given matrix is not in P-Matrix shape" cstar_surface(ZZ[-1 -3 2 0 2 ; -1 -3 0 3 0 ; 1 -4 1 1 1])
    @test_throws "given matrix is not in P-Matrix shape" cstar_surface(ZZ[-1 -3 2 0 0 ; -1 -3 0 3 0 ; 1 -4 1 1 2])
    @test_throws "given matrix is not in P-Matrix shape" cstar_surface(ZZ[-1 -3 2 0 0 0 0 ; -1 -3 0 3 0 0 0; 1 -4 1 1 1 -1 -1])



end
