import OptimalTransport1D as OT

@testset ExtendedTestSet "W1" begin
    N = 5; a = rand(N); b = rand(N)


    @test OT.W1_1D(a,b) ≈ sum(abs.(cumsum(a) - cumsum(b)))
    @test OT.W1_1D(a,b) == OT.W1_1D_unsafe(a,b)

    pop!(b)
    @test_throws DimensionMismatch OT.W1_1D(a, b)
end

@testset ExtendedTestSet "solve_convex_OT" begin
    # Test with full vectors
    N = 5; a = rand(N); b = rand(N)
    a ./= sum(a); b ./= sum(b)
    x = collect(1:N); y = collect(1:N)
    c(x,y) = abs(x-y)
    @test OT.W1_1D(a,b) ≈ OT.solve_convex_OT_cost(a, b, x, y, c)


    # Test with non-full vectors
    N = 10; a = rand(N); b = rand(N)
    x = collect(1:N); y = collect(1:N)
    I = [3,6,8]; J = [2,4] # some index we will delete
    a[I] .= 0; b[J] .= 0
    a ./= sum(a); b ./= sum(b)
    a2 = a[a.>0]; x2 = x[a.>0]
    b2 = b[b.>0]; y2 = y[b.>0]

    @test OT.W1_1D(a,b) ≈ OT.solve_convex_OT_cost(a, b, x, y, c)

    @test OT.W1_1D(a,b) ≈ OT.solve_convex_OT_cost(a2, b2, x2, y2, c)

    # It is symmetric (as long as the cost is)
    @test OT.solve_convex_OT_cost(b2, a2, y2, x2, c) ≈ OT.solve_convex_OT_cost(a2, b2, x2, y2, c)
    
end
