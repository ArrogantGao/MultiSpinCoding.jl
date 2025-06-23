using MultiSpinCoding
using Graphs
using Test

@testset "io" begin
    g, J = load_instance("../example/square16x16.txt")
    @test g isa SimpleGraph{Int}
    @test nv(g) == 256
    @test size(J) == (256, 256)
end