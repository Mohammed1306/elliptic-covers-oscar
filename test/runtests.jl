using Test
using Oscar
using EllipticCoversOscar

@testset "ExplicitMorphism" begin
    phi = ExplicitMorphism("C", "E", "x^2", "y", ("x", "y"), ("u", "v"))

    @test phi.source == "C"
    @test phi.target == "E"
    @test formulas(phi) == ("x^2", "y")
end