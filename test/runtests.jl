using Test
using Oscar
using EllipticCoversOscar

@testset "ExplicitMorphism" begin
    phi = ExplicitMorphism("C", "E", "x^2", "y", ("x", "y"), ("u", "v"))

    @test phi.source == "C"
    @test phi.target == "E"
    @test formulas(phi) == ("x^2", "y")
end

@testset "Genus 2 from one point" begin
    E = elliptic_curve(QQ, [0, -8, 0, 8, 0])
    P = [QQ(1), QQ(1)]

    C, F, morphisms = genus2_cover_from_point(E, P)

    @test cover_genus(C) == 2
    @test cover_genus(F) == 1
    @test haskey(morphisms, "C_to_E")
    @test haskey(morphisms, "C_to_F")
    @test haskey(morphisms, "parameters")
end

@testset "Genus 2 from two elliptic curves" begin
    E = elliptic_curve(QQ, [0, -3, 0, 2, 0])
    F = elliptic_curve(QQ, [0, -7, 0, 10, 0])

    alpha_roots = [QQ(0), QQ(1), QQ(2)]
    beta_roots = [QQ(0), QQ(2), QQ(5)]

    C, morphisms = genus2_cover_from_two_elliptic_curves(E, F, alpha_roots, beta_roots)

    @test cover_genus(C) == 2
    @test haskey(morphisms, "C_to_E")
    @test haskey(morphisms, "C_to_F")
    @test haskey(morphisms, "parameters")
end

@testset "Genus 3 from three elliptic curves" begin
    E1 = elliptic_curve(QQ, [0, 3, 0, 2, 0])
    E2 = elliptic_curve(QQ, [0, 5, 0, 6, 0])
    E3 = elliptic_curve(QQ, [0, 7, 0, 12, 0])

    C, morphisms = genus3_cover_from_three_elliptic_curves(
        E1,
        E2,
        E3,
        QQ(0),
        QQ(0),
        QQ(0),
    )

    parameters = morphisms["parameters"]

    @test parameters["A"] == (QQ(3), QQ(5), QQ(7))
    @test parameters["B"] == (QQ(2), QQ(6), QQ(12))
    @test parameters["Delta"] == (QQ(1), QQ(1), QQ(1))
    @test parameters["testing_factor"] == QQ(-128)
    @test parameters["case"] == "plane_quartic"

    @test haskey(morphisms, "C_to_F1")
    @test haskey(morphisms, "C_to_F2")
    @test haskey(morphisms, "C_to_F3")

    @test morphisms["C_to_F1"] isa ExplicitMorphism
    @test morphisms["C_to_F2"] isa ExplicitMorphism
    @test morphisms["C_to_F3"] isa ExplicitMorphism
end
