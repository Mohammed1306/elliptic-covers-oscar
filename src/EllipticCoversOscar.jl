module EllipticCoversOscar

using Oscar

export ExplicitMorphism
export formulas, as_dict
export HyperellipticModel, ProjectiveCoverModel, cover_genus
export genus2_cover_from_two_elliptic_curves
export genus2_cover_from_point
export genus2_cover_from_two_points
export genus3_cover_from_three_elliptic_curves

struct ExplicitMorphism
    source::Any
    target::Any
    x_map::Any
    y_map::Any
    source_coordinates::Any
    target_coordinates::Any
end

function formulas(phi::ExplicitMorphism)
    return phi.x_map, phi.y_map
end

function as_dict(phi::ExplicitMorphism)
    return Dict(
        "source" => phi.source,
        "target" => phi.target,
        "x" => phi.x_map,
        "y" => phi.y_map,
        "source_coordinates" => phi.source_coordinates,
        "target_coordinates" => phi.target_coordinates,
    )
end

struct HyperellipticModel
    base_field::Any
    polynomial::Any
end

function cover_genus(C::HyperellipticModel)
    return div(degree(C.polynomial) - 1, 2)
end

struct ProjectiveCoverModel
    base_field::Any
    ambient_space::Any
    scheme::Any
    equations::Any
    coordinates::Any
end

function _same_base_field(E, F)
    if base_field(E) != base_field(F)
        error("the elliptic curves must be defined over the same base field")
    end
end

function _make_projective_model(equation_builder, K, variable_names)
    P = projective_space(K, variable_names)
    S = homogeneous_coordinate_ring(P)
    vars = Tuple(S[i] for i in 1:length(variable_names))
    equations = equation_builder(vars)
    X = subscheme(P, equations)
    return ProjectiveCoverModel(K, P, X, equations, vars)
end

function _sqrt_in_base_field(value, name)
    if iszero(value)
        return value
    end

    if !is_square(value)
        error("$name is not a square in the current base field")
    end

    return sqrt(value)
end

function _curve_polynomial(E, K)
    R, x = polynomial_ring(K, :x)
    a1, a2, a3, a4, a6 = [K(c) for c in a_invariants(E)]

    if !iszero(a1) || !iszero(a3)
        error("the elliptic curve must be in the form y^2 = f(x)")
    end

    return R, x, x^3 + a2*x^2 + a4*x + a6
end

function genus2_cover_from_two_elliptic_curves(E, F, alpha_roots, beta_roots)
    _same_base_field(E, F)

    K = base_field(E)

    if characteristic(K) == 2
        error("the base field must have characteristic different from 2")
    end

    if length(alpha_roots) != 3 || length(beta_roots) != 3
        error("alpha_roots and beta_roots must each contain exactly 3 roots")
    end

    R, x, f = _curve_polynomial(E, K)
    _, _, g = _curve_polynomial(F, K)

    if !is_separable(f)
        error("f must be separable")
    end

    if !is_separable(g)
        error("g must be separable")
    end

    alpha_roots = [K(a) for a in alpha_roots]
    beta_roots = [K(b) for b in beta_roots]

    for alpha in alpha_roots
        if !iszero(f(alpha))
            error("alpha_roots must contain roots of f")
        end
    end

    for beta in beta_roots
        if !iszero(g(beta))
            error("beta_roots must contain roots of g")
        end
    end

    alpha1, alpha2, alpha3 = alpha_roots
    beta1, beta2, beta3 = beta_roots

    a1 = (
        (alpha3 - alpha2)^2 / (beta3 - beta2)
        + (alpha2 - alpha1)^2 / (beta2 - beta1)
        + (alpha1 - alpha3)^2 / (beta1 - beta3)
    )

    b1 = (
        (beta3 - beta2)^2 / (alpha3 - alpha2)
        + (beta2 - beta1)^2 / (alpha2 - alpha1)
        + (beta1 - beta3)^2 / (alpha1 - alpha3)
    )

    a2 = (
        alpha1*(beta3 - beta2)
        + alpha2*(beta1 - beta3)
        + alpha3*(beta2 - beta1)
    )

    b2 = (
        beta1*(alpha3 - alpha2)
        + beta2*(alpha1 - alpha3)
        + beta3*(alpha2 - alpha1)
    )

    if iszero(a1) || iszero(b1) || iszero(a2) || iszero(b2)
        error("a1, b1, a2 and b2 must be nonzero")
    end

    Delta_f = discriminant(f)
    Delta_g = discriminant(g)

    A = Delta_g * a1 / a2
    B = Delta_f * b1 / b2

    if iszero(A) || iszero(B)
        error("A and B must be nonzero")
    end

    factor1 = (
        A*(alpha2 - alpha1)*(alpha1 - alpha3)*x^2
        + B*(beta2 - beta1)*(beta1 - beta3)
    )

    factor2 = (
        A*(alpha3 - alpha2)*(alpha2 - alpha1)*x^2
        + B*(beta3 - beta2)*(beta2 - beta1)
    )

    factor3 = (
        A*(alpha1 - alpha3)*(alpha3 - alpha2)*x^2
        + B*(beta1 - beta3)*(beta3 - beta2)
    )

    h = -(factor1 * factor2 * factor3)

    if degree(h) != 6
        error("the resulting polynomial h must have degree 6")
    end

    if !is_separable(h)
        error("the resulting polynomial h must be separable")
    end

    C = HyperellipticModel(K, h)

    t1 = -(A / B) * (b2 / b1)

    t2 = (
        beta1*(beta3 - beta2)^2 / (alpha3 - alpha2)
        + beta2*(beta1 - beta3)^2 / (alpha1 - alpha3)
        + beta3*(beta2 - beta1)^2 / (alpha2 - alpha1)
    ) / b1

    s1 = -(B / A) * (a2 / a1)

    s2 = (
        alpha1*(alpha3 - alpha2)^2 / (beta3 - beta2)
        + alpha2*(alpha1 - alpha3)^2 / (beta1 - beta3)
        + alpha3*(alpha2 - alpha1)^2 / (beta2 - beta1)
    ) / a1

    S, (xC, yC) = polynomial_ring(K, [:xC, :yC])
    L = fraction_field(S)
    xC = L(xC)
    yC = L(yC)

    C_to_F = ExplicitMorphism(
        C,
        F,
        t1*xC^2 + t2,
        (Delta_f / B^3)*yC,
        ("xC", "yC"),
        ("xF", "yF"),
    )

    C_to_E = ExplicitMorphism(
        C,
        E,
        s1 / xC^2 + s2,
        (Delta_g / A^3) * (yC / xC^3),
        ("xC", "yC"),
        ("xE", "yE"),
    )

    morphisms = Dict(
        "C_to_E" => C_to_E,
        "C_to_F" => C_to_F,
        "parameters" => Dict(
            "A" => A,
            "B" => B,
            "t1" => t1,
            "t2" => t2,
            "s1" => s1,
            "s2" => s2,
            "cover_polynomial" => h,
        ),
    )

    return C, morphisms
end

function _point_xy(P)
    if P isa Tuple || P isa AbstractVector
        return P[1], P[2]
    end

    error("pass points as [x, y] or (x, y) in this OSCAR port")
end

function genus2_cover_from_point(E, P)
    K = base_field(E)

    if characteristic(K) == 2
        error("the base field must have characteristic different from 2")
    end

    a1, a2, a3, a4, a6 = [K(c) for c in a_invariants(E)]
    xP_raw, yP_raw = _point_xy(P)
    xP = K(xP_raw)
    yP = K(yP_raw)

    R, x = polynomial_ring(K, :x)

    b2 = a1^2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3^2 + 4*a6

    completed_polynomial = (
        x^3
        + (b2 / 4)*x^2
        + (b4 / 2)*x
        + (b6 / 4)
    )

    XP = xP
    YP = yP + (a1*xP + a3) / 2

    transformed_polynomial = completed_polynomial(x + XP)

    if transformed_polynomial(K(0)) != YP^2
        error("coordinate translation did not move P correctly")
    end

    if iszero(transformed_polynomial(K(0)))
        error("P must not become a 2-torsion point after the coordinate changes")
    end

    c2 = coeff(transformed_polynomial, 2)
    c1 = coeff(transformed_polynomial, 1)
    c0 = coeff(transformed_polynomial, 0)

    H = transformed_polynomial(x^2)

    if degree(H) != 6
        error("the resulting polynomial H must have degree 6")
    end

    if !is_separable(H)
        error("the resulting polynomial H must be separable")
    end

    C = HyperellipticModel(K, H)

    complementary_polynomial = 1 + c2*x + c1*x^2 + c0*x^3

    if degree(complementary_polynomial) != 3
        error("the complementary polynomial must have degree 3")
    end

    if !is_separable(complementary_polynomial)
        error("the complementary polynomial must be separable")
    end

    F = HyperellipticModel(K, complementary_polynomial)

    S, (z, Y) = polynomial_ring(K, [:z, :Y])
    L = fraction_field(S)
    z = L(z)
    Y = L(Y)

    C_to_E = ExplicitMorphism(
        C,
        E,
        z^2 + XP,
        Y - (a1*(z^2 + XP) + a3) / 2,
        ("z", "Y"),
        ("x", "y"),
    )

    C_to_F = ExplicitMorphism(
        C,
        F,
        1 / z^2,
        Y / z^3,
        ("z", "Y"),
        ("w", "v"),
    )

    morphisms = Dict(
        "C_to_E" => C_to_E,
        "C_to_F" => C_to_F,
        "parameters" => Dict(
            "XP" => XP,
            "YP" => YP,
            "a_invariants" => (a1, a2, a3, a4, a6),
            "b_invariants" => (b2, b4, b6),
            "completed_polynomial" => completed_polynomial,
            "transformed_polynomial" => transformed_polynomial,
            "cover_polynomial" => H,
            "complementary_polynomial" => complementary_polynomial,
        ),
    )

    return C, F, morphisms
end

function genus2_cover_from_two_points(E, P, Q)
    a1, _, a3, _, _ = a_invariants(E)

    xP, yP = _point_xy(P)
    xQ, yQ = _point_xy(Q)

    if xQ == xP && yQ == -yP - a1*xP - a3
        C, F, morphisms = genus2_cover_from_point(E, P)
        morphisms["phase_3_translation"] = Dict(
            "original_points" => Dict("P" => P, "Q" => Q),
            "relation" => Dict("opposite_points" => "Q = -P"),
        )
        return C, F, morphisms
    end

    error("automatic half-point translation is not implemented yet in the OSCAR port; pass opposite points for now")
end

function _genus3_normal_form_from_root(E, root, K)
    a1, a2, a3, a4, a6 = [K(c) for c in a_invariants(E)]
    root = K(root)

    if characteristic(K) == 2
        error("the base field must have characteristic different from 2")
    end

    R, x = polynomial_ring(K, :x)

    b2 = a1^2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3^2 + 4*a6

    completed_polynomial = (
        x^3
        + (b2 / 4)*x^2
        + (b4 / 2)*x
        + (b6 / 4)
    )

    if !iszero(completed_polynomial(root))
        error("the chosen root is not a root of the completed cubic")
    end

    translated_polynomial = completed_polynomial(x + root)

    if !iszero(coeff(translated_polynomial, 0))
        error("translation failed: constant term is not zero")
    end

    A = coeff(translated_polynomial, 2)
    B = coeff(translated_polynomial, 1)

    if iszero(B)
        error("B must be nonzero in the model y^2 = x(x^2 + A*x + B)")
    end

    Delta = A^2 - 4*B

    if iszero(Delta)
        error("Delta = A^2 - 4B must be nonzero")
    end

    E_normal = elliptic_curve(K, [K(0), A, K(0), B, K(0)])

    return Dict(
        "original_curve" => E,
        "normal_curve" => E_normal,
        "root" => root,
        "A" => A,
        "B" => B,
        "Delta" => Delta,
        "completed_polynomial" => completed_polynomial,
        "translated_polynomial" => translated_polynomial,
    )
end

function genus3_cover_from_three_elliptic_curves(E1, E2, E3, root1, root2, root3)
    if base_field(E1) != base_field(E2) || base_field(E1) != base_field(E3)
        error("the elliptic curves must be defined over the same base field")
    end

    K = base_field(E1)

    if characteristic(K) == 2
        error("the base field must have characteristic different from 2")
    end

    data1 = _genus3_normal_form_from_root(E1, root1, K)
    data2 = _genus3_normal_form_from_root(E2, root2, K)
    data3 = _genus3_normal_form_from_root(E3, root3, K)

    A1, A2, A3 = data1["A"], data2["A"], data3["A"]
    B1, B2, B3 = data1["B"], data2["B"], data3["B"]
    Delta1, Delta2, Delta3 = data1["Delta"], data2["Delta"], data3["Delta"]

    d1 = _sqrt_in_base_field(Delta1, "Delta1")
    d2 = _sqrt_in_base_field(Delta2, "Delta2")
    d3 = _sqrt_in_base_field(Delta3, "Delta3")

    R_value = d1*d2*d3

    testing_factor = (
        R_value
        * (
            A1^2 / Delta1
            + A2^2 / Delta2
            + A3^2 / Delta3
            - 1
        )
        - 2*A1*A2*A3
    )

    parameters = Dict(
        "normal_forms" => Dict("E1" => data1, "E2" => data2, "E3" => data3),
        "A" => (A1, A2, A3),
        "B" => (B1, B2, B3),
        "Delta" => (Delta1, Delta2, Delta3),
        "d" => (d1, d2, d3),
        "R" => R_value,
        "testing_factor" => testing_factor,
    )

    t_ring, t = polynomial_ring(K, :t)

    if iszero(testing_factor)
        coef_a = (R_value*B1 / 2) * (-B1 / Delta1 + B2 / Delta2 + B3 / Delta3)
        coef_b = (R_value*B2 / 2) * (B1 / Delta1 - B2 / Delta2 + B3 / Delta3)
        coef_c = (R_value*B3 / 2) * (B1 / Delta1 + B2 / Delta2 - B3 / Delta3)

        d_abs = _sqrt_in_base_field(1 / (B2*B3), "1/(B2*B3)")
        e_abs = _sqrt_in_base_field(1 / (B1*B3), "1/(B1*B3)")
        f_abs = _sqrt_in_base_field(1 / (B1*B2), "1/(B1*B2)")

        found_signs = false
        coef_d = d_abs
        coef_e = e_abs
        coef_f = f_abs

        for sd in [K(1), K(-1)]
            for se in [K(1), K(-1)]
                for sf in [K(1), K(-1)]
                    trial_d = sd*d_abs
                    trial_e = se*e_abs
                    trial_f = sf*f_abs

                    if A1 == -coef_a*trial_e*trial_f &&
                       A2 == -coef_b*trial_d*trial_f &&
                       A3 == -coef_c*trial_d*trial_e
                        coef_d = trial_d
                        coef_e = trial_e
                        coef_f = trial_f
                        found_signs = true
                    end
                end
            end
        end

        if !found_signs
            error("could not choose signs for d, e, f in the hyperelliptic case")
        end

        C = _make_projective_model(K, [:W, :X, :Y, :Z]) do vars
            W, X, Y, Z = vars
            equation1 = W^2 * Z^2 - (coef_a*X^4 + coef_b*Y^4 + coef_c*Z^4)
            equation2 = coef_d*X^2 + coef_e*Y^2 + coef_f*Z^2
            return [equation1, equation2]
        end

        q1 = (coef_a*coef_e^2 + coef_b*coef_d^2)*t^4 + 2*coef_a*coef_e*coef_f*t^2 + (coef_a*coef_f^2 + coef_c*coef_d^2)
        q2 = (coef_a*coef_e^2 + coef_b*coef_d^2)*t^4 + 2*coef_b*coef_d*coef_f*t^2 + (coef_b*coef_f^2 + coef_c*coef_e^2)
        q3 = (coef_a*coef_f^2 + coef_c*coef_d^2)*t^4 + 2*coef_c*coef_d*coef_e*t^2 + (coef_b*coef_f^2 + coef_c*coef_e^2)

        F1 = HyperellipticModel(K, q1)
        F2 = HyperellipticModel(K, q2)
        F3 = HyperellipticModel(K, q3)

        S, (W, X, Y, Z) = polynomial_ring(K, [:W, :X, :Y, :Z])
        L = fraction_field(S)
        W, X, Y, Z = L(W), L(X), L(Y), L(Z)

        morphisms = Dict(
            "C_to_F1" => ExplicitMorphism(C, F1, Y / Z, coef_d * W / Z, ("W", "X", "Y", "Z"), ("t", "v")),
            "C_to_F2" => ExplicitMorphism(C, F2, X / Z, coef_e * W / Z, ("W", "X", "Y", "Z"), ("t", "v")),
            "C_to_F3" => ExplicitMorphism(C, F3, X / Y, coef_f * W * Z / Y^2, ("W", "X", "Y", "Z"), ("t", "v")),
        )

        parameters["case"] = "hyperelliptic"
        parameters["curve_type"] = "hyperelliptic genus 3 curve"
        parameters["coefficients"] = Dict("a" => coef_a, "b" => coef_b, "c" => coef_c, "d" => coef_d, "e" => coef_e, "f" => coef_f)
        parameters["equations"] = C.equations
        parameters["quotient_curves"] = Dict("F1" => F1, "F2" => F2, "F3" => F3)
        parameters["quotient_polynomials"] = Dict("F1" => q1, "F2" => q2, "F3" => q3)

        morphisms["parameters"] = parameters

        return C, morphisms
    end

    coef_d = (K(1) / 2) * (-A1*A2 + (A3*R_value) / Delta3)
    coef_e = (K(1) / 2) * (-A1*A3 + (A2*R_value) / Delta2)
    coef_f = (K(1) / 2) * (-A2*A3 + (A1*R_value) / Delta1)

    C = _make_projective_model(K, [:X, :Y, :Z]) do vars
        X, Y, Z = vars
        quartic = (
            B1*X^4
            + B2*Y^4
            + B3*Z^4
            + coef_d*X^2*Y^2
            + coef_e*X^2*Z^2
            + coef_f*Y^2*Z^2
        )
        return [quartic]
    end

    q1 = (coef_d^2 - 4*B1*B2)*t^4 + (2*coef_d*coef_e - 4*B1*coef_f)*t^2 + (coef_e^2 - 4*B1*B3)
    q2 = (coef_d^2 - 4*B1*B2)*t^4 + (2*coef_d*coef_f - 4*B2*coef_e)*t^2 + (coef_f^2 - 4*B2*B3)
    q3 = (coef_e^2 - 4*B1*B3)*t^4 + (2*coef_e*coef_f - 4*B3*coef_d)*t^2 + (coef_f^2 - 4*B2*B3)

    F1 = HyperellipticModel(K, q1)
    F2 = HyperellipticModel(K, q2)
    F3 = HyperellipticModel(K, q3)

    S, (X, Y, Z) = polynomial_ring(K, [:X, :Y, :Z])
    L = fraction_field(S)
    X, Y, Z = L(X), L(Y), L(Z)

    C_to_F1 = ExplicitMorphism(
        C,
        F1,
        Y / Z,
        2*B1*(X / Z)^2 + coef_d*(Y / Z)^2 + coef_e,
        ("X", "Y", "Z"),
        ("t", "v"),
    )

    C_to_F2 = ExplicitMorphism(
        C,
        F2,
        X / Z,
        2*B2*(Y / Z)^2 + coef_d*(X / Z)^2 + coef_f,
        ("X", "Y", "Z"),
        ("t", "v"),
    )

    C_to_F3 = ExplicitMorphism(
        C,
        F3,
        X / Y,
        2*B3*(Z / Y)^2 + coef_e*(X / Y)^2 + coef_f,
        ("X", "Y", "Z"),
        ("t", "v"),
    )

    quotient_elliptic_models = Dict(
        "E1_twisted_jacobian_model" => elliptic_curve(K, [
            K(0),
            2*B1*coef_f - coef_d*coef_e,
            K(0),
            B1*(B1*coef_f^2 + B2*coef_e^2 + B3*coef_d^2 - coef_d*coef_e*coef_f - 4*B1*B2*B3),
            K(0),
        ]),
        "E2_twisted_jacobian_model" => elliptic_curve(K, [
            K(0),
            2*B2*coef_e - coef_d*coef_f,
            K(0),
            B2*(B2*coef_e^2 + B1*coef_f^2 + B3*coef_d^2 - coef_d*coef_e*coef_f - 4*B1*B2*B3),
            K(0),
        ]),
        "E3_twisted_jacobian_model" => elliptic_curve(K, [
            K(0),
            2*B3*coef_d - coef_e*coef_f,
            K(0),
            B3*(B3*coef_d^2 + B1*coef_f^2 + B2*coef_e^2 - coef_d*coef_e*coef_f - 4*B1*B2*B3),
            K(0),
        ]),
    )

    parameters["case"] = "plane_quartic"
    parameters["curve_type"] = "plane quartic genus 3 curve"
    parameters["coefficients"] = Dict("d" => coef_d, "e" => coef_e, "f" => coef_f)
    parameters["equation"] = C.equations[1]
    parameters["quotient_curves"] = Dict("F1" => F1, "F2" => F2, "F3" => F3)
    parameters["quotient_polynomials"] = Dict("F1" => q1, "F2" => q2, "F3" => q3)
    parameters["quotient_elliptic_models"] = quotient_elliptic_models
    parameters["note"] = "In the plane quartic case, the quotient elliptic models may become isomorphic to the input elliptic curves after a quadratic field extension."

    morphisms = Dict(
        "C_to_F1" => C_to_F1,
        "C_to_F2" => C_to_F2,
        "C_to_F3" => C_to_F3,
        "parameters" => parameters,
    )

    return C, morphisms
end

end
