# Elliptic Covers OSCAR

OSCAR/Julia-compatible package for computing covers of elliptic curves.

This repository contains the OSCAR version of the package development work for the project **Computing Covers of Elliptic Curves**.

The package currently focuses on:

- constructing genus 2 covers of elliptic curves
- returning explicit formulas for the corresponding maps
- implementing a genus 3 construction from three elliptic curves
- mirroring the behavior of the SageMath package as closely as possible in OSCAR/Julia

---

## Requirements

This package is meant to be used with **Julia** and **OSCAR**.

Check that Julia is available:

```bash
julia --version
```

Check that OSCAR is available:

```bash
julia -e 'using Oscar; println("OSCAR loaded")'
```

---

## Installation for Development

From the root of the repository, run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This activates the local Julia project and installs the package dependencies.

To enter the Julia REPL with the local project activated, run:

```bash
julia --project=.
```

Then inside Julia:

```julia
using EllipticCoversOscar
```

---

## Package Structure

```text
elliptic-covers-oscar/
├── src/
│   └── EllipticCoversOscar.jl
├── test/
│   └── runtests.jl
├── Project.toml
├── Manifest.toml
├── README.md
└── .gitignore
```

The full implementation is currently in:

```text
src/EllipticCoversOscar.jl
```

The tests are in:

```text
test/runtests.jl
```

---

## Implemented Functions

### 1. Genus 2 construction from two elliptic curves

```julia
genus2_cover_from_two_elliptic_curves(E, F, alpha_roots, beta_roots)
```

This function implements the genus 2 construction from two elliptic curves.

Input:

- `E`: elliptic curve in the form `y^2 = f(x)`
- `F`: elliptic curve in the form `y^2 = g(x)`
- `alpha_roots`: ordered roots of `f`
- `beta_roots`: ordered roots of `g`

The order of the roots determines the matching of the 2-torsion points.

Output:

```julia
C, morphisms
```

where:

- `C` is the constructed genus 2 hyperelliptic model
- `morphisms` contains explicit formulas for the maps:
  - `C_to_E`
  - `C_to_F`

Example:

```julia
using Oscar
using EllipticCoversOscar

E = elliptic_curve(QQ, [0, -3, 0, 2, 0])
F = elliptic_curve(QQ, [0, -7, 0, 10, 0])

alpha_roots = [QQ(0), QQ(1), QQ(2)]
beta_roots = [QQ(0), QQ(2), QQ(5)]

C, morphisms = genus2_cover_from_two_elliptic_curves(
    E,
    F,
    alpha_roots,
    beta_roots,
)

cover_genus(C)
```

---

### 2. Genus 2 construction from one elliptic curve and one point

```julia
genus2_cover_from_point(E, P)
```

This function constructs a genus 2 cover starting from one elliptic curve `E` and one finite point `P` on `E`.

The function performs the following steps:

1. Removes the `a1` and `a3` terms from the elliptic curve equation.
2. Translates the x-coordinate so that the point `P` has x-coordinate `0`.
3. Constructs the genus 2 curve by substituting `x^2` into the transformed cubic.
4. Constructs the complementary genus 1 curve.
5. Returns explicit formulas for the maps.

Output:

```julia
C, F, morphisms
```

where:

- `C` is the genus 2 hyperelliptic model
- `F` is the complementary genus 1 hyperelliptic model
- `morphisms` contains explicit formulas for:
  - `C_to_E`
  - `C_to_F`

Example:

```julia
using Oscar
using EllipticCoversOscar

E = elliptic_curve(QQ, [0, -8, 0, 8, 0])
P = E([QQ(1), QQ(1)])

C, F, morphisms = genus2_cover_from_point(E, P)

cover_genus(C)
cover_genus(F)
```

---

### 3. Genus 2 construction from one elliptic curve and two points

```julia
genus2_cover_from_two_points(E, P, Q)
```

This function starts with one elliptic curve `E` and two finite points `P` and `Q` on `E`.

The goal is to translate the points so that they become opposites.

The function looks for a point `T` such that:

```text
2*T = -(P + Q)
```

Then it defines:

```text
P_new = P + T
Q_new = Q + T
```

so that:

```text
Q_new = -P_new
```

After this, it calls the one-point construction on `P_new`.

Output:

```julia
C, F, morphisms
```

where:

- `C` is the genus 2 hyperelliptic model
- `F` is the complementary genus 1 hyperelliptic model
- `morphisms` contains the formulas inherited from the one-point construction, plus information about the Phase 3 translation

Example:

```julia
using Oscar
using EllipticCoversOscar

E = elliptic_curve(QQ, [0, -8, 0, 8, 0])

P = E([QQ(1), QQ(1)])
Q = -P

C, F, morphisms = genus2_cover_from_two_points(E, P, Q)

cover_genus(C)
cover_genus(F)
```

---

### 4. Genus 3 construction from three elliptic curves

```julia
genus3_cover_from_three_elliptic_curves(E1, E2, E3, root1, root2, root3)
```

This function constructs a genus 3 curve from three elliptic curves.

It takes three elliptic curves and one chosen root for each curve. Each chosen root is translated to `0`, so that each elliptic curve is put into the form:

```text
y^2 = x(x^2 + A_i*x + B_i)
```

The function then computes:

```text
Delta_i = A_i^2 - 4*B_i
```

and a testing factor.

Depending on the testing factor, it returns either:

- the hyperelliptic genus 3 construction, if the testing factor is zero
- the plane quartic genus 3 construction, if the testing factor is nonzero

Output:

```julia
C, morphisms
```

where:

- `C` is the constructed genus 3 model
- `morphisms` contains explicit quotient maps from `C` to three genus-one quotient curves, together with construction parameters

Example:

```julia
using Oscar
using EllipticCoversOscar

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

morphisms["parameters"]["A"]
morphisms["parameters"]["B"]
morphisms["parameters"]["Delta"]
morphisms["parameters"]["testing_factor"]
morphisms["parameters"]["case"]
```

---

## Morphism Output Format

The package returns maps using the `ExplicitMorphism` struct.

This struct stores:

- the source curve
- the target curve
- the x-coordinate formula
- the y-coordinate formula
- the source coordinate names
- the target coordinate names

Example structure:

```julia
morphisms = Dict(
    "C_to_E" => ExplicitMorphism(...),
    "C_to_F" => ExplicitMorphism(...),
    "parameters" => Dict(
        ...
    ),
)
```

The coordinate formulas can be accessed by:

```julia
x_formula, y_formula = formulas(morphisms["C_to_E"])
```

or directly by:

```julia
morphisms["C_to_E"].x
morphisms["C_to_E"].y
```

The morphisms are explicit OSCAR/Julia rational expressions. They are not yet formal OSCAR scheme morphism objects.

For the two-point construction, the dictionary also contains:

```julia
morphisms["phase_3_translation"]
```

This stores the translation point and the translated points used to reduce the two-point case to the one-point case.

---

## Running Tests

Run all tests with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

To enter the test environment manually:

```bash
julia --project=.
```

Then inside Julia:

```julia
using Pkg
Pkg.test()
```

---

## Development Notes

The implemented genus 2 constructions require characteristic different from `2`.

The genus 3 construction also requires characteristic different from `2`.

The code checks that:

- input points lie on the given elliptic curve
- points are finite
- elliptic curves are in the supported Weierstrass form when required
- relevant polynomials are separable
- the constructed genus 2 polynomial has degree `6`
- the complementary polynomial has degree `3`
- the chosen genus 3 roots are valid roots of the completed cubic
- the genus 3 normal forms have nonzero `B_i`
- the values `Delta_i = A_i^2 - 4B_i` are nonzero
- the required square roots exist in the current base field

The two-point genus 2 construction depends on OSCAR being able to find a point `T` satisfying:

```text
2*T = -(P + Q)
```

over the current base field.

The genus 3 construction currently requires the relevant square roots to exist in the current base field. Later versions may support automatic field extensions.

---

## Current Status

Implemented:

- genus 2 construction from two elliptic curves
- genus 2 construction from one elliptic curve and one point
- genus 2 construction from one elliptic curve and two points
- genus 3 construction from three elliptic curves
- explicit morphism formula objects using `ExplicitMorphism`
- explicit quotient map objects for the genus 3 construction
- unit tests for genus 2 and genus 3 functions

Current differences from the Sage version:

- hyperelliptic curves are represented by `HyperellipticModel`
- projective genus 3 curves are represented by `ProjectiveCoverModel`
- genus is accessed with `cover_genus(C)` instead of `C.genus()`
- morphisms are explicit rational expressions, not formal OSCAR scheme morphism objects

Future work:

- convert explicit morphism formulas into formal OSCAR morphism objects
- improve support for field extensions when half-points or square roots are not defined over the current base field
- add more mathematical examples and documentation
- expand tests for more input curves and point configurations
