# EllipticCoversOscar.jl

An OSCAR/Julia package for constructing covers of elliptic curves.

The package currently provides:

- explicit morphism containers
- hyperelliptic and projective cover model containers
- genus 2 constructions from elliptic curve data
- a genus 3 construction from three elliptic curves
- tests on simple examples

## Requirements

This package is intended to be used with Julia and OSCAR.

## Installation

From the root of the repository, run:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Running Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Basic Usage

```julia
using Oscar
using EllipticCoversOscar

E = elliptic_curve(QQ, [0, -8, 0, 8, 0])
P = [QQ(1), QQ(1)]

C, F, morphisms = genus2_cover_from_point(E, P)

cover_genus(C)
cover_genus(F)
formulas(morphisms["C_to_E"])
```

## Genus 3 Example

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

## Current Limitations

The package stores maps as explicit rational expressions using `ExplicitMorphism`.
They are not formal OSCAR scheme morphism objects yet.

The two-point genus 2 construction currently supports the case where the two input points are already opposites.
Automatic half-point translation is not implemented yet.
