# 122 — polynomial gcds, lcms and ideals

This weekly experiment covers PolynomialGCD, PolynomialLCM, PolynomialReduce,
Eliminate, CylindricalDecomposition and GroebnerBasis at stress sizes. It has
40 cases, each with an exact integer check.

```
python3 benchmarks/run_all.py --only 122 --check-labels
```

## Halves

- `polynomial_ideals.m`: Mathilda, and the commercial column `run_all.py`
  times with the same file.
- `polynomial_ideals.mac`: Maxima.
  - `gcd` is core Maxima.
  - LCM is `quotient(a, gcd(a, b)) * b`. The `functs` `lcm` takes 28 s at
    degree 400 and never finishes at degree 1000.
  - Groebner bases, division by a basis and elimination use the `grobner`
    package over its integer ring.
  - The Modulus Groebner cases and the cylindrical decompositions are SKIP:
    the package has no F_p ring, and Maxima has no CAD.
- `polynomial_ideals.gp`: PARI/GP.
  - `gcd` and `lcm`, over Q and with `Mod(1, p)`. The reduction mod p is
    inside the timing, as it is in every other half.
  - Everything else is SKIP.

## Checks

The checks ignore normalisation, so content, sign and monic-ness cannot cause
a disagreement:

- **One polynomial:** `10^40 deg_x(g) + Round[10^6 g(p1)/g(p2)]`, where p1 and
  p2 are two rational points.
- **A list** (a basis or the eliminants): the same, summed over its elements.
- **Under Modulus:** the ratio is taken mod p at two integer points.
- **A decomposition:** membership on a rational grid off every boundary,
  folded into a bitmask.

A reduced Groebner basis is unique for a fixed order. So is the remainder of a
division by one, up to scale. The checks therefore pin the algorithm.

## Results (v0.6080, 2026-10-29; ms, best of three after a warm-up)

| case | Mathilda v0.6076 | Mathilda v0.6080 | PARI/GP | Maxima |
|---|---|---|---|---|
| PolynomialGCD deg-1000, deg-500 gcd | 28.2 | 29.8 | 87.5 | 4180 |
| PolynomialGCD 10 bivariate arguments | 1.27 | 0.239 | 0.488 | 0.960 |
| PolynomialGCD Modulus 32003, deg-1000 | 58.6 | 0.393 | 0.645 | 180 |
| PolynomialGCD Modulus 32003, trivariate | unevaluated | 0.523 | 31.0 | 115 |
| PolynomialLCM deg-1000, deg-500 gcd | 2541 | 44.5 | 96.5 | 4190 |
| PolynomialLCM 10 bivariate arguments | 2789 | 3.8 | 6620 | 31.0 |
| PolynomialReduce by the katsura-4 DRL basis | 681 | 15.8 | — | 195 |
| Eliminate parametric surface | 5.28 | 0.350 | — | 2.3 |
| Eliminate 2 variables from 3 quadratics | 0.199 (wrong) | 0.081 | — | 0.094 |
| CylindricalDecomposition two ellipses | 55.2 | 3.4 | — | — |
| CylindricalDecomposition Whitney umbrella | 516 | 117 | — | — |
| GroebnerBasis cyclic-5 DRL | 29.2 | 3.7 | — | 16.4 |
| GroebnerBasis katsura-6 DRL | > 300 s | 211 | — | 3790 |
| GroebnerBasis cyclic-5 Lex | 114 | 7.6 | — | 24.8 |
| GroebnerBasis katsura-4 Lex | > 300 s | 5.2 | — | 66.6 |
| GroebnerBasis cyclic-6 DRL, Modulus 32003 | 29 280 | 21.7 | — | — |
| GroebnerBasis cyclic-6 DRL | > 300 s | 2830 | — | 3910 |

Mathilda is ahead of Maxima and PARI/GP on every case they run.

Across all four columns, `run_all.py` classifies 37 cases AHEAD, 2 SLOWER and
1 CHECK-FAIL. The three that are not AHEAD:

- **GroebnerBasis cyclic-6 DRL over Q**: 2.83 s, against 0.28 s for the
  fastest system.
  - The same basis mod p takes 22 ms, so the time is the coefficient growth of
    the fraction-free normal forms over Z.
  - The fix is a modular Buchberger: images mod p, CRT, rational
    reconstruction. To certify the result exactly for an inhomogeneous ideal
    like cyclic-6, the system has to be homogenized (a Hilbert-function
    argument), or the S-pairs have to be checked over Q.
- **CylindricalDecomposition, Whitney umbrella**: 117 ms warm (162 ms cold),
  against 26 ms.
  - What remains is the n-variable lift (`cad_leaf`): a `Solve` per fibre and
    evaluator round trips.
  - The fix is a native lift on qqbar and FLINT root isolation.
- **CylindricalDecomposition, quartic curve** (`x^4 + y^4 - 3xy + x < 1`): it
  declines, so the check fails.
  - The section over a root of the degree-12 discriminant needs a degree-48
    algebraic fibre, which is beyond the budget of `rru_algebraic_fiber_roots`.

`results/baseline.json` is the v0.6076 run; `results/final-v0.6080.json` is this
one.
