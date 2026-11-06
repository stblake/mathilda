# 122 — polynomial gcds, lcms and ideals

This weekly experiment covers PolynomialGCD, PolynomialLCM, PolynomialReduce,
Eliminate, CylindricalDecomposition and GroebnerBasis at stress sizes. It has
69 cases, each with an exact integer check (29 of them the GroebnerBasis
family sweep below).

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

## GroebnerBasis over Q: the standard families (2026-11-06)

The last section of `polynomial_ideals.m` times the families Groebner engines
are compared on: katsura-n, cyclic-n, eco-n, noon-n, reimer-n and the dense
random systems of experiment 141 (four and five quadrics, six-digit
coefficients, two quadrics and two cubics), in DRL and in Lex. Maxima and
PARI/GP skip them: the `grobner` package needs 3.8 s on katsura-6 DRL already,
and PARI has no Groebner bases.

ms, best of three after a warm-up (katsura-8 DRL, katsura-7 Lex and cyclic-7
DRL: one run). The commercial column is the same file, run the same day.
`ratio` is Mathilda over the commercial column.

| case | Mathilda | commercial | ratio |
|---|---|---|---|
| 3-variable quadratics, Lex | 0.133 | 0.623 | 0.21 |
| cyclic-5, DRL | 3.77 | 14.8 | 0.25 |
| katsura-4, DRL | 1.67 | 5.85 | 0.28 |
| elimination form, plane curve | 0.35 | 1.18 | 0.30 |
| cyclic-5, Modulus 32003 | 1.26 | 68.7 | 0.02 |
| katsura-5, DRL | 14.7 | 98.2 | 0.15 |
| cyclic-5, Lex | 12.2 | 66.0 | 0.18 |
| cyclic-6, DRL, Modulus 32003 | 5.84 | 228 | 0.03 |
| cyclic-6, DRL | 55.1 | 291 | 0.19 |
| katsura-4, Lex | 2.41 | 29.1 | 0.08 |
| katsura-6, DRL | 48.6 | 691 | 0.07 |
| noon-4, DRL | 2.55 | 7.25 | 0.35 |
| eco-6, DRL | 3.73 | 7.72 | 0.48 |
| reimer-4, DRL | 16.4 | 17.5 | 0.94 |
| dense 4 quadrics, DRL | 4.24 | 8.21 | 0.52 |
| dense 4 quadrics, 6-digit, DRL | 19.0 | 7.99 | 2.38 |
| dense 5 quadrics, DRL | 24.1 | 44.0 | 0.55 |
| dense 2 quadrics + 2 cubics, DRL | 23.6 | 26.2 | 0.90 |
| eco-7, DRL | 25.9 | 49.2 | 0.53 |
| noon-5, DRL | 47.7 | 74.1 | 0.64 |
| katsura-7, DRL | 241 | 3869 | 0.06 |
| eco-8, DRL | 126 | 473 | 0.27 |
| reimer-5, DRL | 262 | 1209 | 0.22 |
| noon-6, DRL | 326 | 877 | 0.37 |
| katsura-5, Lex | 17.3 | 253 | 0.07 |
| reimer-4, Lex | 26.7 | 123 | 0.22 |
| noon-4, Lex | 31.9 | 214 | 0.15 |
| eco-6, Lex | 15.5 | 32.7 | 0.48 |
| dense 4 quadrics, Lex | 5.76 | 25.8 | 0.22 |
| dense 4 quadrics, 6-digit, Lex | 28.4 | 45.1 | 0.63 |
| dense 5 quadrics, Lex | 73.5 | 298 | 0.25 |
| dense 2 quadrics + 2 cubics, Lex | 62.9 | 215 | 0.29 |
| eco-7, Lex | 112 | 214 | 0.53 |
| cyclic-6, Lex | 164 | 1171 | 0.14 |
| katsura-6, Lex | 344 | 2308 | 0.15 |
| noon-5, Lex | 522 | 5149 | 0.10 |
| reimer-5, Lex | 906 | 6274 | 0.14 |
| katsura-8, DRL | 3988 | 23240 | 0.17 |
| katsura-7, Lex | 6885 | 38673 | 0.18 |
| cyclic-7, DRL | 9670 | 440186 | 0.02 |

Every Groebner case is AHEAD except the six-digit dense quadrics in DRL.

### The engine

- **F4 over GF(p)** (`src/poly/gb_f4.c`) computes every modular basis:
  normal-strategy selection under the Gebauer-Moeller criteria, symbolic
  preprocessing on a monomial hash table, row reduction against the reducer
  rows and the S-rows already reduced (one-word branch-free accumulation for
  p < 2^32, two-word 128-bit sums for 62-bit primes), then a heap
  interreduction. Mod 2^31-1, cyclic-7 DRL takes 0.29 s (Buchberger: 2.8 s) and
  katsura-8 0.20 s (1.5 s).
- **Over Q, multi-modular with an exact certificate** (`src/poly/gb_flint.c`):
  - The images mod 62-bit primes are computed in parallel batches. Unlucky
    primes, whose leading monomials differ, are dropped.
  - The images are combined by CRT and rational reconstruction.
  - The candidate is then proved, never trusted. Every Gebauer-Moeller
    S-pair and every input must reduce to 0 over Z. In a degree order, the
    Hilbert series of the leading terms must also equal that of the
    homogenised system mod p, an upper bound on the ideal's Hilbert function.
  - When the homogenised system has components at infinity (cyclic, noon,
    reimer), the leading terms come from a certified basis of the homogenised
    ideal instead.
  - In Lex, every element must reduce to 0 by the certified DRL basis. For a
    zero-dimensional ideal each prime's Lex image comes from FGLM on the F4
    DRL basis.
  - In shape position, h and the rational univariate representation are
    lifted rather than the Lex basis, whose coefficients are D times larger.
  - Small systems are tried directly first, under a 4 ms budget.

### Where it is behind, and why

- **Dense 4 quadrics, 6-digit coefficients, DRL**: 19 ms against 8.
  - The reduced DRL basis has about 1860-bit coefficients, so reconstruction
    needs a 3700-bit modulus: some 60 primes.
  - The route's time is spread out: about 5 ms of F4 images across 16 worker
    threads, 1.3 ms of CRT, 2.5 ms of reconstruction, 4 ms of exact
    verification (26 S-pairs over Z with 1860-bit coefficients), and the 4 ms
    direct attempt that fails first.
  - Exact fraction-free Buchberger takes 18 ms, so it is no faster.
  - The way forward is p-adic lifting of the basis (Hensel), or a cheaper
    exact verification. Either would remove most of the primes or most of the
    S-pair reductions.
- **reimer-4 DRL and dense 2 quadrics + 2 cubics DRL** are level with the
  commercial column (0.9x). Each is a few ms of fixed cost: the Hilbert check,
  for reimer-4 the certified homogenised basis, and the verification.

The multi-modular routes use worker threads, so CPU time can be several times
the wall time shown (cyclic-7 DRL: 9.7 s wall, about 44 s CPU).
