# The ISSAC '97 Systems Challenge (experiments 126–135)

Ten problems posed to the computer algebra systems at ISSAC '97, one experiment
per problem, each timing **every method** of the 1997 solution notebook at
machine precision, at arbitrary precision (the notebook's precision and a
higher one for the cost-vs-digits curve), and exactly where the problem is
exact. The cases are the notebook's code, modernised only where 1997 syntax no
longer runs (string method names, `RandomReal`, a 1 × 1 matrix `Jacobian`,
nodes for `GaussianQuadratureWeights` computed in-file, no plots).

```bash
python3 benchmarks/run_all.py --only 126,127,128,129,130,131,132,133,134,135 -v
cd benchmarks/129-issac97-p4-polynomial-coefficient && ../../Mathilda -file polynomial_coefficient.m
```

Each file sets its own `(* bench-timeout: N *)` floor (the runner uses
`max(--timeout, N)`), and each case is capped at `$cap` = 120 s by
`TimeConstrained`. A capped or unevaluated case emits only its CHECK
(`"timeout"` / `"unevaluated"` / `"failed"`), so it classifies INCOMPLETE and
never pairs a timing with a non-answer. Checks are leading significant digits
and the decimal exponent, `sig[v, d] = {leading d digits, exponent}` (the
answers run from 1e-1 to 1e36300), from `benchmarks/issac97_common.m`.

**Target: Mathilda clearly faster than the baseline column on every case**, at
machine and at arbitrary precision. Correctness comes first: a fast wrong
answer is a CHECK-FAIL and its timing is discarded.

## The problems and their answers

| exp | problem | answer |
|---|---|---|
| 126 | 2-norm condition number of the 256 × 256 Hilbert matrix | 1.7659543288e389 |
| 127 | `Integrate[x^x^x, {x, 1, 6}]` | 1.1026649993619406556e36300 |
| 128 | `Sum[1/(n^Pi + n^2 + n^Sqrt[2] + 1)^(1/3), {n, 1, Infinity}]` | 21.193240377711540943 |
| 129 | coefficient of x^3000 in (x+1)^2000 (x^2+x+1)^1000 (x^4+x^3+x^2+x+1)^500 | 3.97394226558004e1426 (exact, 1427 digits) |
| 130 | largest zero of `LaguerreL[1000, x]` | 3943.2473948452711821852759 |
| 131 | lex Gröbner basis of four dense quadrics in w, x, y, z | shape basis: degree 16 in w, x, y, z linear |
| 132 | `Limit[f[g[x]]/g[f[x]], x -> 0]`, f, g twelve-term nested Sin/Tan/Sinh/Tanh sums | 3060352455.6720 (exact rational) |
| 133 | `Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}]` | 0.21088285956559351177 (closed form) |
| 134 | largest eigenvalue of the integral operator with kernel E^(x+y+x^2+y^2+xy+x^2y^2) on [0, 1] | 37.529145560335 |
| 135 | y'' = x^3 + y^3 + y'^3, y(0) = y'(0) = 0: the blow-up point r and y(r) | r = 1.6443766903388407202, y(r) = 0.93193876511028105335 |

The machine contour row of 129 has no machine-precision answer in either
system: the integrand peaks at ~1e1430, so a system that does not escalate
returns 0 and one that does but keeps its default precision goal returns a
rounding-dominated value. It stays CHECK-FAIL until both columns escalate
fully; the true value is the one above.

## First run (2026-11-06, v0.6225, loaded host — ratios indicative)

Mathilda column at v0.6225, before any fix; the baseline column on the same host
at the same time. Three rows were measured with a since-corrected case (P5
Newton: scalar Jacobian; P7 ratio: recomputed series; P8 machine: a 12-digit
check) and read as CHECK-FAIL/SLOWER here only for that reason.

| case | Mathilda ms | baseline ms | ratio | status |
|---|---|---|---|---|
| p1 m1 eigenvalues machine n=10 | 0.1 | 0.1 | 1.23 | AHEAD |
| p1 m1 eigenvalues wp200 n=64 | 79.0 | 414.1 | 0.19 | AHEAD |
| p1 m1 eigenvalues wp800 n=256 | 16084 | 164314 | 0.10 | AHEAD |
| p1 m2 power machine n=10 | 0.8 | 1.1 | 0.75 | AHEAD |
| p1 m2 power wp800 n=256 | 11141 | 7905 | 1.41 | AHEAD (within 1.5) |
| p1 m3 LUDecomposition only wp500 n=256 | 3518 | 10738 | 0.33 | AHEAD |
| p1 m3 LUBackSubstitution wp500 n=256 |  | 12170 |  | ABSENT |
| p1 m3 LinearSolve machine n=10 | 0.3 | 0.7 | 0.50 | AHEAD |
| p1 m3 LinearSolve wp500 n=256 | 36371 | 11113 | 3.27 | CHECK-FAIL |
| p2 m1 NIntegrate machine |  | 6.7 |  | INCOMPLETE |
| p2 m1 NIntegrate DE wp30 | 0.6 | 12.0 |  | CHECK-FAIL |
| p2 m1 NIntegrate wp50 | 0.7 | 31.0 |  | CHECK-FAIL |
| p2 m1 NIntegrate wp100 | 2.0 | 117.8 |  | CHECK-FAIL |
| p2 m2 ProductLog substitution machine | 0.2 | 23.6 |  | CHECK-FAIL |
| p2 m2 ProductLog substitution wp30 | 2.3 | 88.2 |  | CHECK-FAIL |
| p3 NSum machine | 54.4 | 53.3 |  | CHECK-FAIL |
| p3 NSum wp40 | 345.4 | timeout (120 s) |  | AHEAD |
| p3 NSum Method->NIntegrate wp120 | 1872 | 2213 | 0.85 | AHEAD |
| p4 exact Coefficient[Expand] | 149.7 | 16344 | 0.01 | AHEAD |
| p4 m4 exact CoefficientList square | 30.2 | 1425 | 0.02 | AHEAD |
| p4 m3 exact Series to order 3000 | 29188 | 1631 | 17.9 | SLOWER |
| p4 m2 exact Residue | timeout | 3999 |  | INCOMPLETE |
| p4 m1 contour NIntegrate machine | 0.3 | 1222 |  | CHECK-FAIL (both wrong) |
| p4 m1 contour Trapezoidal wp400 | 2.3 | 1659 |  | CHECK-FAIL |
| p5 m1 Newton wp100 | 74.3 | (stale) |  |  |
| p5 m1 secant wp100 |  | 1108 |  | INCOMPLETE |
| p5 m1 Newton wp200 | 82.8 | (stale) |  |  |
| p5 machine FindRoot L1000 |  | (stale) |  | INCOMPLETE |
| p5 machine FindRoot L100 | 0.1 | (stale) |  |  |
| p5 m1 RootCount contour wp100 | timeout | unevaluated |  | INCOMPLETE |
| p5 m2 all 1000 roots in sequence wp100 |  | 30247 |  | INCOMPLETE |
| p5 NSolve L1000 wp50 | timeout | timeout |  | INCOMPLETE |
| p6 GroebnerBasis exact lex | 18.0 | 26.2 | 0.69 | AHEAD |
| p6 GroebnerBasis InexactNumbers machine | 18.5 | 11.2 | 1.66 | SLOWER |
| p6 NSolve machine | 16.2 | 102.0 | 0.16 | AHEAD |
| p7 m1 Series f order 25 | 2462 | 1864 | 1.32 | AHEAD (within 1.5) |
| p7 m1 Series g order 25 | 1817 | 624.4 | 2.91 | SLOWER |
| p7 m2 numeric approach wp200i | 21.5 | 43.6 | 0.49 | AHEAD |
| p7 direct Limit | timeout | timeout |  | INCOMPLETE |
| p8 exact Integrate direct / both substitutions |  | 5833–7258 |  | INCOMPLETE |
| p8 exact indefinite (both) |  | 3218 |  | INCOMPLETE |
| p8 NIntegrate wp30 | 167.0 | 155.6 | 1.07 | AHEAD (within 1.5) |
| p8 NIntegrate wp60 | 1607 | 1090 | 1.47 | AHEAD (within 1.5) |
| p9 m1 Gauss nodes Golub-Welsch wp20 | 15.3 | 74.4 | 0.21 | AHEAD |
| p9 m1 Nystrom machine | 17.8 | 16.9 | 1.06 | AHEAD (within 1.5) |
| p9 m1 Nystrom wp20 | 640.1 | 322.7 | 1.98 | SLOWER |
| p9 m1 Nystrom wp40 | 665.9 | 391.0 | 1.70 | SLOWER |
| p9 m2 operator iteration FunctionInterpolation | SKIP | timeout |  | ABSENT |
| p10 naive NDSolve machine domain end |  | 1.8 |  | INCOMPLETE |
| p10 m1 inverse function machine / wp40 |  | 1.4 / 77.2 |  | INCOMPLETE |
| p10 m2 derive polar ODE (Together/Factor) | 4.2 | 2.8 | 1.52 | SLOWER |
| p10 m2 polar machine / wp30 |  | 6.8 / 431.2 |  | INCOMPLETE |

## Deficiency register

Status: **open**, **fixed vX**. Ordered by fix priority: wrong answers, then
argument/option handling, then absences, then bottlenecks.

### Wrong answers
| # | deficiency | status |
|---|---|---|
| W1 | NIntegrate at WorkingPrecision 30–100 on x^x^x over [1, 6]: 6.60e25591 for 1.10e36300 (P2) | **fixed v0.6235** |
| W2 | NIntegrate of the ProductLog-substituted P2 integrand: 0.0 (machine), 7.64e35391 (WP 30) | **fixed v0.6235** |
| W3 | Trapezoidal NIntegrate at WP 400 on the periodic P4 contour: 1.5457e1428 for 3.9739e1426, after 2 ms (stops early) | **fixed v0.6236** |
| W4 | machine NIntegrate whose integrand leaves the double range returns unevaluated (P2) or a silent 0.0 (P4) instead of escalating | **fixed v0.6235 (escalation to MPFR)** |
| W5 | machine NSum on a series with n^-1.047 terms: 5e-5 relative error, no warning (P3) | **fixed v0.6226** |
| W7 | machine complex powers past the double range give `Overflow[]` (`(1.0 + 2.0 I)^1000`) where machine reals switch to MPFR; NIntegrate now escalates around it (P4) | open |
| W8 | the half-line Fourier probe rewrote the integrand in place (`Exp[-x^2]` -> `Power[E, -1, x^2]`): every WP integral over [0, oo) wrong (found via P2) | **fixed v0.6234** |
| W9 | MPFR endpoints taken as doubles (`2 Pi` off by 2.4e-16); one-sided tail cut before a climbing side's bulk (`(1 + Cos[x])^2000` lost half) (found via P4) | **fixed v0.6235** |
| W10 | `Limit[2 ArcTanh[1 - 2 z] - Log[1 - z], z -> 1, Direction -> "FromBelow"]` = 2 I Pi (true 0); two-sided gives Indeterminate (found via P8) | **fixed v0.6249** |
| W11 | Equal/Less of arbitrary-precision reals use a machine-sized tolerance: `N[1 - 10^-20, 30] < 1` is False and `== 1` True (found via holdout 145) | **fixed v0.6286-v0.6288** (lowest-accuracy rule; exact cancellation keeps accuracy; FindInstance witness check) |
| W12 | `N[expr, n]` with deep cancellation returned 0 or wrong digits silently (h7.4/h7.5/h7.8; two exact zeros counted as agreement) | **fixed v0.6283-v0.6284** (certified digits; N::meprec at the default cap, correct at $MaxExtraPrecision = 1000) |
| W13 | `N[0, 30]` (and inexact zero parts such as Re of `N[Erf[I], 30]`) is an inexact zero `0.e-30` where an exact 0 is expected; `N[10^-7, 12]` stores no guard bits and displays its binary value `9.99999999999e-08` (found by the printer fix) | **fixed v0.6285** (exact zeros stay exact; N[10^-7, 12] is the correctly rounded 40-bit value, displayed correctly) |
| W14 | printer follow-ups to v0.6264: a MachinePrecision number past the double range prints 15 digits with trailing zeros (`1.00000000000000e+400`, numloop_tests); a non-whole precision one digit short; MPFR -0 loses its sign | **fixed v0.6289** (digit count floor(L) - floor(L - p); a 53-bit number past the double range shows its digits; inexact zero unsigned) |
| W15 | `N[x, n]` stores exactly ceil(n log2 10) bits with no guard bits, so the last displayed digit can be wrong (`N[Sqrt[2], 50]` ends ...770, true ...769); also a timeout during MPFR's constant computation poisoned the Pi cache (fixed v0.6291) | **fixed v0.6292-v0.6293** (16 guard bits; precision = width - 16; wrong last digits 68/360 -> 0; Precision[N[Pi, 30]] = 30.) |
| W16 | `1/3^700 == 0` and `PossibleZeroQ[1/3^700]` True (exact rationals beyond the double range read as 0; made h4.9 look wrong) | **fixed v0.6296** |
| W17 | `PossibleZeroQ[Pi/3^700]` and `Pi/3^700 == 0` True: the zero test's numeric stage measures magnitudes in doubles (1e-334 underflows) | **fixed v0.6320-v0.6321** (log2 magnitudes throughout the zero test; Equal/Less of exact quantities beyond the double; also Log[3] 10^309, Exp[1 + 10^-20] - E, Pi + 10^-20 == Pi) |
| W18 | NIntegrate: a weak second peak next to a strong one (1e-6 of the mass) and an interior kink `Abs[x - 1/3]` on [0, 1] are wrong at ~1e-11 (with a convergence warning, not silent) | open |
| W6 | LinearSolve on N[Hilbert[256], 500] reports `LinearSolve::sing1` and returns garbage (P1) | **fixed v0.6241** |

### Argument and option handling
| # | deficiency | status |
|---|---|---|
| A1 | NDSolve carries `HoldAll`, so `NDSolve[ode, ...]` with the equations in a variable fails with `NDSolve::nodeqs` (P10) | **fixed v0.6232** |
| A2 | FindMinimum `{t, t0, t1}` (two starting points) read as a constraint: `FindMinimum::nimpl` (P10) | **fixed v0.6230** |
| A3 | FindRoot/FindMinimum reject Real-valued WorkingPrecision/AccuracyGoal/MaxIterations (`N[FindRoot[...]]` makes them Real): `badopt` (P5, P10) | **fixed v0.6228** |
| A4 | FindRoot rejects a 1 × 1 matrix `Jacobian -> {{g[x]}}` ("non-real derivative") (P5) | **fixed v0.6228** |
| A5 | FindRoot secant at WP 100 diverges where |f| ~ 1e867 (P5 m1, and so all of m2) | **fixed v0.6228** |
| A6 | machine FindRoot on L_1000 (|f| ~ 1e867) fails after 3.3 s instead of escalating (P5) | **fixed v0.6228** |
| A7 | FindMinimum at WorkingPrecision 30 returns a machine-precision minimum value; MPFR bracket-finding fails on an InterpolatingFunction (P10) | **fixed v0.6233 (value precision: the number carries 30 digits; the printer drops trailing zeros, see P1)** |
| A8 | `InterpolatingFunction[...]["Domain"]` (and the other property queries) unevaluated (P10) | **fixed v0.6231** |
| A10 | the printer drops trailing zeros of arbitrary-precision reals: `N[1, 30]` prints `1.0` (P10) | **fixed v0.6264** |
| A9 | `Precision[matrix]` maps over the elements instead of returning the minimum (P1) | **fixed v0.6240** |

### Absences
| # | deficiency | status |
|---|---|---|
| X1 | `Integrate` of polynomial × PolyLog[3, rational function of x], definite and indefinite (P8) | **fixed v0.6237** (0.05-0.28 s per row; v0.6238 complex values, v0.6239 FunctionExpand) |
| X2 | `Limit` of P7 (needs a series-based quotient limit with deep cancellation; the baseline also fails) | **fixed v0.6248 (~0.3 s)** |
| X3 | `LUBackSubstitution` (P1 m3) | **fixed v0.6241** |
| X4 | `GaussianQuadratureWeights` (P9 m1, verbatim) | **fixed v0.6243** |
| X5 | `FunctionInterpolation` (P9 m2) | **fixed v0.6244/v0.6268** (m2 correct; 22 s, SLOWER-by-default: the baseline times out at 120 s) |

### Bottlenecks
| # | deficiency | Mathilda | baseline | status |
|---|---|---|---|---|
| B1 | `Residue[p/x^3001, {x, 0}]` of a polynomial (P4) | 176 s | 4.0 s | **fixed v0.6247 (176 s -> ~0.24 s)** |
| B2 | `Series` of a product of polynomial powers to order 3000 (P4) | 29.2 s | 1.6 s | **fixed v0.6246 (29 s -> ~0.15 s)** |
| B3 | `Series` of nested Sin/Tan/Sinh/Tanh compositions to order 25 (P7) | 1.8–2.5 s | 0.6–1.9 s | **fixed v0.6246 (~2 s -> ~30 ms)** |
| B4 | 800-digit 256 × 256 Dot (power method) (P1 m2) | 11.1 s | 7.9 s | **fixed v0.6242 (800-digit power method 11.1 s -> ~3 s)** |
| B5 | MPFR Eigenvalues at 20–40 digits, 64 × 64 nonsymmetric (P9 Nyström) | 640 ms | 323 ms | **fixed v0.6250-v0.6253** (Nystrom wp20 3.9 s -> 0.32 s same host; wp40 0.50 s vs 0.39 s still SLOWER) |
| B6 | contour NIntegrate of L'_1000/L_1000 at 100 digits (P5 RootCount): the integrand is evaluated symbolically once, which expands LaguerreL[1000, rho E^(I phi)] into its degree-1000 polynomial (~93 ms and ~250 digits of cancellation per sample); needs sampling of the unexpanded form or precision-loss detection | > 120 s | unevaluated (same expansion) | open |
| B7 | `NSolve[LaguerreL[1000, x] == 0, x, WorkingPrecision -> 50]` | > 600 s | > 120 s | **fixed v0.6254 (> 600 s -> 4-9 s)** |
| B8 | GroebnerBasis with CoefficientDomain -> InexactNumbers (P6) | 18.5 ms | 11.2 ms | **fixed v0.6255 (70 -> 18 ms same load)** |
| B9 | NIntegrate at 60 digits on x^2 PolyLog[3, 1/(x+1)] (P8): PolyLog[3, z] near z = 1 costs ~20 ms per call at 80 digits | 1.6 s | 1.1 s | **fixed v0.6267** (WP 60 11.1 s -> 0.44 s on the same host) |
| B11 | P2 after the W1/W4 fixes: machine 71 ms, WP 30 117 ms, WP 50 805 ms, WP 100 3.1 s, ProductLog machine/WP 30 677/1406 ms (load ~100) | | 6.7 / 12 / 31 / 118 / 24 / 88 ms | **partly fixed v0.6265-v0.6271** (machine 8.6 / WP 30 9.1 ms ahead; WP 50 81, WP 100 316, ProductLog 45 / 120 ms still SLOWER) |
| B10 | the polar ODE derivation: Together/Numerator/Factor of a trig rational (P10) | 4.2 ms | 2.8 ms | open |

## Holdout (experiments 136–145) — failures still open at v0.6242

Target: every one of the 240 holdout cases correct (against the established
true answer) **and** faster in Mathilda than in the baseline column — a SLOWER
case, even within 1.5×, is open work.

The holdout ([`ISSAC97_HOLDOUT.md`](ISSAC97_HOLDOUT.md)) at v0.6225: 127 AHEAD,
49 SLOWER, 18 INCOMPLETE, 46 CHECK-FAIL (Mathilda wrong in 31). At v0.6242
(after W1–W9, A1–A9, X1, X3, B4): 152 / 42 / 16 / 30. Remaining, by area
(H = holdout case. The classifier's "who" column names the system that MATCHES the true
answer: a CHECK-FAIL marked "Mathilda" is the baseline being wrong. Rows below are
the ones where Mathilda is wrong, incomplete or slower; "neither" = both wrong):

| area | cases | status |
|---|---|---|
| NIntegrate | h2.3 machine + wp40 (70 s), h2.7 wp40 timeout, h2.8 (neither) (h2.9 was the baseline wrong) | **fixed v0.6268-v0.6269** (137: 14 AHEAD, 3 SLOWER, 0 wrong; SLOWER: h2.5 wp40, h2.7 machine and wp40) |
| NSum | h3.2, h3.4, h3.6, h3.7, h3.8 machine wrong; h3.2/3.4/3.7/3.8 wp30 timeouts | **fixed v0.6261-v0.6263** (all 20 rows correct; h3.6 wp30 still SLOWER, 3.3x by wall under load) |
| SeriesCoefficient / Residue | h4.4, h4.8, h4.9 (Tan, k = 401) wrong; h4.5 Residue timeout; h4.10 Residue wrong | **fixed v0.6245–v0.6247** (all correct; h4.8 283 ms vs 106 ms still SLOWER) |
| FindRoot / special functions | h5.2 LegendreP[1000] wp50 wrong (36 s); h5.3 JacobiP unevaluated; h5.4 Chebyshev machine; h5.10 RiemannSiegelZ absent | **fixed v0.6272-v0.6277** (140: 20/20 AHEAD; new JacobiP, RiemannSiegelZ/Theta) |
| NSolve | h6.5, h6.7 machine: solutions missing | **fixed v0.6256** (all ten systems complete; h6.5 91 ms vs 98.5) |
| Limit | h7.2, h7.10 (neither), h7.3 unevaluated, h7.3 machine point wrong | **fixed v0.6248** (h7.2 4096/729, h7.3 9765625/1048576, h7.10 1/64; the h7.3 machine point was the baseline wrong) |
| N[expr, 30] with cancellation | h7.4, h7.5, h7.8 (0 or wrong digits; neither) | **fixed v0.6283** (rows now ask with $MaxExtraPrecision = 1000; correct) |
| Integrate | h8.5, h8.8 unevaluated | **fixed v0.6278-v0.6281** (143: 30/30 correct and AHEAD; h8.5 55 ms vs 1.5 s, h8.8 99 ms vs 13.8 s) |
| eigenproblem | h9.9 sine Galerkin wp40 n = 200 (Mathilda wrong: MPFR pencils fell to LAPACK doubles) | **fixed v0.6308** (correct, 1.5 s vs 26 s) |
| NDSolve WorkingPrecision 30 | h10.1, h10.8, h10.9 were the baseline wrong (Mathilda correct); h10.7 singular endpoint (Mathilda 2.96e113) and h10.3/h10.4 Lane-Emden wp30 SLOWER | **fixed v0.6258-v0.6260** (h10.7 correct to 23 digits; Lane-Emden 1.7 s -> 30 ms, 4.3 s -> 0.11 s); h10.7's start data corrected in the benchmark |

## Run at v0.6293 (2026-11-06 evening, load ~3-5, both columns fresh)

**ISSAC '97 (126-135):** 47 rows AHEAD, 6 SLOWER, 0 Mathilda wrong. Six rows
disagree and in every one Mathilda is right (P3 wp40: baseline timeout; P4 machine
contour: baseline wrong; P5 machine L1000 and L100: baseline returns the start point;
P5 NSolve: baseline timeout; P7 direct Limit: baseline timeout). One row is
INCOMPLETE in both (P5 RootCount, B6). SLOWER: P1 machine rows at clock resolution
(1.1-1.2x), P4 Trapezoidal wp400 (454 vs 182 ms), P7 numeric approach (17 vs 6 ms,
now timed as one run), P10 Together/Factor (3.3 vs 2.8 ms).

**Holdout A (136-145):** 188 AHEAD / 38 SLOWER / 0 INCOMPLETE / 14 CHECK-FAIL
(v0.6225: 127 / 49 / 18 / 46). Of the 14, Mathilda matches the true answer in 12;
Mathilda is wrong in 2: h4.9 `SeriesCoefficient[Tan[x], {x, 0, 401}]` (correct at
v0.6249, a regression since) and h9.9 sine Galerkin at 40 digits, n = 200. The h7
"point N[.., 30]" SLOWER rows were timing the baseline's cache of N on repeats; they
are now single runs.

Wave 3 (every SLOWER case, by area): arbitrary-precision SVD/eigen (h1.2 9.8x, h1.4
5.8x) and small machine matrices; NIntegrate (h2.3/2.5/2.6/2.7, h8.10, the contour
Trapezoidal rows h4.3/4.4/4.8/4.9/4.10 and P4); exact coefficients (h4.7, h4.8) and
NSum h3.6; exact lex GroebnerBasis (h6.3/6.7/6.8/6.10, up to 6.2x); N at a point
(h7.x, re-measured as single runs), FindRoot h5.9, sine Galerkin h9.9, P7, P10.

Wave 3 so far: h1.2 / h1.4 singular values at 100 / 60 digits now 74 / 652 ms vs
91 / 695 (v0.6294); h4.7 Coefficient[Expand] 5 vs 26 ms, h4.8 12 vs 35 ms, h3.6 wp30
25 vs 33 ms (v0.6297-v0.6299); h4.9 was never wrong (W16).
NIntegrate / contour wave 3 (v0.6300-v0.6307): 14 of 15 rows now ahead (h2.3 machine
0.11 vs 9.9 ms, h2.7 0.13 / 8.9 vs 3.9 / 9.0, h4.4 machine 8.5 vs 24.5, h4.4 wp40 613 vs
1168, P4 Trapezoidal wp400 114 vs 182). Still SLOWER: h2.5 wp40 (~110 vs 15 ms; needs an
adaptive Gauss-Kronrod rule on [0, oo) at high precision); h2.6 wp40 at parity.
N at a point / eigen wave 3 (v0.6308-v0.6316): h7.3 N 286 -> 5.6 ms (vs 12), h7.4/h7.5/h7.9/h7.10
ahead; still behind at sub-ms: h7.1 N 0.22 vs 0.17, h7.6 0.17 vs 0.09; P7 m2 11.4 vs 8.0 and h9.9 machine
8.7 vs 4.6 (evaluator cost: agent); P10 Together 3.4 vs 2.8 (agent); h2.5 wp40 (agent).
Evaluator wave 3 (v0.6317-v0.6319): exact Table entries 1.5-6x faster (Kahan 6.2 -> 1.05 us); h1.4 machine
3.3 vs 4.0 ms (ahead), h1.7 even; h1.2 / h1.10 machine still behind (machine N of E^(-k/25) goes through MPFR).
Narrow features (v0.6324): h2.5 wp40 2.6 vs 15.3 ms, h2.6 2.0 vs 10.7 ms; holdout 137 all AHEAD or baseline-wrong.
Together kernels (v0.6322-v0.6323): P10 polar ODE derivation 1.7 vs 2.8 ms.
Groebner engine (v0.6326-v0.6334): F4 mod p, certified multi-modular DRL/DegLex/Lex over Q, shape-position Lex,
Lex mod p by FGLM. Holdout 141 and ISSAC 131 all AHEAD (h6.8 27-36 vs 40 ms, h6.3 56 vs 260, P6 5 vs 23); benchmark
122 families 39 of 40 ahead (cyclic-7 DRL 9.7 s vs 440 s, katsura-8 DRL 4.0 vs 23.2 s, katsura-7 Lex 6.9 vs 38.7 s);
behind: dense 4 quadrics with 6-digit coefficients in DRL (19 vs 8 ms).

## Run at v0.6334 (2026-11-07, fresh baseline for 132/142 as single runs)

**ISSAC '97 (126-135):** 51 AHEAD, 2 SLOWER, 0 Mathilda wrong; 6 rows differ, all with
the baseline wrong or timed out; 1 INCOMPLETE in both (P5 RootCount, B6). SLOWER: P1 m2
power machine 0.615 vs 0.581 ms; P7 m2 numeric approach 13.7 vs 9.5 ms.

**Holdout A (136-145):** 226 AHEAD / 2 SLOWER / 0 INCOMPLETE / 12 CHECK-FAIL, and in all 12
Mathilda matches the true answer (the baseline is wrong): Mathilda is correct on all 240
cases. SLOWER: h6.6 NSolve machine 49 vs 29 ms; h7.6 N 0.2 vs 0.1 ms.
(v0.6225: 127 / 49 / 18 / 46, Mathilda wrong in 31.)

## Sealed holdout B (146-155), run once at v0.6340 (2026-11-07)

258 new cases that no fix agent saw, classified against the true answers
([`ISSAC97_HOLDOUT_B_RESULTS_v0.6340.md`](ISSAC97_HOLDOUT_B_RESULTS_v0.6340.md)):
Mathilda correct 212/258 (baseline 242/258); where both are correct (205), Mathilda faster
in 192. Generalised well: linear algebra 24/24, polynomial systems 30/30 (baseline 28/30),
closed-form integrals 29/30. Gaps: cancellation limits 14/32 (Limit unevaluated on most new
limit types); absent ZetaZero, AiryAiZero, MathieuCharacteristicA, CountRoots (wrong on one);
13 wrong answers (NSum b3.8/b3.9/b3.10, NIntegrate b2.6/b2.8 at 30 digits, contour b4.3,
Limit b7.1 -> ComplexInfinity, NDSolve b10.6 at 50 digits, b10.8 shooting -- the last also
wrong in the baseline). Verdict: speed generalises; correctness partly overfit to the
families worked on. Next: wave 4 on these gaps (B then becomes training data) and a fresh
sealed holdout C to re-measure.
