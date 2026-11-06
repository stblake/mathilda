# ISSAC '97 holdout B (experiments 146–155) — sealed

A second, **sealed** holdout for the ISSAC '97 campaign. Experiments 126–135
time the ten challenge problems and 136–145 (holdout A) a hundred more; holdout
A's failures have been sent to fix workers, so it no longer measures
generalisation. Holdout B is 100 new problems, ten families of ten mirroring the
ten challenge problems, chosen to be **different in kind** from the originals
(other matrix families, other special functions, other kernels, other
singularity types, other closed-form constants), not other parameters of the
same formula. It is run **once, on the final code**: nobody fixing Mathilda
looks at these files or their results before then.

```bash
python3 benchmarks/issac97_holdoutb_table.py --run --save /tmp/hb   # run Mathilda on 146-155 and classify
python3 benchmarks/issac97_holdoutb_table.py --from /tmp/hb          # re-classify saved outputs
python3 benchmarks/run_all.py --only 146,147,148,149,150,151,152,153,154,155 --system mathilda   # group J in the weekly job
```

The baseline column was recorded once and is stored in
`issac97b_baseline/<experiment>.base` (the raw tagged output of each file), so
the final run needs only Mathilda. The classifier compares every CHECK with the
**true** answer, not with the other column: CORRECT / WRONG / NO ANSWER / SKIP
per system, and for cases both systems get right, whether Mathilda is faster or
slower (with the ratio), with a per-family summary.

## Format

Same harness as 126–135 (`Get["../harness.m"]; Get["../issac97_common.m"];`,
`run[label, expr, chk, reps, cap]`, `sig[v, d]` keeping the value's own
precision, a per-case `$cap` of 120 s, a per-file `(* bench-timeout: N *)`).
Labels are `b<family>.<k> <method> <machine|wpNN|exact>`. Every problem has a
machine-precision case and an arbitrary-precision case; the exact families
(1, 4, 6, 7, 8 and the exact sub-questions of 5, 10) have an exact case too.
Where a double cannot carry the full problem, the machine case is a stated
reduced instance (named in the label, e.g. `machine n=100`). Digit counts `d`
sit below what the method guarantees — machine 3–10 (each machine linear-algebra
and limit case was also run in IEEE double to confirm the margin), arbitrary
precision at most half the working precision — and one digit lower again
wherever the true value lies within 0.05 units of a rounding edge.
One-equation Jacobians would be 1 × 1 matrices; FindRoot targets built from
polynomial recurrences are `f[x_?NumericQ]` definitions.

| exp | family | cases | machine | wp | exact |
|---|---|---|---|---|---|
| 146 | 1 ill-conditioned linear algebra | 24 | 10 | 10 | 4 |
| 147 | 2 hard definite integrals | 21 | 11 | 10 | 0 |
| 148 | 3 slowly convergent series | 20 | 10 | 10 | 0 |
| 149 | 4 exact coefficients and residues | 30 | 10 | 10 | 10 |
| 150 | 5 zeros of special polynomials / functions | 25 | 10 | 13 | 2 |
| 151 | 6 polynomial systems / Groebner | 30 | 10 | 10 | 10 |
| 152 | 7 limits with deep cancellation | 32 | 10 | 10 | 12 |
| 153 | 8 closed-form definite integrals | 30 | 10 | 10 | 10 |
| 154 | 9 operator and discretised eigenproblems | 22 | 11 | 11 | 0 |
| 155 | 10 ODE singularities, high-precision ODE values | 24 | 13 | 10 | 1 |
| | **total** | **258** | **105** | **104** | **49** |

## True answers and how they were established

Neither benchmarked system was used to establish any answer. Each value was
obtained twice, by independent routes: a closed form derived by hand, exact
integer/rational arithmetic in Python (`fractions`, big integers), mpmath with
a different algorithm from the one the file times, sympy (exact Groebner bases,
Sturm counts, real-root isolation), or a published value. The file headers give
every derivation; the table below gives the answer and the two routes.

### 146 — ill-conditioned linear algebra (vs P1, Hilbert's condition number)

| # | problem | true answer | established by |
|---|---|---|---|
| 1.1 | kappa_2 of symmetric Pascal P_60 (P_10 machine) | 1.05270543880198426471473e69 (4155205697.177144) | eigsy at 200 digits; power iteration with the reciprocal-spectrum property |
| 1.2 | kappa_2 of the Lotkin matrix, n = 40 (8) | 2.370457593474166087334342e59 (23839682200.153) | mpmath SVD; eigenvalues of A^T A |
| 1.3 | top coefficient of the Vandermonde solve for Runge, n = 30 (10) | exact rational, 121163222.3328001737527216 (151.31895461171818) | Fraction elimination; mpmath LU at 300 digits |
| 1.4 | sigma_min of the Kahan matrix s = 3/5, n = 100 (12) | 1.026723597609272540268764e-47 (9.956864726388686e-6) | mpmath SVD; exact inverse + power iteration |
| 1.5 | smallest eigenvalue of the Frank matrix, n = 40 (12) | 0.007509065981719966651839415 (0.031028060644010015) | mpmath eig at 250 digits; exact characteristic polynomial roots + exact sign change |
| 1.6 | least-squares degree-20 fit to Runge at 101 points (8 / 21) | exact rational, 4806.764449233890721337973 (11.906072974689439) | exact normal equations; mpmath QR |
| 1.7 | det of the Catalan Hankel matrix, n = 40 / 25 / 7 | 1 | theorem; integer Bareiss |
| 1.8 | smallest eigenvalue of the prolate matrix w = 1/4, n = 40 (12) | 1.24903564200568087964963e-29 (1.778272859770121e-8) | eigsy; Sylvester-inertia bisection |
| 1.9 | gap of the top pair of Wilkinson W_41^+ (W_11^+) | 1.272072120490893356864979e-37 (7.428822929311607e-5) | eigsy; Sturm bisection |
| 1.10 | det of the Cauchy matrix 1/(i^2 + j), n = 30 (6) | exact rational, 2.857745163831741209353925e-655 (8.541575038817563e-21) | Cauchy determinant formula; Bareiss |

### 147 — hard definite integrals (vs P2, x^x^x)

| # | integral | true answer | established by |
|---|---|---|---|
| 2.1 | int_0^1 sin(1/x) | sin 1 - Ci(1) = 0.504067061906928371989856 | derivation; quadosc |
| 2.2 | int_0^1 x^(-999/1000) cos x | 999.760310630203964532180 | termwise series; quadrature after x = u^1000 |
| 2.3 | int_{-1}^1 dx/(10^-10 + x^2) | 2 10^5 arctan 10^5 = 314157.265358979390512931 | closed form; split quadrature |
| 2.4 | int_0^oo sin(x^3) | Gamma(4/3)/2 = 0.446489755784624605609282 | Mellin transform; quad + quadosc |
| 2.5 | int_0^50 exp(x^2) (25 machine) | 5.448684490011405885572774e1083 (5.437876228800106e269) | sqrt(pi)/2 erfi; asymptotic series |
| 2.6 | int_0^1 (e^x - 1 - x)/x^2 | 0.599620322995358659499721 | termwise series; quadrature at 120 digits |
| 2.7 | Watson's integral int_[0,pi]^3 1/(3 - cos x - cos y - cos z) | 15.67249523473857324485605 | Glasser-Zucker Gamma product; pi^3 int e^(-3t) I_0(t)^3 |
| 2.8 | int_0^4 log abs(x^2 - 2) | 3.601330233221994217662537 | antiderivative; split quadrature |
| 2.9 | int_0^1 cos(10^4 x^2) | 0.006251292347636025417761906 | Fresnel C; quadrature over 3000 half-periods |
| 2.10 | int_0^oo J_0(x) cos(x/2) | 2/sqrt 3 = 1.154700538379251529018298 | the cosine transform of J_0; quadosc |

### 148 — slowly convergent series (vs P3)

| # | series | true answer | established by |
|---|---|---|---|
| 3.1 | sum (-1)^(n+1)/log(n+1) | 0.924299897222938855959570 | Cohen-Villegas-Zagier; mpmath nsum |
| 3.2 | sum_{n>=2} 1/(n log^2 n) | 2.109742801236891974479257 | Euler-Maclaurin by hand; published (OEIS A115563) |
| 3.3 | sum sin(n)/n | (pi - 1)/2 | sawtooth Fourier series; Im(-log(1 - e^i)); Riesz means |
| 3.4 | sum log(n)/n^2 | -zeta'(2) = 0.937548254315843753702574 | zeta derivative; Euler-Maclaurin |
| 3.5 | sum H_n/n^2 | 2 zeta(3) | Euler's identity; Euler-Maclaurin |
| 3.6 | sum (-1)^(k+1)/sqrt k | (1 - sqrt2) zeta(1/2) = 0.604898643421630370247266 | eta(1/2); CVZ |
| 3.7 | sum (-1)^n log(n)/n | gamma log2 - log^2 2/2 = 0.159868903742430971756948 | eta'(1); CVZ |
| 3.8 | sum d(n)/n^2 | zeta(2)^2 = pi^4/36 | Dirichlet series; sieve to 1e6 + tail |
| 3.9 | sum phi(n)/n^3 | zeta(2)/zeta(3) = 1.368432777620205875736766 | Dirichlet series; sieve to 1e6 + tail |
| 3.10 | sum J_0(n) | 1/2 | Poisson summation; Riesz means to 1e6 |

### 149 — exact coefficients and residues (vs P4)

| # | quantity | true answer | established by |
|---|---|---|---|
| 4.1 | [x^1000] prod 1/(1 - x^k) = p(1000) (p(100) machine) | 24061467864032622473692149727991 | pentagonal recurrence; Rademacher |
| 4.2 | [x^400] sec x (order 20 machine) | abs(E_400)/400! = 4.53899084277566e-79 | Entringer integers; mpmath eulernum |
| 4.3 | Res_{z=-20} Gamma(z)^3 | (1/20!^3)((9/2)(H_20 - gamma)^2 + pi^2/4 + (3/2)H_20^(2)) = 3.18864502745108196e-54 | derivation; trapezoidal contour at 80 digits |
| 4.4 | [x^10000] 1/prod(1 - x^c), c in {1,5,10,25,50,100} (x^100 machine) | 139946140451 (293) | two DPs; contour at 80 digits |
| 4.5 | Motzkin M_1000 (M_30 machine) | 6.11327659767719e472 (473 digits) | recurrence; Catalan sum |
| 4.6 | Ramanujan tau(1000) | -30328412970240000 | power series; Hecke relations |
| 4.7 | Res_{z=0} z^20 exp(z + 1/z) (z^3 machine) | I_21(2) = 2.04822323834388e-20 (I_4(2)) | Bessel generating function; series |
| 4.8 | [x^1000] T_2000 (T_40 machine) | 714-digit integer, 7.00127388815131e713 | closed form; integer recurrence |
| 4.9 | Bell number B_500 (B_30 machine) | 1.60607260103999e843 | Bell triangle; Dobinski |
| 4.10 | Res_{z=0} pi cot(pi z)/z^40 (/z^10 machine) | -2 zeta(40) | Laurent series; contour |

### 150 — zeros of special polynomials and functions (vs P5)

| # | zero | true answer | established by |
|---|---|---|---|
| 5.1 | largest zero of H_1000 (H_100 machine) | 44.20915249799639770159036338329 | Jacobi-matrix Sturm bisection; Newton |
| 5.2 | smallest positive zero of P_1000 | 0.001570010480083193829005023042 | Sturm bisection; Newton |
| 5.3 | j_{0,1000} | 3140.807295225078628895545453471 | besseljzero; McMahon |
| 5.4 | gamma_1000 | 1419.4224809459956864659890380799 | zetazero; Odlyzko's table |
| 5.5 | a_500 (Airy) | -177.00925104647920682546615344202 | airyaizero; asymptotic series |
| 5.6 | Gamma(x) = 10^100 | 70.957574457353546136215496627741 | Newton on LogGamma; bisection |
| 5.7 | Wilkinson's perturbed polynomial: real-root count, largest real root | 10; 20.846908101482256914928772892631 | sympy Sturm/isolation; polyroots |
| 5.8 | largest zero of P^(3/2,-1/2)_400 | 0.99993721886612628233934595246627 | Sturm bisection; Newton |
| 5.9 | Mignotte x^40 - 2(10x - 1)^2: roots in (0, 1/5), their gap (n = 8 machine) | 2; 1.41421356237309504880e-21 (1.4142137179366e-5) | fixed-point form x = (1 +- x^20/sqrt2)/10; Newton; sympy isolation (n = 8) |
| 5.10 | first zero of Erf in Q1 | 1.45061616324367559142 + 1.88094300015331537194 I | secant; Newton from another start; Fettis-Caslin-Cramer table |

### 151 — polynomial systems (vs P6)

Katsura-5 (32 solutions), three dense cubics (27), Lagrange critical points
(14, max multiplier 3/2 confirmed by direct maximisation), Arnborg-Lazard (20),
an inconsistent system (basis {1}), implicitisation of (t^3 - 3t, t^4 - t^2 + t)
(equal to the resultant) and its intersection with a circle (8), 20-digit
coefficients (36; the eliminant divides the resultant), power sums (24
solutions, eliminant of degree 4 = (t + 1)(t^3 - 2t^2 + t - 1), by hand), the
non-radical Cox-Little-O'Shea system (8 solutions: 2 with y = 0 plus 6 from
(3z^3 + 40)^2 = 96 z^5, by hand), Trinks (10). Exact invariants from sympy's
Buchberger/FGLM lex bases; solution statistics from mpmath polyroots of the
eliminant at 60 digits.

### 152 — limits with deep cancellation (vs P7)

zeta(4)/4 from the LogGamma series; Gruntz's -e^2; Arnold's 1; the Stieltjes
constant -gamma_1 = 0.07281584548367672486; 1/3; 1 (x^(x^x) at 0+); -1/12;
-13/5760 (two Bessel compositions agreeing to x^5); -15/8 (erfc asymptotics);
30^(1/3) Var(log 2, log 3, log 5)/2 = 0.2183579121076573342. Each derived by
hand and confirmed by evaluation at 10^-10 ... 10^-100 (or 10^10 ... 10^30) in
mpmath at up to 9000 digits, and by sympy series where sympy can expand.

### 153 — closed-form definite integrals (vs P8)

pi log2/8 (Serret), -pi^4/45, 2G, -sqrt(pi)(gamma + 2 log 2)/4, -pi^4/72,
-pi^2/18, -3 zeta(3)/4, pi G/2 - 7 zeta(3)/8, pi^2/4, (gamma + log 2pi)/(2pi)
(Kummer). Each derived by hand (derivations in the file header) and confirmed
by tanh-sinh quadrature at 40 digits.

### 154 — operator and discretised eigenproblems (vs P9)

e^(-abs(x-y)) on [0, 1] (secular equation, 0.738810809416454997194); the Mehler
kernel on R (1/32); the Poisson kernel on the circle (0.729 up to the exact
circulant correction); min(i, j) (closed form); the Kac matrix (exactly n);
Mathieu a_0(25) = -40.25677954656678727580 (Fourier matrix; continued
fraction); the quartic oscillator E_0 = 1.06036209048418289965 (published;
harmonic basis); the sinc kernel's lambda_3(4) = 0.11021098701480231197; the
2-D Laplacian and a P1 finite-element pencil (closed forms).

### 155 — ODE singularities and high-precision values (vs P10)

Riccati pole sqrt(2 j_{-1/4,1}) = 2.00314735942688470800; pendulum period
4 K(sin^2(179 deg/2)) = 24.51111529810688044083; Lorenz x(20) = 13.79319959512861883160 and x(5) (Taylor
integration at two precisions and two methods); Lane-Emden xi_1(3) =
6.89684861937696037545 (published); Blasius f''(0) = 0.33205733621519629894
(published); Ai(8) and Ai(4) by unstable forward integration; the double pole
of y'' = y^2 at sqrt(3/2) B(1/6, 1/2)/3 = 2.97447742540217556147; the
Thomas-Fermi slope -1.58807102261137531272 (published); the log singularity
pi/sqrt2 of y'' = e^y; Kepler's equation for the e = 0.99 orbit.

## Baseline column

Recorded 2026-11-06, each file run sequentially on an otherwise lightly loaded
host; raw outputs in `issac97b_baseline/`. Against the true answers the
baseline column is CORRECT on **243 of 258** cases:

| family | cases | baseline correct |
|---|---|---|
| 1 ill-conditioned linear algebra | 24 | 24 |
| 2 hard definite integrals | 21 | 16 |
| 3 slowly convergent series | 20 | 16 |
| 4 exact coefficients and residues | 30 | 27 |
| 5 zeros of special functions | 25 | 25 |
| 6 polynomial systems | 30 | 28 |
| 7 cancellation limits | 32 | 32 |
| 8 closed-form definite integrals | 30 | 29 |
| 9 operator eigenvalues | 22 | 22 |
| 10 ODE singularities | 24 | 24 |

The 15 it does not answer correctly (each true answer re-checked by both of
its independent routes):

| case | baseline | true |
|---|---|---|
| b2.2 NIntegrate machine / wp30 | 539.2 / 999.7555578 | 999.7603106 |
| b2.9 NIntegrate machine / wp30 | 0.0048082 / 0.00480819728 | 0.00625129235 |
| b2.7 NIntegrate wp20 | timeout | 15.67249523 |
| b3.8 NSum machine / wp30 (DivisorSigma terms) | 2.48120 / 2.48120186455 | 2.70580808428 |
| b3.9 NSum machine / wp30 (EulerPhi terms) | 1.34158 / 1.34157618212 | 1.36843277762 |
| b4.1 contour wp60, b4.6 contour wp40 | timeout | p(1000), tau(1000) |
| b4.7 Residue exact (essential singularity) | unevaluated | I_21(2) |
| b6.1 NSolve machine (Katsura-5) | 30 solutions | 32 |
| b6.10 NSolve machine (Trinks) | 9 solutions | 10 |
| b8.10 Integrate exact | unevaluated | (gamma + log 2pi)/(2 pi) |

b10.8 (truth corrected 2026-11-07): the case's true answer was first given as the
published slope -1.588071023, which the shooting method reaches only from exact
start values.  The program starts at x0 = 1e-4 from the two-term series
y = 1 - B x + (4/3) x^(3/2), y' = -B + 2 x^(1/2), which omits the next terms
(2/5) b x^(5/2) of y and b x^(3/2) of y' (b = -B): the start slope is off by
B x0^(3/2) = 1.6e-6, and the criterion's root moves by the same amount.  A
Taylor-series integrator in s = sqrt(x) (Y' = 2 s P, P' = 2 Y^(3/2), order 25,
30 digits, independent of both systems) running the program's own 34-step
bisection gives, with the program's start values, -1.58807260308 (X = 40) and
-1.58807260058 (X = 60); with the four-term start series, -1.58807134 (X = 20),
-1.588071025 (X = 40) and -1.5880710227 (X = 60), the published value, which
also reproduces the IEEE-double DOP853 figures above.  So the true answer of the
case as programmed is -1.58807260 (9 digits); both systems return it, the
baseline at every X (also X = 100) because the bias is in the start values, not
in X.  The machine row (X = 40, 6 digits) is unaffected.
