# ISSAC '97 holdout C (experiments 156-165) -- SEALED

The third holdout for the ISSAC '97 Systems Challenge (126-135).  Holdouts A
(136-145) and B (146-155) were both used for fixing; C was authored afterwards,
by someone who had not read either, and is meant to be run **once**, on the
final code, as the fresh measurement of overfitting.  Do not fix against it.

* 100 problems, ten families of ten, one experiment per family, each family the
  counterpart of one challenge problem but different in kind from the challenge
  instances, from the textbook examples, and from the structures named in
  holdout B's description (Hilbert/Pascal/Lotkin/Kahan/Frank/prolate/Wilkinson/
  Cauchy/Catalan-Hankel matrices, Riccati/Lorenz/Blasius/Thomas-Fermi/Kepler
  ODEs, Mathieu/quartic/FEM eigenproblems, ...).
* 249 cases: **100 machine**, **100 arbitrary precision**, **49 exact**.  Every
  problem has a machine row and an arbitrary-precision row; the exact families
  (c1, c4, c6, c7, c8) have an exact row too.  Where a double cannot carry the
  problem the machine row is a stated reduced instance.
* Every true answer was established independently of both benchmarked systems,
  by two methods (A and B below) wherever possible: a closed form derived by
  hand, exact integer/Fraction arithmetic in Python, mpmath with a different
  algorithm, PARI/GP, Maxima's Buchberger, scipy at rtol 1e-13, or a published
  value.  The other CAS's run was used only as a check on the harness, never as
  a truth.  Digits checked are at or below what the stated method guarantees:
  machine 5-8, arbitrary precision 15-25 (never more than half the working
  precision), exact 15 leading digits plus residues mod 10^9+7.

Run and classify:

    python3 benchmarks/issac97_holdoutc_table.py --run --mathilda ./Mathilda \
        --save /tmp/hc --md /tmp/hc/HOLDOUT_C.md        # run Mathilda, classify
    python3 benchmarks/issac97_holdoutc_table.py --from /tmp/hc -v
    python3 benchmarks/issac97_holdoutc_table.py --baseline-only

The classifier checks each CHECK value against the true answer (one unit of
tolerance in the last checked digit; exact parts must match), reports
correct / wrong / no answer per system and family, and compares timings only
where both columns are correct.  The baseline column is
`issac97c_baseline/r<NNN>.base`, one run of the other CAS on this machine.
`run_all.py` files the ten experiments under group **K ISSAC'97 holdout C
(sealed)**.

The full statement of every problem, its derivation and both methods are in the
header of each `.m` file; the tables below are the summary.

## c1 -- ill-conditioned linear algebra (156, counterpart of P1)

| # | problem | exact | machine | arbitrary precision |
|---|---|---|---|---|
| 1 | Vandermonde i^(j-1) | Det n=30 = prod k! (2.1555e383) | kappa n=8 = 952111739.18 | kappa n=60 @300 = 6.50027973410168833640e123 |
| 2 | Gaussian Toeplitz q^((i-j)^2) | Det n=40, q=9/10 = prod (1-q^2k)^(n-k) | lambda_min n=30 = 4.6374958257e-9 | lambda_min n=120, q=19/20 @60 = 3.67554667929975235737e-20 |
| 3 | Moler L L^T | (M^-1)_11 n=30 = 1+(4^29-1)/3 | kappa n=12 = 87522898.094 | kappa n=150 @250 = 2.02335096218540765444e93 |
| 4 | grounded path Laplacian, conductances 2^-k | LinearSolve n=60: x_j = 2^(j+1)-2 | kappa n=20 = 2682440.2337 | kappa n=60 @60 = 2949381270885395289.11 |
| 5 | circulant (1+eps) I - P | -- | kappa n=32 = 2000001 | kappa n=48, eps=1e-30 @60 = 2e30+1 |
| 6 | KMS r^abs(i-j) | Det n=50, r=999/1000 = (1-r^2)^49 | lambda_min = 5.0074416717e-4 | lambda_min r=1-1e-12 @60 = 5.00493805098963979720e-13 |
| 7 | Hankel of Bernoulli numbers | Det n=14 (Al-Salam-Carlitz) = -3.1050410217561563 | kappa n=8 = 7885.9266359 | kappa n=20 @60 = 224185454513517586077.46 |
| 8 | Lauchli least squares | x_20, eps=1e-10 (Sherman-Morrison) | x_20, eps=1e-6 = 9500000.050000525 | eps=1e-30 @80 = 9.5e30 + 0.05 |
| 9 | Hankel of factorials | Det n=30 = prod k!^2 | kappa n=6 = 350717286.23 | kappa n=60 @400 = 3.14466193056846718046e205 |
| 10 | Dorr matrix, D x = 1 | x_25 n=50 = 5868.3713440195491651 | x_30 n=60 = 8081580.4333 | x_50 n=100 @60 = 2261277808846905453744539.12 |

A: mpmath svd_r / eigsy / lu_solve at >= twice the digits lost, or Fraction
arithmetic.  B: the closed form, or PARI qfjacobi of the exact matrix; for the
three largest cases qfjacobi of the matrix and of its exact inverse.

## c2 -- hard definite integrals, numerically (157, counterpart of P2)

| # | integral | true value | difficulty |
|---|---|---|---|
| 1 | int_0^oo x J0(x) J0(2x)^2 | 2/(Pi Sqrt 15) = 0.164374518416399948 | Sonine-Gegenbauer; decays like x^-1/2 |
| 2 | int_0^1 Sin[1/x]/x | Pi/2 - Si(1) = 0.624713256427713604 | infinite oscillation at 0 |
| 3 | int_-1^1 Sqrt[1-x^2]/(x-z), z = 1/2 + i delta | -Pi (z - Sqrt[z-1] Sqrt[z+1]) | pole delta = 1e-3 / 1e-10 from the path |
| 4 | int_0^oo x^(a-1) E^-x/(1+x), a = 1/1000 | Gamma(a) E Gamma(1-a,1) = 998.828034469465152 | x^-0.999 endpoint |
| 5 | int_0^oo Cos[t^3/3 - 2t] | Pi Ai(-2) = 0.714421505810163771 | accelerating phase |
| 6 | int_0^2 1/Sqrt[abs(x^2-2)] | Pi/2 + Log[1+Sqrt 2] = 2.452169913814439644 | unannounced interior singularity |
| 7 | int_0^oo J0(x) J0(2x) | K(1/4)/Pi = 0.536591003574682188 | Weber-Schafheitlin, ~1/x decay |
| 8 | int_0^Pi Log[1 - 2r Cos t + r^2] | 2 Pi Log r | r = 1+1e-4 / 1+1e-12: O(1) values cancel to 1e-12 |
| 9 | int int_[0,1]^2 -Log[xy]/(1-xy) | 2 Zeta(3) = 2.404113806319188571 | 2-D corner singularity |
| 10 | int_0^oo prod_{k=0}^7 Sinc[x/(2k+1)] | Pi/2 - 6879714958723010531/935615849440640907310521750000 Pi | Borwein: differs from Pi/2 at 2.3e-11 |

A: the closed form at 60 digits; B: an independent mpmath quadrature
(quadosc, split tanh-sinh, or a Levin-accelerated sum of period integrals),
agreeing to >= 25 digits (17 for c2.10, whose A is the published rational).

## c3 -- slowly convergent series (158, counterpart of P3)

| # | series | true value | B |
|---|---|---|---|
| 1 | sum over primes 1/p^2 | 0.452247420041065498506543 | PARI Moebius log-zeta |
| 2 | sum LiouvilleLambda[n]/n^2 | Pi^2/15 | PARI partial sums to 1e6 |
| 3 | sum (-1)^(n+1) H_n/n | Pi^2/12 - Log[2]^2/2 | mpmath nsum |
| 4 | Lambert sum q^n/(1-q^n), q = 0.999 | 7480.97786909167252765 | PARI suminf |
| 5 | sum_{n>=2} 1/(n Log[n]^2) | 2.10974280123689197448 | explicit Euler-Maclaurin vs OEIS A115563 |
| 6 | sum_{n>=2} (-1)^n/Log[n] | 0.924299897222938855960 | PARI sumalt |
| 7 | sum Cos[n]/Sqrt[n] | Re Li_1/2(E^I) = -0.194108935092182049739 | PARI, Jonquiere via Hurwitz zeta |
| 8 | sum H_n^2/n^2 | 17 Pi^4/360 | partial sum + asymptotic tail |
| 9 | sum n^(-1.01) | Zeta(1.01) = 100.577943338496872490 | PARI zeta |
| 10 | sum Log[n]/n^2 | -Zeta'(2) = 0.937548254315843753703 | PARI zeta' |

Machine rows are NSum (or, for the arithmetic series, a direct sum whose
truncation error is below the checked digits); wp rows are NSum at 40 digits,
or the identity that makes the series computable.

## c4 -- exact coefficients and residues (159, counterpart of P4)

| # | quantity | exact true value | machine | wp40 |
|---|---|---|---|---|
| 1 | overpartitions pbar(300) | 175155936072122455560 | contour r = 0.913 | trapezoid |
| 2 | large Schroeder S(1000) | 9.10543459812736067919e760 | S(150) (reduced) | S(1000) |
| 3 | central Delannoy D(2000) | 1.62223628660764197568e1529 | D(150) | D(500) |
| 4 | Res Gamma(z)^2 at -10 | 2 PolyGamma(11)/(10!)^2 = 3.57186931397414815388e-13 | circle 1/2 | circle 1/2 |
| 5 | Res Zeta(s)^3 at 1 | 3 gamma^2 - 3 gamma_1 = 1.21798130787418619754 | circle 1/2 | circle 1/2 |
| 6 | r_8(1000) = [q^1000] theta_3^8 | 18365675328 (divisor formula) | r = 0.99 | r = 0.99 |
| 7 | 300! [x^300] Exp[x/(1-x)] | 7.99540643814440858751e626 | n = 60 | n = 300 |
| 8 | Res 1/(z^15 Sin z Sinh z) at 0 | 13739/1136785104000 | circle 5/2 | circle 5/2 |
| 9 | 500! [x^500] Exp[x + x^2/2] | 4.23107565308608549952e575 | n = 100 | n = 500 |
| 10 | Res z^100/((z-1)^30 (z+1)^20) at 1 | 7693711204564513690663977952285/35184372088832 | circle 1 | circle 1 |

A: exact Python (series products, P-recursive recurrences, binomial sums,
Fractions); B: PARI polcoeff, a second closed form, and the mpmath trapezoid
rule on the stated contour (which also measured each contour's cancellation).

## c5 -- zeros of special functions and polynomials (160, counterpart of P5)

| # | equation | root |
|---|---|---|
| 1 | w E^w = -1/5, branch -1 | -2.54264135777352642429 |
| 2 | w E^w = 1, branch 2 | -2.40158510486800288417 + 10.77629951611507089850 I |
| 3 | StruveH[0, x] = 0 | 4.33323782040642167053 |
| 4 | 3F2(-30, 31, 1/2; 1, 1; x) = 2F1(-15, 31/2; 1; x)^2 (Clausen): double roots | 0.00620338140950646343366 |
| 5 | L(1/2 + I t, chi_-4) = 0 | t = 6.02094890469759665490 |
| 6 | PolyGamma[0, x] = 0 in (-3, -2) | -2.61072086844414465000 |
| 7 | KelvinBer[0, x] = 0 | 2.84891782079513985275 |
| 8 | EllipticK[m] = 2 EllipticE[m] | 0.826114765984970336177 |
| 9 | sum_{k<=40} z^k/k! = 0, least modulus | -12.07594063041221103760 + 0.71848257151412454217 I |
| 10 | Zeta'(s) = 0 in (-3, -2) | -2.71726282920457410157 |

A: mpmath findroot / polyroots; B: PARI (lambertw with branch, solve,
lfunzeros, ellK/ellE, polroots) or, for Struve and Kelvin, the defining power
series summed directly.

## c6 -- polynomial systems (161, counterpart of P6)

| # | system | GB length / eliminant degree | solutions | max Re (last variable) |
|---|---|---|---|---|
| 1 | 3-player 2-strategy Nash equilibria | 3 / 2 | 2 | 2/3 |
| 2 | mass-action equilibrium with conservation | 4 / 4 | 4 | 1.65840319545635272232 |
| 3 | planar 3-RPR forward kinematics | 4 / 4 | 4 | 0.984693406500657076633 |
| 4 | repressilator steady states | 3 / 9 | 9 | 1.51598022769282058968 |
| 5 | two dense quartics | 2 / 16 | 16 | 1.79759338472934391992 |
| 6 | Lagrange, max x+2y+3z on x^4+y^4+z^4 = 1 | 4 / 36 | 36 | (1+2^(4/3)+3^(4/3))^(3/4)/4 = 1.17206247758678450467 |
| 7 | Noonburg network (not shape position) | 4 / 14 | 21 | 1.68372096585234179108 |
| 8 | normals from (3,3) to a cubic | 2 / 7 | 7 | 3.23552961063904388209 |
| 9 | two coupled Brusselator cells | 4 / 5 | 5 | 6.27469440017220052664 |
| 10 | 4-player 2-strategy Nash equilibria | 4 / 9 | 9 | 3.93580934808760382950 |

A: sympy lex Groebner basis, eliminant roots at 80 digits, every solution
back-substituted (residual <= 2e-74); B: Maxima's poly_reduced_grobner (same
length, same eliminant up to a scalar).  c6.7's count is the quotient-ring
dimension, radical by a generic-form eliminant of degree 21, squarefree.

## c7 -- limits with deep cancellation (162, counterpart of P7)

| # | expression | limit |
|---|---|---|
| 1 | x^2 (Gamma[x+1/2]/(Sqrt[x] Gamma[x]) - 1 + 1/(8x)), x -> oo | 1/128 |
| 2 | (Zeta[s] - 1/(s-1) - EulerGamma)/(s-1), s -> 1 | -StieltjesGamma[1] |
| 3 | x^2 Log[x] (Log[Log[x+1]] - Log[Log[x]] - 1/(x Log[x])), x -> oo | -1/2 (correction 1/(2 Log x)) |
| 4 | Sqrt[x] (Sqrt[x + Sqrt[x + Sqrt[x]]] - Sqrt[x] - 1/2), x -> oo | 1/8 |
| 5 | (Erf[Sin[x]] - Sin[Erf[x]])/x^3, x -> 0 | (4 - Pi)/(3 Pi^(3/2)) |
| 6 | n^2 (H_n - Log[n] - EulerGamma - 1/(2n)), n -> oo | -1/12 |
| 7 | (x^(x^(1/x)) - x)/Log[x]^2, x -> oo | 1 |
| 8 | (Gamma[x] - 1/x + EulerGamma)/x, x -> 0 | (EulerGamma^2 + Pi^2/6)/2 |
| 9 | (x^(x^x) - x)/(x^2 Log[x]^2), x -> 0+ | 1 |
| 10 | (PolyLog[2,x] - Pi^2/6)/(1-x) - Log[1-x], x -> 1- | -1 |

Rows: the exact Limit; N at a machine point P1 (cancellation <= 1e7); N[., 30]
at a point P2 where 20-60 digits cancel (checked against the value AT P2).
A: the asymptotic expansion by hand / sympy series; B: mpmath at 300 digits
30-60 orders closer.  Point values: mpmath at 150 and 250 digits.

## c8 -- closed-form definite integrals (163, counterpart of P8)

| # | integral | closed form |
|---|---|---|
| 1 | int_0^(Pi/3) Log[2 Sin[x/2]] | -Cl_2(Pi/3) |
| 2 | int_0^(1/2) LogGamma[x] | (5/24) Log 2 + (3/2) Log A + (1/4) Log Pi (A = Glaisher) |
| 3 | int_0^(Pi/2) x Log[Sin x] | 7 Zeta(3)/16 - (Pi^2/8) Log 2 |
| 4 | int_0^oo ArcTan[t]/(E^(2 Pi t) - 1) | 1/2 - Log[2 Pi]/4 (Binet) |
| 5 | int_0^oo J0(x) Sin[x/2]/x | Pi/6 |
| 6 | int_(Pi/4)^(Pi/2) Log[Log[Tan x]] | (Pi/2) Log[Gamma(3/4) Sqrt(2 Pi)/Gamma(1/4)] (Vardi) |
| 7 | int_0^1 x LogGamma[x] | Log[2 Pi]/4 - Log A (identified by PSLQ) |
| 8 | int_0^1 Log[x] Log[1-x]^2/x | -Pi^4/180 |
| 9 | int_0^1 x^(-1/3) (1-x)^(-2/3) Log[x] | Pi^2/3 - Sqrt[3] Pi Log 3 |
| 10 | int_0^1 ArcTan[x] Log[x]/x | -Pi^3/32 |

A: the closed form at 50 digits; B: mpmath tanh-sinh at 45 digits (quadosc
for c8.5), agreeing to >= 1e-31.

## c9 -- eigenproblems (164, counterpart of P9)

| # | problem | eigenvalue |
|---|---|---|
| 1 | tridiag(eps, 1, 1), non-normal | 1 + 2 Sqrt[eps] Cos[Pi/(n+1)] (n=60, eps=1e-8: 1.000199734817976966) |
| 2 | Bessel operator of order 1/3, singular at 0 | BesselJZero[1/3,1]^2 = 8.425006929499198575 |
| 3 | Love's kernel 1/(Pi (1+(x-y)^2)) on [-1,1] | 0.453572864122606855682 |
| 4 | PT-symmetric p^2 + i x^3 | 1.156267071988113256 (Bender-Boettcher 1.1562670719881132) |
| 5 | lollipop graph K_10 + P_10, algebraic connectivity | 0.0379599847847906607432 |
| 6 | -u'' + abs(x) u on the line | -a'_1 = 1.018792971647471089 |
| 7 | e^(xy) on [0,1] (degenerate kernel) | 1.353030164745734853309 |
| 8 | Toeplitz T_(i,i-1) = 1, T_(i,i+2) = 1/4, non-normal | rho(n=12) = 1.124913232291133132, rho(n=40) = 1.183677632402317698 |
| 9 | Jordan block closed by eps | eps^(1/n) |
| 10 | radial problem on the annulus 1 < r < 2 | k^2, J0(k)Y0(2k) = J0(2k)Y0(k): 9.753322124750714911 |

A/B per problem: closed forms; mpmath vs PARI (besseljzero vs solve,
eigsy vs polroots of the exact characteristic polynomial, airyaizero vs
derivnum); mpmath Nystrom at n = 40/50/60; oscillator basis at N = 80/100/120
plus the published value; degenerate kernel vs Nystrom.

## c10 -- ODE singularities and high-precision values (165, counterpart of P10)

| # | problem | true value |
|---|---|---|
| 1 | Painleve I y'' = 6y^2 + x, y(0)=y'(0)=0: first pole | 2.61557120988237384515 |
| 2 | Painleve II y'' = 2y^3 + xy, y(0)=1, y'(0)=0: first pole | 1.26277886619324886643 |
| 3 | Emden-Fowler y'' = x y^2, y(0)=1, y'(0)=0: blow-up | 2.84089482764194948767 |
| 4 | Van der Pol mu = 50 from (2, 0): first zero | 40.94520292840095473813 |
| 5 | eps y'' + (1+x) y' + y = 0, y(0)=0, y(1)=1: y(1/2) | eps=1/100: 1.33597477915866113865; eps=1/1000: 1.33359306918159219891 |
| 6 | Rossler (a=b=1/5, c=57/10) from (1,1,1) | x(50) = -7.06153226337260248572; x(100) = 9.65823343087484059219 |
| 7 | Ueda x'' + x'/20 + x^3 = (15/2) Cos t from (3, 4) | x(20) = -0.72512057794546881170; x(40) = 0.59549188064240961447 |
| 8 | Kummer from next to its regular singular point | M(1/3, 1/2, 5) = 77.41577886044270771806 |
| 9 | Painleve IV (alpha = beta = 0), y(0)=1, y'(0)=0 | y(1/2) = 1.33910862429232758493 |
| 10 | Hastings-McLeod y'' = xy + 2y^3 from the Airy asymptote | y(0) = 0.36706155154807842775 |

Poles are located from NDSolve's domain end with the Laurent correction
p = x1 + k y/y' (k = 2 double pole, 1 simple pole).  A: mpmath odefun (Taylor
series, arbitrary precision) at two precisions or two offsets; B: scipy
DOP853/Radau/LSODA at rtol 1e-12..1e-13, a closed form (c10.5 Erfi, c10.8
hypergeometric), or the published value (c10.10).
