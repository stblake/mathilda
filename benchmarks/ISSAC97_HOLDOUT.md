# The ISSAC '97 holdout (experiments 136–145)

An anti-overfitting holdout for the ISSAC '97 Systems Challenge (experiments
126–135, [`ISSAC97.md`](ISSAC97.md)). Fixes are made against those ten problems;
these 100 problems probe the same capabilities on **different instances**, so a
fix tuned to a challenge input — rather than to the capability — shows up here
as a failure. Ten families of ten, one experiment per family, each family
mirroring one challenge problem:

| exp | family | mirrors | what varies |
|---|---|---|---|
| 136 | ill-conditioned matrices | 126 Hilbert condition number | Pascal, Vandermonde, Cauchy, Kahan, Frank, Wilkinson, prolate, Moler, KMS, Gaussian-kernel matrices |
| 137 | hard definite integrals | 127 ∫ x^x^x | towers, E^E^E^x, Sin[1/x], x^-x, tiny/huge values, boundary layers, unbounded oscillation |
| 138 | slowly convergent series | 128 terms ~ n^-1.047 | p barely > 1, log factors, alternating and trigonometric conditionally convergent sums |
| 139 | exact coefficients | 129 middle coefficient | central trinomial, product powers, rational and algebraic series, high-order residues, Tan |
| 140 | special-function roots | 130 largest zero of L₁₀₀₀ | Hermite, Legendre, Jacobi, Chebyshev, Laguerre, Bessel, Airy, Riemann–Siegel Z |
| 141 | polynomial systems | 131 lex Gröbner basis | new dense systems (up to 5 variables / 56 solutions), grevlex, Katsura-5, cyclic-5 |
| 142 | cancellation limits | 132 nested-trig limit | inverse-function compositions, deeper nestings, exponential cancellation at ∞ |
| 143 | polylog integrals | 133 x² PolyLog[3, 1/(1+x)] | Log^k, PolyLog[2,3], ArcTan × Log, Catalan, Li₄(1/2) |
| 144 | operator eigenvalues | 134 integral-operator eigenvalue | six smooth kernels, a rank-2 kernel, Mathieu, anharmonic oscillator, Sturm–Liouville |
| 145 | ODE singularities | 135 ODE blow-up | Riccati blow-ups, Lane–Emden, Bratu, Blasius, Gauss at x = 1, a square-root branch point |

```bash
python3 benchmarks/run_all.py --only 136,137,138,139,140,141,142,143,144,145 -v
cd benchmarks/139-issac97h-coefficients && ../../Mathilda -file coefficients.m
python3 benchmarks/issac97_holdout_table.py DIR -v   # DIR/rN.mth, DIR/rN.base: classify vs the truth
```

**Rules.** Every problem has a machine-precision case and an arbitrary-precision
case; 139, 141, 142 and 143 also have an exact case. Where a double cannot carry
the full problem (an answer or integrand beyond 1e308, a condition number beyond
1e16, cancellation deeper than 16 digits), the machine case is a **reduced
instance** a double resolves, stated per problem. Each problem is at least as
hard as its challenge counterpart in size, conditioning, digits, degree or
cancellation depth. Checks are `sig[v, d]` (leading digits plus exponent, from
`issac97_common.m`) with `d` at or below what the method guarantees (machine
default-goal methods 4–9 digits; arbitrary precision at or below half the
working precision). Every case is capped at 120 s.

**Truth, not agreement.** Each true answer was established independently of both
benchmarked systems, by two methods wherever possible (columns *provenance*):
closed forms derived by hand, pure integer arithmetic in Python, high-precision
recomputation in mpmath by a different algorithm, Sturm-sequence counts, or a
published value. `issac97_holdout_table.py` compares each system with the true
answer, so a CHECK-FAIL names the system that is wrong. One answer first
obtained by a root finder (two Legendre/Jacobi zeros) was wrong in the 8th digit
— the finder had stalled — and was caught by the Sturm count; both are corrected
below, which is why every value here has a second method.

## Problems and true answers

### 136 — ill-conditioned matrices (machine 4 digits, arbitrary precision 12)
Two algorithms at two precisions each: eigenvalues (or singular values) at 400
digits; singular values / eigenvalues of the inverse or of AᵀA / roots of the
exact characteristic polynomial at 600 digits. All agree to ≥ 25 digits.

| # | problem | true value | wp case | machine instance |
|---|---|---|---|---|
| h1.1 | κ₂ of the 100 × 100 symmetric Pascal matrix | 9.191581826323644019e116 | 250 digits | n = 12: 8.763949109e11 |
| h1.2 | σ_min, 50 × 50 Vandermonde on (i−1)/49 | 1.214052974061843950e-42 | 100 | n = 12: 5.670141115e-9 |
| h1.3 | λ_min, 60 × 60 Cauchy 1/(i² + j²) | 1.580934514544763949e-63 | 120 | n = 8: 1.147829916e-8 |
| h1.4 | σ_min, 120 × 120 Kahan, θ = 6/5 | 4.488013531293049260e-20 | 60 | n = 40: 6.928154746e-7 |
| h1.5 | λ_min, 40 × 40 Frank | 0.007509065981719966652 | 100 | n = 12: 0.03102806064 |
| h1.6 | λ₁ − λ₂ of the Wilkinson W₄₁⁺ | 1.272072120490893357e-37 | 60 | n = 13: 2.021759535e-6 |
| h1.7 | λ_min, 64 × 64 prolate, w = 1/4 | 6.738206080100924396e-48 | 80 | n = 12: 1.778272860e-8 |
| h1.8 | λ_min, 80 × 80 Moler | 6.158049892052418769e-48 | 80 | n = 12: 5.364405880e-7 |
| h1.9 | κ₂, 60 × 60 KMS, ρ = 1 − 10⁻³⁰ | 1.199177720852744324e32 | 60 | ρ = 1 − 10⁻⁴: 1.196723905e6 |
| h1.10 | λ_min, 40 × 40 Gaussian kernel E^−((i−j)/5)² | 5.423169387591532737e-21 | 50 | n = 10: 8.651257513e-10 |

### 137 — hard definite integrals (machine 6 digits, WorkingPrecision 40: 15)
Provenance: cf = closed form; q = 50-digit adaptive Gauss–Legendre in mpmath with
breakpoints accumulating at the layer.

| # | integral | true value | provenance | machine instance |
|---|---|---|---|---|
| h2.1 | ∫₁⁵ x^x^x | 1.390195540965462193e2180 | q | [1, 4]: 1.476250472508459e151 (q) |
| h2.2 | ∫₀² E^E^E^x | 4.880478191980581312e698 | q | [0, 3/2]: 6.207592136617010e35 (q) |
| h2.3 | ∫₀¹ Sin[1/x] | Sin[1] − Ci(1) = 0.5040670619069283720 | cf, q | same |
| h2.4 | ∫₀¹ x^−x | Σ n^−n = 1.291285997062663540 | cf, q | same |
| h2.5 | ∫₀^∞ E^(−x² − 10⁶/x²) | √π E^−2000/2 = 2.283395465013097565e-869 | cf | 10⁶ → 100: 1.826649837699538e-9 (cf, q) |
| h2.6 | ∫₀¹ 1/(10⁻²⁰ + x²) | 10¹⁰ ArcTan[10¹⁰] = 15707963266.94896619 | cf, q | same |
| h2.7 | ∫₀^∞ Sin[x³] | Γ(4/3)/2 = 0.4464897557846246056 | cf | same |
| h2.8 | ∫₀¹ x^(−999/1000) E^−x | γ(1/1000, 1) = 999.2042906713881633 | cf, series | same |
| h2.9 | ∫₀^∞ x¹⁰⁰⁰ E^−x | 1000! = 4.023872600770937735e2567 | cf | x¹⁰⁰: 100! = 9.332621544394415e157 |
| h2.10 | ∫₀^π Cos[1000 Sin[x]] | π J₀(1000) = 0.07786967112327907847 | cf, q | same |

### 138 — slowly convergent series (NSum: machine 6 digits, WorkingPrecision 30: 14)
Provenance: cf; em = Euler–Maclaurin (direct head + closed-form tail integral +
Bernoulli corrections) at two cut-offs; hz = summand expanded in powers of n for
n ≥ N and summed as Hurwitz zeta values, at two N. A double carries every answer,
so the machine case is the full problem.

| # | series | true value | provenance |
|---|---|---|---|
| h3.1 | Σ n^(−101/100) | ζ(101/100) = 100.5779433384968725 | cf, em |
| h3.2 | Σ_{n≥2} 1/(n Log²n) | 2.109742801236891974 | em (N = 10⁴, 2·10⁴) |
| h3.3 | Σ (−1)^(n+1)/√n | (1 − √2) ζ(1/2) = 0.6048986434216303702 | cf, η |
| h3.4 | Σ (−1)ⁿ Log[n]/n | γ Log 2 − Log²2/2 = 0.1598689037424309718 | cf, nsum |
| h3.5 | Σ 1/(n^(11/10) + √n + 1) | 8.957773692316167741 | hz (N = 2000, 4000) |
| h3.6 | Σ Sin[n]/n | (π − 1)/2 = 1.070796326794896619 | cf, nsum |
| h3.7 | Σ Cos[n]/√n | Re Li_{1/2}(E^I) = −0.1941089350921820497 | cf (Hurwitz reflection), PolyLog |
| h3.8 | Σ_{n≥3} 1/(n Log n (Log Log n)²) | 38.40676809282178632 | em (N = 10⁴, 4·10⁴) |
| h3.9 | Σ (−1)^(n+1)/n^(1/10) | (1 − 2^(9/10)) ζ(1/10) = 0.5222702824645705055 | cf, η |
| h3.10 | Σ 1/(n^(21/20) + 1) | 1/2 + Σ_k (−1)^k (ζ(21(k+1)/20) − 1) = 19.64669717802233418 | cf, hz |

### 139 — exact coefficients (exact; contour machine 5 digits, trapezoidal wp 14)
Exact values by pure integer arithmetic in Python (power-by-squaring convolution,
recurrences, closed sums — two routes where one exists). The check is {15 leading
digits, exponent} plus numerator and denominator mod 10⁹. Machine case: the Cauchy
integral on a reduced instance; wp case: the trapezoidal rule on the full one.

| # | coefficient | true value | machine instance |
|---|---|---|---|
| h4.1 | [x²⁰⁰⁰] (1+x+x²)²⁰⁰⁰ | 1.909454965696130777e952 (953 digits) | n = 200 |
| h4.2 | [x²²⁵⁰] (1+x)¹⁵⁰⁰ (1+x+x²+x³)¹⁰⁰⁰ | 3.984929649324310828e1051 | 150/100, k = 225 |
| h4.3 | [x²²⁰⁰] (x+2)¹⁵⁰⁰ (x³−x+1)⁷⁰⁰ (x²+3x+1)⁴⁰⁰ | 2.837757356724861352e992 | 150/70/40, k = 220 |
| h4.4 | [x²⁰⁰⁰⁰] 1/((1−x)(1−x²)(1−x³)(1−x⁵)(1−x⁷)(1−x¹¹)) | 11585914270681804 | k = 1000: 3876087157; wp contour at k = 2000: 119680802813 |
| h4.5 | Res_{x=1} 1/((x−1)²⁰⁰ (x+2)³⁰⁰) | −Binomial[498, 199]/3⁴⁹⁹ = −1.002978091496547155e-94 | 20/30 |
| h4.6 | [x¹⁰⁰⁰] (1−x)¹⁰⁰⁰ (1+x)²⁰⁰⁰ | −6.656284598533242651e449 | 100/200, k = 100 |
| h4.7 | [x²⁵²⁵] Π_{k≤100} (1+x^k) | 1.731024005948725017e27 | Π to 30, k = 232: 8679280 |
| h4.8 | [x¹⁵⁰⁰] 1/(√(1−4x)(1−x)) | Σ_{k≤1500} C(2k,k) = 2.389556237309696426e901 | k = 150 |
| h4.9 | [x⁴⁰¹] Tan[x] | T₄₀₁/401! = 2.889611317106376218e-79 | k = 41 |
| h4.10 | Res_{x=0} 1/Sin[x]²⁰¹ | C(200,100)/2²⁰⁰ = 0.05634847900925642225 | Sin²¹: 0.176197052001953125 |

### 140 — special-function roots (FindRoot: machine 9 digits, WorkingPrecision 50: 20)
Orthogonal polynomials: a 60-digit Newton iteration on the three-term recurrence
(mpmath) and a Sturm-sequence count at r ∓ 10⁻³⁰ that brackets the root and proves
it extreme (h5.2, h5.3 also a 150-step bisection on that count). Functions are
wrapped `f[x_?NumericQ]` so they are never expanded symbolically.

| # | root | true value | provenance | machine instance |
|---|---|---|---|---|
| h5.1 | largest zero of H₅₀₀ (|H| ~ 1e1100 there) | 31.05074638009002022 | Newton, Sturm | H₁₀₀: 13.40648733814491014 |
| h5.2 | 1 − largest zero of P₁₀₀₀ | 2.888701924489430124e-6 | Newton, Sturm, bisection | same |
| h5.3 | 1 − largest zero of P₃₀₀^(3, 1/2) | 2.227861165971667687e-4 | Newton, bisection, exact CountRoots | same |
| h5.4 | 1 − largest zero of T₂₀₀₁ | 1 − Cos[π/4002] = 3.081169277385704946e-7 | cf | same |
| h5.5 | 100th zero of J₀ | 313.3742660775278447 | library zero routine, Newton | same |
| h5.6 | first zero of J₁₀₀ | 108.8361658984097744 | library zero routine | same |
| h5.7 | largest zero of L₈₀₀^(10) | 3167.337473886621269 | Newton, Sturm | L₈₀^(10): 316.1495015360838567 |
| h5.8 | smallest zero of L₂₀₀₀ | 7.227175802102313270e-4 | Newton, Sturm | same |
| h5.9 | 50th zero of Ai | −38.02100867725525443 | library zero routine | same |
| h5.10 | 10th zero of Z(t) (Riemann–Siegel) | 49.77383247767230218 | published (Odlyzko), zeta-zero routine | same |

### 141 — polynomial systems (exact GroebnerBasis; NSolve machine 8 digits, wp40 18)
Random coefficients from Python `random.seed(97)`, fixed in the file. The ratio
(constant/leading coefficient of the univariate basis member) was checked against
(−1)^d Π vᵢ over a 60-digit NSolve's solutions; counts are Bézout numbers (32 and
70 for Katsura-5 and cyclic-5); Max Re agrees between a 60-digit and a machine
NSolve (h6.5: between two independent 40-digit solvers).

| # | system | basis | ratio | solutions | Max Re (check variable) |
|---|---|---|---|---|---|
| h6.1 | 4 dense quadrics, x y z w | shape, deg 16 | 0.8696578864525617678 | 16 | w 0.4997973735980418593 |
| h6.2 | 3 dense cubics, x y z | shape, deg 27 | 0.01501362839974060569 | 27 | z 3.134381580032942796 |
| h6.3 | 5 dense quadrics, v x y z w | shape, deg 32 | 192.0781013427937613 | 32 | w 3.922013766919696123 |
| h6.4 | 4 dense quadrics, grevlex | 12 elements | — | 16 | w 0.9016989787499557851 |
| h6.5 | Katsura-5 | shape, deg 32 | 9.976845909679557194e-11 | 32 | u1 0.3100445581289618 |
| h6.6 | cyclic-5 | 11 elements, deg 15 | −1 | 70 | e 2.118033988749894848 |
| h6.7 | degrees 7 and 8 in x y | shape, deg 56 | 1.393462280185198192 | 56 | y 1.599325999854060862 |
| h6.8 | 4 quadrics, 6-digit coefficients | shape, deg 16 | 36.78270028974930448 | 16 | w 3.359152893765975268 |
| h6.9 | 3 quadrics + 1 cubic | shape, deg 24 | 7.172513660038497679 | 24 | w 0.9978641849262317102 |
| h6.10 | 2 quadrics + 2 cubics | shape, deg 36 | −32.06374098385247742 | 36 | w 3.388175712657317654 |

### 142 — cancellation limits (exact Limit 20 digits; point machine 5; N[e(x_w), 30] 20)
Limits from the leading Taylor terms (f ~ a x^p, g ~ b x^q ⇒ f[g]/g[f] → a^(1−q) b^(p−1))
or asymptotics at ∞ (h7.4, h7.5 are Gruntz's thesis examples), confirmed by mpmath
at 3000 and 6000 digits at points approaching the limit. Machine case: the value
at x_m = 1/10 (or 10); for h7.2, h7.3, h7.10 the reduced quotient f/g, because a
double cannot evaluate f[g[x]] anywhere near the limit. wp case: N[·, 30] at an
exact deep point x_w (100–450 digits cancel).

| # | expression | limit | x_w |
|---|---|---|---|
| h7.1 | (Sin[Tan] − Tan[Sin])/(ArcSin[ArcTan] − ArcTan[ArcSin]) | 1 | 10⁻⁸ |
| h7.2 | f₂[g₂]/g₂[f₂], f₂ = Tan[Sin] − Sin[Tan], g₂ = Tan[Tanh] − Tanh[Tan] | 4096/729 | 10⁻⁶ |
| h7.3 | f₃[g₃]/g₃[f₃], the 132 sums with ArcSin/ArcTan/ArcSinh/ArcTanh (~x¹¹) | (5/4)¹⁰ | 10⁻⁴ |
| h7.4 | E^x (E^(1/x − E^−x) − E^(1/x)), x → ∞ | −1 | 1000 |
| h7.5 | (E^(x E^−x/(E^−x + E^(−2x²/(x+1)))) − E^x)/x, x → ∞ | −E² | 1000 |
| h7.6 | (1 + 1/x)^(x²) E^−x, x → ∞ | E^(−1/2) | 10³⁰ |
| h7.7 | (Tan∘Tan∘Tan − Sin∘Sin∘Sin)/(Tan∘Tan − Sin∘Sin) | 3/2 | 10⁻¹² |
| h7.8 | (Log[Cos x] + x²/2 + x⁴/12 + x⁶/45)/x⁸ | −17/2520 | 10⁻¹⁰ |
| h7.9 | (Sinh[Tan] − Tan[Sinh])/(ArcSinh[ArcTan] − ArcTan[ArcSinh]) | 1 | 10⁻⁸ |
| h7.10 | f₁₀[g₁₀]/g₁₀[f₁₀], f₁₀ = Tan∘Sin∘Tan∘Sin − Sin∘Tan∘Sin∘Tan, g₁₀ = Sinh[Tanh] − Tanh[Sinh] | 1/64 | 10⁻⁵ |

The point values are in the file header.

### 143 — polylog integrals (exact 20 digits; NIntegrate machine 8, wp30 14)
Closed forms derived by termwise integration (∫₀¹ x^(n−1) Log^k x = (−1)^k k!/n^(k+1))
and partial fractions in n, or by parts; each checked against a 40-digit mpmath
quadrature (≥ 38 digits).

| # | integral | closed form | value |
|---|---|---|---|
| h8.1 | ∫₀¹ Log²x Log(1−x)/x | −π⁴/45 | −2.164646467422276383 |
| h8.2 | ∫₀¹ x³ Li₃(x) | ζ(3)/4 − π²/96 + 25/768 | 0.2302579299452177524 |
| h8.3 | ∫₀¹ Log³(1+x)/x | π⁴/15 + π² Log²2/4 − Log⁴2/4 − 6 Li₄(1/2) − 21 ζ(3) Log 2/4 | 0.1425141979357109157 |
| h8.4 | ∫₀¹ ArcTan x Log x/x | −π³/32 | −0.9689461462593693805 |
| h8.5 | ∫₀¹ Log(1+x²)/(1+x²) | π Log 2/2 − G | 0.1728274509745820502 |
| h8.6 | ∫₀¹ x Li₂(1/(1+x)) | π²/12 + 1/4 − Log 2 | 0.3793198528641679088 |
| h8.7 | ∫₀¹ Log³x/(1−x) | −π⁴/15 | −6.493939402266829149 |
| h8.8 | ∫₀¹ Log(1−x) Log(1+x)/x | −5 ζ(3)/8 | −0.7512855644747464284 |
| h8.9 | ∫₀^∞ Log²x/(1+x²) | π³/8 | 3.875784585037477522 |
| h8.10 | ∫₀¹ x Li₃(−x) | −3ζ(3)/8 + π²/48 − 1/16 | −0.3076545803288195525 |

### 144 — operator eigenvalues (machine 9 digits, WorkingPrecision 40: 18)
Nystrom (composite 32-point Gauss–Legendre by Golub–Welsch, panels doubled to a
tolerance) for the kernels; Nystrom at 60 and 90 nodes and 80 digits agree to ≥ 33
digits. Second methods as named.

| # | operator | true value | second method |
|---|---|---|---|
| h9.1 | λ₁, E^(xy) on [0,1] | 1.353030164745734853 | node count |
| h9.2 | λ₁, 1/(1+x+y) on [0,1] | 0.5362060346834763502 | node count |
| h9.3 | λ₁, Cos[xy] on [0,π] | 1.253314101714016907 | node count |
| h9.4 | λ₂, E^−(x−y)² on [0,2] | 0.5359573175891977563 | node count |
| h9.5 | λ₁, 1/(1+(x−y)²) on [−1,1] | 1.424941177795263197 | node count |
| h9.6 | λ₁, J₀(xy) on [0,4] | 1.665886980965023224 | node count |
| h9.7 | Mathieu a₀(5) (Fourier matrix; machine 20, wp 40 modes) | −5.800046020851508419 | the special function |
| h9.8 | ground state of −y'' + x⁴y (Hermite basis; 40 / 100 functions) | 1.060362090484182900 | published (Hioe–Montroll), 200 functions |
| h9.9 | λ₁ of −y'' = λ(1+x)y, Dirichlet (sine Galerkin; 20 / 200 modes; 7 / 14 digits) | 6.548395306000592515 | Airy characteristic equation |
| h9.10 | λ₁, xy + E^(x+y) on [0,1] (rank 2) | (−1 + 3E² + √(169 − 30E² + 9E⁴))/12 = 3.509384443552449795 | closed form |

### 145 — ODE singularities (machine 6 digits, WorkingPrecision 30: 15)
True values by closed form and/or a 35–40-digit Taylor-series ODE integration in
mpmath (`odefun`).

| # | problem | method in the file | true value | provenance |
|---|---|---|---|---|
| h10.1 | y' = x² + y², y(0) = 0: blow-up | continue with u = 1/y, FindRoot u = 0 | √(2j), J₋₁/₄(j) = 0: 2.003147359426884708 | cf, odefun |
| h10.2 | y' = 1 + y³, y(0) = 0: blow-up | x(s), s = y/(1+y), to s = 1 | 2π/(3√3) = 1.209199576156145234 | cf, quadrature |
| h10.3 | Lane–Emden n = 3: ξ₁ | series start at 10⁻⁶, FindRoot | 6.896848619376960375 | published (Horedt), odefun |
| h10.4 | Lane–Emden n = 3/2: ξ₁ | same | 3.653753736219122425 | published (Horedt), odefun |
| h10.5 | Bratu, λ = 1, lower branch: y'(0) | shooting (secant) | θ Tanh[θ/4], θ = √2 Cosh[θ/4]: 0.5493527287752708190 | cf, odefun |
| h10.6 | Blasius f''(0) | machine shooting at 10; wp scaling at 12 | 0.3320573362151962989 | published (Boyd), odefun |
| h10.7 | Gauss equation, a, b, c = 1/3, 1/2, 2: y(1) | series start, NDSolve to the singular point | Γ(c)Γ(c−a−b)/(Γ(c−a)Γ(c−b)) = 1.159595266963928366 | cf |
| h10.8 | y' = x + y², y(0) = 0: blow-up | machine: naive domain end (5 digits); wp: reciprocal | zero of Bi'(0)Ai(−x) − Ai'(0)Bi(−x): 1.986352707430472813 | cf, odefun |
| h10.9 | y'' = x + y² + y'³: branch point r and y(r) | y(x) to 3/2, then inverse function X(t) | r = 1.617137534968669344, y(r) = 1.149567134607272215 | odefun, two switch points |
| h10.10 | x y'' + y' + x y = 0, y(0) = 1: y(10) | series start at 10⁻⁶ | J₀(10) = −0.2459357644513483352 | cf |

## Baseline run

First run, 2026-11-06: Mathilda v0.6225 (commit 8edb9c4d, before any ISSAC '97
fix), the baseline column on the same host at the same time. The host was
heavily loaded (concurrent builds), so timings are indicative; classes are not
affected. Mathilda files ran in parallel, baseline files sequentially. 240 cases
(100 problems: 100 machine, 100 arbitrary-precision, 40 exact). Classes are from
`issac97_holdout_table.py`, i.e. **against the true answer**: CHECK-FAIL means
at least one system answered wrongly, and the last column says which system(s)
match the truth.

| exp | family | AHEAD | SLOWER | INCOMPLETE | CHECK-FAIL |
|---|---|---|---|---|---|
| 136 | ill-conditioned matrices | 14 | 6 | 0 | 0 |
| 137 | hard definite integrals | 8 | 2 | 0 | 10 |
| 138 | slow series | 9 | 0 | 3 | 8 |
| 139 | exact coefficients | 13 | 3 | 1 | 13 |
| 140 | special-function roots | 11 | 3 | 4 | 2 |
| 141 | polynomial systems | 16 | 12 | 0 | 2 |
| 142 | cancellation limits | 13 | 7 | 4 | 6 |
| 143 | polylog integrals | 20 | 5 | 5 | 0 |
| 144 | operator eigenvalues | 14 | 5 | 0 | 1 |
| 145 | ODE singularities | 9 | 6 | 1 | 4 |
| | **total** | 127 | 49 | 18 | 46 |

### Mathilda wrong answers (silent or with a warning, but returned as a number)

| case | Mathilda | true | repro |
|---|---|---|---|
| h2.1 wp | 8.55e1673 | 1.390e2180 | `NIntegrate[x^x^x, {x, 1, 5}, WorkingPrecision -> 40]` |
| h2.2 wp | 9.06e687 | 4.880e698 | `NIntegrate[E^E^E^x, {x, 0, 2}, WorkingPrecision -> 40]` |
| h2.3 | 0.504166 / 0.50406721577 | 0.5040670619 | `NIntegrate[Sin[1/x], {x, 0, 1}]` (and at WP 40) |
| h2.5 wp | 2.13e-888 | 2.283e-869 | `NIntegrate[E^(-x^2 - 10^6/x^2), {x, 0, Infinity}, WorkingPrecision -> 40]` |
| h2.7 wp | -1.6e619 | 0.44649 | `NIntegrate[Sin[x^3], {x, 0, Infinity}, WorkingPrecision -> 40]` |
| h2.8 | 131.05 / 274.00 | 999.2043 | `NIntegrate[x^(-999/1000) E^-x, {x, 0, 1}]` (and at WP 40) |
| h2.9 wp | 6.96e799 | 1000! = 4.024e2567 | `NIntegrate[x^1000 E^-x, {x, 0, Infinity}, WorkingPrecision -> 40]` |
| h2.10 wp | last of 15 digits | | `NIntegrate[Cos[1000 Sin[x]], {x, 0, Pi}, WorkingPrecision -> 40]` |
| h3.1 | 100.487 | 100.578 | `NSum[n^(-101/100), {n, 1, Infinity}]` |
| h3.2 | 2.10831 | 2.109743 | `NSum[1/(n Log[n]^2), {n, 2, Infinity}]` |
| h3.4 | 0.159883 | 0.159869 | `NSum[(-1)^n Log[n]/n, {n, 1, Infinity}]` |
| h3.6 | 1.11244 (both precisions) | 1.0707963 | `NSum[Sin[n]/n, {n, 1, Infinity}]` |
| h3.7 | -0.391483 / 1.8e511 | -0.1941089 | `NSum[Cos[n]/Sqrt[n], {n, 1, Infinity}]` (and at WP 30) |
| h3.8 | 38.2541 | 38.40677 | `NSum[1/(n Log[n] Log[Log[n]]^2), {n, 3, Infinity}]` |
| h4.4 | 0 (k = 20000); 1705835779/26400 at k = 100 (true 70055) | 11585914270681804 | `SeriesCoefficient[1/((1-x)(1-x^2)(1-x^3)(1-x^5)(1-x^7)(1-x^11)), {x, 0, 100}]` |
| h4.8 | 0 (already at k = 200) | Sum[Binomial[2k, k], {k, 0, 1500}] | `SeriesCoefficient[1/(Sqrt[1 - 4 x] (1 - x)), {x, 0, 200}]` |
| h4.9 | a negative rational (-9.2e-252 at k = 201) | T_401/401! = 2.89e-79 | `SeriesCoefficient[Tan[x], {x, 0, 201}]` |
| h4.10 | a rational equal to 15+ digits, wrong numerator/denominator (denominator divisible by 10^9) | Binomial[200, 100]/2^200 | `Residue[1/Sin[x]^201, {x, 0}]` (exact at Sin^41) |
| h4.1-h4.4, h4.6-h4.10 wp | the trapezoidal rule stops early: off by 1-15000 orders of magnitude | | `NIntegrate[..., Method -> "Trapezoidal", WorkingPrecision -> 30]` on a periodic contour (deficiency W3 of ISSAC97.md, on new inputs) |
| h5.2 wp | the start point | 1 - x = 2.888701924e-6 | `f[x_?NumericQ] := LegendreP[1000, x]; FindRoot[f[x], {x, 999997/1000000}, WorkingPrecision -> 50]` |
| h5.4 machine | 1 - x = 3.07990e-7 | 3.08117e-7 | `f[x_?NumericQ] := ChebyshevT[2001, x]; FindRoot[f[x], {x, 0.9999997}]` |
| h6.5, h6.7 machine | 24 of 32, 46 of 56 solutions | | `NSolve[katsura5, vars]`; `NSolve[{deg-7, deg-8 bivariate}, {x, y}]` (systems in the file) |
| h7.2 | 0 | 4096/729 | `Limit[f2[g2[x]]/g2[f2[x]], x -> 0]`, f2 = Tan[Sin[x]] - Sin[Tan[x]], g2 = Tan[Tanh[x]] - Tanh[Tan[x]] |
| h7.10 | 0 | 1/64 | `Limit[f10[g10[x]]/g10[f10[x]], x -> 0]` (definitions in the file) |
| h7.4, h7.5 wp | 0 | -1.0010005..., -7.374307... | `N[E^1000 (E^(1/1000 - E^-1000) - E^(1/1000)), 30]` |
| h7.8 wp | wrong from the 18th digit | -0.006746031746031746031768 | `N[(Log[Cos[x]] + x^2/2 + x^4/12 + x^6/45)/x^8 /. x -> 10^-10, 30]` |
| h9.9 wp | 6.5483953059745 | 6.5483953060006 | `Eigenvalues[{A, B}]` (generalized, 200 x 200, 40 digits; the sine-Galerkin pencil of the file) |
| h10.1 wp | last of 15 digits | 2.003147359426885 | the reciprocal continuation at WorkingPrecision 30 |
| h10.7 wp | 1.15959526487621 | 1.15959526696393 | NDSolve of the Gauss equation to its singular point x = 1 at WorkingPrecision 30 |

### Mathilda absences and incompletes

| case | what | repro |
|---|---|---|
| h3.2, h3.4, h3.8 wp | NSum at WorkingPrecision 30 times out (> 120 s) | `NSum[1/(n Log[n]^2), {n, 2, Infinity}, WorkingPrecision -> 30]` |
| h4.5 exact | Residue at a 200-fold pole times out | `Residue[1/((x - 1)^200 (x + 2)^300), {x, 1}]` |
| h5.3 | JacobiP with a non-integer parameter does not evaluate numerically | `JacobiP[300, 3, 1/2, 0.9998]` |
| h5.10 | RiemannSiegelZ is not implemented | `RiemannSiegelZ[49.77]` |
| h7.3 | Limit times out | `Limit[f3[g3[x]]/g3[f3[x]], x -> 0]` (the ArcSin/ArcTan/ArcSinh/ArcTanh sums) |
| h7.2, h7.3, h7.10 wp | N[·, 30] returns unevaluated (both systems) | `N[f2[g2[x]]/g2[f2[x]] /. x -> 10^-6, 30]` |
| h8.2, h8.5, h8.6, h8.8, h8.10 | definite Integrate unevaluated | `Integrate[x^3 PolyLog[3, x], {x, 0, 1}]`, `Integrate[Log[1 + x^2]/(1 + x^2), {x, 0, 1}]`, `Integrate[x PolyLog[2, 1/(1 + x)], {x, 0, 1}]`, `Integrate[Log[1 - x] Log[1 + x]/x, {x, 0, 1}]`, `Integrate[x PolyLog[3, -x], {x, 0, 1}]` |
| h10.8 machine | `InterpolatingFunction[...]["Domain"]` unevaluated (deficiency A8) | `(y /. First[NDSolve[{y'[x] == x + y[x]^2, y[0] == 0}, y, {x, 0, 3}]])["Domain"]` |

### Baseline wrong answers

h2.8 (machine 538.67, WP 40 999.20462 — true 999.2042907), h2.9 WP 40 (wrong
from the 12th digit), h3.7 WP 30 (5.9e14), h3.8 (38.3285 and 38.33569 — true
38.40677), h7.2 Limit (0), h7.4, h7.5, h7.8 N[·, 30] (0), h9.9 WP 40 (times
out), h10.1/h10.8/h10.9 WorkingPrecision 30 (wrong from the 12th digit), h10.7
WP 30 (0). Neither system is right on h2.8, h3.7 WP 30, h3.8, h7.2 (Limit and
N), h7.3 (Limit and N, both time out), h7.4/h7.5 N, h7.10 N, h10.1 WP 30 and
h10.7 WP 30. (Before the h4.4 trapezoidal case was moved from k = 20000 to
k = 2000, both systems' rule was fooled by aliasing there.)

### All cases

| case | Mathilda ms | baseline ms | ratio | class | matches truth |
|---|---|---|---|---|---|
| h1.1 eigenvalues machine n=12 | 0.7 | 0.5 | 1.35 | AHEAD |  |
| h1.1 eigenvalues wp250 n=100 | 2994.6 | 20809.5 | 0.14 | AHEAD |  |
| h1.2 singular values machine n=12 | 1.4 | 0.8 | 1.77 | SLOWER |  |
| h1.2 singular values wp100 n=50 | 9954.1 | 1138.6 | 8.74 | SLOWER |  |
| h1.3 eigenvalues machine n=8 | 0.4 | 0.3 | 1.16 | AHEAD |  |
| h1.3 eigenvalues wp120 n=60 | 508.2 | 3343.4 | 0.15 | AHEAD |  |
| h1.4 singular values machine n=40 | 72.1 | 47.9 | 1.51 | SLOWER |  |
| h1.4 singular values wp60 n=120 | 51062.2 | 9420.4 | 5.42 | SLOWER |  |
| h1.5 eigenvalues machine n=12 | 0.5 | 0.8 | 0.69 | AHEAD |  |
| h1.5 eigenvalues wp100 n=40 | 672.1 | 2321.5 | 0.29 | AHEAD |  |
| h1.6 eigenvalue gap machine n=13 | 0.8 | 0.7 | 1.07 | AHEAD |  |
| h1.6 eigenvalue gap wp60 n=41 | 87.6 | 1107.8 | 0.08 | AHEAD |  |
| h1.7 eigenvalues machine n=12 | 3.7 | 2.0 | 1.84 | SLOWER |  |
| h1.7 eigenvalues wp80 n=64 | 1010.2 | 3473.6 | 0.29 | AHEAD |  |
| h1.8 eigenvalues machine n=12 | 0.5 | 0.6 | 0.74 | AHEAD |  |
| h1.8 eigenvalues wp80 n=80 | 999.1 | 9418.1 | 0.11 | AHEAD |  |
| h1.9 eigenvalues machine rho=1-1e-4 | 155.8 | 863.7 | 0.18 | AHEAD |  |
| h1.9 eigenvalues wp60 rho=1-1e-30 | 1570.2 | 3909.5 | 0.40 | AHEAD |  |
| h1.10 eigenvalues machine n=10 | 3.1 | 0.6 | 5.19 | SLOWER |  |
| h1.10 eigenvalues wp50 n=40 | 303.9 | 1033.3 | 0.29 | AHEAD |  |
| h2.1 NIntegrate machine [1,4] | 0.1 | 9.4 | 0.01 | AHEAD |  |
| h2.1 NIntegrate wp40 [1,5] | 3.0 | 379.5 |  | CHECK-FAIL | baseline |
| h2.2 NIntegrate machine [0,3/2] | 0.1 | 29.3 | 0.00 | AHEAD |  |
| h2.2 NIntegrate wp40 [0,2] | 2.4 | 130.2 |  | CHECK-FAIL | baseline |
| h2.3 NIntegrate machine | 2.3 | 74.0 |  | CHECK-FAIL | baseline |
| h2.3 NIntegrate wp40 | 167028.0 | 1326.5 |  | CHECK-FAIL | baseline |
| h2.4 NIntegrate machine | 0.1 | 6.8 | 0.02 | AHEAD |  |
| h2.4 NIntegrate wp40 | 192.0 | 153.2 | 1.25 | AHEAD |  |
| h2.5 NIntegrate machine b=100 | 0.2 | 9.4 | 0.02 | AHEAD |  |
| h2.5 NIntegrate wp40 b=10^6 | 33.0 | 114.6 |  | CHECK-FAIL | baseline |
| h2.6 NIntegrate machine | 27.1 | 22.3 | 1.22 | AHEAD |  |
| h2.6 NIntegrate wp40 | 1900.9 | 91.5 | 20.78 | SLOWER |  |
| h2.7 NIntegrate machine | 1670.3 | 22.7 | 73.50 | SLOWER |  |
| h2.7 NIntegrate wp40 | 21.8 | 73.7 |  | CHECK-FAIL | baseline |
| h2.8 NIntegrate machine | 310.9 | 35.7 |  | CHECK-FAIL | neither |
| h2.8 NIntegrate wp40 | 374101.0 | 197.1 |  | CHECK-FAIL | neither |
| h2.9 NIntegrate machine x^100 | 3.1 | 78.6 | 0.04 | AHEAD |  |
| h2.9 NIntegrate wp40 x^1000 | 0.5 | 534.3 |  | CHECK-FAIL | neither |
| h2.10 NIntegrate machine | 2.0 | 309.3 | 0.01 | AHEAD |  |
| h2.10 NIntegrate wp40 | 3334.4 | 2251.3 |  | CHECK-FAIL | baseline |
| h3.1 NSum machine | 363.3 | 39.9 |  | CHECK-FAIL | baseline |
| h3.1 NSum wp30 | 270.0 | 200.3 | 1.35 | AHEAD |  |
| h3.2 NSum machine | 358.8 | 52.3 |  | CHECK-FAIL | baseline |
| h3.2 NSum wp30 | timeout | 197.8 |  | INCOMPLETE | baseline correct |
| h3.3 NSum machine | 0.1 | 204.6 | 0.00 | AHEAD |  |
| h3.3 NSum wp30 | 1.7 | 229.2 | 0.01 | AHEAD |  |
| h3.4 NSum machine | 3273.6 | 261.7 |  | CHECK-FAIL | baseline |
| h3.4 NSum wp30 | timeout | 312.3 |  | INCOMPLETE | baseline correct |
| h3.5 NSum machine | 23.7 | 55.1 | 0.43 | AHEAD |  |
| h3.5 NSum wp30 | 1095.0 | 837.3 | 1.31 | AHEAD |  |
| h3.6 NSum machine | 0.1 | 123.4 |  | CHECK-FAIL | baseline |
| h3.6 NSum wp30 | 33.0 | 277.4 |  | CHECK-FAIL | baseline |
| h3.7 NSum machine | 244.1 | 242.2 |  | CHECK-FAIL | baseline |
| h3.7 NSum wp30 | 22.0 | 992.7 |  | CHECK-FAIL | neither |
| h3.8 NSum machine | 368.9 | 184.3 |  | CHECK-FAIL | neither |
| h3.8 NSum wp30 | timeout | 685.9 |  | INCOMPLETE | baseline wrong |
| h3.9 NSum machine | 0.1 | 199.1 | 0.00 | AHEAD |  |
| h3.9 NSum wp30 | 1.6 | 216.1 | 0.01 | AHEAD |  |
| h3.10 NSum machine | 5.7 | 44.5 | 0.13 | AHEAD |  |
| h3.10 NSum wp30 | 129.5 | 195.1 | 0.66 | AHEAD |  |
| h4.1 exact Coefficient[Expand] | 392.0 | 2244.1 | 0.17 | AHEAD |  |
| h4.1 contour NIntegrate machine n=200 | 0.4 | 102.3 | 0.00 | AHEAD |  |
| h4.1 contour Trapezoidal wp30 | 1.1 | 271.3 |  | CHECK-FAIL | baseline |
| h4.2 exact Coefficient[Expand] | 524.3 | 3662.6 | 0.14 | AHEAD |  |
| h4.2 contour NIntegrate machine | 0.6 | 61.9 | 0.01 | AHEAD |  |
| h4.2 contour Trapezoidal wp30 | 1.2 | 414.8 |  | CHECK-FAIL | baseline |
| h4.3 exact Coefficient[Expand] | 675.5 | 4247.6 | 0.16 | AHEAD |  |
| h4.3 contour NIntegrate machine | 1.2 | 153.8 | 0.01 | AHEAD |  |
| h4.3 contour Trapezoidal wp350 | 14.7 | 1512.1 |  | CHECK-FAIL | baseline |
| h4.4 exact SeriesCoefficient k=20000 | 2537.5 | 4714.3 |  | CHECK-FAIL | baseline |
| h4.4 contour NIntegrate machine k=1000 | 2566.8 | 100.8 | 25.48 | SLOWER |  |
| h4.4 contour Trapezoidal wp40 k=2000 | 63042.2 | 5803.3 |  | CHECK-FAIL | baseline |
| h4.5 exact Residue | timeout | 31.4 |  | INCOMPLETE | baseline correct |
| h4.5 contour NIntegrate machine 20/30 | 0.6 | 52.6 | 0.01 | AHEAD |  |
| h4.5 contour Trapezoidal wp30 | 5394.3 | 95.2 | 56.68 | SLOWER |  |
| h4.6 exact Coefficient[Expand] | 249.3 | 1854.6 | 0.13 | AHEAD |  |
| h4.6 contour NIntegrate machine | 10.4 | 199.3 | 0.05 | AHEAD |  |
| h4.6 contour Trapezoidal wp70 | 1.5 | 1106.0 |  | CHECK-FAIL | baseline |
| h4.7 exact Coefficient[Expand] | 412.6 | 77.2 | 5.35 | SLOWER |  |
| h4.7 contour NIntegrate machine m=30 | 29.7 | 267.5 | 0.11 | AHEAD |  |
| h4.7 contour Trapezoidal wp30 | 103064.0 | 25040.0 |  | CHECK-FAIL | baseline |
| h4.8 exact SeriesCoefficient k=1500 | 144.4 | 106.4 |  | CHECK-FAIL | baseline |
| h4.8 contour NIntegrate machine k=150 | 4.4 | 51.9 | 0.09 | AHEAD |  |
| h4.8 contour Trapezoidal wp60 k=1500 | 9718.6 | 641.9 |  | CHECK-FAIL | baseline |
| h4.9 exact SeriesCoefficient Tan k=401 | 5117.7 | 1891.6 |  | CHECK-FAIL | baseline |
| h4.9 contour NIntegrate machine k=41 | 1.2 | 62.6 | 0.02 | AHEAD |  |
| h4.9 contour Trapezoidal wp40 k=401 | 1.1 | 185.6 |  | CHECK-FAIL | baseline |
| h4.10 exact Residue | 42372.6 | 4478.4 |  | CHECK-FAIL | baseline |
| h4.10 contour NIntegrate machine Sin^21 | 0.6 | 58.4 | 0.01 | AHEAD |  |
| h4.10 contour Trapezoidal wp50 | 2.1 | 61.9 |  | CHECK-FAIL | baseline |
| h5.1 FindRoot machine H100 | 0.1 | 1.7 | 0.05 | AHEAD |  |
| h5.1 FindRoot wp50 H500 | 33.2 | 66.4 | 0.50 | AHEAD |  |
| h5.2 FindRoot machine | 0.1 | 0.6 | 0.20 | AHEAD |  |
| h5.2 FindRoot wp50 | 57823.4 | 428.0 |  | CHECK-FAIL | baseline |
| h5.3 FindRoot machine | unevaluated | 26.4 |  | INCOMPLETE | baseline correct |
| h5.3 FindRoot wp50 | unevaluated | 196.4 |  | INCOMPLETE | baseline correct |
| h5.4 FindRoot machine | 1.0 | 0.7 |  | CHECK-FAIL | baseline |
| h5.4 FindRoot wp50 | 35.1 | 1.9 | 18.85 | SLOWER |  |
| h5.5 FindRoot machine | 0.1 | 1.0 | 0.06 | AHEAD |  |
| h5.5 FindRoot wp50 | 1.7 | 7.6 | 0.22 | AHEAD |  |
| h5.6 FindRoot machine | 0.1 | 6.6 | 0.01 | AHEAD |  |
| h5.6 FindRoot wp50 | 2.4 | 120.8 | 0.02 | AHEAD |  |
| h5.7 FindRoot machine L80 | 0.1 | 1.4 | 0.07 | AHEAD |  |
| h5.7 FindRoot wp50 L800 | 108.6 | 116.7 | 0.93 | AHEAD |  |
| h5.8 FindRoot machine | 81.4 | 24.4 | 3.33 | SLOWER |  |
| h5.8 FindRoot wp50 | 142.3 | 301.1 | 0.47 | AHEAD |  |
| h5.9 FindRoot machine | 0.1 | 0.4 | 0.20 | AHEAD |  |
| h5.9 FindRoot wp50 | 20.8 | 2.9 | 7.28 | SLOWER |  |
| h5.10 FindRoot machine | unevaluated | 2.7 |  | INCOMPLETE | baseline correct |
| h5.10 FindRoot wp50 | unevaluated | 106.8 |  | INCOMPLETE | baseline correct |
| h6.1 GroebnerBasis exact lex | 151.7 | 46.9 | 3.23 | SLOWER |  |
| h6.1 NSolve machine | 122.2 | 171.6 | 0.71 | AHEAD |  |
| h6.1 NSolve wp40 | 177.5 | 230.2 | 0.77 | AHEAD |  |
| h6.2 GroebnerBasis exact lex | 296.2 | 65.0 | 4.55 | SLOWER |  |
| h6.2 NSolve machine | 253.2 | 350.5 | 0.72 | AHEAD |  |
| h6.2 NSolve wp40 | 517.8 | 591.4 | 0.88 | AHEAD |  |
| h6.3 GroebnerBasis exact lex | 4575.6 | 489.0 | 9.36 | SLOWER |  |
| h6.3 NSolve machine | 1432.2 | 825.7 | 1.73 | SLOWER |  |
| h6.3 NSolve wp40 | 1797.4 | 1183.0 | 1.52 | SLOWER |  |
| h6.4 GroebnerBasis exact grevlex | 28.4 | 13.3 | 2.14 | SLOWER |  |
| h6.4 NSolve machine | 101.9 | 159.2 | 0.64 | AHEAD |  |
| h6.4 NSolve wp40 | 138.0 | 222.6 | 0.62 | AHEAD |  |
| h6.5 GroebnerBasis exact lex | 428.6 | 398.1 | 1.08 | AHEAD |  |
| h6.5 NSolve machine | 566.4 | 98.5 |  | CHECK-FAIL | baseline |
| h6.5 NSolve wp40 | 944.7 | 1124.7 | 0.84 | AHEAD |  |
| h6.6 GroebnerBasis exact lex | 49.7 | 109.1 | 0.46 | AHEAD |  |
| h6.6 NSolve machine | 2151.5 | 54.7 | 39.31 | SLOWER |  |
| h6.6 NSolve wp40 | 4910.2 | 4622.6 | 1.06 | AHEAD |  |
| h6.7 GroebnerBasis exact lex | 648.3 | 101.8 | 6.37 | SLOWER |  |
| h6.7 NSolve machine | 996.7 | 1021.3 |  | CHECK-FAIL | baseline |
| h6.7 NSolve wp40 | 2198.3 | 3198.5 | 0.69 | AHEAD |  |
| h6.8 GroebnerBasis exact lex | 1875.8 | 81.0 | 23.15 | SLOWER |  |
| h6.8 NSolve machine | 336.7 | 179.7 | 1.87 | SLOWER |  |
| h6.8 NSolve wp40 | 378.2 | 265.1 | 1.43 | AHEAD |  |
| h6.9 GroebnerBasis exact lex | 548.5 | 102.6 | 5.35 | SLOWER |  |
| h6.9 NSolve machine | 267.6 | 400.2 | 0.67 | AHEAD |  |
| h6.9 NSolve wp40 | 366.4 | 491.5 | 0.75 | AHEAD |  |
| h6.10 GroebnerBasis exact lex | 2773.4 | 337.6 | 8.22 | SLOWER |  |
| h6.10 NSolve machine | 914.1 | 805.4 | 1.13 | AHEAD |  |
| h6.10 NSolve wp40 | 1308.1 | 1139.5 | 1.15 | AHEAD |  |
| h7.1 Limit exact | 244.3 | 300.7 | 0.81 | AHEAD |  |
| h7.1 point machine | 0.2 | 0.0 | 7.91 | SLOWER |  |
| h7.1 point N[.., 30] | 1.5 | 0.0 | 42.43 | SLOWER |  |
| h7.2 Limit exact | 27380.7 | 107907.0 |  | CHECK-FAIL | neither |
| h7.2 point machine reduced f/g | 0.0 | 0.0 | 1.48 | AHEAD |  |
| h7.2 point N[.., 30] | unevaluated | unevaluated |  | INCOMPLETE | baseline also none |
| h7.3 Limit exact | timeout | timeout |  | INCOMPLETE | baseline also none |
| h7.3 point machine reduced f/g | 1.3 | 0.2 |  | CHECK-FAIL | Mathilda |
| h7.3 point N[.., 30] | unevaluated | unevaluated |  | INCOMPLETE | baseline also none |
| h7.4 Limit exact | 19.9 | 145.1 | 0.14 | AHEAD |  |
| h7.4 point machine | 0.0 | 0.0 | 0.83 | AHEAD |  |
| h7.4 point N[.., 30] | 0.3 | 0.1 |  | CHECK-FAIL | neither |
| h7.5 Limit exact | 5768.4 | 553.8 | 10.42 | SLOWER |  |
| h7.5 point machine | 0.0 | 0.1 | 0.45 | AHEAD |  |
| h7.5 point N[.., 30] | 0.7 | 0.3 |  | CHECK-FAIL | neither |
| h7.6 Limit exact | 7.8 | 83.3 | 0.09 | AHEAD |  |
| h7.6 point machine | 0.0 | 0.0 | 0.71 | AHEAD |  |
| h7.6 point N[.., 30] | 0.5 | 0.1 | 5.32 | SLOWER |  |
| h7.7 Limit exact | 180.5 | 1072.4 | 0.17 | AHEAD |  |
| h7.7 point machine | 0.0 | 0.0 | 0.67 | AHEAD |  |
| h7.7 point N[.., 30] | 0.5 | 0.1 | 3.99 | SLOWER |  |
| h7.8 Limit exact | 5.5 | 103.3 | 0.05 | AHEAD |  |
| h7.8 point machine | 0.0 | 0.0 | 0.85 | AHEAD |  |
| h7.8 point N[.., 30] | 0.3 | 0.2 |  | CHECK-FAIL | neither |
| h7.9 Limit exact | 177.1 | 570.6 | 0.31 | AHEAD |  |
| h7.9 point machine | 0.1 | 0.0 | 1.70 | SLOWER |  |
| h7.9 point N[.., 30] | 1.5 | 0.1 | 15.68 | SLOWER |  |
| h7.10 Limit exact | 49096.1 | timeout |  | CHECK-FAIL | neither (baseline none) |
| h7.10 point machine reduced f/g | 0.0 | 0.1 | 0.64 | AHEAD |  |
| h7.10 point N[.., 30] | unevaluated | unevaluated |  | INCOMPLETE | baseline also none |
| h8.1 Integrate exact | 413.6 | 5003.4 | 0.08 | AHEAD |  |
| h8.1 NIntegrate machine | 0.3 | 8.7 | 0.03 | AHEAD |  |
| h8.1 NIntegrate wp30 | 23.2 | 36.8 | 0.63 | AHEAD |  |
| h8.2 Integrate exact | unevaluated | 3156.4 |  | INCOMPLETE | baseline correct |
| h8.2 NIntegrate machine | 0.1 | 7.9 | 0.01 | AHEAD |  |
| h8.2 NIntegrate wp30 | 1104.5 | 685.1 | 1.61 | SLOWER |  |
| h8.3 Integrate exact | 160.4 | 3839.1 | 0.04 | AHEAD |  |
| h8.3 NIntegrate machine | 0.1 | 7.3 | 0.01 | AHEAD |  |
| h8.3 NIntegrate wp30 | 27.3 | 10.9 | 2.51 | SLOWER |  |
| h8.4 Integrate exact | 3594.6 | 3189.8 | 1.13 | AHEAD |  |
| h8.4 NIntegrate machine | 0.2 | 7.6 | 0.03 | AHEAD |  |
| h8.4 NIntegrate wp30 | 102.9 | 24.3 | 4.24 | SLOWER |  |
| h8.5 Integrate exact | unevaluated | 1496.7 |  | INCOMPLETE | baseline correct |
| h8.5 NIntegrate machine | 0.1 | 5.4 | 0.01 | AHEAD |  |
| h8.5 NIntegrate wp30 | 42.0 | 9.6 | 4.38 | SLOWER |  |
| h8.6 Integrate exact | unevaluated | 3825.4 |  | INCOMPLETE | baseline correct |
| h8.6 NIntegrate machine | 0.1 | 4.9 | 0.02 | AHEAD |  |
| h8.6 NIntegrate wp30 | 640.7 | 522.5 | 1.23 | AHEAD |  |
| h8.7 Integrate exact | 0.9 | 2567.6 | 0.00 | AHEAD |  |
| h8.7 NIntegrate machine | 0.2 | 8.2 | 0.03 | AHEAD |  |
| h8.7 NIntegrate wp30 | 22.8 | 23.5 | 0.97 | AHEAD |  |
| h8.8 Integrate exact | unevaluated | 13849.2 |  | INCOMPLETE | baseline correct |
| h8.8 NIntegrate machine | 0.2 | 8.6 | 0.02 | AHEAD |  |
| h8.8 NIntegrate wp30 | 25.2 | 22.2 | 1.14 | AHEAD |  |
| h8.9 Integrate exact | 5.6 | 8773.2 | 0.00 | AHEAD |  |
| h8.9 NIntegrate machine | 0.1 | 9.7 | 0.01 | AHEAD |  |
| h8.9 NIntegrate wp30 | 16.8 | 47.3 | 0.36 | AHEAD |  |
| h8.10 Integrate exact | unevaluated | 1416.9 |  | INCOMPLETE | baseline correct |
| h8.10 NIntegrate machine | 0.1 | 6.7 | 0.01 | AHEAD |  |
| h8.10 NIntegrate wp30 | 16182.6 | 58.7 | 275.48 | SLOWER |  |
| h9.1 Nystrom machine | 44.1 | 21.2 | 2.08 | SLOWER |  |
| h9.1 Nystrom wp40 | 607.7 | 1388.8 | 0.44 | AHEAD |  |
| h9.2 Nystrom machine | 56.1 | 29.0 | 1.94 | SLOWER |  |
| h9.2 Nystrom wp40 | 560.5 | 1167.5 | 0.48 | AHEAD |  |
| h9.3 Nystrom machine | 59.6 | 26.1 | 2.28 | SLOWER |  |
| h9.3 Nystrom wp40 | 589.2 | 1628.8 | 0.36 | AHEAD |  |
| h9.4 Nystrom machine | 81.0 | 58.9 | 1.37 | AHEAD |  |
| h9.4 Nystrom wp40 | 808.3 | 1891.9 | 0.43 | AHEAD |  |
| h9.5 Nystrom machine | 88.4 | 62.1 | 1.42 | AHEAD |  |
| h9.5 Nystrom wp40 | 714.5 | 1832.1 | 0.39 | AHEAD |  |
| h9.6 Nystrom machine | 47.8 | 382.5 | 0.13 | AHEAD |  |
| h9.6 Nystrom wp40 | 707.2 | 3191.8 | 0.22 | AHEAD |  |
| h9.7 Fourier matrix machine n=20 | 2.0 | 1.8 | 1.12 | AHEAD |  |
| h9.7 Fourier matrix wp40 n=40 | 31.1 | 321.5 | 0.10 | AHEAD |  |
| h9.8 Hermite basis machine n=40 | 38.5 | 75.2 | 0.51 | AHEAD |  |
| h9.8 Hermite basis wp40 n=100 | 36258.7 | 2866.0 | 12.65 | SLOWER |  |
| h9.9 sine Galerkin machine n=20 | 40.4 | 24.2 | 1.67 | SLOWER |  |
| h9.9 sine Galerkin wp40 n=200 | 4955.6 | timeout |  | CHECK-FAIL | neither (baseline none) |
| h9.10 Nystrom machine | 55.3 | 45.9 | 1.20 | AHEAD |  |
| h9.10 Nystrom wp40 | 613.0 | 3030.8 | 0.20 | AHEAD |  |
| h10.1 reciprocal machine | 0.7 | 4.1 | 0.16 | AHEAD |  |
| h10.1 reciprocal wp30 | 8298.1 | 459.6 |  | CHECK-FAIL | neither |
| h10.8 naive NDSolve domain end machine | unevaluated | 6.8 |  | INCOMPLETE | baseline correct |
| h10.8 reciprocal wp30 | 3553.4 | 274.5 |  | CHECK-FAIL | Mathilda |
| h10.2 substitution machine | 0.3 | 1.7 | 0.18 | AHEAD |  |
| h10.2 substitution wp30 | 516.9 | 199.4 | 2.59 | SLOWER |  |
| h10.3 Lane-Emden n=3 machine | 0.6 | 4.2 | 0.15 | AHEAD |  |
| h10.3 Lane-Emden n=3 wp30 | 14466.4 | 295.8 | 48.91 | SLOWER |  |
| h10.4 Lane-Emden n=3/2 machine | 0.6 | 6.5 | 0.10 | AHEAD |  |
| h10.4 Lane-Emden n=3/2 wp30 | 11287.9 | 597.7 | 18.88 | SLOWER |  |
| h10.5 Bratu shooting machine | 1.0 | 11.4 | 0.09 | AHEAD |  |
| h10.5 Bratu shooting wp30 | 6458.7 | 547.6 | 11.80 | SLOWER |  |
| h10.6 Blasius shooting machine | 3.4 | 25.1 | 0.13 | AHEAD |  |
| h10.6 Blasius scaling wp30 | 1849.7 | 441.5 | 4.19 | SLOWER |  |
| h10.7 Gauss value at x=1 machine | 1.2 | 7.1 | 0.16 | AHEAD |  |
| h10.7 Gauss value at x=1 wp30 | 4128.4 | 4761.6 |  | CHECK-FAIL | neither |
| h10.9 inverse function machine | 0.8 | 4.1 | 0.19 | AHEAD |  |
| h10.9 inverse function wp30 | 2487.2 | 458.3 |  | CHECK-FAIL | Mathilda |
| h10.10 Bessel y(10) machine | 0.6 | 3.0 | 0.20 | AHEAD |  |
| h10.10 Bessel y(10) wp30 | 2356.6 | 231.6 | 10.17 | SLOWER |  |
