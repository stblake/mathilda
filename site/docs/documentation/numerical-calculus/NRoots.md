# NRoots

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NRoots[lhs == rhs, var]`**

yields a disjunction of equations var==r1 || var==r2 || ... giving numerical approximations to the roots of the polynomial equation in var. Roots of multiplicity k appear as k identical equations; a single root yields a bare equation. Real and complex coefficients are handled at machine and arbitrary precision. Method -\> Automatic uses the Aberth-Ehrlich simultaneous iteration; "CompanionMatrix" uses companion-matrix eigenvalues (real QR directly, complex via a real 2n embedding); "JenkinsTraub" uses the three-stage Jenkins-Traub algorithm.

<details>
<summary>Notes</summary>

Options: Method (Automatic | "Aberth" | "CompanionMatrix" | "JenkinsTraub"), PrecisionGoal (Automatic = machine; a digit count selects arbitrary precision), AccuracyGoal (default MachinePrecision; a root whose residual exceeds the goal triggers an NRoots::accgl warning), MaxIterations, StepMonitor.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= NRoots[1 + 2 x + 3 x^2 + 4 x^3 == 0, x]
Out[1]= x == -0.60583 || x == -0.0720852 - 0.638327*I || x == -0.0720852 + 0.638327*I

In[2]:= NRoots[x^2 - 2 == 0, x]
Out[2]= x == -1.41421 || x == 1.41421

In[3]:= NRoots[x^2 + 1 == 0, x]
Out[3]= x == 0.0 - 1.0*I || x == 0.0 + 1.0*I

In[4]:= NRoots[(x - 1)^3 == 0, x]
Out[4]= x == 1.0 || x == 1.0 || x == 1.0

In[5]:= NRoots[x^2 - (3 + 4 I) == 0, x]
Out[5]= x == -2.0 - 1.0*I || x == 2.0 + 1.0*I
```

### Options (1)

```mathematica
In[6]:= NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 30]
Out[6]= x == -1.414213562373095048801688724209 || x == 1.414213562373095048801688724209
```

### Worked examples (2)

```mathematica
In[7]:= NRoots[1==0, x]
Out[7]= False

In[8]:= NRoots[1==1, x]
Out[8]= True
```

### Applications (4)

The two real square roots of 2

```mathematica
In[9]:= NRoots[x^2 - 2 == 0, x]
Out[9]= x == -1.41421 || x == 1.41421
```

A complex-conjugate pair

```mathematica
In[10]:= NRoots[x^2 + 1 == 0, x]
Out[10]= x == 0.0 - 1.0*I || x == 0.0 + 1.0*I
```

The three cube roots of unity

```mathematica
In[11]:= NRoots[x^3 - 1 == 0, x]
Out[11]= x == 1.0 || x == -0.5 - 0.866025*I || x == -0.5 + 0.866025*I
```

A digit count selects arbitrary precision

```mathematica
In[12]:= NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 20]
Out[12]= x == -1.4142135623730950488 || x == 1.4142135623730950488
```

## Algorithm

```text
nroots.c — NRoots[lhs == rhs, var, opts]   (see nroots.h)
```

Numerically finds every root of a univariate polynomial equation and returns

```text
a disjunction of equations  var==r1 || var==r2 || ...  (a bare equation when
```

there is a single root), repeating identical equations for roots of multiplicity > 1.

Pipeline:

```text
  1. parse  Equal[lhs, rhs]  and the variable; form poly = lhs - rhs.
  2. Expand, validate polynomial-in-var, extract coefficients.
  3. numericalise each coefficient to a complex MPFR value (working prec).
  4. strip a trailing x^m factor (exact zero roots), then dispatch the
     reduced polynomial to the selected engine (Aberth / CompanionMatrix /
     JenkinsTraub), which returns all roots with multiplicity.
  5. chop noise, sort canonically, cluster multiple roots to identical
     values, round to the target precision, assemble the disjunction.
```

All numeric work is MPFR; without USE_MPFR, NRoots returns unevaluated.

Memory contract: never frees `res`; returns a fresh Expr* or NULL.

## Implementation notes

**Algorithm.** `builtin_nroots` parses `Equal[lhs, rhs]` and the variable, forms
`poly = lhs − rhs`, expands it, extracts coefficients (each numericalised to a
complex MPFR value at working precision), strips a trailing `x^m` factor (exact
zero roots), and dispatches the reduced polynomial. The output is a disjunction
`x == r1 || … || x == rd` (a bare `Equal` for one root), with multiple roots
emitted as repeated identical equations. Three engines:

- **Aberth–Ehrlich** (`nroots_aberth.c`), the default: all `d` roots iterated
  together from **Bini's convex-hull (Newton-polygon) initial placement** — the
  upper hull of `(k, log|a_k|)` seeds starting points on circles — using the
  Aberth correction `w_i = N_i / (1 − N_i S_i)` with `N_i = p(z_i)/p'(z_i)` and
  `S_i = Σ_{j≠i} 1/(z_i − z_j)`, applied Gauss–Seidel in place. Cubically
  convergent at simple roots; a multiplicity-`m` root emerges as a tight
  cluster. Every root is Newton-polished against the original polynomial.
- **CompanionMatrix** (`"CompanionMatrix"`, also the Automatic machine fast
  path): the Frobenius companion's eigenvalues — real coefficients via a real
  QR (LAPACK `dgeev` on the machine path, the MPFR real-QR kernel otherwise,
  "`numpy.roots` exactly"), complex coefficients via the `ℂ^{n×n} → ℝ^{2n×2n}`
  embedding `[[Re, −Im], [Im, Re]]` whose spectrum doubles the roots, filtered
  by residual and read for multiplicity by synthetic division.
- **JenkinsTraub** (`nroots_jt.c`): the three-stage shifted-deflation CPOLY —
  Stage 1 no-shift, Stage 2 fixed-shift near a Cauchy-circle zero (sweeping up
  to 12 angles), Stage 3 variable-shift (cubically convergent) — the found zero
  polished against the original polynomial, then deflated; a multiplicity-`m`
  zero is found `m` times.

An exact-integer polynomial of degree ≥ 2 is first **Yun squarefree-
decomposed** so a high multiplicity like `(x²−2)^30` stays well conditioned.

**Data structures.** `NrPoly { int deg; ncpx* c; int is_real; mpfr_prec_t prec }`,
where `ncpx` is the MPFR complex type (a pair of `mpfr_t`). An Automatic machine
request works in `double`/`double _Complex` (real/interleaved companion arrays
for LAPACK), rounding the final roots back to `double`/`Complex`; a
`PrecisionGoal` digit count runs the *same* Aberth loop entirely in MPFR complex
arithmetic at `target_bits + max(48, target_bits/2)` bits (the convergence
tolerance `2^{−(wp−8)}` scales with it). Post-processing chops noise, single-
linkage clusters near-equal roots to a shared value (so multiplicities become
numerically identical), symmetrises conjugate pairs for real coefficients, and
canonically sorts. The whole numeric core is `#ifdef USE_MPFR`; without MPFR,
`NRoots` emits `nompfr` and stays unevaluated.

**Complexity / limits.** Aberth is `O(d²)` per sweep over `MaxIterations`
(default `100 + 20d`); companion is the `O(d³)` QR; Jenkins–Traub is one root at
a time with fixed inner caps (`MaxIterations` is ignored by companion and JT).
`PrecisionGoal` (`Automatic`/`Infinity` → machine; a digit count → MPFR) selects
precision; `AccuracyGoal` (default `MachinePrecision`) drives only the post-solve
residual check — a root whose Newton correction exceeds the combined tolerance
triggers `NRoots::accgl` (clustered/stationary roots are skipped in that test).
Declines (stay unevaluated) on a non-equation, non-polynomial, non-numeric
coefficient, non-variable, unknown method, or engine non-convergence (`conv`);
all diagnostics route through `mth_message`. `NRoots` is `Protected`, not
`Listable`.

**Attributes:** `Protected`.

## References

**See also:** [Expand](../../algebra/Expand/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- D. A. Bini, *Numerical computation of polynomial zeros by means of Aberth's method*, Numer. Algorithms **13** (1996) 179–200.
- O. Aberth, *Iteration methods for finding all zeros of a polynomial simultaneously*, Math. Comp. **27** (1973) 339–344.
- M. A. Jenkins and J. F. Traub, *Algorithm 419: Zeros of a complex polynomial*, Comm. ACM **15** (1972) 97–99 (CPOLY, ACM TOMS 419).
- Source: [`src/numerical_roots/nroots.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_roots/nroots.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_nroots.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nroots.c)

## Notes & additional examples

### Notes

`NRoots[lhs == rhs, x]` returns a disjunction `x == r1 || x == r2 || ...` of
numerical roots of a polynomial equation. A root of multiplicity `k` appears as
`k` identical equations, and a single root yields a bare equation. Real and
complex coefficients are handled.

`Method -> Automatic` uses the **Aberth-Ehrlich** simultaneous iteration (all
roots refined at once from Bini's convex-hull initial placement); an Automatic
machine request first tries LAPACK companion-matrix eigenvalues (`numpy.roots`
exactly) and falls back to Aberth. `Method -> "CompanionMatrix"` forces the
companion eigenvalue route (real QR directly, complex via a real `2n x 2n`
embedding), and `Method -> "JenkinsTraub"` uses the three-stage shifted-deflation
algorithm (CPOLY, ACM TOMS 419). An exact-integer polynomial is first squarefree-
decomposed (Yun), so high multiplicities like `(x^2 - 2)^30` stay well
conditioned.

`PrecisionGoal` selects precision: `Automatic` (or `Infinity`) gives machine
`double`/`Complex`; a positive digit count runs the whole solve in MPFR complex
arithmetic. `AccuracyGoal` (default `MachinePrecision`) drives only the post-solve
residual check — a root whose Newton correction exceeds the goal triggers an
`NRoots::accgl` warning. `NRoots` requires the MPFR build.
