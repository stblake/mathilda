# NSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NSum[f, {i, imin, imax}]`**

gives a numerical approximation to the sum of f for i from imin to imax.

**`NSum[f, {i, imin, imax, di}] uses step di. imax may be Infinity. NSum[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional sum (an inner bound may depend on an outer index). The index is localised (HoldAll). Method -> Automatic picks Euler-Maclaurin for a monotone series whose summand extends off the integers, the Cohen-Villegas-Zagier method for alternating series, and Wynn's epsilon (partial-sum acceleration) otherwise, with Levin's u-transform as a last resort; an integer-only summand such as 1/Prime[n] cannot use Euler-Maclaurin (it has no continuous tail integral) and is extrapolated instead. Large finite sums use the difference of two infinite tails. Method -> "Levin" forces Levin's transformation ("LevinU" | "LevinT" | "LevinV" select the u/t/v variant). Machine or arbitrary precision via WorkingPrecision.`**

<details>
<summary>Notes</summary>

Options: Method (Automatic | EulerMaclaurin | AlternatingSigns | WynnEpsilon | "Levin"), WorkingPrecision (default MachinePrecision), NSumTerms (head terms summed explicitly, default 15), NSumExtraTerms, WynnDegree, VerifyConvergence (default True; a divergent sum gives ComplexInfinity), AccuracyGoal, PrecisionGoal.

</details>

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= NSum[1/i^2, {i, 1, Infinity}] - Pi^2/6 // N
Out[1]= 2.22045e-16

In[2]:= NSum[1/2^i, {i, 0, Infinity, 2}]
Out[2]= 1.33333

In[3]:= NSum[Log[x]/x^(2 + 2 I), {x, 1, Infinity}]
Out[3]= -0.182175 - 0.136618*I

In[4]:= NSum[1/i^2, {i, 100, 10^6}]
Out[4]= 0.0100492

In[5]:= NSum[(-1)^n (2/n)^k/k^2, {n, 2, Infinity}, {k, 1, n}]
Out[5]= 0.770188

In[6]:= NSum[2^i, {i, 0, Infinity}] NSum::div: the sum does not appear to converge
Out[6]= Optional[ComplexInfinity NSum::div, appear converge does not sum the to]
```

### Options (4)

```mathematica
In[7]:= NSum[(-5)^i/i!, {i, 0, Infinity}, NSumTerms -> 25] - Exp[-5]
Out[7]= -2.4182e-15

In[8]:= NSum[1/n^(11/10), {n, 1, Infinity}, WorkingPrecision -> 40] - Zeta[11/10]
Out[8]= -2.9387358770557187699218413430556141945467e-39

In[9]:= NSum[(-1)^x/(1 + (x - 12)^2), {x, 0, Infinity}, Method -> "AlternatingSigns", WorkingPrecision -> 30]
Out[9]= 0.2751938594139530395689715615907

In[10]:= NSum[1/n^2, {n, 1, Infinity}, Method -> "Levin", WorkingPrecision -> 30]
Out[10]= 1.644934066848226436472415166646
```

### Applications (5)

```mathematica
In[11]:= NSum[1/n^2, {n, 1, Infinity}]
Out[11]= 1.64493

In[12]:= NSum[(-1)^(n+1)/n, {n, 1, Infinity}]
Out[12]= 0.693147

In[13]:= NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[13]= 1.644934066848226436472415166646

In[14]:= NSum[1/n^4, {n, 1, Infinity}, WorkingPrecision -> 30]
Out[14]= 1.082323233711138191516003696546

In[15]:= NSum[1/n^2, {n, 1, Infinity}, Method -> "Levin"]
Out[15]= 1.64493
```

## Algorithm

```text
nsum.c — NSum[f, {i, imin, imax (, di)}, opts]   (see nsum.h)
```

Strategy -------- NSum holds its arguments, evaluates the iterator bounds, then Block-localises

```text
the index and evaluates the summand once per term.  Terms are reindexed to
```

k = 0, 1, 2, … with the actual index value x_k = imin + k·di, so a step di is handled uniformly and multidimensional sums fall out by making the summand of the outer sum an inner NSum[...] (HoldAll + localisation lets a dependent inner bound such as {k,1,n} see the bound outer index).

Methods are layered: this file currently provides Direct (small finite sums) and WynnEpsilon (partial-sum extrapolation, shared seqaccel kernels), machine

```text
and MPFR, real and complex.  Euler–Maclaurin and Cohen–Villegas–Zagier are
```

added on top of the same term machinery.

Memory: receives `res` owned by the evaluator; returns a fresh Expr* on

```text
success or NULL (unevaluated).  Never frees `res`.  Every temporary index
```

binding is removed on all return paths.

## Implementation notes

**Algorithm.** `builtin_nsum` (HoldAll) parses the iterator `{i, imin, imax,
di}`, reindexes terms to `k = 0, 1, 2, …` with `x_k = imin + k·di` built and
evaluated under a Block-style index binding, and sums the first `NSumTerms`
(default 15) explicitly before accelerating the tail. `ns_choose_method` picks
from a sampled profile:

- **Euler–Maclaurin** (`ns_em_*`) for a monotone summand defined off the
  integers: `Σ ≈ (1/di)∫_N^∞ f + f(N)/2 − Σ_j B_{2j}/(2j)! · di^{2j-1}
  f^{(2j-1)}(N)`. The tail integral is a dedicated double-exponential (exp-sinh)
  quadrature (`dequad_halfline_*`), **not** `NIntegrate`; the derivative
  corrections are a hybrid — symbolic `D` while the derivative tree stays small,
  switching to Taylor coefficients recovered from a circle DFT (Cauchy's
  formula, as `NSeries` does) once it balloons — and the correction series is
  truncated asymptotically at its smallest term.
- **Cohen–Villegas–Zagier** (`ns_cvz_*`) for a strictly alternating real series:
  Chebyshev weights `d_n = ((3+√8)^n + (3+√8)^{-n})/2` in a single linear pass
  giving ≈ 2.54 n bits.
- **Wynn's epsilon** (`seqaccel.c`, the iterated Shanks transform) otherwise,
  with **Levin's u/t/v transform** as a last resort for logarithmically/
  algebraically convergent tails.

An integer-only summand such as `1/Prime[n]` is detected *behaviourally* — the
summand is probed at two non-integer points, and if neither is finite,
Euler–Maclaurin is forbidden (no continuous tail) and the series is extrapolated
instead (with a larger head-term count). Multidimensional sums nest an inner
`NSum` as the summand (dependent inner bounds see the outer index via HoldAll).
A large finite sum of decaying terms is computed as the difference of two
infinite tails. With `VerifyConvergence -> True` (default), a summand whose
magnitude is still rising far into the sampled tail emits `NSum::div` and
returns `ComplexInfinity`.

`NProduct` (`nprod.c`) is evaluated as **`Exp[NSum[Log[f], …]]`** (Keiper 1992):
it maps `NProductFactors -> NSumTerms` (so the default factor count is also 15),
`NProductExtraFactors -> NSumExtraTerms`, runs the inner `NSum` at ten guard
digits above the request (because `Exp` turns the exponent's absolute error into
the product's relative error), and rounds the final `Exp` back. A divergent
log-sum propagates as `ComplexInfinity`.

**Data structures.** Two parallel implementations gated on `USE_MPFR`. The
machine path is `double _Complex` throughout — partial-sum sequence `P[]`, the
Wynn ε-table `(terms+1)²`, Levin's `a`/`omega` arrays. The MPFR path keeps split
`mpfr_t` real/imag buffers at `2·target + 32` internal bits (absorbing
cancellation in near-1 Euler–Maclaurin samples), rounded back to the target.
`WorkingPrecision` (`MachinePrecision` or a digit count) selects between them;
the extrapolation sequence length scales with the bit count. A machine request
draws every term through one of three cached compiled programs (real, at
precision, and complex-input for the contour), so a single compile serves all
methods; the MPFR path stays on the interpreter per term.

**Complexity / limits.** Linear in the head terms plus `O(seq²)` for the
ε-table (`NS_MAX_SEQ = 64`). Options: `Method` (`Automatic` | `EulerMaclaurin` |
`AlternatingSigns` | `WynnEpsilon` | `"Levin"`/`"LevinU"`/`"LevinT"`/`"LevinV"`),
`WorkingPrecision`, `NSumTerms`, `NSumExtraTerms`, `WynnDegree`,
`VerifyConvergence`, `AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal`
(default Automatic); a shortfall against the combined tolerance warns
`NSum::ncvg`. NSum declines (stays unevaluated) on a non-numeric finite bound, a
summand that never numericalises, or a non-decaying finite sum beyond
2 000 000 terms. Per-term internal probes are muted via
`arith_warnings_mute`, and all diagnostics route through `mth_message`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Block](../../scoping-constructs/Block/), [Chop](../../elementary-functions/Chop/), [Integrate](../../calculus/Integrate/), [D](../../calculus/D/), [BernoulliB](../../special-functions/BernoulliB/), [NLimit](../../numerical-calculus/NLimit/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- H. Cohen, F. Rodriguez Villegas and D. Zagier, *Convergence acceleration of alternating series*, Experiment. Math. **9** (2000) 3–12.
- D. Levin, *Development of non-linear transformations for improving convergence of sequences*, Internat. J. Comput. Math. **3** (1973) 371–388.
- P. Wynn, *On a device for computing the e_m(S_n) transformation*, MTAC **10** (1956) 91–96 — the epsilon algorithm.
- J. B. Keiper, *Numerical computation of infinite products*, Wolfram Research tech. report (1992) — the Exp[NSum[Log]] reduction used by NProduct.
- Source: [`src/numerical_calculus/nsum.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/nsum.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_nprod.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nprod.c)
- Tests: [`tests/test_nsum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsum.c)

## Notes & additional examples

### Notes

`NSum[f, {i, imin, imax}]` numerically sums a series, with `imax` allowed to be
`Infinity`. The first two cases recover the Basel sum `Pi^2/6 = 1.64493...` and
the alternating harmonic sum `Log[2] = 0.693147...`. With `WorkingPrecision -> 30`
the Basel sum is computed to 30 digits, and `Sum[1/n^4]` returns
`Pi^4/90 = 1.082323233711...`. `Method -> Automatic` chooses Euler–Maclaurin for
monotone series, the Cohen–Villegas–Zagier method for alternating series, and
Wynn's epsilon otherwise, with Levin's u-transform as a last resort. Any
accelerator can be forced: `Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`)
selects Levin's transformation, which reaches full `WorkingPrecision` on smooth
series. With `VerifyConvergence -> True` (default) a divergent sum gives
`ComplexInfinity`.
