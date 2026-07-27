# NumericFunction heads without a machine-precision fast path

_Generated 2026-07-27 from the build. Regenerate with:_

```bash
cd tests/build && COVERAGE_REPORT=1 ./compile_coverage_tests
```

## Why this matters

`compile_expr` returns `NULL` when it meets a head it cannot lower, the caller
quietly interprets, and the answer is still **correct** — just 10-40x slower. A
missing fast path is therefore invisible at the call site: it looks exactly like
a working one.

It is also not proportional. The compilable subset is a **cliff, not a slope**:
one unsupported head anywhere in a body costs the *entire* body, including every
other head in it that would have compiled. `Compile[{{z, _Complex}}, Sin[z] +
Gamma[z]]` does not lose the `Gamma`; it loses the `Sin` too, and the loop around
them.

So the question this document answers is not "which functions have a kernel" but
**"which bodies silently drop to the interpreter, and is that the right
answer"** — because sometimes it is. The compiled path must never answer where
the interpreter declines, nor differently; when the interpreter itself returns
`ComplexInfinity`, `Indeterminate`, or nothing at all, bailing is correct
behaviour rather than a gap.

## Summary

| | count | of 103 |
|---|---|---|
| Full fast path (real **and** complex) | 63 | 61% |
| Real only — **complex bails** | 30 | 29% |
| No fast path at all | 10 | 10% |

_(52 / 41 / 10 when first measured. Closed since: `Sign`, `FractionalPart`,
`Rescale`, `Gamma`, `LogGamma`, and the exponential-integral family —
`ExpIntegralEi`, `LogIntegral`, `SinIntegral`, `CosIntegral`, `SinhIntegral`,
`CoshIntegral`.)_

**The headline is the middle row, and it is the one least visible today.** Real
coverage went 55 → 93 of 103 during M5 and is effectively finished: all ten
remaining are deliberate. Complex coverage was never separately measured, and it
sits at 52 of 103 — so *the single largest remaining source of silent
interpreter fallback in the whole engine is complex arguments*, not exotic
functions.

---

## Class A — no fast path at all (10 heads)

All ten are deliberate. None is pending work.

### A1. The interpreter itself declines on machine reals (5)

A kernel here would make the compiled path **answer where the interpreter does
not**, which is the one divergence the engine forbids and a far worse defect
than a bail.

| head | interpreter on reals |
|---|---|
| `BesselJZero` | `BesselJZero[2., 1.]` unevaluated |
| `BarnesG` | `BarnesG[3.5]` unevaluated |
| `Hyperfactorial` | `Hyperfactorial[3.5]` unevaluated |
| `Factorial2` | `Factorial2[5.5]` unevaluated |
| `FactorialPower` | `FactorialPower[5.5, 2.]` unevaluated |

Closing any of these means **first** teaching the interpreter to evaluate them at
real (and complex) argument, and only then writing the kernel. Check this before
writing numerics for anything here: `Binomial[5.5, 2.]` *does* evaluate, to
12.375, which is why `Binomial` is covered and these are not.

### A2. Not a single machine number (5)

| head | why |
|---|---|
| `GCD`, `LCM`, `DigitSum` | exact-integer semantics a double cannot represent faithfully |
| `ReIm`, `QuotientRemainder` | return a two-element list, not one machine number |

`ReIm` and `QuotientRemainder` could in principle get bespoke lowerings that
produce two registers, but nothing in the type system expresses a multi-value
result today, and no measured workload wants it.

---

## Class B — real works, complex bails (41 heads)

**This is the real gap.** Every one of these compiles for `_Real` and drops the
whole body for `_Complex`.

```
AiryAi  AiryAiPrime  AiryBi  AiryBiPrime  Ceiling  Clip  CosIntegral
CoshIntegral  Erf  Erfc  Erfi  ExpIntegralEi  Factorial  Fibonacci  Floor
FractionalPart  FresnelC  FresnelS  Gamma  HarmonicNumber  Hypergeometric1F1
Hypergeometric2F1  HypergeometricPFQ  IntegerPart  InverseErf  InverseErfc
LerchPhi  LogGamma  LogIntegral  LucasL  Mod  ProductLog  Quotient  Rescale
Round  Sign  SinIntegral  Sinc  SinhIntegral  UnitStep  Zeta
```

The cause is uniform and structural: `sf_machine.c` kernels were registered
`real_closed` with a `real` function only and no `cplx`, so `try_kernel` has
nothing to emit for a complex operand. It is one missing function pointer per
head, not a missing algorithm — but the algorithms themselves differ, so this is
real numerical work, not a mechanical sweep.

### B1. Genuine gaps — the interpreter answers, the compiler bails (35)

Verified against the interpreter directly:

| head | `f[1.0 + 1.0 I]` (interpreter) |
|---|---|
| `Gamma` | `0.165915 + 0.149463 I` |
| `LogGamma` | `-1.4992 + 0.733281 I` |
| `Erf` | `1.31615 + 0.190453 I` |
| `Erfc` | `-0.316151 - 0.190453 I` |
| `Erfi` | `0.190453 + 1.31615 I` |
| `AiryAi` | `0.0604583 - 0.15189 I` |
| `AiryBi` | `0.716658 + 0.619889 I` |
| `ExpIntegralEi` | `1.76463 + 2.38777 I` |
| `LogIntegral` | `1.41126 + 1.22471 I` (at `2 + I`) |
| `SinIntegral` | `1.10422 + 0.882454 I` |
| `CosIntegral` | `0.882172 + 0.287249 I` |
| `SinhIntegral` | `0.882454 + 1.10422 I` |
| `CoshIntegral` | `0.882172 + 1.28355 I` |
| `FresnelC` | `2.55579 + 2.55579 I` |
| `FresnelS` | `-2.06189 + 2.06189 I` |
| `ProductLog` | `0.656966 + 0.32545 I` |
| `Fibonacci` | `-0.345569 - 1.79838 I` (at `1.5 + I`) |
| `LucasL` | `4.42158 + 5.92656 I` (at `1.5 + I`) |
| `HarmonicNumber` | `1.38699 + 0.457248 I` (at `1.5 + I`) |
| `Hypergeometric1F1` | `1.37802 + 0.909331 I` |
| `Hypergeometric2F1` | `1.12876 + 0.22088 I` |
| `HypergeometricPFQ` | `1.43609 + 1.1295 I` |
| `LerchPhi` | `1.08511 + 0.0337871 I` |
| `Sinc` | `0.966711 - 0.331747 I` |
| `Floor` | `1 + 2 I` (at `1.5 + 2.5 I`) |
| `Ceiling` | `2 + 3 I` |
| `Round` | `2 + 2 I` |
| `IntegerPart` | `1 + 2 I` |
| `FractionalPart` | `0.5 + 0.5 I` |
| `Sign` | `0.707107 + 0.707107 I` |
| `Quotient` | `2` (at `5.5 + I, 3.`) |
| `Rescale` | `0.5 + 0.5 I` |

Plus `Zeta`, `AiryAiPrime`, `AiryBiPrime` by the same argument.

Three tiers, and the first two are **not** the same size they first appeared:

- **~~Cheap and mechanical~~ — DONE (3):** `Sign`, `FractionalPart`, `Rescale`
  are closed. `Sign`'s complex kernel had been in the registry the whole time;
  the compiler's own inline lowering for `Sign` was shadowing it and bailing
  before `try_kernel` was ever consulted.

- **TYPE-BLOCKED, not mechanical (4):** `Floor`, `Ceiling`, `Round`,
  `IntegerPart` on a complex return **`Complex[Integer, Integer]`** —
  `Head /@ ReIm[Floor[1.5 + 2.5 I]]` is `{Integer, Integer}`. The compile
  engine's type lattice has `CT_INT` and `CT_COMPLEX` (a `double _Complex`) but
  nothing for a complex *integer*, so answering with `CT_COMPLEX` would produce
  `Complex[Real, Real]` and differ from the interpreter in HEAD. These need a
  new type, not a new function pointer. (This is why the tier was originally
  listed as eight cheap heads and is not: worth checking result *types*, not
  just values, before calling anything mechanical.)

- **~~`Gamma`/`LogGamma`~~ — DONE.** Both already had a `double complex` Lanczos
  series inside the *interpreter* (`gamma.c`, `loggamma.c`); it is now exposed
  and shared rather than reimplemented, so the compiled and interpreted paths
  agree bit for bit by construction. Fixing the branch cut was the actual work —
  see below. `Beta`, `Binomial`, `Pochhammer` and `Hypergeometric0F1` turned out
  to have complex kernels already.

- **~~The exponential-integral family~~ — DONE (6).** `ExpIntegralEi`,
  `LogIntegral`, `SinIntegral`, `CosIntegral`, `SinhIntegral`, `CoshIntegral`.
  Again no new numerics: each already had a `double complex` ascending series in
  its module, dead behind `#ifndef USE_MPFR`. The work was a **cancellation
  gate** — the series converges everywhere but is only usable where the terms
  do not dwarf the value they sum to. Peak term over result is measured, and the
  kernel declines above a 1e3 budget (~10 bits), which holds the error to
  ~1e-13 while still covering `|z| <= 12`–32 depending on the function.
  `LogIntegral` is `Ei[Log z]` with a principal log, matching the symbolic path.

- **Genuine complex numerics (~20 left):** the rest of the special functions.
  The remaining exponential-integral work
  needs branch-cut handling; the hypergeometrics need the complex series (which
  the pFq machine kernel is already close to, since it sums with complex terms
  internally). Each needs a parity test against the MPFR path, and each must
  respect the same branch cuts the interpreter uses — a kernel that picks a
  different branch is worse than no kernel.

`Quotient` came off this list for a different reason — see the interpreter bug
note below.

### B2. Correct declines — the interpreter declines too (6)

Not gaps. Leave them.

| head | interpreter at complex |
|---|---|
| `Mod[5.5 + I, 3.]` | unevaluated |
| `UnitStep[1. + I]` | unevaluated |
| `Clip[1. + I, {0., 2.}]` | unevaluated |
| `Factorial[1.5 + I]` | unevaluated (note `Gamma` **does** evaluate — an interpreter inconsistency, not a compiler one) |
| `InverseErf[0.3 + 0.1 I]` | unevaluated |
| `InverseErfc[0.3 + 0.1 I]` | unevaluated |

---

## Class C — a fast path that covers only part of its domain

These compile, and then decline at run time for some inputs, which sends that
*call* to the interpreter. Cheaper than a bail (the program still exists) but the
same invisibility.

### C1. ~~Airy, `2.5 < |x| < 8`~~ — CLOSED 2026-07-27

The ascending series and the asymptotic expansion do not meet in double
precision (the series lost roughly `2ζ/ln 10` digits computing a decaying `Ai`
as a difference of growing terms; the asymptotic bottomed out at ~1e-5 at
`|x| = 3.5`). The band is now covered by a **third method: Taylor marching of
`y'' = x y`**, seeded from whichever expansion is exact at the nearer end.
Errors across the former band are 1e-16 to 2.6e-15 — better than the ascending
series achieves at `|x| = 2.5`. See the changelog for why the marching
*direction* is the whole problem.

No Airy argument on the real axis declines any more.

### C2. Real input, complex output (correct as long as the declared type is Real)

A `real_closed` kernel must decline where the true value is not real, because
the compiled program's declared result type is `Real` and returning the real part
would be silently wrong:

| call | interpreter |
|---|---|
| `CosIntegral[-1.]` | `0.337404 + 3.14159 I` |
| `CoshIntegral[-1.]` | `0.837867 + 3.14159 I` |
| `LogIntegral[-1.]` | `0.0736679 + 3.42273 I` |
| `ProductLog[-0.5]` | `-0.794024 + 0.770112 I` |
| `LogGamma[-1.5]` | `0.860047 - 6.28319 I` |
| `PolyLog[2, 2.]` | `2.4674 - 2.17759 I` |

These become non-declines only once Class B is done **and** `infer_type` widens
the result to `CT_COMPLEX` for a real argument in those regions — which it
cannot do statically, since the region depends on the runtime value. The honest
fix is a complex-typed body; the real-typed one is right to decline.

### C3. Correct declines: the answer is not a machine number

| call | interpreter |
|---|---|
| `LogIntegral[1.]` | `-Infinity` |
| `LogIntegral[0.]` | `Indeterminate` |
| `ExpIntegralEi[0.]` | `-Infinity` |
| `Gamma[-3.]` | `ComplexInfinity` |
| `Zeta[1.]` | `ComplexInfinity` |
| `Beta[0.5, -1.]` | `ComplexInfinity` |

Nothing to do — a double cannot carry these, and the engine correctly routes the
call back to the interpreter which can.

---

## Recommended order

1. ~~**Airy `2.5 < |x| < 8`**~~ — **done**, see C1.
2. ~~**The mechanical complex kernels**~~ — **done** for the three that are
   cleanly typed (`Sign`, `FractionalPart`, `Rescale`). The other four
   (`Floor`, `Ceiling`, `Round`, `IntegerPart`) turned out to be type-blocked,
   not mechanical: they return `Complex[Integer, Integer]`, which the type
   lattice cannot express. Closing them means adding a complex-integer type —
   a real design change, and worth doing only if something wants it.
3. ~~**Complex `Gamma`/`LogGamma`**~~ — **done**. Look for an existing
   `static double complex` implementation inside the interpreter before writing
   a kernel: both were already there, and so was `Sign`'s. The work that
   remained was a branch-cut bug, not numerics.
4. ~~**Complex exponential-integral family**~~ — **done**, gated to the region
   where a double can carry the answer. **Extending that region is the natural
   follow-on**: the real `sf_machine_ei` already switches to a continued
   fraction past `|x| = 40`, and the same
   `E1(z) = e^-z / (z + 1 - 1^2/(z + 3 - 2^2/(z + 5 - ...)))` converges for
   complex `z` with `|arg z| < Pi`. With `Ei(z) = -E1(-z)`,
   `Chi +/- Shi = Ei(+/-z)` and `Ci`/`Si` from `E1(+/-iz)`, one continued
   fraction would lift the whole family's declines — with branch care at each
   relation. Complex hypergeometrics are the other bulk item.
5. **Interpreter first, then kernel** for Class A1, if those functions are ever
   wanted at machine precision.

## Invariants any new kernel must hold

- **Never answer where the interpreter declines.** Check with the interpreter
  before writing numerics — that is how `Factorial2`, `FactorialPower`,
  `BesselJZero`, `BarnesG` and `Hyperfactorial` turned out to be exclusions.
- **Never answer differently.** Same branch cuts, same principal values.
- **Decline out of domain rather than approximate.** Shipping a fast path that is
  quietly wrong to five digits is far worse than being slower over one interval.
- **Parity-test against the MPFR path**, a few hundred points, and check *which
  side is wrong* when it fails — these tests have now found three interpreter
  bugs (`ProductLog` on `[0.35, 1/e]`, `Zeta` at 0, `HypergeometricPFQ` at
  negative real `z`) and no kernel bugs.
