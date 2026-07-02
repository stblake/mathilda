# PossibleZeroQ — Flaky / Failing Cases

**Date:** 2026-06-25
**Component:** `src/zero_test.c` (`PossibleZeroQ`, hybrid symbolic–numeric zero recognition)

This file collects inputs on which `PossibleZeroQ` is **flaky** — it returns
`True` on some process runs and `False` on others for the *same* input — or
**consistently wrong** (returns `False` on an expression that is identically
zero). These undermine any caller that trusts the verdict (e.g. the
`Integrate` verification gate, `integrate_jeffrey_tests`).

## How the flakiness arises

`PossibleZeroQ` falls back to a **Schwartz–Zippel random sampler**
(`decide_schwartz_zippel`) for expressions with free symbols. It draws random
(optionally complex) sample points for each free symbol, substitutes, and
numerically tests whether the result is zero:

- The global PRNG (`random.c`) is seeded **non-deterministically per process**,
  so the same input can sample different points on different runs → different
  verdicts. (This is why `integrate_jeffrey_tests` reports 0 or 1 failures on
  repeated identical runs.)
- The builtin collapses `UNKNOWN → True` and returns `False` **only** when the
  sampler decisively classifies a point as non-zero. So for a genuine zero,
  **any `False` is a sampling false-negative.**

Two distinct false-negative mechanisms appear in the cases below:

1. **Catastrophic cancellation.** When a true-zero expression is a difference
   of two large nearly-equal quantities (e.g. `Gamma[x+1]` vs `x*Gamma[x]` at a
   complex sample where `Gamma` is large), the machine/MPFR residual at the
   sampled point can exceed the zero threshold, so the point reads as non-zero.
   The precision ladder (`PRECISION_LADDER = {53,200,500,1000}`) is supposed to
   shrink a genuine zero geometrically, but the verdict still depends on the
   sampled magnitude.
2. **Branch-cut / saturation.** Complex samples (imaginary part up to
   `2^ZT_IMAG_NUMERATOR_BITS`, i.e. `±8`) push transcendental heads
   (`Tan`, `Log`, `Tanh`, antiderivatives with `ArcTan`/`Log`/`Floor`) across
   branch cuts where the *symbolic* identity legitimately does not hold, so the
   numeric residual is genuinely non-zero at that point.

## Measurement method

Each case is `PossibleZeroQ[<expr>]` where `<expr>` is **identically zero**
(unless flagged otherwise). The battery was run **40 times in separate
processes** (fresh seed each) and the `True`/`False` verdicts tallied. A case is
**FLAKY** if both `True` and `False` occur; **ALWAYS-FALSE** if every run says
`False`.

Reproduce:

```bash
M=./Mathilda
for i in $(seq 1 40); do printf 'PossibleZeroQ[Gamma[x + 1] - x Gamma[x]]\n' | "$M" 2>/dev/null | grep '^Out'; done | sort | uniq -c
```

## Confirmed flaky cases (genuine zeros)

| # | Expression (identically zero) | True/False over 40 runs | Mechanism |
|---|-------------------------------|-------------------------|-----------|
| 1 | `Gamma[x + 1] - x Gamma[x]` | 20 / 20 | cancellation (Gamma large at complex samples) |
| 2 | `D[Integrate[1/(5 + 3 Cosh[x]), x] /. Floor[_] -> 0, x] - 1/(5 + 3 Cosh[x])` | 19 / 21 | branch-cut (antiderivative has `ArcTan`/`Log`) — the `integrate_jeffrey_tests` flake |
| 3 | `D[Integrate[Cosh[x] Cosh[2 x], x, Method -> "Weierstrass"] /. Floor[_] -> 0, x] - Cosh[x] Cosh[2 x]` | 21 / 19 | branch-cut — the other `integrate_jeffrey_tests` flake |

## Branch-dependent (not a clean identity — `False` is defensible but still non-deterministic)

| # | Expression | True/False over 40 runs | Note |
|---|------------|-------------------------|------|
| A | `Log[Exp[x]] - x` | 23 / 17 | only holds on the principal strip; complex samples off it are genuinely non-zero |

## Stable controls (no flakiness observed — for regression reference)

All `True` 40/40: `Sin[x]^2+Cos[x]^2-1`, `Sin[2x]-2 Sin[x] Cos[x]`,
`Cos[2x]-1+2 Sin[x]^2`, `Cosh[x]^2-Sinh[x]^2-1`, `Cosh[2x]-1-2 Sinh[x]^2`,
`Tan[x] Cos[x]-Sin[x]`, `Sin[x+y]-Sin[x]Cos[y]-Cos[x]Sin[y]`,
`Exp[x+y]-Exp[x]Exp[y]`, `Sin[3x]-3Sin[x]+4Sin[x]^3`,
`Cos[3x]-4Cos[x]^3+3Cos[x]`, `Tan[2x]-2Tan[x]/(1-Tan[x]^2)`,
`Cot[x]-Cos[x]/Sin[x]`, `Tan[x]+Cot[x]-1/(Sin[x]Cos[x])`,
`Sinh[2x]-2Sinh[x]Cosh[x]`, `Tanh[x]-Sinh[x]/Cosh[x]`, `Sqrt[x]^2-x`,
`D[Integrate[1/(2+Sin[x]),x] /. Floor[_]->0, x] - 1/(2+Sin[x])`.

Notably the `Tan`-heavy identities (#11, #13) and one integral roundtrip
(#21 = `1/(2+Sin[x])`) were **stable** here, so flakiness is input-specific,
not a blanket property of `Tan`/integral roundtrips.

---

## 2026-07-01 — cases blocking `Integrate\`GoursatAlgebraic` verification

The square-root Goursat descent verifies its candidate antiderivative with
`PossibleZeroQ[D[result, x] - f]` (this is the correct design — verification
should use `PossibleZeroQ`, never `Simplify`). Two cases were flagged. **Resolved
2026-07-02**: B1 turned out **not** to be a `PossibleZeroQ` defect, and B2 was a
performance issue that is now fixed. The descent keeps its positive-real numeric
gate (`diff_back_ok`) because it encodes domain knowledge `PossibleZeroQ` cannot
and should not assume (see B1).

### Case B1 — RESOLVED: not a bug (the integral is only valid on the positive real axis)

```
f  = (t - (-1)^(2/3)) / ((t + 2 (-1)^(2/3)) Sqrt[t^3 - 1])
r  = Integrate[f, t, Method -> "GoursatAlgebraic"]   (* period-3 reduction, 139 leaves *)
PossibleZeroQ[D[r, t] - f]   ==>  False              (* CORRECT, not a false negative *)
```

`D[r,t] - f` is **not** identically zero on `ℂ`. The antiderivative `r` is only
valid on the branch it was built on — the positive real axis, where `Sqrt[t^3-1]`
is on its principal sheet (`t > 1`). Off that branch it is a genuine non-zero:

| point | `Abs[D[r,t] - f]` |
|-------|-------------------|
| `t = 3.4`  (on branch)  | `1.9e-14`  (zero) |
| `t = -3.4` (off branch) | `0.20`     (**non-zero**) |

This is the **identical numeric signature** of `Sqrt[x^2] - x` and
`Log[x^2] - 2 Log[x]` — zero on the positive axis, non-zero on the negative — both
of which the test suite *requires* to return `False`. No sampling policy (real,
complex, or exact-algebraic point evaluation) can return `True` for B1 while
keeping those `False`: they are the same equivalence class. Only genuine
algebraic-function-field zero testing (radicals as field generators over
`Q(t)(ζ_n)(√…)`) could distinguish them, and that is out of scope.

So `PossibleZeroQ` returning `False` is correct. The Goursat descent verifies its
antiderivative with `diff_back_ok` (positive-real fixed samples) precisely because
it *knows* the answer is a positive-real-axis identity — domain knowledge a general
zero-test must not assume. That gate is the right design and is retained.

### Case B2 — FIXED (performance): large parametric outputs now settle in ~1.9 s

```
f = (1 + k x) / ((-1 + k x) Sqrt[(1 - x) x (1 - k^2 x)])
r = Integrate[f, x, Method -> "GoursatAlgebraic"]   (* parametric V4, 9495 leaves *)
PossibleZeroQ[D[r, x] - f]   ==>  True   (correct)   4.8 s  →  1.9 s
```

`PossibleZeroQ` was always **correct** here; the cost was climbing the full
`{53,200,500,1000}`-bit precision ladder on the ~60 700-leaf residual tree at every
confirm sample. Fixed in `src/zero_test.c` by (a) a **deep-zero early-exit**: once a
residual has both shrunk geometrically and fallen far below any plausible
cancellation floor, the identity is settled `True` without the 500/1000-bit rungs;
and (b) reusing the precision-independent operand scale across rungs instead of
recomputing it (a redundant second numericalisation pass per rung). Every verdict
is unchanged; the cost drops to ~1.9 s, comfortably inside the descent's 6 s budget.

**Reproduce** (both): set `Method -> "GoursatAlgebraic"`, take the printed `r`,
and evaluate `Timing[PossibleZeroQ[D[r, var] - f]]`.
