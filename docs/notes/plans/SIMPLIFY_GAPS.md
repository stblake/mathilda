# Simplify gaps surfacing in the Risch integrator

> **RESOLVED (2026-07-15) — Families 1, 2, 3.** A fast-path exact trig/exp-kernel zero
> test (`src/simp/simp_trigexp_zero.c`, FP2) now proves these diff-back identities by exact
> rational point-evaluation on a Nullstellensatz grid (no numeric sampling): kernelize
> `Sin[k x]/Cos[k x]/…` to `t = E^(i x)`, treat `Log[...]` subterms and parameters as
> independent generators, and certify a pole-free zero grid. Wired into both `Simplify`
> (top-level fast path, `simp_builtins.c`) and `PossibleZeroQ` (a symbolic `zero_test.c`
> stage). `Simplify[D[G3] − Sec^3] → 0` (was non-terminating > 40 s); the symbolic
> `∫1/(a+b Sin x)` diff-back no longer hangs (Family 3). Family 2's `½Log[1+Tan²]` collapse
> is handled by the FP1 trig-log canonicalization (`transform_trig_log_canon`) up to the
> principal-branch `−Log[Cos]` step, which remains assumption-gated. See the 2026-07-13
> changelog. Tests: `tests/test_trigexp_zero.c`.

**Date:** 2026-07-15
**Scope:** Cases where the transcendental Risch integrator produces a *correct* closed form
but `Simplify`/`FullSimplify` cannot reduce it to the expected/real shape (or to prove a
diff-back). These are **`Simplify` deficiencies, not integrator bugs** — the antiderivatives
are correct (numerically verified). Three families are documented: (1) the `Sec`/`Csc`-power
exp-route forms, (2) the hypertangent `½Log[1+Tan²] = −Log[Cos]` collapse, and (3) the
symbolic-parameter I-laden diff-back blow-up that makes the integrator *hang*.

---

## Family 3 — `Simplify[D[G] − f]` blows up (hang) on symbolic-parameter I-laden forms

`∫ 1/(a + b Sin[x]) dx` with SYMBOLIC parameters `a, b`, via
`Method->"RischTranscendental"`, runs **> 90 s** (CPU-bound). The **numeric-coefficient**
sibling `∫ 1/(2 + 3 Sin[x]) dx` returns in milliseconds. The hang is NOT in the integration
— `TrigToExp`, the kernelization, the Rothstein–Trager resultant / `TranscendentalLogPart`,
and the LogToReal reconstruction all finish in **< 0.2 s** and produce the *correct*
antiderivative

```
G = (Sqrt[-16 a^2 + 16 b^2] · ArcTan[(b^3 Sqrt[-16 a^2+16 b^2] - a^2 b Sqrt[-16 a^2+16 b^2]
        + 4 a^2 b^2 E^(I x) - 4 b^4 E^(I x)) / (4 a^3 b - 4 a b^3)]) / (4 a^2 - 4 b^2)
```

The hang is **entirely** in the correctness gate `rt_verify_antideriv`, i.e.
`Simplify[D[G, x] − 1/(a + b Sin[x])] === 0`. This difference **is** identically `0`, but
`Simplify` does not terminate in reasonable time on the symbolic `a, b` form.

**Root cause.** To reduce `D[G] − f` to `0`, `Simplify` must (i) differentiate the `ArcTan`
(giving a rational function of `E^(I x)`, `a`, `b`), (ii) recognize `E^(2 I x) = (E^(I x))²`
to combine the integrand's kernel with `G`'s, and (iii) cancel a large rational function
over the symbolic parameter field `Q(a, b)[√(b²−a²)]`. Its transformation search explodes
on the symbolic parameters — the numeric case collapses immediately because the coefficients
are concrete. A naive fast substitute (`Together[(D[G] − f) /. E^(I x) → w]`) does **not**
work: `E^(2 I x)` is a distinct kernel that the substitution misses, so the rational identity
is not exposed — which is precisely the recognition step `Simplify` performs (correctly, but
slowly).

**Fix directions.**
1. **A fast exact diff-back for a single-exponential-kernel result** (in the integrator):
   kernelize *both* `D[G]` and `f` with `rt_exp_kernelize` (`E^(k x) → t^k`, so
   `E^(2 I x) → t²`), then test the resulting rational function of `{x, t, a, b}` for a zero
   numerator via `Together` — exact and O(polynomial), no `Simplify` search. This is the
   same exact-coordinate idea the tower cases already use; applying it to `rt_frac_lrt`'s
   gate removes the hang and keeps the (correct) closed form.
2. **A `Simplify` improvement** for symbolic-parameter rational-of-exponential + `ArcTan`
   differences (multiple-angle exponential recognition + parameter-field cancellation),
   which would also make the generic diff-back fast.

Until one lands, the integrator **hangs** rather than declining — a robustness bug in the
verification gate (the integration itself is correct and fast).

## Family 2 — hypertangent `½Log[1+Tan[u]²]` not collapsing to `−Log[Cos[u]]`

`∫ Tan[Log[x]]/x dx` (and the direct `∫Tan[x]`-family through the §5.10 hypertangent driver)
comes back as the **raw hypertangent form**

```
Integrate[Tan[Log[x]]/x, x, Method->"RischTranscendental"]  ==>  (1/2) Log[1 + Tan[Log[x]]^2]
```

which is **correct** — real, and `Simplify[D[½Log[1+Tan[Log[x]]²]] − Tan[Log[x]]/x] === 0`
(the diff-back of the *actual* form reduces, so the integrator's own gate and the test's
`assert_tan_real` both pass). But it does **not** collapse to the canonical `−Log[Cos[Log[x]]]`:

- **`Simplify` FAILS:** `Simplify[½Log[1 + Tan[Log[x]]²] − (−Log[Cos[Log[x]]])]` does not
  return `0`, though it is identically `0` (`1+Tan² = Sec²`, `½Log[Sec²] = Log[Sec] =
  −Log[Cos]`).

**Root cause.** `Simplify`/`TrigReduce` lacks the chain `1+Tan[u]² → Sec[u]²`, then
`½Log[Sec[u]²] → Log[Sec[u]] = −Log[Cos[u]]` (a Pythagorean-identity + `Log`-of-power +
`Sec→1/Cos` rewrite) — the same multiple-angle/trig-log weakness as Family 1.

**Fix directions (in `Simplify`, not the integrator).** (a) A `Simplify`/`TrigReduce`
rule `½Log[1+Tan[u]²] → −Log[Cos[u]]` (via `1+Tan²→Sec²`, `Log[z²]→2Log[z]`, `Sec→1/Cos`);
or (b) integrator-side, apply that collapse to the hypertangent driver's `c·Log(t²+1)`
output before returning. The integrator is already correct; this is purely a canonical-form
improvement. (The `test_real_hypertangent` suite asserts correctness via `assert_tan_real`,
not the specific `−Log[Cos]` string, for exactly this reason.)

---

## Family 1 — `∫ Sec[x]^n`, `∫ Csc[x]^n` (n ≥ 3)

**Context:** The transcendental Risch integrator's real/exp-route reconstruction is gated by
the **exact symbolic diff-back** `Simplify[D[G, x] − f] === 0` (`rt_verify_antideriv`). The
former **numeric** sampling gate (`rt_realify_numverify`, 4 fixed points) was removed —
a Risch decision procedure must never certify by numeric sampling. The consequence: any
integrand whose Risch antiderivative `G` is *correct* but whose diff-back `Simplify` cannot
reduce to `0` now **declines cleanly** (returns unevaluated) instead of shipping a
numerically-rubber-stamped answer.

These declines are **not integrator bugs** — the integrator builds a correct closed form
(verified numerically below) — they are **`Simplify` deficiencies**. Closing them is a
`Simplify`/`TrigReduce` capability improvement; until then the integrator correctly declines
(never wrong). Companion: `SIMPLIFY_DEFICIENCIES.md` (the I-laden `Tan`→real gap).

The exp route (`rt_trig_frontend` → `rt_exp_ratreduce_case`, kernel `t = E^(I x)`) produces
a correct antiderivative in a **multiple-angle** basis — `Sin[3x]`, `Sin[5x]`, `Cos[4x]`,
`Cos[6x]`, `Log[2 ± 2 Sin[x]]`, … The diff-back of that form is mathematically `0` but
`Simplify`/`FullSimplify` cannot reduce it (it does not return `0`; on the `Sec[x]^3` form
below it runs **> 40 s without terminating**).

### `∫ Sec[x]^3 dx` — the form the integrator builds

```
G3 = (8 Sin[x] + 4 Sin[5 x] + 12 Sin[3 x]
      + (-10 - 15 Cos[2 x] - 6 Cos[4 x] - Cos[6 x]) Log[2 - 2 Sin[x]]
      + ( 10 + Cos[6 x] + 6 Cos[4 x] + 15 Cos[2 x]) Log[2 + 2 Sin[x]])
     / (40 + 4 Cos[6 x] + 24 Cos[4 x] + 60 Cos[2 x])
```

- **Correct (numeric diff-back = 0):** `Chop[N[(D[G3,x] − Sec[x]^3) /. x -> 11/10]] === 0`,
  and `=== 0` at `x -> 21/10`, `7/10`.
- **`Simplify` FAILS:** `Simplify[D[G3, x] − Sec[x]^3]` does **not** return `0` (runs
  > 40 s). `FullSimplify` likewise.
- **A clean equivalent DOES reduce:** `Simplify[D[(Sec[x] Tan[x] + Log[Sec[x] + Tan[x]])/2,
  x] − Sec[x]^3] === 0`. So the failure is specific to the **multiple-angle basis**, not
  to `Sec[x]^3` itself.

### `∫ Csc[x]^3 dx` — the analogous form

```
G3c = (-8 Cos[x] - 4 Cos[5 x] + 12 Cos[3 x]
       + (-10 + Cos[6 x] - 6 Cos[4 x] + 15 Cos[2 x]) Log[2 + 2 Cos[x]]
       + ( 10 - 15 Cos[2 x] - Cos[6 x] + 6 Cos[4 x]) Log[2 - 2 Cos[x]])
      / (40 - 60 Cos[2 x] - 4 Cos[6 x] + 24 Cos[4 x])
```

Same profile: numerically a correct antiderivative of `Csc[x]^3`; `Simplify` of the
diff-back does not reduce to `0`. `Sec[x]^4` and higher `Sec`/`Csc` powers are in the same
family (each raises the multiple-angle degree).

Reachable regression (both close on clean HEAD via the old numeric gate; both now decline):
`Integrate[Sec[x]^3, x, Method->"RischTranscendental"]`,
`Integrate[Csc[x]^3, x, Method->"RischTranscendental"]`.

---

## Root cause

`D[G3, x] − Sec[x]^3` is a rational function of `{Sin[k x], Cos[k x], Log[2 ± 2 Sin[x]]}`
whose value is identically `0`. To *prove* that, `Simplify` must:

1. Differentiate the `Log[2 ± 2 Sin[x]]` terms (→ `Cos[x]/(1 ± Sin[x])` weights) and the
   multiple-angle rational envelope;
2. Collapse the whole **multiple-angle** combination (`Sin[3x]`, `Sin[5x]`, `Cos[6x]`, …)
   back to a common single-angle (`Sin[x]`/`Cos[x]`) power basis — i.e. a **`TrigReduce`
   power-reduction** across a large mixed numerator/denominator — and cancel.

Mathilda's `Simplify`/`TrigReduce` does not perform step 2 at this expression size: the
multiple-angle → power collapse either is not attempted on this shape or the intermediate
swell defeats the complexity-guided search (cf. `[[project_simplify_multigenerator_explosion]]`,
`[[project_trigrat_radical_generators]]`). So the diff-back never certifies.

---

## Fix directions (in `Simplify`, not the integrator)

1. **Multiple-angle power reduction in the diff-back path.** Teach `TrigReduce` (or a
   `Simplify` transform) to rewrite `Sin[k x]`/`Cos[k x]` combinations over a common angle
   into a `Sin[x]`/`Cos[x]` power basis (Chebyshev expansion), then cancel — the standard
   route for `D[trig-antiderivative] − trig-integrand === 0`.
2. **Half-angle / Weierstrass normal form for the zero test.** Substituting
   `s = Tan[x/2]` turns the whole diff-back into a *rational function of `s`* whose
   `Together`-numerator is exactly `0`; a `Simplify` strategy that recognizes a pure-trig
   rational identity and routes it through the half-angle rationalization would decide these
   in one step. (This must stay a **symbolic** decision — no numeric sampling.)
3. **Integrator-side (secondary):** emit the answer in the reduced single-angle power basis
   (apply `TrigReduce`/power-collapse to `G` *before* the gate), so the diff-back `Simplify`
   starts from a shape it can reduce. This narrows the gap without a general `Simplify` win.

Any of these restores `Sec^n`/`Csc^n` (n ≥ 3) to the integrator behind a genuine **exact**
certificate — the correct closure, replacing the removed numeric hack.

---

## Not in this bucket (genuinely non-elementary / different cause)

- `∫ Tan[Log[x]] dx` — **non-elementary** (`u = Log[x] ⇒ ∫ e^u Tan[u] du`, no elementary
  form). Declining is *correct*, not a `Simplify` gap.
- `∫ Tan[Log[x]]/x dx = −Log[Cos[Log[x]]]` — elementary and simple, but declines through the
  *nested-tangent* dispatch (a pre-existing masked failure in
  `integrate_risch_transcendental_tests` on clean HEAD, unrelated to the numeric gate).
  Tracked separately as a hypertangent-dispatch gap, not a `Simplify` deficiency.
