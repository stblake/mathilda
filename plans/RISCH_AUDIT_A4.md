# Risch Transcendental — Independent Further Audit (A4)

**Date:** 2026-07-15
**Reference:** Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed. (2004).
**Subject:** `src/calculus/integrate_risch_transcendental.c` (5237 LoC),
`src/calculus/integrate_risch_rde.c` (1002 LoC, extracted since the A1–A3 audits),
the `risch_{field,canonical,structure,hermite,hypertangent,coupled}` modules, and the
`Risch\`` decision/DE builtins.
**Method:** independent code review of the paths that changed since A1–A3 (recursive
RDE stack, RT_TAN tower monomial, the P3 decision procedure, the residue split + FTC
rule) **plus empirical black-box probing of the built binary** — diff-back of every
emission at multiple real points, and adversarial construction of wrong-decision /
wrong-answer inputs. All findings below are reproduced on the prebuilt `./Mathilda`.

> **Headline:** `RISCH_AUDIT_FINDINGS.md` states *"No CRITICAL findings and no soundness
> defects anywhere … a wrong antiderivative cannot be produced."* This audit **falsifies
> that claim.** Two genuine soundness defects exist — one **wrong answer** (public
> builtins) and one **wrong decision** (`ElementaryIntegralQ` → `False` and an
> `Integrate::nonelem` message on integrands that *do* have elementary antiderivatives).
> Both post-date A1–A3 (the RDE extraction and the RT_TAN/P3 work), which is why the
> earlier audits missed them. The **emitted antiderivatives of the elementary integrator
> remain diff-back-correct** across a wide battery — the "never ship a wrong closed form"
> property holds for the *tower success paths*; the defects are in the **base-field DE
> solver** and the **decision/verdict layer** built on top of it.

---

## Severity-ranked findings

| # | Type | Sev | Site | One-line repro | Status |
|---|------|-----|------|----------------|--------|
| 1 | WRONG-ANSWER | **High** | `integrate_risch_rde.c` b=0 gcd/quot path; unguarded `rde_base`; root cause `PolynomialGCD[·,0]→1` | `Risch\`RischDE[-2/x, x, x]` → `0` (should decline) | **FIXED 2026-07-15** |
| 2 | WRONG-DECISION (false "non-elementary") | **High** | `integrate_risch_transcendental.c:3882` | `Risch\`ElementaryIntegralQ[Sec[x]^2 E^Tan[x], x]` → `False` | **FIXED 2026-07-15** |
| 3 | Non-termination | Med | RT_TAN / trig field | `Integrate[1/(a+b Sin[x]), x, Method->"RischTranscendental"]` runs >90 s | open |
| 4 | ~~Faithfulness regression~~ | — | residue-split partial emission | `Integrate[…]` returns `logs + Integrate[remainder, x]` | **not a defect** — intended Mathilda behavior (per maintainer) |
| 5 | Latent diff-back hole | Med (latent) | FTC rule `deriv.c:820` + `rt_verify_antideriv` | diff-back certificate is structurally foolable by any result containing `Integrate[…]` | open |
| 6 | Hackish numeric gate | Med | `rt_realify_numverify` (4 fixed positive points) | numeric sampling used as a soundness certificate in a decision procedure | **FIXED 2026-07-15** |
| 7 | Under-decision | Low | tangent `Dc≠0` certificate discarded | `ElementaryIntegralQ[x Tan[x], x]` → unevaluated (safe) | open |
| 8 | Cosmetic (pre-existing) | Low | `Risch\`PolyDivide` (`risch_field.c`) | `Risch\`PolyDivide[t^3+t+1, t^2+1, t]` → remainder `1+t+t^3-t(1+t^2)` unexpanded (= 1); soft-FAIL in `risch_field_tests` on clean HEAD too | open |

---

## Finding 1 — WRONG-ANSWER: base RDE/SPDE fabricates spurious solutions

**The base-field Risch-DE solver returns a `q` that is not a solution, and the base path
has no correct-by-construction gate to catch it.** Reproduced through public builtins:

```
Risch`RischDE[-2/x, x, x]   ==>  0        D[0] - (2/x)·0 - x = -x  ≠ 0   (WRONG)
Risch`RischDE[-1/x, 1, x]   ==>  0        residual -1 ≠ 0                (WRONG)
Risch`SPDE[x, 0, 1, x, 0]   ==>  {0,0,0,0,0}   claims x·q' = 1 solved by q=0 (WRONG)
Risch`RischDE[1, x, x]      ==>  -1 + x   (correct — the solver works normally)
```

`y' − (2/x) y = x` **is** solvable, but only by the *non-rational* `y = x² Log[x]`
(independently verified: `Simplify[D[x^2 Log[x],x] − (2/x)(x^2 Log[x]) − x] === 0`). The
correct return is therefore "no rational solution" (decline / unevaluated); returning `0`
is a plain wrong answer.

**Root cause — `PolynomialGCD` with a zero operand:**

```
PolynomialGCD[x, 0]   ==>  1     (mathematically gcd(a,0)=a, i.e. x)
PolynomialGCD[x^2, 0] ==>  1
```

`rde_spde` / `rde_spde_field` compute `g = gcd(a, b)` and decline iff `g ∤ c`
(Bronstein Thm 6.4.1). In the `b = 0` sub-problem the correct `g = gcd(a,0) = a`, and
`a ∤ c` would correctly decline. But `rde_gcd` = `PolynomialGCD` returns `1`, so the
divisibility test passes spuriously; the subsequent "exact" quotient
`z = (c − b·r)/a` via `PolynomialQuotient` **silently truncates** the nonzero remainder
(`PolynomialQuotient[1, x, x] = 0`), manufacturing a `c = 0` subproblem that the `n<0 ∧
c=0` terminal accepts, reconstructing `q = 0`.

**Why it escapes:** `rde_base` / `rt_solve_rde` return the reconstructed `q` with **no
`D[q] + f q == g` verification** — unlike `rde_tower`, which gates every result with the
exact identity check at `integrate_risch_rde.c:984-1000`. `builtin_rischtranscendental`
applies no diff-back either, trusting "correct by construction." The base RDE path is thus
the one place a reconstruction slip escapes as a wrong closed form.

**Blast radius:** directly wrong via the public `Risch\`RischDE` / `Risch\`SPDE` builtins.
A wrong closed form was *not* surfaced through top-level `Integrate` in the audit window
(flat exp integrands like `E^x/x` produce coefficient `f = +1`, not a pole `−n/x`, and
correctly yield `ExpIntegralEi`), but the defect is proven at the `rde_base`/`rde_spde`
layer that the whole *"authoritative NULL"* claim rests on. Second-order consequence
(hypothesized, not demonstrated): a fabricated success returned instead of NULL can
**mask a non-elementary term as elementary**, i.e. a `False`→`True` flip of
`ElementaryIntegralQ` in the opposite direction.

**Fix direction:** either fix `PolynomialGCD[a,0]→a` (a general CAS correctness bug), or
make `rde_gcd`/`rde_spde` defensive about a zero operand (Bronstein assumes `gcd(a,0)=a`),
and add the `D[q]+fq==g` self-verify to the base path so a slip declines instead of
shipping.

**✅ FIXED (2026-07-15).** Two changes in `integrate_risch_rde.c`: (a) `rde_gcd` now
returns the nonzero operand when one operand is zero (root cause — the field `rc_gcd`
already handled this via its Euclidean loop and only delegates to `rde_gcd` in the base
case, so the single guard corrects every caller); (b) a `D[y]+f y−g == 0` self-verify was
added to `rde_base`'s Expr fallback (the one path lacking the gate `rde_tower` already
has), so any residual reconstruction slip declines instead of shipping. Verified:
`Risch\`RischDE[-2/x,x,x]` and `Risch\`SPDE[x,0,1,x,0]` now decline; the genuine
`Risch\`RischDE[1,x,x] → -1+x` and all `risch_rde_tower_tests` / `risch_field_tests` are
unaffected; valgrind unchanged (13,440 B / 420 blocks, the module baseline).

---

## Finding 2 — WRONG-DECISION: `ElementaryIntegralQ` calls elementary integrands "non-elementary"

`ElementaryIntegralQ[f,x] → False` is documented as *"the recursive Risch decision
procedure PROVES no elementary antiderivative exists."* It is wrong for a whole,
easily-generated family:

```
Risch`ElementaryIntegralQ[Sec[x]^2 E^Tan[x], x]     ==>  False
Integrate[Sec[x]^2 E^Tan[x], x]                     ==>  E^Tan[x]      (Mathilda's own)
Simplify[D[E^Tan[x],x] - Sec[x]^2 E^Tan[x]]         ==>  0            (⇒ elementary!)

Risch`ElementaryIntegralQ[Sech[x]^2 E^Tanh[x], x]   ==>  False   (antideriv E^Tanh[x])
Risch`ElementaryIntegralQ[Tan[x] Sec[x]^2 E^Tan[x], x] ==> False (antideriv Tan[x]E^Tan[x]-E^Tan[x])
Risch`ElementaryIntegralQ[Sec[x]^2 E^(2 Tan[x]), x] ==>  False   (antideriv E^(2Tan[x])/2)
```

It is also **user-facing** via a false message:

```
Integrate[Sec[x]^2 E^Tan[x], x, Method->"RischTranscendental"]
  Integrate::nonelem: The integrand Sec[x]^2 E^Tan[x] has no antiderivative elementary in x.
```

**Root cause — `rt_field_rde`, `integrate_risch_transcendental.c:3866–3882`:**

```c
if (L >= 1 && (T->kind[L-1] == RT_LOG || T->kind[L-1] == RT_EXP)) {
    ... rde_tower(...) ...          /* recursive solve only for LOG/EXP lower tops */
}
rt_dec_nonelem();                   /* line 3882 — ALSO fires for the RT_TAN case */
return NULL;
```

For `Sec[x]^2 E^Tan[x]` the tower is `{Tan[x] (RT_TAN), E^Tan[x] (RT_EXP)}`; integrating
the `E^Tan[x]` level needs the field RDE over `k = C(x, Tan[x])` whose coefficient
`w' = 1+Tan[x]²` is *not* free of the tower variable, so the base branch is skipped and,
because the lower top is `RT_TAN`, the `rde_tower` dispatch guard is *also* skipped —
control falls straight to the **unconditional** `rt_dec_nonelem()`. The surrounding
comment (lines 3875–3881) itself admits the tangent case is *"Gap 3 … out of current
scope"*, yet the nonelem flag is raised *before* that acknowledgement. **An out-of-scope
decline must reach `RT_DEC_UNKNOWN` (→ unevaluated), never `RT_DEC_NONELEMENTARY`
(→ `False`).**

The bug is **not uniform** — `Csc[x]^2 E^Cot[x]` (a `Cot` spelling of the same integrand)
returns *unevaluated*, the safe outcome — which is itself evidence this is an accidental
scope leak rather than a theorem. Related non-authoritative NULLs that flow into the same
line-3882 verdict (harder to trigger, same class): `rde_tower` returns NULL for scope, not
theorem, reasons at `integrate_risch_rde.c:923, 937, 967` (its own header correctly says
NULL *"= no solution / out of the current increment's scope"*, but `rt_field_rde` treats
every such NULL as authoritative).

**Fix direction:** line 3882 must return NULL leaving `g_rt_decision` UNKNOWN when the
lower top is `RT_TAN` (or the `b = Dz/z` limited-integration branch); only LOG/EXP
lower-tops that actually *ran* `rde_tower` to a no-solution conclusion are authoritative.

**✅ FIXED (2026-07-15).** In `rt_field_rde` the `rt_dec_nonelem()` was moved *inside* the
`RT_LOG || RT_EXP` branch — it now fires only after `rde_tower` was actually run and
returned NULL (an authoritative "no rational solution in K_L"). The `RT_TAN` /
out-of-scope fall-through returns NULL **without** raising the certificate, so
`ElementaryIntegralQ` stays `undec` rather than a false `False`. Verified:
`Sec[x]^2 E^Tan[x]`, `Sec[x]^2 E^(2 Tan[x])`, `Sech[x]^2 E^Tanh[x]` now return
unevaluated and no longer print `Integrate::nonelem`; every legitimate `False`
(`E^(x^2)`, `E^x/x`, `E^(E^x)/(1+E^(E^x))`, `1/Log[x]`) and `True` (`x E^x`) is preserved;
`risch_elementaryq_tests` green. *Residual (unchanged, not reached by any constructed
elementary integrand):* the `rc_ispoly` gate at `rde_tower` (`integrate_risch_rde.c:937`)
and the cancellation branch can still return a scope-NULL that the RT_LOG/RT_EXP branch
treats as authoritative; making `rde_tower` signal scope-vs-theorem distinctly is the
principled follow-up, but no elementary integrand was found to reach it.

---

## Finding 3 — Non-termination on tangent/trig fields

Two elementary or symbolic-parameter integrands run the RischTranscendental engine
CPU-bound without returning:

```
Integrate[1/(a+b Sin[x]), x, Method->"RischTranscendental"]      >90 s  (killed)
Risch`ElementaryIntegralQ[Log[Tan[x]] Sec[x]^2, x]               >30 s  (killed)
```

The second is elementary (`Integrate[Log[Tan[x]] Sec[x]^2, x] = Tan[x] Log[Tan[x]] −
Tan[x]` returns instantly), so the decision procedure hangs where the integrator does not.
Both involve the tangent/hypertangent field. Under the pipe protocol a hang also
manifests as a *silent no-response* (clean exit, no result line) when a subsequent
expression's `quit` reaps the process — worth a guard so the engine always terminates
(or bounds) rather than spinning.

---

## Finding 4 — Faithfulness regression: `Integrate` returns a partial non-answer

```
Integrate[1/(x Log[x]) + 1/(x (Log[x]^2 - x)), x]
  ==>  Log[ … Log[x] … ] + Integrate[1/(-x^2 + x Log[x]^2), x]
```

The residue-split feature (Thm 5.6.1 partial log part) integrates the constant-residue
part and hands the non-constant-residue part back as an **inert inner `Integrate[…]`**,
surfaced as a *final* answer. Mathematica returns this **fully unevaluated**. The emitted
logs are mathematically correct, so this is faithfulness, not correctness — but a user
cannot easily distinguish an "answer" containing `Integrate[…]` from a failure, and the
lone non-elementary term (without the elementary sibling) correctly stays unevaluated, so
the behavior is inconsistent.

---

## Finding 5 — Latent diff-back hole from the FTC rule

The rule `D[Integrate[f,x],x] → f` (`deriv.c:820`) is applied **unconditionally** for any
2-arg `Integrate` whose variable matches, and `rt_verify_antideriv`
(`integrate_risch_transcendental.c:~1960`) does `Simplify[D[result] − f] === 0` with **no
guard rejecting an unevaluated `Integrate` head in `result`**. Demonstrated:

```
D[Integrate[NonElem[x], x], x]                                    ==>  NonElem[x]
Simplify[D[99 x^2 + Integrate[Sin[x] - D[99 x^2, x], x], x] - Sin[x]]  ==>  0
```

So the universal diff-back certificate returns `True` for arbitrary garbage plus a
self-referential remainder. **Not currently exploited** — both live partial-emission paths
build the remainder as an *independent* partial-fraction split (never `f − D[logs]`), and
the coupled tower path verifies an *exact* `D_tower[Q] == F` tower-coordinate identity
before back-substitution — but the safety rests on an **unenforced invariant** with no
assertion. Any future path (or refactor) leaving an unevaluated `Integrate` in a candidate
result would be rubber-stamped.

**Fix direction:** have `rt_verify_antideriv` reject any result containing an unevaluated
`Integrate[…,x]`, and gate partial results through the exact-tower identity exclusively.

---

## Finding 6 — hackish numeric gate in a decision procedure

`rt_realify_numverify` checked `Σ|N[(D[g,x]−f)/.x→pₖ]| < 1e-9` at exactly **four fixed
positive rational points** `{0.7, 1.3, 2.3, 3.7}`. It was the *sole* soundness guard for
two paths where "correct by construction" no longer holds: the transcendental-argument
tangent extension (`rt_hypertan_family`) and the real reconstruction of the I-laden exp
route (`rt_realify` / `cx_reim`, `rt_exp_ratreduce_case`). A numeric sample has no place
certifying a Risch decision procedure — a finite positive-only derivative sample can pass
a form whose `D[g]−f` merely vanishes at those points.

**✅ FIXED (2026-07-15).** `rt_realify_numverify` was **deleted**; all four call sites now
use the exact symbolic gate `rt_verify_antideriv` (`Simplify[D[g,x]−f] === 0`). Strictly
more conservative — only a *proven-zero* form ships, so no wrong result can pass — and no
sampling anywhere in the integrator. Coverage cost: forms whose diff-back `Simplify` cannot
reduce now **decline** (the honest outcome). Most rational-trig is unaffected (`Tan`, `Sec`,
`Csc`, `1/(2+Cos)`, `1/(5+4Cos)`, `Sec²`, `1/(3+Tan²)` still ship real forms); the
higher-power `Sec³`/`Csc³`/`Sec⁴` families now decline — a **`Simplify` deficiency**, not an
integrator bug, written up in `SIMPLIFY_GAPS.md`. Verified: zero new test failures
(`integrate_risch_transcendental_tests` FAIL set byte-identical to clean HEAD), valgrind
unchanged (13,440 B / 420 blocks).

---

## Finding 7 — Tangent `Dc≠0` certificate discarded (under-decision, safe)

The genuine `Dc≠0 ⇒ non-elementary` certificate (`risch_hypertangent.c:267-270`) and the
§5.6 `beta=False` residue are computed but thrown away at the tower dispatch boundary
`rt_int_hypertangent_field` (line ~3530 requires `beta=True`), so the tangent field never
surfaces a `False`: `ElementaryIntegralQ[x Tan[x], x]` and `[Tan[x]/x, x]` return
unevaluated though both are genuinely non-elementary. This is the *opposite* of Finding 2
— an under-decision, always safe — but combined with Finding 2 it shows the tangent
decision layer is simply unfinished and inconsistent, not merely conservative.

---

## What is sound (independently confirmed)

- **Emitted antiderivatives diff-back-correct.** A ~30-integrand direct-engine battery
  (`Method->"RischTranscendental"`, complex-aware `Abs` diff-back at 4 points): tan/cot/
  tanh/sec/csc, `E^x/(1+E^(2x))→ArcTan[E^x]`, nested logs, `E^x E^(E^x)/(1+E^(2E^x))→
  ArcTan[E^(E^x)]`, `Sec^3`, Hermite repeated poles, etc. — every emitted form verifies.
  No wrong *closed form* from a tower success path was found.
- **The `rde_tower` exact-identity gate** (`integrate_risch_rde.c:984-1000`) holds; the
  recursive tower path emits only behind a `Together`-zero tower-variable identity.
- **Ch. 6 base transcription** (`rde_spde`, `RdeNormalDenominator` `en|dn·h²` guard,
  `rde_polyrischde_nocancel1/2`, degree bounds, weak normalizer) is faithful to the book
  boxes *except* for the `gcd(·,0)` slip of Finding 1.
- **Genuine non-elementary verdicts are correct:** `E^(x^2)`, `E^x/x`, `1/Log[x]`,
  `E^(E^x)/(1+E^(E^x))`, `E^(x^2)/(x-1)` → `False`; scope cases (`Sqrt[x]`,
  `1/Sqrt[1-x^2]`, multivariate) → unevaluated, never guessed.

---

## Disposition

The subsystem's *elementary integrator* is as strong as advertised — its emissions are
verified and its tower success paths are gated by exact identities. The two soundness
defects are localized:

- **Finding 1** is a base-DE reconstruction escaping because the base path uniquely lacks
  the self-verify that the tower path has — plus a general `PolynomialGCD[·,0]` bug.
- **Finding 2** is a single mis-placed `rt_dec_nonelem()` that turns a `Gap 3` scope
  decline into a false theorem.

Both are small, well-localized fixes with outsized correctness impact, and both **directly
contradict the `RISCH_AUDIT_FINDINGS.md` "no soundness defects" headline** — which should
be revised. Findings 3–5 are robustness/faithfulness items; 6–7 are hardening. None
requires touching `src/external/`.
