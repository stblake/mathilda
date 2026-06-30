# Plan: fast, fully-general `Together` / `Cancel` / `Simplify` over algebraic-number extensions

This is the implementation plan for the deficiencies catalogued in
[`SIMPLIFY_ALGEBRA_SLOW.md`](SIMPLIFY_ALGEBRA_SLOW.md). The goal is to make
`Together`, `Cancel`, `Simplify`, and `PossibleZeroQ` over an algebraic-number
field `Q(α)` (cyclotomic and radical constants) run at the same order of cost as
over `Q(x)` — **microseconds-to-milliseconds, not tens of seconds** — and to do
so in full generality (arbitrary cyclotomic order, arbitrary radical degree,
towers of several generators, and combinations).

---

## 0. Reframing: the substrate already exists

The investigation that preceded this plan established a crucial fact that
changes the strategy. Mathilda **already has** a complete, correct
arithmetic substrate for `Q(α)[x]`:

| Layer | File | What it provides |
|-------|------|------------------|
| Field element `Q(α)` | `src/poly/qa.{c,h}` — `QAExt`, `QANum` | `qa_add/sub/mul/inverse/div`, `reduce_mod_palpha`; element is a dense `mpq_t` vector mod the minimal polynomial `P_α` |
| Univariate poly over `Q(α)` | `src/poly/qaupoly.{c,h}` — `QAUPoly` | `qaupoly_divrem`, **`qaupoly_gcd`** (Euclidean, monic), `qaupoly_eval`, `qaupoly_shift` |
| Surface ↔ substrate bridge | `src/poly/qafactor.{c,h}` | `qa_resolve_extension` (α-form → `QAExt`), `qa_expr_to_qaupoly` (lift), `qaupoly_to_expr_alpha` (render), Trager norm/factor, tower/primitive-element machinery |
| Builtin wiring | `src/poly/poly.c` | `PolynomialGCD/LCM/Quotient/Remainder` all detect extensions via `extension_autodetect_args` and dispatch into the substrate |

Because `QANum` is a fixed `Q`-basis vector reduced modulo `P_α`, the identities
that `Simplify` currently has to *rediscover* after the fact — `α^3 = 1`,
`1 + α + α^2 = 0` for `α = Exp[2πI/3]` — are **structurally built in**:
`α^2` is just the stored vector `(-1, -1)`. Once an expression's constants live
in this representation, the super-polynomial blowup cannot occur, because every
field element is bounded to `deg(P_α)` rationals.

**Therefore the blowup is a *reachability* problem, not a *capability* problem.**
The expressions in `SIMPLIFY_ALGEBRA_SLOW.md` blow up because their constants
are **never lifted into the substrate** — they stay as opaque `Power[...]` /
`Exp[...]` / `Complex[...]` trees and get combined by the rational-over-`Q`
path (`cancel_recursive` / `together_recursive` in `src/rat.c`), which treats
each distinct algebraic atom as an independent transcendental symbol and forms
giant common denominators across all of them.

The two concrete reachability gaps:

1. **`extension_autodetect` (`src/poly/qafactor.c`, walker `autodetect_walk`
   ~`:3374`) does not recognize cyclotomic / root-of-unity generators.** It
   matches only `Sqrt[...]` and `Power[integer, p/q]`. So `Exp[2πI/3]`,
   `(-1)^(1/3)`, `Exp[I π/5]`, etc. are invisible to it and fall through to the
   `Q` path. *This is the single highest-leverage fix in this plan.*

2. **For root-of-unity bases, the "natural generator `yⁿ − c`" minimal
   polynomial used by `qa_resolve_extension` is reducible** (e.g. `(-1)^(1/3)`
   → `y^3 + 1 = (y+1)(y^2−y+1)`), which the substrate cannot use as a field
   modulus. The correct minimal polynomial is the **cyclotomic polynomial
   `Φ_n`** (degree `φ(n)`), which `MinimalPolynomial`'s helper
   `mp_root_of_unity_order` (`src/poly/minpoly.c:231`) already knows how to
   identify — but that machinery is not wired into the extension path.

Everything else in this plan (zero-test ordering, reduce-before-combine,
real-subfield output) is secondary optimization on top of closing gap 1+2.

---

## 1. Guiding principles

- **Full generality, no special cases.** No hard-coded `Exp[2πI/3]`. The
  recognizer computes the order `n` of any root of unity and the cyclotomic
  polynomial `Φ_n` of any degree; radical generators of any degree; and uses
  the *existing* primitive-element tower machinery when several distinct
  generators co-occur.
- **Reuse the substrate; do not reinvent it.** Land changes by *routing into*
  `qa.c`/`qaupoly.c`/`qafactor.c`, not by writing parallel arithmetic.
- **Correctness is non-negotiable.** Every phase ships with the
  `SIMPLIFY_ALGEBRA_SLOW.md` reproducers as regression tests plus a numeric
  cross-check (`PossibleZeroQ[old − new]` / sampled diff-back) so a fast wrong
  answer can never slip through.
- **Bound intermediate size.** Reduce each term into the `Q`-basis *before*
  combining (the `radrat` "reduce-before-rationalise" lesson), so no giant
  common denominator over conjugates is ever formed.
- **Measure against the documented numbers.** Definition of done is tied to the
  §6 reproducers in `SIMPLIFY_ALGEBRA_SLOW.md`.

---

## 2. Phase 1 — Recognize cyclotomic / root-of-unity generators (core unblock)

**Goal.** Make `extension_autodetect` (and `qa_resolve_extension`) recognize any
root of unity as an extension generator with its *cyclotomic* minimal
polynomial, so the existing `Q(α)` substrate handles cyclotomic-coefficient
rationals automatically.

### 2.1 Detection (`src/poly/qafactor.c`)

Add a root-of-unity classifier and call it from `autodetect_walk` and from
`qa_resolve_extension`. It must recognize, for rational `r = p/q`:

- `Exp[2 π I r]`            → `ζ = Exp[2πI/n]`, primitive `n`-th root, `n = q/gcd(...)`
- `Exp[I π r]`             → `Exp[2πI (r/2)]`, reduce to the above
- `Power[E, Complex[0, c]]` where `c` is a rational multiple of `π`
- `(-1)^r = Power[-1, p/q]` → `Exp[I π p/q]`, reduce to the above
- `I`                       → already handled (`n = 4`); keep as a fast path
- products `Times[radical, root-of-unity]` and `Times[I, Sqrt[c]]` (already
  partially handled) — route the root-of-unity factor through the new path

Reuse the order/normalization logic already proven in
`src/poly/minpoly.c`:
- `mp_root_of_unity_order` (`minpoly.c:231`) — extract the order `n` from any of
  the above surface forms.
- `mp_make_radical` / `mp_imag` (`minpoly.c:203` / `:185`) — surface-form
  helpers.

Factor these out of `minpoly.c` into a shared header
(`src/poly/rootofunity.{c,h}`, or extend `numbertheory_internal.h`) so both
`MinimalPolynomial` and the extension recognizer call one implementation. **Do
not duplicate** the order logic.

### 2.2 Cyclotomic minimal polynomial (`src/poly/`)

The substrate needs `Φ_n(x) ∈ Z[x]`, degree `φ(n)`, as the `QAExt` modulus.
Add `qaext_cyclotomic(unsigned long n)` (sibling to `qaext_sqrt_si` /
`qaext_root_si` in `src/poly/qa.c:71`). Compute `Φ_n` by the standard divisor
recursion:

```
x^n − 1 = ∏_{d | n} Φ_d(x)      ⟹      Φ_n = (x^n − 1) / ∏_{d | n, d < n} Φ_d
```

- Memoize computed `Φ_d` in a small static cache (orders are tiny in practice).
- All arithmetic is exact polynomial division over `Z` (use the existing
  integer/`mpq` polynomial divrem; the quotient is exact by construction).
- `φ(n)` and the divisors of `n` come from the existing number-theory module
  (`numbertheory/` — `EulerPhi`, `Divisors`); reuse, don't reimplement.

Wire it into `qa_resolve_extension` (`qafactor.c:601`): when the surface α is a
root of unity of order `n`, build `QAExt` from `Φ_n` (not from `yⁿ − c`), and
record the surface render (`Exp[2πI/n]`) for `qaupoly_to_expr_alpha`. This
single change is what makes `(-1)^(1/3)` resolve to a **degree-2** field instead
of a reducible cubic.

### 2.3 Generality: towers and mixed generators

When an expression contains several distinct generators (e.g. a cyclotomic
constant *and* a `Sqrt`, or two different roots of unity), the existing
primitive-element tower path already composes them:
`qa_resolve_extension_tower` / `qa_tower_extend` / `find_primitive_shift`
(`qafactor.c:1547` / `:1425` / `:1334`). Once 2.1–2.2 make each *individual*
cyclotomic generator resolvable, the tower machinery composes them with no
further work. Add tests for the mixed case; expect the recognizer to feed the
tower builder, not a special path.

### 2.4 Files touched

- `src/poly/qafactor.c` — `autodetect_walk`, `qa_resolve_extension`, recognizer.
- `src/poly/qa.c` / `qa.h` — `qaext_cyclotomic`.
- `src/poly/rootofunity.{c,h}` (new) or `minpoly.c` refactor — shared order logic.
- `tests/CMakeLists.txt` — add any new `src/poly/*.c` to `COMMON_SRC`
  (see memory note *Tests COMMON_SRC list*).

### 2.5 Tests (`tests/test_qafactor*.c` or new `test_cyclotomic_extension.c`)

- `qaext_cyclotomic`: `Φ_1=x−1`, `Φ_2=x+1`, `Φ_3=x²+x+1`, `Φ_4=x²+1`,
  `Φ_6=x²−x+1`, `Φ_12=x⁴−x²+1`, `Φ_5`, `Φ_8`, a prime `Φ_p` and a
  prime-power `Φ_{p^2}`.
- Recognizer: `Exp[2πI/3] → (n=3)`, `(-1)^(1/3) → Exp[Iπ/3] → (n=6)`,
  `Exp[Iπ/5] → (n=10)`, `I → (n=4)`.
- End-to-end the §6 *primitive-level* reproducers from
  `SIMPLIFY_ALGEBRA_SLOW.md`:
  ```
  al = Exp[2 Pi I/3];
  Together[F + (F/.t->Sa) + (F/.t->Sb) + (F/.t->Sc), Extension->Automatic]
  ```
  must finish in **< 100 ms** and equal the `Q`-path answer where one exists.
- Numeric cross-check `PossibleZeroQ[old − new]` for every algebraic case.

**Expected impact.** This alone is projected to clear the §1 / §6 blowups,
because the constants now live in a bounded `Q`-basis. Verify before proceeding.

---

## 3. Phase 2 — Zero-test ordering guard (`src/zero_test.c`)

**Goal.** Stop `PossibleZeroQ` from running the `Together`∘`Cancel` Stage 1 on
algebraic-number expressions that still have free symbols — route them straight
to numeric Schwartz–Zippel sampling, which decides true identities without any
symbolic combination.

### 3.1 Change

`zero_test_decide` (`src/zero_test.c:833`) currently calls `decide_rational`
(`:279`) **unconditionally** at `:839`, before the numeric/sampling stages.
`decide_rational` runs `internal_together` ∘ `internal_cancel`
(`:290`–`:292`) — the exact blowup site.

Insert a guard immediately before `:839`:

```
if (has_free_symbols(e) && expr_has_algebraic_constant(e))
    return decide_schwartz_zippel(e, ...);   /* skip Stage 1 entirely */
```

- `has_free_symbols` already exists (`src/zero_test.c:202`).
- Add `expr_has_algebraic_constant(e)`: true if `e` contains `Power[_, Rational]`
  with non-integer exponent, `(-1)^rational`, `Exp[2πI·rational]`, or `I` in a
  position that makes the expression algebraic-over-`Q` with free symbols. Model
  it on the `alg_collect_*` walkers in `src/simp/simp_algebraic.c:93` and the
  Phase 1 recognizer (share the predicate).

### 3.2 Why this is correctness-safe

Stage 1's only trustworthy verdict is **TRUE** (`:840` comment: "never trust
False from Stage 1 alone"). Stage 3 (`decide_schwartz_zippel`, `:788`)
independently and correctly decides true identities by sampling free symbols at
deterministic real rationals (seeded from `expr_hash`, `:800`). Skipping Stage 1
for these inputs therefore loses no decision power — it only avoids the
non-terminating symbolic combination. The numeric ladder (`decide_numeric`,
`:617`) still covers the no-free-symbol case.

### 3.3 Payoff

This generalizes the hand-rolled `sample_clearly_nonzero`
(`src/calculus/integrate_goursat.c:277`) into the language primitive. After it
lands, the Goursat numeric-gate guards become redundant and can be removed in
Phase 5.

### 3.4 Files / tests

- `src/zero_test.c` — guard + `expr_has_algebraic_constant`.
- `tests/test_zero_test*.c` — `PossibleZeroQ[f0+f1+f2+f3]` over `Q(α)` returns in
  **< 100 ms**; plus true-identity and true-nonzero cases over `Q(α)` with free
  symbols; plus the existing rational cases unchanged (no regression).

---

## 4. Phase 3 — Reduce-before-combine and drop the `evaluate` round-trips (`src/rat.c`)

**Goal.** Tighten the extension path in `Together`/`Cancel` so it (a) reduces
each summand into the `Q(α)` basis before combining and (b) calls the C
substrate directly instead of reconstructing option-bearing Exprs and
re-`evaluate`-ing.

### 4.1 Current shape

`cancel_with_extension` (`src/rat.c:812`) and `together_recursive_ext`
(`:651`) implement extension arithmetic by **building Expr trees**
(`PolynomialGCD[…, Extension->α]`, `PolynomialLCM`, `PolynomialQuotient`) and
calling `evaluate` (`:847`, `:673`, `:719`, `:910`). Each round-trip re-parses
options, re-detects the extension, re-lifts to `QAUPoly`. Correct, but heavy,
and it combines before reducing.

### 4.2 Changes

1. **Reduce-before-combine.** In `together_recursive_ext`, lift each summand to
   `QAUPoly` over the (already-resolved) `QAExt` and reduce mod `P_α` *first*,
   then combine numerators over the common `QAUPoly` denominator. Never form a
   common denominator over un-reduced conjugate constants. This is the bound
   that prevents the giant-intermediate explosion (`radrat` lesson).
2. **Call C directly.** Replace the `evaluate(PolynomialGCD[…])` round-trips
   with direct `qa_expr_to_qaupoly` → `qaupoly_gcd` / `qaupoly_divrem` →
   `qaupoly_to_expr_alpha`, reusing the resolved `QAExt` from the single
   `extension_autodetect` call already done at the top of `builtin_cancel_compute`
   / `builtin_together_compute` (`rat.c:1025`). This removes O(depth) redundant
   re-detection/re-lift work.

### 4.3 Risk & sequencing

This is a refactor of working code; it is **lower priority than Phases 1–2** and
should only land *after* Phase 1 proves the substrate handles the cases, with
the cross-checks green. It is a constant-factor / intermediate-size improvement,
not the asymptotic fix. Keep the existing `evaluate`-based path as a fallback
behind the same NULL-means-"can't" contract so any unhandled shape degrades
gracefully rather than regressing.

### 4.4 Files / tests

- `src/rat.c` — `together_recursive_ext`, `cancel_with_extension`.
- Also audit the parallel threading in `src/parfrac.c:35` (Apart) and
  `src/calculus/intrat.c:3531` for the same round-trip pattern; apply the same
  treatment if cheap.
- Tests: timing assertions on multi-term `Q(α)` `Together`/`Cancel`; valgrind
  clean (watch the NULL-out-before-free contract when reusing input subtrees).

---

## 5. Phase 4 — Real-subfield / conjugate-pair output (generality + beauty)

**Goal.** When a cyclotomic computation has a real result (conjugate `α`-term +
`α²`-term ∈ `R`), present it over the real subfield instead of carrying complex
roots of unity — e.g. work in `Q(α)`'s real subfield via the real quadratic
`(t−α)(t−α²) = t²+t+1`. This is the natural route for the Legendre–Clausen
family and yields `ArcTan`/`Log` real closed forms instead of complex ones.

### 5.1 Status

With Phase 1 in place the cyclotomic arithmetic is already **fast and correct**;
Phase 4 is about **output form**, not performance. It is the most research-heavy
phase (detecting that conjugate pairs sum to a real element, and re-expressing
in the maximal real subfield `Q(α)^+`). Treat it as optional / follow-on:

- Detect conjugate-closed inputs (the constant set is stable under
  `α → ᾱ`).
- Compute the real subfield generator (`α + ᾱ`, e.g. `2cos(2π/n)`) and its
  minimal polynomial; reduce the result into `Q(α + ᾱ)`.
- Render via the existing complex→real recombination already used elsewhere
  (`intsimp_finalize`'s `ArcTan` tidy; `radrat`'s conjugate handling — see
  memory notes).

### 5.2 Defer condition

Do **not** start Phase 4 until Phases 1–3 are green and measured. It is the only
phase that may not be strictly necessary for the blocked features if the
complex-form answers from Phase 1 are acceptable (they are mathematically
correct, just less pretty).

---

## 6. Phase 5 — Reconnect the integrators and remove the workarounds

Once Phases 1–3 land:

1. **Remove the numeric-gate workarounds** in
   `src/calculus/integrate_goursat.c` (`sample_clearly_nonzero` and the
   `eigenproj*_raw` numeric gates at `:615`, `:862`, `:1004`, `:1137`,
   `:1349`) — Phase 2 makes `PossibleZeroQ` itself fast on these, so the
   *criterion* no longer needs a hand-rolled sampler. Keep one sampler only if a
   measured case still needs it; otherwise delete.
2. **Unblock Goursat Section 5 / A4.** With cyclotomic `canonic` fast, both
   blockers in `SIMPLIFY_ALGEBRA_SLOW.md` §3 clear:
   - the fix-α / fix-α² basis terms in `period3_reduce` now combine in `Q(ζ_3)`;
   - `goursat_v4`'s descent no longer exhausts the time budget.
   Verify `Integrate[t/((t^3+8) Sqrt[t^3-1]), t]` produces a closed form in
   well under 1 s (the user's hard efficiency requirement), and that the
   single-character basis term
   `Integrate[(t-al)/((t+2 al) Sqrt[t^3-1]), t, Method->"GoursatAlgebraic"]`
   (>60 s today) is now fast.
3. **Drop / raise the `TimeConstrained` budget** in `builtin_gs_run` once the
   underlying operations are fast, so correct-but-slightly-slow cases are not
   spuriously declined.

---

## 7. Sequencing & dependencies

```
Phase 1 (recognize cyclotomic generators)   ──┐  core unblock; do first
Phase 2 (zero-test ordering guard)           ──┤  independent of P1, cheap, do in parallel
Phase 3 (reduce-before-combine in rat.c)     ──┘  after P1 proves substrate
Phase 4 (real-subfield output)               ── optional, after 1–3, output quality only
Phase 5 (reconnect integrators, delete guards) ── after 1–3 green
```

- Phases 1 and 2 are independent and can be developed concurrently (different
  files); 1 is the asymptotic fix, 2 is a cheap orthogonal win that also helps
  every other algebraic-number `PossibleZeroQ` caller.
- Phase 3 depends on Phase 1 (needs a resolved `QAExt` to reduce into).
- Phases 4 and 5 depend on 1–3.

---

## 8. Definition of done (measured against `SIMPLIFY_ALGEBRA_SLOW.md` §6)

A clean `make` build of `Mathilda`, each run **foreground** (never background a
Mathilda process — memory note *No background Mathilda pollers*):

| Reproducer | Today | Target |
|------------|-------|--------|
| `Together[F + F/.t->Sa + F/.t->Sb + F/.t->Sc, Extension->Automatic]` (cyclotomic) | does not finish in seconds | **< 100 ms**, correct |
| `Simplify[F + F/.t->Sa + …]` (cyclotomic) | slow | **< 200 ms**, correct |
| `PossibleZeroQ[f0+f1+f2+f3]` over `Q(α)` | ~20 s+ | **< 100 ms** |
| `Integrate[(t-al)/((t+2 al) Sqrt[t^3-1]), t, Method->"GoursatAlgebraic"]` | > 60 s, killed | **< 1 s**, closed form |
| `Integrate[t/((t^3+8) Sqrt[t^3-1]), t]` (Section 5 / A4) | ~12 s, declines | **< 1 s**, closed form |
| Rational-only `Together`/`Cancel`/`Simplify`/`Factor` | µs | **unchanged** (no regression) |

Plus, for every algebraic case: `PossibleZeroQ[old_answer − new_answer]` is
`True` (or a sampled diff-back is numerically zero), so no fast-but-wrong answer
can ship. All affected `*_tests` binaries pass; `valgrind` clean on the new
paths (diff against a baseline per memory note *macOS valgrind baseline noise*).

---

## 9. Risks and mitigations

- **Wrong minimal polynomial for a root of unity.** Using `yⁿ − c` (reducible)
  instead of `Φ_n` corrupts the field. Mitigation: Phase 1 unit tests pin every
  `Φ_n` against known values, and the recognizer routes *all* root-of-unity
  bases through `qaext_cyclotomic`, never `qaext_root_si`.
- **Mis-recognizing a non-root-of-unity as one** (e.g. `Exp[2 x I]` with free
  `x`). Mitigation: the recognizer fires only when the exponent is a *rational
  multiple of `π I`* with no free symbols; otherwise the constant stays opaque
  and the existing path runs.
- **Tower primitive-element cost.** Composing many distinct generators can be
  expensive even with bounded fields. Mitigation: reuse the existing
  `find_primitive_shift` heuristics and cap generator count exactly as the
  multi-generator `Simplify` guard already does (`simp_search.c:161`).
- **Regression on the fast rational path.** Every change is gated behind "has an
  algebraic constant"; pure-`Q` inputs must take the identical code path they do
  today. Enforced by the "unchanged" row in §8.
- **`src/external/ecm` must not be touched** (it shows dirty in `git status`;
  leave it). All work is under `src/poly/`, `src/rat.c`, `src/zero_test.c`,
  `src/simp/`, `src/calculus/`.

---

## 10. Summary

The headline of `SIMPLIFY_ALGEBRA_SLOW.md` — "the deficiency is arithmetic over
`Q(α)`" — has a sharper diagnosis after mapping the code: **the arithmetic
substrate is already there and already fast; cyclotomic constants simply never
reach it.** Phase 1 (recognize root-of-unity generators and give them their
cyclotomic minimal polynomial) is therefore expected to clear the bulk of the
blowups by itself, with Phase 2 (zero-test ordering) a cheap orthogonal win and
Phase 3 (reduce-before-combine) a constant-factor tightening. Phases 4–5 turn
the now-fast primitives into prettier real output and the unblocked Goursat
Section 5 integrator.
