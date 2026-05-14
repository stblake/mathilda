# Together / Cancel — Layer 4 design: cut internal recursion cost

**Status:** Design only — not implemented. Document scopes the next phase of
work on `Together` / `Cancel` performance after the Layer 0+1+2+3+5
shipped in commits `bb5b6ba`, `8d72345`, `508dad2`, `e41cd67`.

**Goal:** reduce the headline `Simplify[D[Integrate[a x/(x^3+2), x], x]]`
wall-clock from ~134 ms to ≤ 50 ms, and proportionally cut every other
multi-α algebraic-rational Simplify whose tower path is genuinely
needed.

---

## 1. Problem statement

After Layers 0+1+2+3+5 (input-level nested-radical prefilter,
pre-substitution gate, guided fallback, squarefree-prime coalesce,
FactorMemo memoisation), the only remaining slow tail in
`Together[…, Extension -> Automatic]` is on inputs that **legitimately
need** the multi-α tower path: a Plus of fractions over Q(γ)[x]
where qaupoly_gcd over the tower finds non-trivial cancellation.

The headline example is the derivative of an Integrate'd
rational function:

```
Simplify[D[Integrate[a x/(x^3+2), x], x]]
  →  (a x)/(2 + x^3)
```

Current cost breakdown (wall clock, current `main`):

| Step | Cost | Notes |
|------|------|-------|
| `D[Integrate[…]]` arg-evaluation | ~25 ms | ~9 internal `Together` calls |
| `builtin_simplify` alg_top branch | ~85 ms | One `Together[expr, Extension -> Automatic]` |
| Surrounding `simp_bottomup` polish | ~25 ms | Per-subnode work |
| **Simplify total** | **~134 ms** | |

The 85 ms `Together[expr, Extension -> Automatic]` call is the
dominant single line item. **Of that 85 ms**, only ~10 ms is in
the productive work (`extension_autodetect` + `qa_cancel_with_tower`);
the other ~75 ms is in **117 additional internal calls to
`builtin_together_compute`** — the function is reentered ~118 times
during a single user-level `Together[…, Extension -> Automatic]` call.

The Layer-4 goal is to cut that internal recursion.

---

## 2. Profiling evidence

### 2.1 Call counts (current `main`)

For a single user-level `Together[D[Integrate[a x/(x^3+2), x], x],
Extension -> Automatic]`:

- `builtin_together_compute` invocations: **~118**
- `extension_autodetect` invocations: **1** (the top-level one)
- `qa_cancel_with_tower` invocations: **1**
- Total wall: **~85 ms**
- Productive wall (autodetect + qa_cancel direct on the same input
  with the same tower, reused): **~10 ms**
- Overhead wall: **~75 ms** spread across the 117 extra calls

### 2.2 Per-call cost

Most of the 118 calls finish in < 0.5 ms each. They share the
following non-trivial per-call cost:

- `extract_extension_option_full` — option parsing
- `internal_args_contain_inexact` — O(input) tree walk
- `expr_has_nested_radical_radicand` — O(input) tree walk
- For inner calls (no `Extension` option, no inexact): straight to
  `together_recursive(arg)`

`together_recursive` itself dispatches by head:
- `Plus` → fraction-combine via `PolynomialLCM` + `cancel_recursive`
- `List` / relational heads → recurse and rebuild
- Other heads (e.g. `Power`, `Times` at the leaves) → recurse, rebuild
  via `eval_and_free`, then `cancel_recursive` on the result

The cost per inner call is ~0.6 ms on average. Roughly half of that
is the two O(input) prefilter walks; the other half is the actual
recursive descent.

### 2.3 Why so many calls?

The 117 inner calls are NOT direct recursion in `together_recursive`
(which is `static` and never re-enters the `builtin_together`
dispatcher). They come from two sources:

**Source A: arg evaluation in the user-level call.**
`D[Integrate[a x/(x^3+2), x], x]` itself fires ~9 internal `Together`
calls during partial-fraction normalisation and derivative
canonicalisation. These are unavoidable at the call site but
contribute to the 118 total.

**Source B (dominant): `qa_cancel_with_tower`'s Step 4 cascade.**

```c
/* qa_cancel_with_tower, Step 4 (qafactor.c:1850-1854) */
Expr* together_call = expr_new_function(expr_new_symbol("Together"),
    (Expr*[]){arg_internal}, 1);
Expr* arg_combined = evaluate(together_call);
expr_free(together_call);
```

The `arg_internal` is the input after the per-α γ-substitution: a
Plus / Times tree over `Q[γ_internal, free_var_1, ..., free_var_k]`.
The substituted form has many sub-Plus / sub-Times / sub-Power nodes;
each one that `together_recursive` walks through ends up rebuilt via
`eval_and_free(expr_new_function(...))`, and the `evaluate()` step
re-fires the head's builtin — which in turn invokes the Together-on-
sub-result chain implicitly through `cancel_recursive` →
`PolynomialGCD` → polynomial-arithmetic builtins → … → some path
that loops back through `Together`.

**This last step — the exact closure that gets the 117 calls back to
`builtin_together_compute` — is the open question that Layer 4 must
pin down concretely before designing the replacement.** Hypotheses
worth verifying with stack-trace instrumentation:

1. `cancel_recursive`'s `PolynomialGCD` path internally calls
   `Together` to normalise its inputs.
2. The evaluator's infinite-evaluation fixed-point loop re-fires
   `Together` against intermediate results that retain a `Together`
   head.
3. Some builtin called from `cancel_recursive` (e.g. `Cancel`,
   `Apart`, `Factor`) implicitly calls `Together`.

A 30-minute instrumentation experiment (transient `fprintf` in the
wrapper + stack-walking) would settle this before any code is
written.

---

## 3. Design alternatives

### Option A — Native QAUPoly combine, replacing Step 4 entirely

Replace the `evaluate(Together[arg_internal])` call in
`qa_cancel_with_tower` with a direct routine
`qaupoly_combine_sum(arg_internal, free_var, tower)`:

1. Walk `arg_internal`'s top-level Plus.
2. For each summand, extract `(num_expr, den_expr)` via the existing
   `extract_num_den` helper.
3. Lift each pair to `(QAUPoly num_i, QAUPoly den_i)` over
   Q(γ)[free_var] via `qa_expr_to_qaupoly`.
4. Compute the LCM denominator `D = lcm(den_1, …, den_n)` via
   `qaupoly_lcm` (a new wrapper around `qaupoly_gcd` + `qaupoly_mul`
   + `qaupoly_divrem`).
5. Compute combined numerator `N = Σ num_i · (D / den_i)` via
   `qaupoly_mul` + `qaupoly_add`.
6. Cancel `gcd(N, D)`.
7. Return `(N / D)` as either a `(QAUPoly, QAUPoly)` pair (feeding
   directly into `qa_cancel_with_tower`'s existing Step 6/7) or as a
   rebuilt Expr.

**Pros:**
- Bypasses the entire 117-call Expr round-trip.
- Stays in the qaupoly substrate — no allocator overhead from
  Expr trees, no symbol-table lookups, no fixed-point evaluation.
- Naturally composes with Layer 2's input pre-decomposition (the
  Sqrt[15] → Sqrt[3]·Sqrt[5] rewrite still happens at Step 0 entry).

**Cons:**
- Adding a new substrate routine (`qaupoly_lcm` + the combine
  driver) is ~300-500 lines of careful code.
- Edge cases: summands that aren't structurally
  `num · Power[den, -1]` (e.g. `Power[E, x] + 1/(x-1)`) need
  fallback to the Expr path.
- Round-trip into `qa_expr_to_qaupoly` for every summand still
  has cost; not zero work.

**Risk:** medium. Subtle correctness bugs in QAUPoly arithmetic
would manifest as silently-wrong results. Needs extensive
differential testing against the current path.

**Estimate:** 1–2 weeks.

### Option B — Make the inner Together calls cheap

If profiling shows the 117 extra calls are mostly trivial (Together
on already-canonical Times / Power inputs), a smart fast-path in
`builtin_together_compute` could skip the full `together_recursive`
walk for those.

**Attempted variant in this session:** an O(input) "trivial form"
predicate that returned `expr_copy(arg)` when `arg` had no Plus,
no negative-exponent Power, and no Exp. Result: **regressed** the
headline by ~20 ms because the predicate walk costs more on
non-trivial inputs than it saves on trivial ones (most non-trivial
inputs have Plus at the top, so the walk returns false in O(1) —
fine — but the cumulative cost across many calls was still net
negative).

The predicate has to be either:
1. O(1) on non-trivial inputs (the variant above almost achieves
   this — Plus at the top → return false in O(1)) AND
2. catch enough trivial inputs to net-win.

A variant gated on the input being an **atom or a `Power[atom,
non-negative integer]`** is genuinely O(1) and catches simple cases.
But it likely doesn't catch enough — many of the 117 calls have
small but compound inputs (e.g. `Times[2, Sqrt[3], Sqrt[5]]`).

**Pros:** small, local change in `rat.c`. Low risk.

**Cons:** modest payoff. The 117-call problem is fundamentally
structural; trimming per-call cost only saves
`117 × small_constant`.

**Estimate:** 1 day (including the differential-test sweep).

### Option C — Memoize at the Expr level (extend Layer 5)

The current `FactorMemo` caches `Together[res]` and `Cancel[res]`
keyed by the full input Expr. But each of the 117 inner calls
passes a freshly-constructed wrapper expression with a different
internal Expr pointer; structural `expr_eq` should still match if
the wrapped content is identical.

Hypothesis worth checking: are the 117 inner calls actually on
**distinct** sub-inputs (i.e. each call is unique), or are there
duplicates that would hit the cache?

If duplicates exist, instrument the cache to count hit/miss ratio.
A 50%+ hit rate means a tighter cache (e.g. unhashed direct
pointer compare for hot inner loops) could shave more time. A 5%
hit rate confirms Option C alone is insufficient.

**Pros:** orthogonal to Option A — composes additively.

**Cons:** wins are bounded by the hit rate. Likely modest.

**Estimate:** 1–2 days (including instrumentation).

### Option D — Reduce inner Together calls at the source

Find the exact path that closes the loop back to
`builtin_together_compute` (the "Source B" hypothesis verification
in §2.3) and break it. For example: if `cancel_recursive`'s
`PolynomialGCD` path is calling `Together` internally,
parameterise that call to skip when we're already inside a
Together cascade.

**Pros:** potentially the biggest win if the closure is a single
identifiable path. Could cut 117 → 1 call.

**Cons:** requires understanding the closure before designing the
fix. Risk: the closure is needed for correctness in non-Together
callers, and breaking it requires careful flag-passing.

**Estimate:** 1 day to investigate + 2-3 days to implement
depending on what's found.

---

## 4. Recommended approach

**Phase 1 (investigation, 0.5 day):** instrument
`builtin_together_compute` to print a stack snapshot for each of the
117 inner calls during a single headline `Together[…, Extension ->
Automatic]`. Group the call sites by frequency. Confirm or refute the
"Source B" hypotheses in §2.3.

**Phase 2 (one of the following, based on Phase 1):**

- **If the inner calls come from a single bounded closure** (e.g.
  `cancel_recursive` → `PolynomialGCD` → `Together`): go with
  **Option D**. Add a thread-local "inside-Together" flag that
  short-circuits `cancel_recursive`'s `Together` call when set.
  Estimate: 2 days.

- **If the inner calls come from the evaluator's fixed-point loop**:
  add a "no-op fixed point" check at `builtin_together`'s entry —
  if the input is already structurally a fraction, return it
  verbatim. Estimate: 1 day.

- **If the inner calls are genuinely diverse** (no single closure):
  go with **Option A** (native QAUPoly combine). Estimate: 1–2
  weeks.

In all cases, **Option C (extended memoisation) is a free add-on**
and can land in parallel.

---

## 5. Test plan (applies to any option)

### 5.1 Existing coverage to keep green

- `extension_auto_builtins_tests` — Together / Cancel / GCD / LCM /
  Apart with explicit α vs Automatic
- `extension_options_tests` — option parsing & propagation
- `simp_algebraic_cuberoot_tests` — includes the headline
  `test_simplify_headline_d_integrate`
- `qafactor_tests`, `simp_*_tests`, `simplify_tests`, `poly_tests`,
  `parfrac_tests`, `intrat_tests`, `mvfactor_tests`, `factor_*_tests`

### 5.2 New tests

- Performance regression test: `Simplify[D[Integrate[a x/(x^3+2),
  x], x]]` completes within a wall-clock budget (e.g.
  `TimeConstrained[Simplify[…], 0.08]` ≠ `$Aborted`).
- For Option A: differential test asserting `qa_cancel_with_tower`
  with the native combine and with the fallback Expr path produce
  expression-equal results on a corpus of multi-α inputs.

---

## 6. Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Subtle correctness bug in QAUPoly arithmetic (Option A) | Differential test against the current Expr path on a 50+ example corpus; fall back to the Expr path on any divergence. |
| Breaking the inner Together closure breaks other callers (Option D) | Use a thread-local entry-guard, not a global flag; only set when entering `builtin_together` from the top level and clear before any non-Together evaluation. |
| Phase-1 instrumentation is misleading on a single example | Profile 3–5 distinct inputs covering: pure rational, single-α algebraic, multi-α algebraic, multi-α with free vars, multi-α with no free vars. Patterns should agree. |
| Performance variance under `-O3` | Run each timing 5× and take median; gate the perf test on a 2× headroom relative to current `main`. |

---

## 7. Open questions

1. **Where exactly do the 117 inner Together calls come from?** This
   is the gating Phase-1 question. Answer determines whether to go
   with Option A, D, or both.
2. **Are duplicates among the 117?** Determines whether Option C
   (better memoisation) is worth pursuing.
3. **Does `cancel_recursive` call `Together` internally, directly or
   via PolynomialGCD?** Inspect `cancel_recursive` and
   `builtin_polynomial_gcd`.
4. **Does the evaluator's fixed-point loop re-fire `Together`?** Add
   a one-line counter at `evaluate_step`'s Together dispatch.

---

## 8. Out of scope

- Layer 4 does NOT touch the autodetect cost. That's now bounded by
  the Layer-0 input prefilter (which skips autodetect entirely on
  nested-radical inputs) and by the Layer-2 squarefree-prime
  coalesce (which shrinks the tower for dependent-Sqrt cases).
- Layer 4 does NOT alter the public API of `Together`, `Cancel`,
  `PolynomialGCD`, etc. All wins are below the API surface.
- Layer 4 does NOT target the `D[Integrate[…]]` arg-evaluation
  cost (~25 ms in the headline) — that's separate work in the
  `intrat` / `intsimp` modules.

---

## 9. References

- Commits this design supersedes:
  - `bb5b6ba` — Layers 0+1+3+5
  - `8d72345` — alg_top prefilter
  - `508dad2` — simp_algebraic_impl prefilter
  - `e41cd67` — Layer 2 squarefree-prime coalesce
- `docs/spec/changelog/2026-05.md` for shipped-state documentation.
- `src/qafactor.c:qa_cancel_with_tower` — Step 4 is the target.
- `src/rat.c:builtin_together_compute`,
  `src/rat.c:together_recursive`, `src/rat.c:cancel_recursive` —
  the call chain to instrument in Phase 1.
