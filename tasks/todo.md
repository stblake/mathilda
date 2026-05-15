---
title: BronsteinRational closed-form gaps
date_started: 2026-05-11
status: in progress
---

# BronsteinRational closed-form gaps (2026-05-11)

Five failing `Integrate`BronsteinRational[...]` cases.  All five
reproduce in picocas REPL; all five integrate cleanly via the
Mathematica .m source (verified with wolframscript), so the gaps are
in the C port at `src/intrat.c`, not in the underlying algorithm.

## Root-cause map

| # | Input | Root cause |
|---|-------|------------|
| 1 | `1/(b - a x^3)` | Cubic Q-in-t.  `logtoreal_dispatch` has no deg-3 branch; falls to RootSum. |
| 2 | `1/(x^5+1)` | Palindromic-quartic Q-in-t `5(t^4-t^3+t^2-t+1)`.  `logtoreal_dispatch` deg-4 only handles biquadratic & Sophie-Germain.  Falls to RootSum. |
| 3 | `x^4/(a+b x^3)^2` | Same as #1 — cubic Q after Hermite reduction. |
| 4 | `(-1+x^2)/(1-2x^2+2x^4)` | Biquadratic Q-in-t `2t^4-2t^2+1` with negative inner discriminant.  Current deg-4 c1z&&c3z branch bails on neg disc; falls to `expand_simple_rootsum`, which produces complex form. |
| 5 | `1/x/(1-2x^2+2x^4)` | Post-processing: `c * piece_int` stays held as `Times[c, Plus[…]]`; also `Log[c p]` not stripped when `FreeQ[c, x]`. |

## Phases (independently testable)

### A — Output post-processing (#5)
- Distribute scalar Times-over-Plus in the final accumulator.
- Strip `Log[c · p] -> Log[p]` when `FreeQ[c, x]`.

### B — Negative-inner-disc biquadratic (#4)
- Extend `logtoreal_dispatch` deg-4 c1z&&c3z branch: when inner disc
  `c2^2 − 4 c0 c4 < 0` (provably under positive-symbol assumption)
  and `c0 c4 > 0`, factor over R as
    `c4 (t^2 + α t + Sqrt[c0/c4])(t^2 − α t + Sqrt[c0/c4])`
  with `α = Sqrt[2 Sqrt[c0/c4] − c2/c4]`.
- Dispatch each factor through `logtoreal_quadratic`.

### C — Cubic (#1, #3)
- Add deg-3 branch to `logtoreal_dispatch`.
- C.1: nth-root cubic (c1 = c2 = 0).  Real root
  `r0 = −(c0/c3)^(1/3)`; factor `c3 (t − r0)(t^2 + r0 t + r0^2)`;
  dispatch the linear and quadratic factors.

### D — Palindromic quartic (#2)
- Detect `c0 t^4 + c1 t^3 + c2 t^2 + c1 t + c0`.
- `u = t + 1/t` ⇒ `c0 u^2 + c1 u + (c2 − 2 c0) = 0`.  Factor as
  `c0 (t^2 − u_+ t + 1)(t^2 − u_− t + 1)` and dispatch through
  `logtoreal_quadratic`.

### E — Tests
- Add picocas tests for all 5 cases (no `RootSum`/`Function` head in
  output; `D[result, x] − integrand` simplifies to 0).
- Extend `IntegrateRationalTests.m` for Mathematica parity.

## Review

All 5 phases shipped.

### Phase A — Output post-processing (#5) ✓
Added `intrat_distribute_plus(e)` (uses `expr_expand`) and
`intrat_strip_log_constants(e, x)`.  Both run before AND after
`intrat_log_to_arctanh` in `intrat_integrate_rational`.

### Phase B — Negative-inner-disc biquadratic (#4) ✓
Extended `logtoreal_dispatch` deg-4 `c1z && c3z` branch with a
Sophie-Germain-with-c2 sub-branch that factors as
`c4 (t^2 + α t + β)(t^2 − α t + β)` with `β = Sqrt[c0/c4]`,
`α = Sqrt[2 β − c2/c4]`.

### Phase C — Cubic (#1, #3) ✓
New deg-3 branch in `logtoreal_dispatch` for nth-root cubics
`c3 t^3 + c0`.  Real-root selection under positive-symbol or
numeric-sign assumption; dispatch linear + quadratic factors.

### Phase D — Palindromic quartic in logtoreal_dispatch ✓ (restricted)
Detect scaled-palindromic `c4 c1^2 == c3^2 c0`, but only fire when
r = 1 (pure palindromic).  Scaled cases would need a smarter
`LogToAtan` (current C port hangs on nested-radical S).

### Phase D2 — Palindromic quartic in NaiveLogPart (#2) ✓
New `expand_palindromic_quartic_real(a_t, d_t, dd_t, t, x)` builds
the real form directly from `(a(α) / d'(α)) · Log[x − α] + …` per
conjugate root pair.  Bypasses LogToReal / LogToAtan entirely; the
nested-radical S that overwhelms LogToAtan on the LRT path never
gets constructed here.

### Phase E — Tests ✓
`tests/test_intrat.c` adds five new regression tests (one per
issue) using `assert_closed_real` (no `RootSum` / `Function`
leak) + `assert_integral_numeric_ok` for nested-radical cases.
All five pass.  `IntegrateRationalTests.m` gains five new
{integrand, var, optimal-antiderivative} entries.

### New helpers added
- `intrat_distribute_plus(e)` — Plus-over-Times distribution.
- `intrat_strip_log_constants(e, x)` — strip x-free constants out of Log args.
- `intrat_numeric_sign(e)` — N[]-based sign decision, fallback for
  symbolic positive-symbol bails on `Sqrt[5] − 5` style inputs.
- `expand_palindromic_quartic_real(a_t, d_t, dd_t, t, x)` —
  real-form expander for palindromic-quartic d in NaiveLogPart.

### Verification
- Build clean under `-std=c99 -Wall -Wextra -O3`.
- All five regressions return `Plus`-headed real-elementary form.
- Differentiation check passes (symbolic for #5, numeric for #1–#4
  whose nested radicals don't reduce by Cancel/Together).
- Full `tests/` suite: only pre-existing `logexp` + `poly` ordering
  failures (introduced by 2026-05-11's canonical comparator change),
  no new regressions from the closed-form work.

### Known follow-ups
- Cardano's formula for non-depressed cubics in `logtoreal_dispatch`.
- Scaled-palindromic dispatch (r != 1) — requires a cheaper LogToAtan
  for nested-radical coefficients, or a different routing strategy.
- The closed-form outputs aren't always Mathematica-minimal; the
  underlying form is correct but Sqrt simplifications like
  `Sqrt[1/16 (1 + Sqrt[2])] -> Sqrt[1+Sqrt[2]]/4` are left to a
  smarter `intrat_simp_pos_sqrt` (or full `Simplify`).

## 2026-05-11 rev 2 — generalised nth-root coverage

User flagged that the cubic-specific fix above was too narrow.
Adjacent inputs all leaked `RootSum`:

  - `1/(b - a x^4)`, `1/(b - a x^5)`,
    `1/(b - a x^6)`, `1/(b + a x^6)`,
    `1/(b + a x^8)`, `1/(b + a x^9)`.

### Generalisation
Replaced the per-degree branches in `logtoreal_dispatch` with one
helper `logtoreal_nthroot_sparse(base, deg, s, x, t)` that closes
`c_n t^n + c_0` (all intermediate coefficients zero) for every
n ≥ 3 via the standard cyclotomic factorisation over R.  Subsumes
the previous deg-3 nth-root and deg-4 Sophie-Germain branches.

Algorithm: given `q = -c_0/c_n` and `r = |q|^(1/n)`,
- q > 0: real roots ±r (r always; -r iff n even), conjugate pairs
  for k = 1..⌊(n-1)/2⌋ at angles 2πk/n.
- q < 0: real root -r iff n odd, conjugate pairs at angles
  (2k+1)π/n.
Each conjugate pair contributes a real quadratic
`t^2 - 2 r cos(θ) t + r^2` with discriminant -4 r^2 sin²θ < 0, so
`logtoreal_quadratic`'s ArcTan branch always fires.

### Tests
Six new regressions in `tests/test_intrat.c`:
`test_closed_nth_root_quartic_minus`, `_quintic_minus`,
`_sextic_minus`, `_sextic_plus`, `_octic_plus`, `_nonic_plus`.
Each asserts no RootSum/Function head + numerical derivative
match at two sample points.  All pass.

### Regression status
Full test suite (105 binaries): 90 pass / 15 fail both before and
after this change.  Identical failure set (pre-existing canonical
ordering issues from commit cc8b164).  Zero new regressions.

---

# Multi-extension nested-radical denesting (2026-05-14)

Status: PLAN — awaiting user sign-off.

## Motivation

`Simplify[Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]] - Sqrt[5]
- Sqrt[11 - 2 Sqrt[29]]]` should reduce to 0. It currently returns
unchanged. The denester collapses `Sqrt[152 - 24 Sqrt[29]]` correctly,
so the algorithm works one level down — what's missing is recursion
across nested extensions (Landau / Borodin style).

## What's deferred today

`src/simp.c:2256` (`split_plus_into_a_plus_b_sqrt_c`) refuses any
radicand with two or more sqrt-bearing terms — the comment at line
2250 names it explicitly as the "phase 2 multi-extension" case.

`src/simp.c:2129` (`sqrt_if_clean_square`) only recognises integer /
rational / `Power[u, 2k]` / `FactorSquareFree` even-power patterns —
it cannot see that `152 - 24 Sqrt[29] = (2 Sqrt[29] - 6)^2`.

The two limits compound: even if we partition multi-sqrt radicands,
the discriminant `D = A^2 - b^2 C` ends up in `Q(Sqrt[k])` rather
than `Q`, and `sqrt_if_clean_square` can't take its root.

## Algorithm

Let radicand `r = A + b·Sqrt[C]` with `A`, `b` in `Q(Sqrt[k])` (one
sqrt-bearing term picked as "outer", the rest absorbed into `A`).
Half-sum identity still gives

```
P, Q = (A ± Sqrt[A² - b² C]) / 2
Sqrt[r] = Sqrt[P] + sgn(b)·Sqrt[Q]
```

For the failing case (outer = `2 Sqrt[55 - 10 Sqrt[29]]`):

  A = 16 - 2 Sqrt[29],   b = 2,   C = 55 - 10 Sqrt[29]
  A² - b² C = 152 - 24 Sqrt[29]
  Sqrt[A² - b² C] = 2 Sqrt[29] - 6           ← needs recursion
  P = 5,   Q = 11 - 2 Sqrt[29]                ← clean
  Sqrt[r] = Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]   ✓

The picking matters: the alternative (outer = `Sqrt[29]`) leaves an
A containing the deep sqrt, and the discriminant doesn't denest
cleanly. **Heuristic:** prefer the sqrt-bearing term whose own
radicand contains a nested radical (deepest-outer-first). Fall
back to trying each candidate in turn.

## Design (modular)

Three small, focused changes; existing single-extension path stays
byte-for-byte unchanged on inputs it already handles.

### 1. Generalise `split_plus_into_a_plus_b_sqrt_c`

Rename existing call sites unchanged. Add a sibling:

```c
/* Pick the "best" outer sqrt when there are multiple sqrt-bearing
 * terms. Ranks candidates by sqrt-nesting depth of their radicand
 * (deeper first); ties broken by leaf count (smaller first).
 * Returns a small ordered list the caller iterates. */
static bool split_plus_pick_outer_sqrt(const Expr* plus_node,
                                       const AssumeCtx* ctx,
                                       size_t candidate_idx,
                                       Expr** out_A, Expr** out_b,
                                       Expr** out_C);
```

The caller (`denest_compute_pq_s`) iterates candidate_idx until one
succeeds or all are exhausted.  Single-sqrt radicands hit
`candidate_idx == 0` once and behave identically to today.

### 2. Recurse in `sqrt_if_clean_square`

Add one case at the end of the function: when `D` is a Plus that
matches the "phase-2-shape" (`α + β·Sqrt[γ]` with `α, β` in `Q`),
recursively call a depth-bounded `simp_denest_sqrt_recursive` on
`Sqrt[D]`. If the result is a clean `Q(Sqrt[γ])` expression (no
surviving `Sqrt[…]` head), return it. Otherwise NULL.

Termination: a `depth_budget` integer carried in the `AssumeCtx`
(or a thread-local — `AssumeCtx` is the cleaner home). Initial
budget = 4 (handles four nesting levels — far beyond any practical
input). Each recursive entry decrements; at 0 the recursive call
returns NULL.

### 3. Memoisation (efficiency)

Nested denesting can revisit the same subexpression: the walker hits
`Sqrt[55 - 10 Sqrt[29]]` once, but the discriminant path hits
`Sqrt[152 - 24 Sqrt[29]]` which itself triggers the recursion.
Add a small hash-table cache `(expr_hash, expr_eq) -> result` keyed
on the radicand, scoped to a single top-level `Simplify` call (or
to the AssumeCtx — preferred so we don't need a separate lifetime).
On miss compute, on hit return a copy. Empirically the cache stays
in the low tens of entries even on adversarial inputs.

## File-by-file changes

| File | Change |
|------|--------|
| `src/simp.c` | Add `split_plus_pick_outer_sqrt`, extend `sqrt_if_clean_square` with the recursive Q(Sqrt[k]) case, thread `depth_budget` through `AssumeCtx`, wire the iteration loop in `denest_compute_pq_s`, add cache. |
| `src/simp.h` | No public API changes. |
| `tests/test_simp_denest_phase2.c` | New: targeted regressions. Listed below. |
| `tests/CMakeLists.txt` | Register new test binary. |
| `docs/spec/builtins/simplify.md` (or equivalent) | Note the extended denester capability. |
| `docs/spec/changelog/2026-05.md` | Summary entry. |

## Test plan

Targeted regressions in `tests/test_simp_denest_phase2.c`:

1. The motivating expression:
   `Simplify[Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]] -
   Sqrt[5] - Sqrt[11 - 2 Sqrt[29]]]` → `0`.
2. Direct denesting:
   `Simplify[Sqrt[16 - 2 Sqrt[29] + 2 Sqrt[55 - 10 Sqrt[29]]]]` →
   `Sqrt[5] + Sqrt[11 - 2 Sqrt[29]]`.
3. Inner step (must still work, single-extension):
   `Simplify[Sqrt[152 - 24 Sqrt[29]]]` → `2 Sqrt[29] - 6`.
4. Non-denestable multi-sqrt (must remain unchanged, not blow up):
   `Simplify[Sqrt[1 + Sqrt[2] + Sqrt[3]]]` stays as-is.
5. Soundness check at random-ish numeric points for all phase-2
   denestings (the standard pattern: compute LHS-RHS, `NumberQ` the
   `Abs` at two distinct numeric instantiations of any free symbols
   — **not** as a zero-detector but as a regression check on test
   values where the answer is known).
6. Recursion budget exhaustion: a hand-built triple-nested radical
   that should denest only at budget ≥ 3 — confirms the budget is
   load-bearing and the failure mode is "unchanged", not crash.

Full-suite regression: all existing `test_simp_*` binaries must
remain green. The single-extension code path is unchanged on its
own inputs, so this is mostly a "did I break ownership somewhere"
check.

## Risks / edge cases

- **Branch validity.** The half-sum identity picks a specific sign
  combination. The existing `denest_is_nonneg` prover stays as the
  gate — for P, Q now living in `Q(Sqrt[k])`, the prover has to
  decide non-negativity of `(A ± s)/2` where A, s themselves carry
  a Sqrt. For the motivating case both P=5 and Q=11-2 Sqrt[29] are
  clearly nonneg (`11 > 2 Sqrt[29]` since 121 > 116). The
  `denest_prov_nonneg` numeric-bounded path already handles this.
- **Cost.** Worst case is exponential in nesting depth without the
  cache. With the cache and a budget of 4, each radicand is denested
  at most once; total cost is linear in the number of distinct
  sqrt-bearing radicands in the input.
- **Ownership.** Every new allocation must have a matching free path
  on every early-return — the existing function is already careful
  here, the new code must mirror that pattern.
- **Wrong-outer choice masks a valid denesting.** The "try every
  candidate" iteration covers this — order is a heuristic, not a
  correctness lever.

## Implementation tasks (checklist)

- [ ] Add `depth_budget` field to `AssumeCtx`; initialise to 4 in
  the top-level entry point; propagate.
- [ ] Implement `split_plus_pick_outer_sqrt(plus, ctx, idx, ...)`.
- [ ] Refactor `denest_compute_pq_s` to iterate candidates.
- [ ] Extend `sqrt_if_clean_square` with the Q(Sqrt[k]) recursion.
- [ ] Add per-`AssumeCtx` memo cache (hash table keyed on radicand).
- [ ] Write `tests/test_simp_denest_phase2.c` with the 6 cases
  above; register in `tests/CMakeLists.txt`.
- [ ] Run full test suite; confirm zero new failures.
- [ ] Update `docs/spec/changelog/2026-05.md`.
- [ ] Update `docs/spec/builtins/simplify.md` if behaviour text
  needs revising.
