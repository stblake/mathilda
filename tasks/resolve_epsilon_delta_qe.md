# Scope: automatic ε–δ limit proofs via `Resolve`/`ForAll`/`Exists`

## Goal
Make the natural ε–δ limit sentence decide automatically for **semialgebraic**
`f` (polynomial / rational / algebraic):

```
Resolve[ForAll[eps, eps>0,
   Exists[del, del>0,
     ForAll[x, 0 < Abs[x-a] < del, Abs[f - L] < eps]]], Reals]   ->  True | False
```

Keep the hard soundness invariant: on anything out of scope, **decline**
(return unevaluated) — never a wrong verdict. Transcendental `f` (e.g.
`Sin[x]/x`) stays permanently out of scope (no polynomial CAD) and must keep
declining, as `test_quantifiers_decline` already asserts.

## Verified state at v0.205 (empirical)
- The whole statement **declines**. The decision *core* (CAD, real-closed-field)
  is complete; the gaps are all in the QE **front-end** (`src/solve/reduce_qe.c`),
  proven by composing the levels by hand with `Abs`/chains pre-removed:
  level1 `Reduce[ForAll[x, x>2-d&&x<2+d&&x!=2, 3x>6-e&&3x<6+e],{d,e},Reals]` →
  `d<0||d==0||d>0&&e>=3d`; level2 `Reduce[Exists[del,del>0&&(…)],{e},Reals]` →
  `e>0`; level3 `Resolve[ForAll[e,e>0,e>0],Reals]` → `True`. Wrong limit (L=6)
  composes to `e>1` then `False`. So the target is reachable via front-end
  plumbing only.

## Root causes (each independently blocks the shape)

### Gap 1 — parametric QE never runs real-function preprocessing (DOMINANT)
`reduce.c:539-555` already eliminates `Abs`/`Min`/`Max`/`Piecewise`/radicals for
a direct `Reduce[expr,{vars},Reals]`: univariate via `reduce_realfn_preprocess`,
**multivariate via the public `reduce_piecewise_preprocess(e,vars,nv,changed)`**
(`reduce_realfn.h:66`). `eliminate_abs` is variable-agnostic (sign-splits on the
`Abs` argument regardless of which vars appear — `reduce_realfn.c:612`).
But the parametric path `qe_parametric` (`reduce_qe.c:197-215`) calls
`reduce_form_from_expr` + `reduce_cad_qe` **directly**, bypassing that pass. So a
δ,ε-parametric body containing `Abs` (every ε–δ inner `∀x`) is fed raw to the
CAD, which rejects non-polynomial atoms → decline.
- Evidence: `Reduce[ForAll[x, Abs[x]<d, x^2<9], {d}, Reals]` declines; poly twin
  `Reduce[ForAll[x, x^2<d, x^2<9], {d}]` → `d<=9`.
- The fully-quantified **decision** path (`qe_decide`→`qe_call_reduce`→
  `Reduce[...]`) does *not* have this gap (it re-enters `builtin_reduce`, hence
  the preprocessing) — that is why single-block decisions with `Abs` already work.

### Gap 2 — a bounded (3-arg) quantifier wrapping a different-kind quantifier
`qe_normalize` (`reduce_qe.c:121-154`), on a 3-argument quantifier, folds the
condition into `And`/`Or` and sets `cur = body_owned`, ending the peel loop
(`:142`). The subsequent `is_quantifier(cur)` (`:148`) then sees the `And`/`Or`,
not the inner quantifier it now buries → `alternating` stays false → the inner
quantifier reaches the CAD as a raw head → decline. The ε–δ prefix is exactly a
chain of *bounded* alternating quantifiers.
- Evidence: `ForAll[eps,eps>0,Exists[del,del>0,del<eps]]` declines, but the
  **unbounded** `ForAll[eps,Exists[del,del<eps]]` → `True`.
- Bonus correctness bug this also fixes: a same-kind chain of >1 bounded
  quantifiers (`ForAll[a,ca,ForAll[b,cb,M]]`) is currently mis-handled for the
  same reason (only the first condition folds, then the loop ends burying the 2nd
  `ForAll`).

### Gap 3 — chained `Inequality` containing `Abs` in a bounded condition
`0 < Abs[x-a] < del` parses to `Inequality[0,Less,Abs[..],Less,del]`. Once folded
under `Not` by the 3-arg path it fails to eliminate, though the `And`-spelled
condition succeeds.
- Evidence: `Resolve[ForAll[x, 0<Abs[x]<1, x^2<4], Reals]` declines (P6);
  `…, 0<Abs[x]&&Abs[x]<1, …` → `True` (P7); `…, 0<x<1, …` (no `Abs`) → `True`.
- Note: a direct `Reduce[Not[0<Abs[x]<1]||x^2<4, x, Reals]` → `True`, so the
  interaction is inside the QE fold, not the base engine. **Confirm the exact
  trigger with a printf/debugger before fixing**; the fix (normalize chained
  `Inequality` → `And` of binary relations in the condition) is validated by P7
  regardless of mechanism.

## Design — three localized changes, all in `src/solve/reduce_qe.c`

### Phase 1 — Gap 1 (reuse existing multivariate preprocessing)  [~10 LOC]
In `qe_parametric`, before `reduce_form_from_expr`, mirror `reduce.c:545-555`:
```c
Expr* pre = NULL; bool ch = false;
if (reduce_stmt_has_piecewise(body, vall, nvall)
 || reduce_stmt_has_radical (body, vall, nvall))
    pre = reduce_piecewise_preprocess(body, vall, nvall, &ch);
const Expr* use = pre ? pre : body;
RForm* F = reduce_form_from_expr(use, vall, nvall, &ok);
...
expr_free(pre);   /* free before returning; NULL-safe */
```
Include `reduce_realfn.h`. No change to `reduce_realfn.{c,h}`. This alone makes
`P3`/`P4` and every ε–δ inner `∀x` produce a δ,ε formula. Inherited for free by
the alternating recursion (it re-enters `qe_parametric`/`qe_decide`).

### Phase 2 — Gap 2 (thread a block side-condition; stop burying inner quantifiers)  [~40-60 LOC]
Rewrite `qe_normalize` + the two dispatch branches so the condition-fold no
longer terminates the peel:
- `qe_normalize` gains an out-param `Expr** out_cond` (owned, or NULL). On a
  3-arg quantifier: record bound vars, **conjoin** its condition into the running
  block condition, and continue into `args[ac-1]` (same as the 2-arg case) — do
  **not** build `And`/`Or` or set `cur=body_owned`. Stop when `cur` is not a
  same-kind quantifier; set `alternating = is_quantifier(cur)` as today.
- `reduce_qe_dispatch` applies the condition when forming the block matrix `M`:
  - ForAll block: `M := C ? Or[Not[C], M] : M`
  - Exists block: `M := C ? And[C, M] : M`
  - **non-alternating** branch: `M := apply(C, cur)` then `qe_decide`/`qe_parametric`.
  - **alternating** branch: eliminate `cur` recursively → `psi`; `M := apply(C, psi)`;
    then `qe_rebuild_quant(quant,B,nb,M)` and re-dispatch (now a clean 2-arg block).
- Backward-compat: `C==NULL` (all 2-arg) reproduces today's behavior byte-for-byte;
  a single 3-arg block reproduces today's `Or[Not[C],M]`/`And[C,M]` fold. So the
  existing `test_quantifiers_*` stay green; only genuinely-bounded-alternating
  inputs change from *decline* to *decided*.
- Memory: `C` and each intermediate are owned; free on every path; **valgrind
  gate** this file's paths (the module is ownership-careful — see
  `[[feedback_builtin_res_ownership]]`).

### Phase 3 — Gap 3 (expand chained `Inequality` in conditions)  [~15 LOC]
Add a small `expand_inequality_chains(e)` that rewrites
`Inequality[a,op1,b,op2,c,…]` → `And[a op1 b, b op2 c, …]`, and apply it to the
block condition `C` in Phase 2 (before `Not`-wrapping). P7 proves the `And` form
composes end-to-end. Re-confirm the exact trigger first (see Gap 3 note).

## Non-goals / explicitly out of scope (document, don't build now)
- **Quantifiers nested inside boolean connectives between blocks** (`ForAll[eps,
  Implies[eps>0, Exists[...]]]`, I4/I5). The *natural* ε–δ form nests quantifiers
  directly via the 3-arg bounded spelling, so Phases 1-3 cover it. General prenex
  extraction / mini-scoping is a larger separate effort — note as future work.
- **Transcendental `f`** — no polynomial CAD; must keep declining (soundness net).
- **CAD cost** — doubly-exponential; fine for linear/low-degree `f` with a couple
  of `Abs`, but a high-degree `f` or many nested `Abs` may be slow. Existing
  `TimeConstrained`-style guards apply; not a correctness issue.

## Verification
- New `tests/test_reduce.c::test_epsilon_delta` (via `run_test`/`run_contains`):
  - linear `3x-1` at 2: `=5` → `True`, `=6` → `False`;
  - quadratic `x^2`: at 0 `=0` → `True`, at 3 `=9` → `True`, `=8` → `False`;
  - a same-kind bounded chain (Gap-2 bonus) e.g. `ForAll[{a,b}, a>0&&b>0, …]`;
  - **soundness**: `Resolve[…Sin[x]…]` and a high-degree transcendental stay
    unevaluated (`run_contains` "Resolve[").
- Regression: full `test_reduce` suite + the reduce corpus qe-* rows stay green.
- `valgrind --leak-check=full` over the new `reduce_qe.c` paths.
- `make check-c99` (no new POSIX symbols expected). No packed/NDArray/Compile
  surfaces — this is symbolic/structural QE returning `True`/`False`/formulas, so
  the numeric-fastpath audits do not apply.

## Files
- `src/solve/reduce_qe.c` — Phases 1-3 (only file with logic changes).
- `tests/test_reduce.c` — new test group.
- `docs/spec/builtins/solutions-of-equations.md` — note Resolve/ForAll/Exists now
  decide bounded-alternating semialgebraic sentences (with the ε–δ example).
- `docs/spec/changelog/2026-09-21.md` — change summary (Monday of this ISO week).
- `src/version.h` — bump per substantive commit (+0.001), one git tag each.

## Effort & sequencing
Phase 1 first (biggest payoff, ~0.5 day incl. tests) → Phase 2 (~1-1.5 day,
the careful one) → Phase 3 (~0.25 day). Land as small tagged commits so each is
independently green. Total ≈ 2-2.5 days.

## Review section

### Phase 1 — DONE (v0.206), verified
- **Change** (`src/solve/reduce_qe.c`, `qe_parametric`): before `reduce_form_from_expr`,
  run `reduce_piecewise_preprocess(body, vall, nvall, &changed)` guarded by
  `reduce_stmt_has_piecewise || reduce_stmt_has_radical` — the *same* multivariate
  Reals policy `reduce.c:545-555` uses. General, not Abs-specific: `Abs`/`Min`/`Max`/
  `Piecewise`/`Sign`/`UnitStep`/… and √-radicals all case-split uniformly. Ownership
  airtight (`pre` owned + freed on every path; `body` borrowed). `+#include "reduce_realfn.h"`.
- **Verified**: `P3 -> d<=3`, `P4 -> e>=1`, two-Abs inner level ->
  `d<0||d==0||d>0&&e>=3d` (matches hand-derived polynomial), `Min` case works
  (generality). Regressions green: `reduce_tests`, `reduce_corpus_tests` (170/170),
  `solve_tests`; transcendental still declines (soundness). New tests added to
  `test_quantifiers_parametric` incl. form-agnostic parametric≡decision guards.
  valgrind: no new leak (420-block baseline unchanged on trivial input; no leak
  stack implicates the path at 30 callers). Docs + changelog + `$VersionNumber`
  bumped to 0.206.
- **Not yet committed/tagged** — awaiting user go-ahead on git.

### Phase 2 — DONE (v0.207), verified
- **qe_normalize** now accumulates each block's 3-arg restriction as a
  side-condition (`out_cond`) instead of folding it into the body and stopping the
  peel; keeps descending so an inner different-kind quantifier is no longer buried.
  **qe_apply_condition** combines the condition with the block matrix in both
  dispatch branches (`ForAll: !c||M`, `Exists: c&&M`). Strict generalisation:
  `cond==NULL` (2-arg chain) and a single 3-arg block reproduce the old matrix
  exactly — all prior tests unchanged. Also fixed a latent same-kind-chain bug
  (2nd+ restriction was dropped).
- **qe_decide** made robust: `ForAll[{v},g]` now decided by `Reduce[!g] === False`
  (¬g unsatisfiable), not `Reduce[g] === True`, so a tautological *region* that
  `Reduce` doesn't collapse to `True` (`a<=0||b<=0||b>-a`) is no longer misread as
  False. Latent weakness exposed once bounded multi-var ForAll decisions became
  reachable; `Exists` direction unchanged (already robust).
- **Verified**: full ε–δ now decides — `lim(3x-1)@2=5 -> True`, `=6 -> False`;
  `lim x^2@3=9 -> True`, `=8 -> False`; `lim x^2@0=0 -> True`. Bounded alternation
  `ForAll[eps,eps>0,Exists[del,del>0,del<eps]] -> True`; same-kind bounded chain
  `ForAll[a,a>0,ForAll[b,b>0,a+b>0]] -> True`; unbounded `ForAll[a,ForAll[b,a+b>0]]
  -> False`. Transcendental `Sin[x]/x` still declines. Regression green
  (reduce/corpus 170/170/solve); valgrind no new leak. New test groups
  `test_quantifiers_bounded_alternation`, `test_epsilon_delta`. Docs + changelog +
  `$VersionNumber` 0.207.
- **Not yet committed/tagged** — awaiting user go-ahead on git.

### Phase 3 — CLOSED by Phases 1+2 (no code needed), verified
- Investigated first (per "no hacks"): at v0.207 the chained `Inequality`
  (`0<Abs[x-a]<del`) spelling **already decides on every path**. A 10-case
  adversarial battery — decision/parametric × ForAll/Exists, chained+Abs, mixed
  `<`/`<=`, raw `Not[chained+Abs]` multivar, full nested ε–δ (both limits) — all
  correct. So the originally-scoped `Inequality`->`And` normalisation is NOT
  needed; adding it would be a redundant patch.
- **Why it works without new code**: (decision) qe_decide (v0.207) tests the
  NEGATION, so the chained `Inequality` appears positively; (parametric) Phase 1's
  `reduce_piecewise_preprocess` descends into `Not`→`Inequality`→`Abs` and splits
  the `Abs`, leaving a *polynomial* chained `Inequality` that `reduce_form` already
  handled (even under `Not`) at v0.205 (P6b). The v0.205 decline was Gap 1 (no
  parametric Abs pass), not a distinct `Inequality` gap.
- **Locked in**: `test_chained_inequality_abs` + chained ε–δ cases added to
  `test_epsilon_delta`; docs/changelog corrected to drop the "pending" caveat. No
  version bump (behaviour already shipped in 0.207; tests+docs are
  contributor-facing).

## Outcome
`Resolve`/`ForAll`/`Exists` now prove semialgebraic ε–δ limits automatically, in
the natural chained spelling, sound in both directions (True/False) and declining
on transcendentals. Gaps 1 (v0.206) and 2 (v0.207) closed; Gap 3 dissolved.
