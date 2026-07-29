# Project: FLINT-backed multivariate zero test (`is_zero_poly`)

## Goal

Replace the exponential-cost core of Mathilda's exact multivariate polynomial
zero test (`is_zero_poly`, `src/poly/poly.c`) with a FLINT `fmpq_mpoly`-backed
implementation, used **consistently across all of Mathilda**. The public
contract and every verdict must be preserved bit-for-bit; only the cost changes.

## Why (root-cause, already diagnosed)

`is_zero_poly(e)` decides whether `e` is the zero polynomial in its *opaque
generators* (the kernels `collect_variables` returns: bare symbols **and**
maximal non-polynomial subexpressions like `Log[x]`, `Sin[x]`, `x^x`, treated as
independent indeterminates). The current algorithm, `is_zero_poly_depth`
(poly.c ~L1114), at **every recursion level** does:

1. `expr_expand(e)`  — full distributed expansion (Expr schoolbook multiply),
2. `collect_variables`,
3. `is_polynomial`,
4. `internal_coefficientlist(expanded, vars[0])`  — split by one variable,
5. recurse into every coefficient (which re-expands, re-collects, …).

For a polynomial in *k* generators this fans out multiplicatively and each node
re-expands a large multivariate Expr. Profiling the Gruntz Gamma hang (thesis
8.31) showed **~all** wall-time here: `is_zero_poly_depth` self-time dominated
by `expr_expand`/`estimate_terms`/`multiply_two` (expand.c) + `is_polynomial`
(poly.c) + `CoefficientList`. Instrumentation showed **few** top-level calls,
each catastrophically expensive on large mrv-rewritten exp-log polynomials —
i.e. the cost is **per-call on big multi-generator inputs**, not repetition
(so memoization does not help) and not the non-zero-common-case (so numeric
pre-screening does not help — and a sound magnitude-aware numeric screen is
exactly the hard part the existing Schwartz–Zippel stage already solves).

`fmpq_mpoly` does the multiply+collect in packed FLINT arithmetic (already the
engine behind `Expand` via `flint_expand_polynomial`), turning the whole test
into one linear-ish conversion + `fmpq_mpoly_is_zero`.

## Blast radius (why "consistent across all of Mathilda" matters)

`is_zero_poly` has **~100 call sites** — every one benefits, and every one must
keep identical verdicts:

- `src/zero_test.c` — Stage 1 of `PossibleZeroQ` (`decide_rational`): the raw
  Expand path and the Together∘Cancel path both end in `is_zero_poly`.
- `src/rat.c`, `src/calculus/intsimp.c` — Together/Cancel numerator zero-tests.
- `src/calculus/intrat.c`, `intrischnorman.c` — rational + Risch integration.
- `src/poly/*` — `poly.c` internal (GCD/divide/PRS loops), `facpoly.c`,
  `squarefreeq.c`, `subresultants.c`/`subresultantpoly.c`, `zupoly.c`.
- `src/linalg/*` — `inv.c`, `linsolve.c`, `ludecomp.c`, `qrdecomp.c`,
  `svdecomp.c`, `eigen*.c`, `matrank.c`, `nullspace.c`, `util.c` (every
  symbolic pivot/rank/singularity decision).

Because the change is *purely* a faster implementation of the same predicate,
all callers are updated for free — no call-site edits.

## Existing FLINT infrastructure to build on (`src/poly/flint_bridge.{c,h}`)

- `flint_bridge_available()` — runtime gate (also `USE_FLINT` compile gate).
- `flint_expand_polynomial(e)` — expr → `fmpq_mpoly` → expr. Its private
  `collect_vars(e, VarSet*)` (flint_bridge.c ~L99) collects generators, BUT only
  **bare symbols** — it *declines* `Sqrt[..]`, `f[..]`, `Sin`/`Log`/… (comment:
  "those are extension generators, M2/M3 not M1"). This is the one gap to close.
- `flint_is_polynomial_over_q(e)` — predicate matching the above (symbols only).
- Conversion primitives (`fmpq_mpoly_ctx_init`, term build, render-back) already
  present and tested.

## Design

### New primitive

```c
/* flint_bridge.h */
/* Decide whether `e` is the zero polynomial in its opaque generators (the same
 * generator set is_zero_poly/collect_variables uses: symbols AND maximal
 * non-polynomial kernels like Log[x], Sin[x], x^x). Returns:
 *    1  -> provably the zero polynomial,
 *    0  -> provably non-zero,
 *   -1  -> declined (FLINT absent, or a generator/coefficient FLINT can't model
 *          — e.g. an irrational/symbolic exponent that isn't a clean kernel);
 *          caller falls back to the classical is_zero_poly_depth. */
int flint_mpoly_is_zero(const Expr* e);
```

### Kernel-aware conversion (the crux)

`is_zero_poly` treats `Log[x]`, `x^x`, etc. as indeterminates; flint_bridge's
`collect_vars` does not. Bridge them by **kernel substitution**, reusing the
kernel-aware collector that already exists in poly.c:

1. `collect_variables(e, &gens)` (poly.c) → the generator list `{g_1..g_n}`
   (symbols + opaque kernels), exactly what the classical path uses.
2. Substitute each `g_i` → a fresh symbol `$z i$` (structural `xreplace`, no
   eval), yielding a polynomial in bare symbols that flint_bridge already
   accepts. (Kernels are atomic under whole-node replacement — same treatment
   the classical CoefficientList recursion gives them, so **semantics are
   identical**.)
3. Convert to `fmpq_mpoly` via the existing term-build path; if any coefficient
   or exponent is not representable (non-integer exponent, inexact real,
   residual non-poly head after substitution), **return -1 (decline)**.
4. Result is `fmpq_mpoly_is_zero(p, ctx)` → 1/0.

Determinism + purity: no randomness, pure function of the input (matches the
classical predicate's contract).

### Wiring

`is_zero_poly` (poly.c) becomes:

```c
bool is_zero_poly(Expr* e) {
    int f = flint_mpoly_is_zero(e);          /* -1 if declined */
    if (f >= 0) return f == 1;
    return is_zero_poly_depth(e, 0);         /* classical fallback, unchanged */
}
```

The classical `is_zero_poly_depth` stays as the correctness reference and the
fallback for the declined cases and `USE_FLINT=0` builds.

## Correctness plan (paramount — this predicate underpins the CAS)

1. **Differential test harness.** Add a hidden debug mode (env-gated) that runs
   BOTH paths on every `is_zero_poly` call and aborts on any disagreement.
   Run the entire test suite + a large random-polynomial fuzzer (dense/sparse,
   many generators, kernel generators, rational coefficients, genuine zeros via
   `p - Expand[p]`, near-zeros) under it. Zero disagreements is the gate.
2. **Full regression** (must stay green, no new FAILs): `zero_test`,
   `trigexp_zero`, `simplify`, `fullsimplify`(+corpus), `rat`, `ratcanon_*`,
   `rationalize`, `intrat`(+corpus), `integrate_*`, and all `linalg`
   (`inv`/`linsolve`/`ludecomp`/`qr`/`svd`/`eigen`/`matrank`/`nullspace`), plus
   `limit`/`series`/`gruntz`.
3. **Decline-path audit.** Every `-1` must be a genuine "FLINT can't model it"
   (verified by the differential harness reaching the classical path with the
   same verdict), never a silent wrong answer.

## Performance plan

- Micro-bench: `is_zero_poly` on (a) dense univariate high-degree, (b) k-generator
  dense, (c) the mrv/Gamma-8.31 expressions, (d) large `Det`/`Inverse` of
  symbolic matrices. Record before/after.
- Macro-bench: `fullsimplify_corpus`, `intrat_corpus`, and a symbolic-linalg
  batch wall-time before/after.
- **Stretch validation:** re-enable `Gamma` isolation in the Gruntz engine
  (`src/calculus/gruntz.c`, `is_semitractable_head` + the `Exp[LogGamma]` branch
  — see GRUNTZ_STATE.md) and confirm thesis 5.5 / 8.31 / the Stirling ratio now
  resolve **without hanging**. This is the acceptance test that closes the loop
  with the work that motivated this project. If they resolve, add them as passing
  gruntz tests and move them out of `test_honest_abstentions`.

## Files

- **New**: `flint_mpoly_is_zero` in `src/poly/flint_bridge.{c,h}` (+ the
  kernel-substitution helper; may lift/generalize poly.c's `collect_variables`
  or call it directly — poly.c and flint_bridge are already linked).
- **Edit**: `src/poly/poly.c` — `is_zero_poly` dispatches to the FLINT path with
  the classical `is_zero_poly_depth` as fallback. Optionally add the env-gated
  differential-check harness here.
- **Tests**: extend `tests/test_zero_test.c` (kernel-generator zeros/non-zeros,
  decline cases) + a new random fuzzer; `tests/CMakeLists.txt` if a new binary.
- **Docs**: `docs/spec/builtins/simplification.md` (PossibleZeroQ internals note),
  weekly changelog; note the FLINT acceleration in `flint_bridge` docs.

## Risks / watch-items

- **Generator identity vs `collect_variables`.** The FLINT path MUST use the
  identical generator set as the classical path or verdicts can differ on
  algebraically-dependent kernels (e.g. `E^(2x)` vs `(E^x)^2`). Reusing poly.c's
  `collect_variables` verbatim guarantees this. (Note `E^x` = `Power[E,x]` is
  dropped by `collect_variables` because base `E` is numeric — preserve that.)
- **Rational coefficients / content.** `fmpq_mpoly` carries rational
  coefficients natively; ensure `Rational[p,q]` and `Complex[..]` coefficients
  either convert (Q, Q(i)) or cleanly decline.
- **`USE_FLINT=0` and macOS stale-object builds** (see the FLINT stale-artifacts
  lesson): the fallback must be a plain `make`-clean no-op path; test both gates.
- **No behavioural change** is the whole point — the differential harness is the
  primary safety net; do not ship without a clean differential run.

## Definition of done

1. `flint_mpoly_is_zero` implemented + wired; classical fallback intact.
2. Differential harness: zero disagreements across suite + fuzzer.
3. All listed regression suites green (no new FAILs); valgrind-clean.
4. Before/after benchmarks recorded (expect large wins on multi-generator inputs).
5. Gamma-Gruntz acceptance: 8.31 / 5.5 / Stirling-ratio resolve without hanging
   (or, if still limited by residual mrv cost, documented with the new numbers).
6. Docs + changelog + memory note updated.
