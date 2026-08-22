---
ticket: DEMO-2
created: 2026-08-21T00:00:00-04:00
researcher: Michael Sollami
source_sha: ea0f8c3c308cc38a7bb41f410066b9ff290f260d
branch: claude/where-are-you-5f1147
repository: mathilda
topic: "NMinimize returns constrained solutions that violate their constraints by ~1e-4"
tags: [research, codebase, nminimize, feasibility, constraints, DEMO-2]
subsystems: [numerical_calculus]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-21
last_updated_by: Michael Sollami
---

# Research: NMinimize constrained-feasibility bug

**Date**: 2026-08-21
**Researcher**: Michael Sollami
**Git Commit**: ea0f8c3c308cc38a7bb41f410066b9ff290f260d
**Branch**: claude/where-are-you-5f1147
**Repository**: mathilda
**Ticket**: DEMO-2

## TL;DR

Asked how to fix a feasibility threshold tested against a squared penalty. Found the
root cause confirmed, but also that the obvious fix — returning a violation norm from
`nm_eval_pen` — is unsafe: six of eight engines combine `f + 1e6*penalty`
arithmetically, so rescaling the penalty reshapes their search dynamics. Fix is to
square the tolerance at the definition site instead. Nothing unresolved.

## Summary

The bug in `NMINIMIZE_FEASIBILITY_BUG.md` is confirmed exactly as documented:
`nm_eval_pen` accumulates `sum of violation^2` and that total is compared against
`NM_FEAS_EPS = 1e-8`, making the effective tolerance on a real constraint violation
`sqrt(1e-8) = 1e-4`.

The important new finding is about the *fix*, not the bug. `nm_eval_pen` has exactly
one caller (`nm_eval`), so the risk surface is fully enumerable — and enumerating it
refuted the initially-chosen approach. Changing what `nm_eval_pen` returns is safe
only for the three engines that use the penalty purely for comparison (DE,
RandomSearch, SHGO). The other five combine it arithmetically as
`f + NM_PENALTY_MU * p` with `NM_PENALTY_MU = 1e6`, feeding NelderMead's simplex
ranking, three Metropolis-style acceptance functions, and DIRECT's convex-hull
geometry. Rescaling the penalty from quadratic to linear in violation changes what
those five algorithms *search*, not merely what they *accept*.

The safe fix leaves the returned quantity alone and squares the tolerance where it is
defined, so the constant's comment becomes true and the relationship cannot drift.

## Open Questions

### Unresolved

_None._

### Resolved

- [x] **Which fix matches the design intent?** — Initially answered "make `nm_eval_pen`
  return a violation norm". Research then refuted that (see Detailed Findings §2), and
  the revised answer is: **keep the squared penalty, define the tolerance on the actual
  violation and square it at the definition site**, so the arithmetic engines are
  bit-identical and the semantics are explicit in code.
- [x] **Does `NM_FEAS_FINAL` get fixed too?** — Yes. It carries the identical
  squared-vs-unsquared ambiguity (`1e-6` on a squared penalty = `1e-3` on a real
  violation) and gets the same treatment, with its violation tolerance chosen
  deliberately rather than inherited.
- [x] **What is the acceptable speed cost?** — Correctness wins, but measure and report
  before accepting. Hard checkpoint: if constrained cases regress worse than ~3x
  against the DEMO-1 baseline (C1 1.35x, C2 1.17x vs scipy SLSQP), stop and re-decide
  rather than absorb it silently.
- [x] **What about problems that flip to `Infinity`?** — Keep the selection threshold
  tight (`1e-8` violation) but choose `NM_FEAS_FINAL`'s violation tolerance loosely
  enough that no currently-passing test flips to `Infinity`. Four tests assert exact
  `Infinity`/`Indeterminate` and two assert a result must stay finite; all six are
  named in Detailed Findings §3 and must be re-run.
- [x] **Is a test asserting feasibility actually missing?** — Yes, structurally. Every
  existing feasibility check is a *one-sided* tolerance band, which a more-feasible
  result trivially still satisfies. That one-sidedness is the exact mechanism by which
  a 1e-4 violation shipped undetected.

## Requires Approval

The revised fix approach was re-confirmed after research contradicted the first answer.
The speed checkpoint (~3x regression on constrained cases) is a human call already
made. No further sign-off needed before implementation.

---

## Research Question

> Ticket DEMO-2. NMinimize returns constrained solutions that violate their constraints
> by ~1e-4. Root cause is documented in NMINIMIZE_FEASIBILITY_BUG.md at the repo root.
> I need to fix it, and fix the reason no test caught it (29 tests assert the objective
> within 1e-3; nothing asserts feasibility of the returned point).

## Detailed Findings

### 1. The bug, confirmed

- `nm_eval_pen` (`src/numerical_calculus/findmin_nm_common.c:230-268`) accumulates
  `term = m * m` per violated constraint into `total`.
- `NM_FEAS_EPS = 1.0e-8` (`findmin_internal.h:444`), commented "penalty <= this =>
  feasible (selection)".
- Deb's rule (`findmin_nm_common.c:287-293`) tests `pa <= NM_FEAS_EPS`.
- Because the penalty is squared, the effective tolerance on the real violation is
  `sqrt(1e-8) = 1e-4`. Observed: `9.99991e-05`.
- DE reuses the same predicate as its convergence break (`nm_de.c:159-163`), so it
  stops tightening the constraint the moment the squared penalty crosses the threshold.
  That is why DE — and therefore `Method -> Automatic` — is the affected path.

### 2. Why the obvious fix is unsafe (the finding that changed the plan)

`nm_eval_pen` has exactly **one** call site: `nm_eval`
(`findmin_nm_common.c:273-284`). Every consumer is reachable from that one value, so
the risk surface is closed and enumerable. Classifying every consumer:

**Comparison-only (safe under any monotone rescale):**
| Engine | Use |
|---|---|
| `nm_de.c` | `nm_better` + `NM_FEAS_EPS` thresholds (159, 163) |
| `nm_randomsearch.c` | `nm_better` only |
| `nm_shgo.c` | `nm_better` only |

**Arithmetic (NOT safe — rescaling changes the search, not just the threshold):**
| Engine | Combination | Consumed by |
|---|---|---|
| `nm_neldermead.c` | `nm_phi = f + MU*p` (`findmin_nm_common.c:304-309`) | simplex reflect/expand/contract/shrink + convergence (107,135,140,150,158) |
| `nm_sa.c` | `phi = fx + MU*px` (107,139,157) | Metropolis `exp(-d/T)` (158-166) |
| `nm_basin_hopping.c` | `Ecur = fcur + MU*pcur` (27,60) | Metropolis `exp(-d*beta)` (63-70) |
| `nm_dual_annealing.c` | `Ecur = fcur + MU*pcur` (127,146,165,198,212) | Tsallis acceptance (166-174) |
| `nm_direct.c` | `dc_val = fval + MU*pval` (102-105) | convex-hull slopes (254-256), child ranking (326-328) |

`NM_PENALTY_MU = 1.0e6` (`findmin_internal.h:446`). Changing `p` from quadratic to
linear in violation changes the relative weight of infeasibility against the objective
inside every one of those expressions.

**Additional problem with the norm framing:** `PenaltyFunction` lets a user replace
`m*m` with an arbitrary function (`nm_apply_penalty_fn`,
`findmin_nm_common.c:139-148`; docstring examples at `:997-1005` include `Sqrt` and
`(10 #) &`). With a custom penalty, `total` is not a sum of squares and `sqrt(total)`
has no principled interpretation.

### 3. The local polish is provably independent

`nm_local_polish` does **not** use `nm_eval_pen` during its iterations. It runs
`fm_run_penalty` (`findmin_penalty.c:12-162`) over `fm_eval_penalty`
(`findmin_common.c:1323-1351`), a separate function that always computes `sum d^2` and
takes no `penalty_fn` parameter — deliberately, per the comment at
`findmin_nm_common.c:1003-1005`: "the final local polish keeps the differentiable
squared penalty its analytic gradient assumes". The BFGS/augmented-Lagrangian machinery
and its mu-schedule (1 -> 1e8) are unaffected by anything in this ticket.

`nm_local_polish` calls `nm_eval` once on exit purely to report `(f, pen)` upward.

### 4. Why no test caught it

`tests/test_nminimize.c` asserts entirely through Mathilda source text: `check_true`
(`:54-66`) parses+evaluates a string and `strcmp`s its FullForm against `"True"`.
There is no C-side numeric extraction anywhere in the file, and no shared feasibility
helper anywhere in `tests/` (`test_utils.h` has none).

The structural gap: **every existing feasibility check is a one-sided band.**
- `:142` — `(x^2+y^2 /. Last[...]) <= 9.001` (slack 1e-3 above the bound)
- `:470-472` — `<= 25.001`
- `:1013-1017` — equality residuals `Abs[...] < 1*^-3`, slacks `<= 100.001`
A result that violates its constraint by 1e-4 passes all of them, because they only
catch violations *larger* than their slack. A fix that makes the point more feasible
also passes. So the suite cannot distinguish the bug from the fix.

**Tests sensitive to an `Infinity` flip (must be re-run):**
- `:193` `check_eq(... {x, x^2+1 <= 0} ...) == "Infinity"` (genuinely infeasible)
- `:212-213` `{x, x > 2 && x < 1}` -> `"Infinity"` / `"Indeterminate"`
- `:187-190` feasible-but-displaced must stay finite (`< 1.*^-2`)
- `:1471` `o < 1000` — fails if the result becomes `Infinity`

### 5. Test harness mechanics

- `tests/CMakeLists.txt:1405-1406` registers `nminimize_tests`; **no `add_test`** —
  run the binary directly, matching the `findmin_tests`/`nminimize_tests` convention.
- Build/run only this target:
  `cd tests/build && make nminimize_tests && ./nminimize_tests`
- No pass/fail tally: `ASSERT_STR_EQ` `exit(1)`s on first failure, otherwise the binary
  prints "All NMinimize tests passed." (`:1595`).
- Section 4 of the file is already titled "Feasibility of the returned point" (`:137`)
  and holds only two loose checks — the natural home for the new assertions.

## Code References

- `src/numerical_calculus/findmin_internal.h:444-446` — the three constants
- `src/numerical_calculus/findmin_nm_common.c:230-268` — `nm_eval_pen`, `term = m*m`
- `src/numerical_calculus/findmin_nm_common.c:273-284` — `nm_eval`, the sole caller
- `src/numerical_calculus/findmin_nm_common.c:287-293` — `nm_better`, Deb's rule
- `src/numerical_calculus/findmin_nm_common.c:304-309` — `nm_phi`, the arithmetic site
- `src/numerical_calculus/nm_de.c:159-163` — DE convergence break on the same predicate
- `src/numerical_calculus/nm_driver.c:453,456` — `NM_FEAS_FINAL` decision
- `src/numerical_calculus/findmin_common.c:1323-1351` — `fm_eval_penalty`, independent
- `tests/test_nminimize.c:137-144` — the existing feasibility section
- `tests/test_nminimize.c:193,212-213` — the exact-`Infinity` assertions

## Architecture Insights

The bug is an instance of a general hazard: **a threshold constant whose units are
implicit in the quantity it is compared against.** `NM_FEAS_EPS`'s comment says
"penalty <= this => feasible", which is literally true and practically misleading,
because the reader must independently know that "penalty" means "squared violation".
The same file defines `NM_PENALTY_MU = 1e6` immediately below it — a weight chosen for
the *squared* scale — so the two constants are coupled through an unstated convention.

The reason the safe fix is "square the tolerance at the definition site" rather than
"normalise the quantity" is that the squared scale is load-bearing for five optimizers.
Making the tolerance's units explicit fixes the readability defect that caused the bug
without disturbing the convention that six other call sites depend on.

## Historical Context (from thoughts/)

- `thoughts/shared/research/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarking.md` —
  the benchmarking research that surfaced this bug.
- `thoughts/shared/plans/2026-08-21-DEMO-1-nminimize-nmaximize-benchmarks.md` — the
  DEMO-1 plan; its C1/C2 rows are the speed baseline this ticket must not regress ~3x.
- `NMINIMIZE_FEASIBILITY_BUG.md` — the root-cause note this ticket implements.
- `docs/spec/changelog/2026-08-17.md` — the DEMO-1 changelog entry.

## Related Research

- `benchmarks/89-nminimize-nmaximize/README.md` — constrained rows C1/C2/C3
- `benchmarks/90-nminimize-testbed/README.md` — F1/F2 feasibility rows, currently `0`
