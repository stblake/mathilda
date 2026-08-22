# NMinimize: constrained solutions are returned infeasible (~1e-4)

> **RESOLVED 2026-08-21 (DEMO-2).** Fixed by splitting the one overloaded
> threshold into two named ones — `NM_FEAS_RANK` (loose, search-time ranking)
> and `NM_FEAS_RETURN` (tight, enforced on the return path) — in
> `src/numerical_calculus/findmin_internal.h`. `NMinimize[{x^2+y^2, x+y >= 2},
> {x,y}]` now returns `{2.0, {x -> 1.0, y -> 1.0}}` with a residual of ~4e-12,
> and a result that cannot meet the return bound is reported as
> `{Infinity, x -> Indeterminate}` rather than returned and called feasible.
> Two-sided feasibility assertions added to `tests/test_nminimize.c`; 83/83 pass.
> See "Resolution" at the end of this document for what was NOT fixed.

**Found**: 2026-08-21, while building `benchmarks/89-nminimize-nmaximize`.
**Affects**: `NMinimize` / `NMaximize` with `Method -> Automatic` (i.e. the default)
or an explicit `DifferentialEvolution`. Other methods are unaffected.
**Severity**: correctness. The returned point does not satisfy the stated constraint.

---

## What's wrong

`NMinimize` returns a point that **violates the constraint by ~1e-4**, and reports it
as the solution without any warning.

```
In[1]:= NMinimize[{x^2 + y^2, x + y >= 2}, {x, y}]
Out[1]= {1.9998, {x -> 0.999971, y -> 0.999929}}
```

`x + y = 1.99990`. The constraint is `x + y >= 2`. It is not satisfied. The true optimum
is `2` at `(1, 1)`. Both SciPy (`SLSQP`) and Mathematica return feasible points to
machine precision.

The same happens for equality constraints:

```
In[2]:= NMinimize[{x^2 + y^2, x + y == 2}, {x, y}]
Out[2]= {1.9998, {x -> 0.999971, y -> 0.999929}}
```

## Root cause

The feasibility test compares a **squared** quantity against a threshold intended for an
unsquared one.

`nm_eval_pen` (`src/numerical_calculus/findmin_nm_common.c:230-268`) accumulates the
total constraint penalty as the **sum of squared violations**:

```c
double m = D->gens[k].equality ? fabs(d) : (d > 0.0 ? d : 0.0);
if (m == 0.0) continue;
double term;
if (!D->penalty_fn || !nm_apply_penalty_fn(D->penalty_fn, m, &term))
    term = m * m;                     /* <-- penalty is violation SQUARED */
total += term;
```

That total is then tested against `NM_FEAS_EPS`
(`src/numerical_calculus/findmin_internal.h:444`):

```c
#define NM_FEAS_EPS       1.0e-8   /* penalty <= this => feasible (selection)  */
```

used by Deb's feasibility rule (`findmin_nm_common.c:288-289`):

```c
bool fa_feas = (pa <= NM_FEAS_EPS);
bool fb_feas = (pb <= NM_FEAS_EPS);
```

Because the penalty is `violation^2`, a threshold of `1e-8` on the penalty is a
threshold of **`sqrt(1e-8) = 1e-4` on the actual constraint violation.** The comment
says "penalty <= this => feasible", which is literally true and practically misleading:
the tolerance a *user* experiences is the square root of the constant.

Observed violation: `9.99991e-05`. Predicted: `1e-4`. That is the whole bug.

## Why only DifferentialEvolution

DE uses the same predicate as its **convergence break** (`src/numerical_calculus/nm_de.c:159-163`):

```c
if (*penbest <= NM_FEAS_EPS) {
    ...
    if (cnt >= NP / 2 && (fmax - fmin) <= tol) break;
}
```

So DE stops improving the constraint the moment the *squared* penalty crosses `1e-8` —
i.e. as soon as the violation reaches 1e-4 — and never tightens further. The other
engines reach their answer through a path that keeps polishing, so they land at ~1e-12.

## Reproduction

Measured on `b614d1ed`, macOS arm64, GCC 16.2.0. Violation is `2 - (x + y)`; positive
means the constraint is broken.

| `Method` | objective | `x + y` | violation |
|---|---:|---:|---:|
| **`Automatic` (default)** | 1.9998 | 1.99990 | **9.99991e-05** |
| **`DifferentialEvolution`** | 1.9998 | 1.99990 | **9.99991e-05** |
| `NelderMead` | 2.0 | 2.0 | 1.00092e-06 |
| `RandomSearch` | 2.0 | 2.0 | 7.44116e-12 |
| `SimulatedAnnealing` | 2.0 | 2.0 | 6.03118e-12 |
| `DIRECT` | 2.0 | 2.0 | 3.98659e-12 |
| `FindMinimum` (local solver) | 2.0 | 2.0 | 3.94640e-12 |

Script:

```wolfram
p[lbl_, e_] := Module[{r, s, sx, sy},
  r = e; s = Last[r]; sx = x /. s; sy = y /. s;
  Print[lbl, "\tf=", InputForm[First[r]],
        "\tx+y=", InputForm[sx + sy], "\tviol=", InputForm[2 - (sx + sy)]]];
SetAttributes[p, HoldRest];
p["default        ", NMinimize[{x^2+y^2, x+y >= 2}, {x,y}]];
p["NelderMead     ", NMinimize[{x^2+y^2, x+y >= 2}, {x,y}, Method -> "NelderMead"]];
p["RandomSearch   ", NMinimize[{x^2+y^2, x+y >= 2}, {x,y}, Method -> "RandomSearch"]];
p["DIRECT         ", NMinimize[{x^2+y^2, x+y >= 2}, {x,y}, Method -> "DIRECT"]];
```

### Things that do NOT change the result

Each of these was tested and made no difference — which is what rules out the obvious
explanations:

- `AccuracyGoal -> 12`, `PrecisionGoal -> 12` — not a precision-goal problem.
- `MaxIterations -> 2000` — not a budget problem.
- `Method -> {"DifferentialEvolution", "PostProcess" -> False}` — **identical output**,
  which rules out the local polish (`findmin_penalty.c`) as the cause and locates the
  fault in the global search's constraint handling.

## What I'd check next

1. **Is `NM_FEAS_EPS` meant to be compared against a squared penalty?** If yes, the
   constant should be `1e-16` to give a 1e-8 violation tolerance, and the comment should
   say so. If no, the comparison should be against `sqrt(total)`. Either fix is one line;
   deciding *which* is the actual design intent needs whoever wrote it.
   `NM_FEAS_FINAL = 1e-6` (`findmin_internal.h:445`) has the same ambiguity — it implies
   a 1e-3 violation is "feasible" for the final Infinity test.

2. **Does `PenaltyFunction` change the exponent and silently rescale the tolerance?**
   `nm_apply_penalty_fn` replaces `m*m` with an arbitrary user function, so a user
   supplying `Abs` rather than a square gets a *different effective feasibility
   tolerance* (1e-8 instead of 1e-4) without being told. That is the same bug with a
   user-visible knob attached.

3. **Why is `NelderMead` at 1e-6 rather than 1e-12?** It is feasible, but it sits six
   orders off the other three. Possibly the same threshold reached through a different
   route; worth confirming it is not a second instance.

4. **Should an infeasible return warn at all?** Right now the point is returned with no
   indication. `nm_driver.c:456` computes `feasible` against `NM_FEAS_FINAL` for the
   `Infinity` path, so the machinery to detect this exists and is simply not surfaced to
   the user in the near-miss case.

5. **Regression test.** `tests/test_nminimize.c` has 29 tests, all passing, including
   constrained cases — they assert the *objective* within `1e-3`, which `1.9998` passes.
   A test asserting **feasibility of the returned point** would have caught this and
   currently does not exist. That is the gap that let it ship.

## Related

A second, independent robustness finding from the same benchmark run, recorded in
`benchmarks/90-nminimize-testbed/`: on the standard hard global-optimization corpus,
Mathilda solves **2 of 7** landscapes where SciPy solves **4 of 7** (misses Schwefel 5-D,
Griewank 5-D, Rastrigin 10-D; wins drop-wave, which SciPy misses). That is a search-budget
tradeoff rather than a defect — it is why Mathilda returns in 0.5 ms where SciPy takes
96 ms — but it is worth knowing that the default budget is tuned toward speed.

Also from that run, and separately actionable: naming `Method` explicitly cuts DE's
generation budget from `150n` to `100` (`nm_de.c:77-86`), which makes 5 of 6 seeds fail
Rastrigin 5-D that the default solves. See `benchmarks/89-nminimize-nmaximize/README.md`.

---

## Resolution (2026-08-21, DEMO-2)

### What was done

**Two named thresholds, because there were two jobs.** The single `NM_FEAS_EPS`
was steering the search *and* deciding what got returned, and it was compared
against a squared quantity so it meant 1e-4 rather than the 1e-8 it read as.

- **`NM_FEAS_RANK`** (violation 1e-4, squared to 1e-8 — bit-identical to the old
  constant): Deb's rule during the search. Loose on purpose. Deb's rule only
  consults the objective when *both* points count as feasible, so a threshold
  tighter than the achievable residual makes the rule degenerate to pure
  violation-minimisation. Measured: setting it to 1e-16 sent the 15-dimensional
  `test_minimax_chebyshev` from 0.125116 to 1.85479.
- **`NM_FEAS_RETURN`** (violation 1e-5, squared): the return path, and enforced.
  Applied at post-polish selection in `nm_de.c`, the driver's cross-attempt best
  and region-expansion break, and the final feasible-vs-`Infinity` decision.

Both are stated as violations and squared at the definition site, so the
comparison's units are visible in the source instead of implied.

**A last-chance continuous refinement** was added to `nm_local_polish`'s
mixed-integer branch: when the answer would otherwise be rejected and continuous
coordinates exist, pin the integers and refine, adopting the result only on a
Deb improvement.

### Verified

| Check | Result |
|---|---|
| Original repro | `{2.0, {x -> 1.0, y -> 1.0}}`, residual ~4e-12 |
| All 8 methods feasible | yes, 7 at ~4e-12; NelderMead at 1.00092e-06 |
| `test_nminimize.c` | 83/83 pass |
| Genuinely infeasible still `Infinity` | yes |
| Feasible-but-displaced still finite | yes |
| Speed (exp-89 C1/C2) | 0.259 / 0.210 ms vs 0.255 / 0.218 baseline — no regression |
| `make check-c99` | clean |

### What was NOT fixed — follow-up

1. **CORRECTION + OPEN REGRESSION — `test_fixed_charge_flow`.** The commit
   message for this fix, and an earlier version of this section, state that the
   solver's best answer on that instance has a flow residual of **20.0**. **That
   number is wrong.** It was measured while `nm_int_descent`'s move acceptance
   was still (incorrectly) using the tight predicate — a change that was reverted
   before shipping. Re-measured against the shipped code, with the return gate
   opened:

   | build | objective | max flow residual |
   |---|---:|---:|
   | pre-fix (`ea0f8c3c`) | 877.38 | **6.59657e-05** |
   | shipped (`d4eb10fa`) | 476.948 | **0.23574** |

   So this is **not** merely a strict gate rejecting a problem the solver never
   solved. The pre-fix code found a near-feasible point (6.6e-5, comfortably
   inside the test's own 1e-2 tolerance); the shipped code finds one three orders
   of magnitude worse. **The fix degraded feasibility on this problem**, and the
   test now asserting `Infinity` documents that regression rather than an
   inherent limitation.

   Suspected mechanism, not yet confirmed: when no candidate meets
   `NM_FEAS_RETURN`, `nm_better_return` ranks purely by violation
   (`pa < pb`), discarding objective guidance entirely — the same degeneration
   that a tight *ranking* threshold caused in `test_minimax_chebyshev`. Selection
   is still ranking, so it may need the loose predicate too, with the tight bound
   applied only at the final gate in `nm_driver.c`. That was not established
   before shipping.

   **Recommended next step:** either confirm and fix that, or revert
   `nm_better_return` at the selection sites (`nm_de.c`, `nm_driver.c:443,449`)
   and re-check whether the original 1e-4 bug returns.

   The underlying accuracy gap is still real — the MINLP path reaches ~1e-3
   per-constraint against ~4e-12 for continuous problems, which is why
   `NM_FEAS_RETURN_VIOL` is 1e-5 rather than tighter.
2. **NelderMead's 1.00092e-06.** Six orders looser than the other seven methods
   and unchanged by this ticket. Inside the guarantee, so a legitimate result,
   but it looks like a second instance of a threshold set against the wrong
   quantity.
3. **`PenaltyFunction` tolerance semantics.** A user-supplied non-squared penalty
   still changes the effective tolerance, since both constants assume `m*m`.
4. **No warning on an infeasible return.** The result is `Infinity` with no
   message explaining that a feasible point was not found.

---

## The test-suite failure class — a second confirmed instance

`test_fixed_charge_flow` is not a one-off. It is the **same failure class as the
bug this document is about**, and finding two instances in one suite is the
reason to audit the rest.

**The pattern:** a test asserts the *objective* tightly and *feasibility* loosely
(or not at all), so a point that is optimal-looking but constraint-violating
passes. The objective is the thing everyone thinks to assert; feasibility is the
thing that actually makes the answer meaningful.

| instance | asserted | missed |
|---|---|---|
| the original bug | `Abs[f - 2] < 1e-3` across 29 tests | that `x + y` was 1.99990, not 2 |
| `test_fixed_charge_flow` | `o < 1000`, residual `< 1e-2` | that a 1e-2 slack on a *flow-conservation equality* accepts a materially infeasible network |

Both are the same mistake at different scales: the tolerance was chosen to make a
known-imperfect result pass, rather than to express what feasibility means for
the problem.

**What to audit in `tests/test_nminimize.c`.** Two mechanical searches find the
candidates:

- **13 one-sided slack assertions** of the shape `<= bound + epsilon` — `9.001`,
  `25.001`, `3.0001`, `100.0001`, `15.001`, `100.001`/`200.001`, `0.027001`,
  `1.0001`, `5.001`. Each is a ceiling only: it catches overshoot and is blind to
  a point that never reaches the constraint boundary, which is the direction the
  original bug went. A two-sided assertion, or an explicit residual bound, tests
  something the one-sided form cannot.
- **53 assertions at 1e-2 or 1e-3 tolerance.** Not all are wrong — a stochastic
  global optimum genuinely needs a loose objective band. But a loose band on a
  *constraint residual* is different in kind from a loose band on an objective,
  and the two are not currently distinguished anywhere in the file.

**The rule worth adopting:** an objective tolerance may be loose, because the
search is stochastic. A feasibility tolerance should be tight, because
feasibility is not an approximation — it is the difference between an answer and
a wrong answer. Where a feasibility tolerance must be loose, that is a finding
about the solver to record, not a number to quietly choose.

### AUDITED — DEMO-3, 2026-08-22

That audit was done. It was worse than two instances.

Measured by mutating `nm_build_result` to return a deliberately wrong answer:

| probe | before | after |
|---|---:|---:|
| returned point 10% wrong, objective right | 26 / 83 caught | **32 / 83** |
| returned point 1% wrong | 24 / 83 | **30 / 83** |
| one integer coordinate flipped | 7 | **8** |
| objective 10% wrong, point right | 60 / 83 | unchanged |

**31 constrained tests never checked their constraints at all.** The suite caught a wrong
*number* more than twice as well as a wrong *answer*. Git archaeology across 37 commits
found no tolerance ever widened to rescue a failing test — the looseness was authored in,
under a stated policy that governs objectives and is silent on feasibility. The defect was
an absent category, not carelessness.

Fixed by writing that missing category into the file header as an explicit feasibility
policy, adding feasibility assertions to the eleven tests with real constraints, and
documenting the eighteen remaining omissions as deliberate (sixteen box-only where the
solver clamps and an assertion could not fail, one that turned out to be unconstrained and
was misfiled, one deferred with its reason). Full audit:
`thoughts/shared/tickets/DEMO-3/research.md`.
