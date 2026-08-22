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

1. **Mixed-integer feasibility accuracy.** The MINLP path reaches ~1e-3
   per-constraint feasibility against ~4e-12 for continuous problems. This is why
   `NM_FEAS_RETURN_VIOL` is 1e-5 rather than tighter. `test_fixed_charge_flow`
   (50-var coupled flow-conservation MINLP) now asserts `Infinity`, because the
   solver's best answer there has a flow residual of **20.0** — the entire demand
   at node 1. It was previously reported as a solution under a 1e-2 tolerance.
   Tightening the return guarantee further requires fixing this path first.
2. **NelderMead's 1.00092e-06.** Six orders looser than the other seven methods
   and unchanged by this ticket. Inside the guarantee, so a legitimate result,
   but it looks like a second instance of a threshold set against the wrong
   quantity.
3. **`PenaltyFunction` tolerance semantics.** A user-supplied non-squared penalty
   still changes the effective tolerance, since both constants assume `m*m`.
4. **No warning on an infeasible return.** The result is `Infinity` with no
   message explaining that a feasible point was not found.
