---
schema_version: 2
title: Machine learning family 3 — Predict, and the trained-model representation
slug: predict
status: implementation
source: direct-user-request
owner: Michael Sollami
issue: pending
pull_request: pending
started: 2026-08-14
last_updated: 2026-08-14
blocked_by: none
goal_lock:
  status: active
  stamped: 2026-08-14 01:45
  scope:
    - "src/ml/**"
    - "src/list/**"
    - "src/linalg/**"
    - "src/eval.c"
    - "src/core.c"
    - "src/sym_names.c"
    - "makefile"
    - "tests/**"
    - "docs/spec/**"
  success_criteria:
    - "AC-1 — A trained model survives assignment and re-application."
    - "AC-2 — Predict recovers an exactly-linear relationship to machine precision."
    - "AC-3 — Predict agrees with the existing independent Fit implementation."
    - "AC-4 — A singular (collinear) system declines rather than returning one of infinitely many answers."
    - "AC-5 — The representation is designed ONCE and is reusable by DimensionReducerFunction and ClassifierFunction."
---

# Machine learning family 3 — Predict

## Feature Definition

- One-line goal: `Implement Predict and LinearModelFit, and design the trained-model representation the remaining families need.`
- Problem: `None of Predict, PredictorFunction, LinearModelFit or NearestNeighbors existed. More importantly, two earlier families deferred the same blocker to here: DimensionReduce cannot return a reusable reducer and Classify cannot return a classifier until there is a way to represent a FITTED MODEL as an Expr that survives storage and re-application.`
- Requested by: `direct user request, continuing the ML families`

## Current State Study

- `EXPR_COMPILED` exists as an opaque reference-counted node type — the codebase's
  precedent for a callable carrying binary state, and the expected answer.
- `eval.c` dispatches COMPOSITE heads through an explicit chain: `Function[...][args]`
  at step 7, `Association[...][key]` at 7a. That is the existing idiom for callable
  objects, and a fitted model is one more branch.
- `Fit` (src/fit.c, registered in info.c) already solves least squares independently —
  usable as a cross-check rather than a thing to reimplement.

## Scope

### In Scope

- The trained-model representation (done, iteration 12)
- `Predict` with `Method -> "LinearRegression"`, `LinearModelFit` (done, iteration 12)
- `Predict` with `Method -> "NearestNeighbors"`
- Retrofitting `DimensionReducerFunction` onto the same representation

### Out Of Scope

- Wolfram's `FittedModel` regression diagnostics (RSquared, standard errors, ANOVA) —
  deliberately not approximated
- Automatic method selection, cross-validation, `PerformanceGoal`

## Success Criteria

- `AC-1` — a model survives assignment and re-application
- `AC-2` — exact recovery of an exactly-linear relationship
- `AC-3` — agreement with `Fit`
- `AC-4` — a singular system declines
- `AC-5` — the representation is reusable by the other two model kinds

## Tests

- `tests/test_ml_predict.c`, registered with BOTH `add_executable` and `add_test`.
  `eval_tests` is part of the regression set for this family because `eval.c` changed.

## Task List

- [x] `Trained-model representation, designed once` | `independent` | `done (iteration 12)`
- [x] `Predict + LinearRegression + LinearModelFit + application via the evaluator chain` | `depends on: representation` | `done (iteration 12)`
- [x] `Predict with Method -> NearestNeighbors` | `depends on: representation` | `done (iteration 13), with a NeighborsNumber sub-option`
- [x] `Retrofit DimensionReducerFunction onto the representation (family 2 open row)` | `depends on: representation` | `done (iteration 14) -- needed NO new evaluation machinery, which was the payoff for designing the representation once`
- [x] `Abbreviated printing for fitted models` | `independent` | `done (iteration 15) -- followed InterpolatingFunction's existing elision in print.c; FAMILY 3 NOW CLOSED`

## Checkpoints

- [x] start | completed: `2026-08-14 01:45`
- [x] spec / plan created | completed: `2026-08-14 01:45`
- [x] implementation started | completed: `2026-08-14 01:45`
- [ ] implementation complete | completed: `pending`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [ ] feature validated | completed: `pending`
- [ ] PR created | completed: `pending`
- [ ] closeout complete | completed: `pending`

## Decisions

- `THE REPRESENTATION: a fitted model is a plain EXPR_FUNCTION whose head is a symbol (PredictorFunction) and whose arguments are ordinary Exprs. NOT a new node type following EXPR_COMPILED, which was the expected answer -- EXPR_COMPILED exists because compiled code is a VM program (bytecode, register file, refcount), whereas a fitted linear model is a short vector of numbers. A new node type would mean new cases in expr_copy, expr_free, expr_eq, expr_hash, expr_compare and print.c -- touching the core -- to store what existing Expr types already hold. Machine precision is not an argument either: a packed List already holds a dense double buffer and is still a List.`
- `NOT an Association payload, though that is the most Wolfram-ish shape: it buys nothing positional arguments do not, while adding an Association dependency to every read. Named properties are still exposed through the application path (p["Coefficients"]), which is where a user reaches for them.`
- `Application reuses eval.c\'s existing composite-head chain rather than introducing a new evaluation concept. A probe function (ml_model_apply_probe) keeps eval.c from having to know the model head names.`
- `ml_ols solves the NORMAL EQUATIONS by Gaussian elimination rather than calling LAPACK dgels. dgels is numerically better (a QR avoids squaring the condition number), but A\'A is (dim+1) square -- a handful of rows -- and the normal equations make the SINGULAR case directly detectable as a zero pivot, which matters because declining on a singular fit is deliberate behaviour here.`

## Risks And Unknowns

- `Squaring the condition number is a real cost of the normal equations. Acceptable at a small feature count; a large one would be a reason to switch to dgels, and the call site is one function.`
- `A fitted model currently prints its full parameter list. Wolfram abbreviates. Harmless for small models, ugly for a large one.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Wrong abstraction, expensive to change later | med | high | The representation is used by three families, so a bad choice compounds. Mitigated by choosing the SIMPLEST thing that meets the constraints and by writing down why the two richer alternatives were rejected, so a future change is an informed one. | Accepted. | Michael Sollami |
| RISK-2 | Numerical conditioning | low | med | Normal equations square the condition number. Bounded by the small feature count; the singular case is detected and declined rather than silently fitted. | Accepted at this size. | Michael Sollami |

## Test Results

- Command: `cd tests/build && ./ml_predict_tests` plus findclusters x3, list, ml_pca, fit, linalg, eval
- Outcome: `all pass; make check-c99 / check-packed-aware / check-array-exactness OK`

## Proofs Of Completion

- Completion timestamp: `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `RESOLVED (iteration 15): models print method + <>, following InterpolatingFunction. FullForm still reveals everything.`
  - `RESOLVED (iteration 13): promoted to src/ml/mlutil.{c,h} -- ml_list_of_reals, ml_list_matrix, ml_read_data, plus ml_sqdist.`
- Existing tech debt noticed:
  - `na_build_matrix returning a visible NDArray remains a trap for new user-facing builtins.`
  - `find_clusters.c still exports one symbol; its machine distance layer will be wanted by NearestNeighbors.`
- Mitigations taken:
  - `Test rows assert storage/re-application and agreement with Fit.`
- Follow-up needed:
  - `NearestNeighbors; retrofit DimensionReducerFunction; abbreviated printing; promote the list builders.`

## Activity Log

- `2026-08-14 02:30` Iteration 15: abbreviated printing, and family 3 is CLOSED. print.c already had the pattern -- InterpolatingFunction prints its domain then `, <>` and says in its own comment that FullForm reveals the rest -- so this was following an in-tree precedent rather than inventing a scheme, and it needed no new mechanism. Added SYM_PredictorFunction / SYM_DimensionReducerFunction to sym_names.c per CLAUDE.md rather than strcmp-ing in the printer. Asserted BOTH halves: the abbreviated form, and that FullForm plus application plus property access are all untouched -- a change that reached into the object rather than its rendering would show up in the second half. One test-expectation error of mine: this harness renders strings WITH quotes where the REPL's Print does not.
- `2026-08-14 02:20` Iteration 14: DimensionReducerFunction, via a new DimensionReduction[data, k] builtin (Wolfram's split: DimensionReduce returns reduced data, DimensionReduction returns a reducer). This closes family 2's last open row. The retrofit needed NO new node type, NO new evaluation concept and no change to eval.c beyond the probe recognising one more head -- which is the payoff for having designed the representation once in iteration 12 rather than per-family. The one real hazard was centring: new rows must be centred on the TRAINING means, and the plausible bug (centring each incoming batch on its own mean) gives all zeros for a single point, which is indistinguishable from correct on data near the origin. Training data with means {104, 210} makes it detectable, and the test asserts both halves -- a point at the training mean projects to 0, a point away from it projects far from 0. Also note the loop was PAUSED for a while here: the permission classifier was unavailable, which blocked ScheduleWakeup and then Bash; no partial writes occurred because the rejected call never executed.
- `2026-08-14 01:55` Iteration 13: Predict with NearestNeighbors, plus a NeighborsNumber sub-option, plus the promotion of the duplicated helpers to src/ml/mlutil.{c,h} (the third copy had arrived, which was the stated trigger). DECLINED to extract find_clusters.c's machine metric layer: it takes that file's FcData and honours a DistanceFunction, so extracting it means introducing a metric enum in src/ml and rewriting a 2600-line file with 22 pinned answers -- real risk to a finished family, bought for a four-line squared-Euclidean. Recorded the actual trigger instead: the moment Predict or Classify gains a DistanceFunction option, that layer has its second REAL consumer. Cross-checked k=1 against the existing Nearest on both sides of a midpoint. One of my tests was wrong for the exact reason I had written down: the k-NN-differs-from-linear test used data that was exactly y=10x, so both methods agreed and it proved nothing -- fixed with y=x^2.
- `2026-08-14 01:45` Created, and iteration 12 landed: the trained-model representation plus Predict, LinearModelFit and application through eval.c\'s composite-head chain. Rejected a new EXPR node type (the expected answer) with reasons recorded. Cross-checked against the existing Fit, which is an independent least-squares implementation. One test of mine was wrong rather than the code: Head of an unevaluated composite application is the whole head PredictorFunction[...], not the symbol, so the assertion became NumberQ.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- `Family 4 LearnDistribution (ml_gmm_fit/ml_gmm_logpdf already reusable), family 5 Classify.`

## Team Addendum

- `none`
