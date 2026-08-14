---
schema_version: 2
title: Machine learning family 5 — Classify
slug: classify
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
  stamped: 2026-08-14 03:25
  scope:
    - "src/ml/**"
    - "src/eval.c"
    - "src/core.c"
    - "src/sym_names.c"
    - "src/sym_names.h"
    - "src/print.c"
    - "makefile"
    - "tests/**"
    - "docs/spec/**"
  success_criteria:
    - "AC-1 — A k=1 nearest-neighbour classifier reproduces EVERY training label exactly."
    - "AC-2 — Class probabilities sum to 1, including on a non-degenerate split."
    - "AC-3 — A class may be a string, a symbol or a number; structural comparison distinguishes them."
    - "AC-4 — Class ordering is deterministic and the prediction does not depend on it."
    - "AC-5 — Malformed input declines rather than guessing."
---

# Machine learning family 5 — Classify

## Feature Definition

- One-line goal: `Implement Classify, and the categorical label encoding every earlier family avoided.`
- Problem: `None of Classify, ClassifierFunction, NaiveBayes, LogisticRegression, DecisionTree or RandomForest existed. The real work is not the learners: families 1-4 all took numeric features and numeric responses, and Classify is where an arbitrary CLASS arrives. Something has to turn arbitrary expressions into dense indices and back.`
- Requested by: `direct user request, the last of the five families`

## Current State Study

- The model representation (src/ml/predict.h) already supports a new head in three steps;
  ClassifierFunction is its FOURTH user and needed no change to it.
- The k-NN machinery in src/ml/predict.c is the piece to reuse — a classifier votes where
  a predictor averages.
- No categorical encoding existed anywhere in src/ml.

## Scope

### In Scope

- The label vocabulary, in its own module (done, iteration 20)
- `Classify` with `"NearestNeighbors"` (done, iteration 20)
- `"NaiveBayes"` (Gaussian per class — the Multinormal density makes it nearly free)
- `"LogisticRegression"`

### Out Of Scope

- `DecisionTree` / `RandomForest` — they reuse nothing built so far. To be judged honestly
  when the earlier methods are in, and recorded as deferred if they do not fit.
- `ContingencyTable`, inherited from family 4 as a follow-up now that the vocabulary exists.

## Success Criteria

- `AC-1` — k=1 reproduces every training label exactly
- `AC-2` — probabilities sum to 1 on a non-degenerate split
- `AC-3` — strings, symbols and numbers all work as classes
- `AC-4` — deterministic ordering; prediction independent of it
- `AC-5` — malformed input declines

## Tests

- `tests/test_ml_classify.c`, registered with BOTH `add_executable` and `add_test`.
  `eval_tests` and `print_tests` are in the regression set because the head was added to
  the evaluator chain and the printer.

## Task List

- [x] `Label vocabulary in src/ml/encode.{c,h}` | `independent` | `done (iteration 20)`
- [x] `ClassifierFunction on the existing model representation + Method -> NearestNeighbors` | `depends on: vocabulary` | `done (iteration 20)`
- [ ] `Method -> "NaiveBayes"` | `depends on: vocabulary` | `pending`
- [ ] `Method -> "LogisticRegression"` | `depends on: vocabulary` | `pending`
- [ ] `DecisionTree / RandomForest, or a recorded decision not to` | `independent` | `pending`
- [ ] `ContingencyTable (inherited from family 4)` | `depends on: vocabulary` | `pending`

## Checkpoints

- [x] start | completed: `2026-08-14 03:25`
- [x] spec / plan created | completed: `2026-08-14 03:25`
- [x] implementation started | completed: `2026-08-14 03:25`
- [ ] implementation complete | completed: `pending`
- [ ] critic pass | completed: `pending`
- [ ] risk-register reviewed | completed: `pending`
- [ ] feature validated | completed: `pending`
- [ ] PR created | completed: `pending`
- [ ] closeout complete | completed: `pending`

## Decisions

- `The label vocabulary gets its own module (src/ml/encode.{c,h}) rather than living inside classify.c, designed once for the same reason the model representation was: a ContingencyTable, a categorical FEATURE encoder and any second classifier all need exactly this, and three copies would be three chances to order the classes differently.`
- `Order is FIRST APPEARANCE, not sorted. expr_compare could give a total order, but first appearance is stable under adding classes later, matches how fc_assign_from_uf numbers clusters, and reads naturally beside the data. What matters is only determinism.`
- `Class comparison is expr_eq, so "a" and a are two classes. That is the same structural distinction the pattern matcher makes, and pretending otherwise would be a silent surprise.`
- `Only the rule form {features -> class} is accepted. A matrix with the class in the last column cannot work because a class need not be numeric -- na_load_matrix would refuse the whole input. Predict accepts a matrix precisely because its response IS numeric.`
- `k defaults to 1 for a CLASSIFIER but 3 for the k-NN PREDICTOR. A regression averages, so smoothing helps; a classifier votes, and k=1 gives the exactly-checkable training-set property.`
- `The payload stores class INDICES with the features and the vocabulary separately, so the numeric part stays numeric and a class is named in exactly one place.`

## Risks And Unknowns

- `A classifier is easy to test badly -- an accuracy figure hides which rows are wrong. Mitigated by asserting the two ABSOLUTE properties instead: exact training-set reproduction at k=1, and probabilities summing to 1.`
- `LogisticRegression needs an iterative solver (IRLS or gradient descent); ml_ols is only the linear piece. That is the one remaining method with real new numerics.`

## Risk Register

| ID | Category | Likelihood | Impact | Mitigation | Residual | Owner |
|---|---|---|---|---|---|---|
| RISK-1 | Weak verification | med | med | An accuracy percentage would hide which rows are misclassified. Two absolute checks instead: every training label reproduced exactly at k=1, and probabilities summing to 1 on a non-degenerate split. | Closed for this method. | Michael Sollami |
| RISK-2 | Encoding divergence | med | high | Three copies of a label vocabulary could order classes differently, making two models disagree about what class 0 means. One module, one ordering rule, tested. | Closed. | Michael Sollami |

## Test Results

- Command: `cd tests/build && ./ml_classify_tests` plus the ten other suites
- Outcome: `all pass (6 test groups); check-c99 / check-packed-aware / check-array-exactness OK`

## Proofs Of Completion

- Completion timestamp: `pending`

## Tech Debt Review

- Potential tech debt introduced:
  - `ml_classifier_apply rebuilds the label vocabulary from the payload on EVERY application, which is O(classes) allocation per call. Fine for interactive use; a batch classifier over many points would want it hoisted.`
  - `Only Euclidean distance. A DistanceFunction option here would be the second real consumer of find_clusters.c\'s private metric layer -- which is exactly the trigger recorded for that extraction.`
- Existing tech debt noticed:
  - `ContingencyTable still unimplemented, now unblocked by the vocabulary.`
- Mitigations taken:
  - `Two absolute assertions rather than an accuracy figure.`
- Follow-up needed:
  - `NaiveBayes, LogisticRegression, a decision on DecisionTree/RandomForest, ContingencyTable.`

## Activity Log

- `2026-08-14 03:25` Created, and iteration 20 landed: the label vocabulary (src/ml/encode.{c,h}) plus Classify with a k-NN majority vote. ClassifierFunction is the FOURTH head on the model representation from iteration 12 and needed no change to it -- adding it was the documented three steps plus the printer entry. The two verifications are ABSOLUTE rather than comparative: k=1 reproduces all six training labels exactly (each point is its own nearest neighbour, so there is no "close enough"), and probabilities sum to 1 including at k=5 where the split is a genuine 0.6/0.4 rather than degenerate. Strings, symbols and integers all work as classes, and structural comparison keeps "a" and a distinct.

## Reflection

- What went well: `pending`
- What went wrong: `pending`
- Gaps to close: `pending`
- Skill / AGENTS.md updates to propose: `pending`

## Follow-Up

- `This is the last family; when it closes, the whole arc should be summarised for the user.`

## Team Addendum

- `none`
