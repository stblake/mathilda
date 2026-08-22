/* predict.h -- Predict, LinearModelFit, and the trained-model representation.
 *
 * THE TRAINED-MODEL REPRESENTATION, designed once here and shared.
 *
 * Three builtins need to hand back a fitted object that survives being stored in a
 * variable and re-applied to new data: Predict (a PredictorFunction), DimensionReduce
 * (a DimensionReducerFunction) and, later, Classify (a ClassifierFunction). This is
 * the one place that decision gets made.
 *
 * The decision: a plain EXPR_FUNCTION whose head is the symbol `PredictorFunction`
 * and whose arguments are ordinary Exprs.
 *
 * Two alternatives were considered and rejected, and the reasons are worth keeping
 * because they are not obvious:
 *
 *   - A NEW NODE TYPE, following EXPR_COMPILED. That is the codebase's existing
 *     precedent for an opaque callable carrying binary state, and it was the expected
 *     answer -- but EXPR_COMPILED exists because compiled code is a VM program:
 *     bytecode, a register file, a reference count. A fitted linear model is a short
 *     vector of numbers. A new node type means new cases in expr_copy, expr_free,
 *     expr_eq, expr_hash, expr_compare and print.c -- touching the core of the system
 *     -- to store something the existing Expr types already hold perfectly. Machine
 *     precision is not an argument for it either: a packed List already holds a dense
 *     double buffer and is still an ordinary List.
 *
 *   - AN ASSOCIATION payload, `PredictorFunction[<|"Method" -> ..., ...|>]`, which is
 *     the most Wolfram-ish shape. Rejected only because it buys nothing here that
 *     positional arguments do not, while adding a dependency on the Association
 *     machinery to every read. Named PROPERTIES are still exposed, through the
 *     application path (`p["Coefficients"]`), which is where a user reaches for them.
 *
 * Application uses the evaluator's existing composite-head chain: eval.c already
 * dispatches `Function[...][args]` and `Association[...][key]` by matching on the
 * head's own head, and a fitted model is one more branch in the same chain. So no new
 * evaluation concept is introduced.
 *
 * Layout (positional arguments of the PredictorFunction head):
 *   0  method name, a String
 *   1  coefficient vector, a List of Reals: intercept first, then one per feature
 *   2  feature count, an Integer -- carried so application can reject a wrongly
 *      shaped input instead of reading past the end of a short vector
 */
#ifndef ML_PREDICT_H
#define ML_PREDICT_H

#include <stddef.h>
#include <stdbool.h>
#include "expr.h"

/* Ordinary least squares with an intercept.
 *
 * `x` is row-major n x dim, `y` has n entries; `coef` receives dim + 1 values with the
 * intercept first.
 *
 * Returns false when the system is singular -- collinear features, or fewer
 * observations than parameters. Declining is deliberate: a pseudo-inverse would return
 * one of infinitely many answers and look like a fit. */
bool ml_ols(const double* x, const double* y, size_t n, size_t dim, double* coef);

/* Apply a fitted model to arguments. Returns NULL when `head` is not a model this
 * module owns, so the evaluator can fall through to its other composite-head cases. */
/* Cheap test for "is this head a fitted model", so eval.c's dispatch chain does not
 * have to know the model head names. Separate from ml_model_apply because the chain
 * needs to decide which BRANCH to take before committing to it. */
bool ml_model_apply_probe(Expr* head);

Expr* ml_model_apply(Expr* head, Expr** args, size_t argc);

void ml_predict_init(void);

#endif /* ML_PREDICT_H */
