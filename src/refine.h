/*
 * refine.h
 *
 * Refine[expr, assum] -- give the form of `expr` that would be obtained if the
 * symbols in it were replaced by explicit values satisfying the assumptions
 * `assum`. Refine[expr] uses the ambient $Assumptions (as set by any enclosing
 * Assuming[] construct).
 *
 * Refine is a thin orchestrator over the shared assumption engine (simp.h /
 * simp_internal.h): it delegates the assumption-driven rewrites to
 * apply_assumption_rules, folds domain/predicate statements via element_decide /
 * zero_test_decide_assuming, and proves the harder inequality/positivity facts
 * through the Reduce/CAD engine (reduce.h). It is one of the transformations
 * tried by Simplify.
 */
#ifndef MATHILDA_REFINE_H
#define MATHILDA_REFINE_H

#include "expr.h"

/* Refine[expr] / Refine[expr, assum] with options Assumptions and
 * TimeConstraint. Follows the builtin ownership contract (returns a new tree or
 * steals from res; never frees res). */
Expr* builtin_refine(Expr* res);

/* Register Refine (Protected), its options and docstring. Called from
 * core_init(). */
void refine_init(void);

#endif /* MATHILDA_REFINE_H */
