/* radicals.h
 *
 * ToRadicals: convert held Root[Function[poly], k] objects into
 * closed-form radical expressions.
 *
 * Closed forms are available unconditionally for polynomials of degree
 * <= 4 (linear, quadratic, Cardano cubic, Ferrari quartic).  Binomial
 * polynomials a x^n + b of any degree are also handled.  Other Root
 * objects of degree >= 5 are left untouched.
 *
 * ToRadicals walks its argument recursively, so List, Equal, Less,
 * And, Or, ... arguments thread automatically: every Root[..] node
 * anywhere in the tree is attempted independently.
 */
#ifndef RADICALS_H
#define RADICALS_H

#include "expr.h"

/* Built-in entry point for `ToRadicals[expr]`.  Takes ownership of the
 * outer res node per the standard contract; returns a freshly owned
 * Expr* on success, or NULL when the input is malformed (wrong argc),
 * in which case the evaluator preserves the unevaluated form. */
Expr* builtin_to_radicals(Expr* res);

/* Register the ToRadicals builtin, attach Protected, and install the
 * docstring placeholder (the full docstring lives in info.c so the
 * existing `?ToRadicals` lookup pipeline finds it). */
void radicals_init(void);

#endif /* RADICALS_H */
