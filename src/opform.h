/*
 * opform.h -- generic operator (curried) forms: h[o1, ..., ok][x].
 *
 * Many Wolfram heads have an "operator form": a call missing its data
 * argument, h[o1, ..., ok], stays inert, and applying it to one expression x
 * supplies that argument:
 *
 *     Select[crit][x]          == Select[x, crit]
 *     Lookup[k][a]             == Lookup[a, k]
 *     Map[f][x]                == Map[f, x]
 *     Insert[e, n][x]          == Insert[x, e, n]
 *
 * Rather than special-casing each head, a head registers ONE row in a table:
 * how many operator arguments it accepts (min_ops..max_ops) and the 0-based
 * position at which x is inserted into the full call.  The evaluator consults
 * the table when it reaches a call whose head is itself a call h[...] of a
 * registered symbol h, with exactly one argument and no other rule for the
 * composite head (see the composite-head chain in eval.c).
 *
 * The rewritten call h[..., x, ...] is evaluated immediately.  If it comes
 * back structurally unchanged (h declined those arguments) the curried form
 * is left unevaluated, which is what Mathematica shows:
 * CountsBy[f][<|...|>] stays CountsBy[f][<|...|>], not CountsBy[<|...|>, f].
 *
 * Only registered (builtin, Protected) heads are affected, so user-defined
 * f[x][y] expressions are untouched.  Heads whose builtin already rewrites
 * its one-argument call into a Function (SortBy, Cases, MapAt, ...) keep
 * working as before; they never reach this table.
 */
#ifndef MATHILDA_OPFORM_H
#define MATHILDA_OPFORM_H

#include "expr.h"
#include <stdbool.h>
#include <stddef.h>

/* Register `head` (any spelling; it is interned) as having operator forms
 * head[o1..ok][x] for min_ops <= k <= max_ops.  x is inserted at 0-based
 * position insert_at of the rewritten call (0 = first argument, 1 = after the
 * first operator argument, ...); insert_at must be <= min_ops.
 * Re-registering a head replaces its row.  Call from a module's *_init(). */
void opform_register(const char* head, size_t min_ops, size_t max_ops, size_t insert_at);

/* True when `head` (the composite head h[o1..ok]) is a registered operator
 * form with an admissible number of operator arguments. */
bool opform_matches(const Expr* head);

/* Apply the operator form `head` = h[o1..ok] to `x`.  Returns the evaluated
 * result (caller owns it), or NULL to leave head[x] unevaluated.  Neither
 * argument is consumed. */
Expr* opform_apply(const Expr* head, const Expr* x);

/* Registers the built-in operator-form table.  Called from core_init(). */
void opform_init(void);

#endif /* MATHILDA_OPFORM_H */
