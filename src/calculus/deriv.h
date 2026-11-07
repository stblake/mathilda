/*
 * deriv.h -- Symbolic differentiation (D, Dt, Derivative)
 *
 * Provides a native C implementation of Mathematica-style differentiation,
 * replacing the earlier rule-based src/internal/deriv.m bootstrap. All
 * elementary differentiation identities (sums, products, powers, exp/log,
 * trig, hyperbolic, inverse trig/hyperbolic, chain rule on unknown
 * functions) are encoded directly in C so that differentiation does not
 * depend on the pattern matcher for its core cases. See deriv.c for the
 * full implementation and a discussion of the evaluation strategy.
 */
#ifndef MATHILDA_DERIV_H
#define MATHILDA_DERIV_H

#include "expr.h"

/* Register D, Dt and Derivative with the symbol table. Safe to call
 * multiple times; idempotent with respect to attributes. */
void deriv_init(void);

/* Numeric values of derivatives that stay symbolically inert.
 *
 * D[RiemannR[x], x] is the inert Derivative[1][RiemannR][x] (as in the
 * reference), yet Derivative[1][RiemannR][2.] should be a number, so that
 * FindRoot and friends can Newton on it.  A special-function module
 * registers a C evaluator for its function name; the evaluator calls it for
 * every application Derivative[n..][f][args..] whose innermost head is f.
 * The evaluator gets the whole application (borrowed) and returns a new
 * value, or NULL to leave it inert (symbolic or exact arguments, an order it
 * does not provide).  This is a system rule, not a user SubValue: it does
 * not appear in SubValues[Derivative].  fname must be interned. */
typedef Expr* (*DerivNumericFn)(const Expr* call);
void deriv_register_numeric(const char* fname, DerivNumericFn fn);
/* The registered evaluator for dop = Derivative[n..][f], or NULL. */
DerivNumericFn deriv_numeric_for(const Expr* dop);
/* Apply the registered evaluator to call = Derivative[n..][f][args..]:
 * a new value, or NULL. */
Expr* deriv_apply_numeric(const Expr* call);

/* Public built-ins. The semantics follow Mathematica:
 *   D[f, x]              -- partial derivative of f with respect to x.
 *   D[f, {x, n}]         -- n-th partial derivative (n >= 0).
 *   D[f, x, y, ...]      -- mixed derivative (chained).
 *   Dt[f]                -- total derivative of f (symbols are dependent
 *                            unless they are numeric or recognized constants).
 *   Dt[f, x]             -- equivalent to D[f, x] (partial form).
 *   Dt[f, {x, n}]        -- equivalent to D[f, {x, n}].
 *   Dt[f, x, y, ...]     -- mixed total derivative.
 *   Derivative[...][f]   -- the higher-order derivative operator. The
 *                            built-in below only provides attribute hooks
 *                            for the head symbol Derivative; the composite
 *                            forms Derivative[n][f][x] are simplified
 *                            implicitly by D itself when they appear inside
 *                            a larger expression.
 */
Expr* builtin_d(Expr* res);
Expr* builtin_dt(Expr* res);
Expr* builtin_derivative(Expr* res);

/* Reduce Derivative[n1,...,nm][Function[...]] to an equivalent Function[body']
 * whose body has been differentiated in place.
 *
 *   deriv_head  -- the Derivative[n1,...,nm] head (must have arity >= 1,
 *                  all arguments must be nonnegative integers)
 *   pure_fn     -- the Function[...] expression being wrapped
 *
 * Supports all Function signatures handled by apply_pure_function:
 *   Function[body]              (slot form,  differentiate wrt Slot[i])
 *   Function[Null, body, attrs] (slot form with attributes)
 *   Function[x, body]           (single named parameter)
 *   Function[{x1,...}, body]    (named parameter list)
 *
 * Returns a fresh Function[...] expression on success or NULL if the shape
 * is not supported. The input expressions are NOT consumed. */
Expr* derivative_of_pure_function(Expr* deriv_head, Expr* pure_fn);

/* Reduce Derivative[n1,...,nm][f] where ``f`` is a plain symbol with one
 * or more DownValues, by constructing a Function[{t1,...,tm}, f[t1,...,tm]]
 * whose body is the DownValue-expanded result and then differentiating
 * that body via derivative_of_pure_function.
 *
 *   deriv_head -- the Derivative[n1,...,nm] head (must have arity >= 1,
 *                 all arguments must be nonnegative integers)
 *   fsym       -- the wrapped target; must be EXPR_SYMBOL
 *
 * Returns a fresh Function[...] expression on success. Returns NULL when
 *   - fsym is not a symbol,
 *   - fsym has no DownValues,
 *   - the constructed call f[t1,...,tm] does not rewrite (e.g. arity or
 *     pattern mismatch), or
 *   - any orders in deriv_head are non-integer or negative.
 * The input expressions are NOT consumed. */
Expr* derivative_of_symbol(Expr* deriv_head, Expr* fsym);

#endif /* MATHILDA_DERIV_H */
