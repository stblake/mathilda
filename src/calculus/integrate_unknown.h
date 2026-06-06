/* integrate_unknown.h
 *
 * Integration of expressions containing *undefined* (unknown) functions
 * and their derivatives, following Kelly Roach, "Indefinite and Definite
 * Integration" (1992), §1.7 "Undefined Functions".
 *
 * Treats each undefined function value u(g) and its derivative tower
 * u'(g), u''(g), ... as differential-field generators and integrates
 * rational expressions in them.  Examples:
 *
 *     Integrate[f'[x], x]                    -> f[x]
 *     Integrate[x f'[x] + f[x], x]           -> x f[x]
 *     Integrate[f'[x] g[x] + f[x] g'[x], x]  -> f[x] g[x]
 *
 * Provides the package builtin `Integrate`Undefined[f, x]` and a direct
 * entry point used by the Integrate dispatcher.
 */
#ifndef INTEGRATE_UNKNOWN_H
#define INTEGRATE_UNKNOWN_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher.  Returns a freshly
 * owned antiderivative of `f` with respect to symbol `x`, or NULL when
 * the integrand contains no undefined-function derivative or the method
 * does not apply.  Does not take ownership of `f` or `x`. */
Expr* integrate_unknown_try(Expr* f, Expr* x);

/* `Integrate`Undefined[f, x]` builtin (strict; returns NULL / unevaluated
 * on failure). */
Expr* builtin_integrate_unknown(Expr* res);

/* Register the package builtin + attributes + docstring. Called from
 * integrate_init(). */
void integrate_unknown_init(void);

#endif /* INTEGRATE_UNKNOWN_H */
