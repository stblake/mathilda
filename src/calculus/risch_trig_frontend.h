/* risch_trig_frontend.h — trig/hyperbolic front-end for the Risch integrator.
 *
 * Integrates trigonometric / hyperbolic integrands by exponentialising them
 * (TrigToExp) to a Laurent-rational function of E^(I x) / E^x, running the
 * exponential machinery, converting back (ExpToTrig), and reconstructing a real
 * closed form (rt_realify, via an exact real/imaginary split, diff-back gated).
 * Also the real hypertangent §5.10 case (rational functions of Tan/Tanh) and
 * the rational-function-of-a-single-exponential closer.  Defined in
 * risch_trig_frontend.c.
 */

#ifndef MATHILDA_RISCH_TRIG_FRONTEND_H
#define MATHILDA_RISCH_TRIG_FRONTEND_H

#include "expr.h"

/* Exponentialise-integrate-realify front-end (NULL == decline). */
Expr* rt_trig_frontend(Expr* f, Expr* x);

/* Real §5.10 hypertangent case: rational function of Tan[u]/Tanh[u]. */
Expr* rt_hypertangent_case(Expr* f, Expr* x);

/* Rational-function-of-a-single-exponential closer (kernelize E^(I x)). */
Expr* rt_exp_ratreduce_case(Expr* f, Expr* x);

#endif /* MATHILDA_RISCH_TRIG_FRONTEND_H */
