/*
 * reduce_algfiber.h
 *
 * Real-algebraic-coefficient fibre isolation for `Reduce`'s CAD (REDUCE_PLAN.md,
 * Phase 6b).  The n-variable CAD lifts a cell by substituting the outer sample
 * point into a projection factor and isolating the resulting univariate's real
 * roots.  When an outer coordinate is an IRRATIONAL algebraic number (a section
 * at a non-innermost level), that substitution produces a univariate polynomial
 * whose coefficients are real-algebraic numbers -- which `Solve[..., Reals]`
 * (factoring over Q) cannot isolate.  This module supplies that missing
 * capability by TOWER-GENERAL iterated-resultant projection back down to Q:
 *
 *   - the outer assignment is a tower Q <= Q(a0) <= Q(a0,a1) <= ... in which each
 *     irrational coordinate a_i is a root of a KNOWN factor `defs[i]` over Q
 *     (identified by the CAD's root provenance);
 *   - substitute the rational coordinates directly, then eliminate each algebraic
 *     tower variable v_i (deepest first) by Resultant[R, defs[i], v_i], leaving a
 *     univariate R(var) over Q whose real roots CONTAIN the fibre's real roots
 *     (plus spurious roots from the conjugates the resultant also eliminates);
 *   - isolate R's real roots with the ordinary integer-coefficient path, then keep
 *     exactly the candidates b for which the ORIGINAL fibre vanishes at the true
 *     assignment -- factor(vals.., b) == 0, decided EXACTLY by the qqbar oracle.
 *
 * Soundness invariant (inherited by the CAD callers): any undecidable qqbar test,
 * an unclean Solve shape, a resource-budget overrun, a degenerate identically-zero
 * resultant, or FLINT being absent makes the routine return false (bail), so
 * `Reduce` stays unevaluated rather than emitting a wrong formula.
 */
#ifndef REDUCE_ALGFIBER_H
#define REDUCE_ALGFIBER_H

#include "expr.h"
#include "reduce_opts.h"
#include <stdbool.h>

/* Isolate & order the distinct real roots of the fibre polynomial `factor` (a
 * polynomial in `var` whose coefficients are polynomials in the outer variables
 * vv[0..nlev-1]) at the outer assignment vals[0..nlev-1].  Each coordinate i is
 * either RATIONAL (defs[i] == NULL: vals[i] is substituted directly) or ALGEBRAIC
 * (defs[i] != NULL: a factor over Q, a polynomial in vv[0..i], whose root in vv[i]
 * at the assignment is vals[i]).
 *
 * On success the surviving real roots (owned Root[]/rational Exprs) are appended
 * to (*arr, *n, *cap); when `prov` is non-NULL it is grown in lockstep and each
 * pushed root records `factor_id` (mirroring rru_collect_roots, for the CAD's
 * symbolic sector emission).  `opts` forwards the Cubics/Quartics radical flags
 * onto the internal integer-coefficient isolation (NULL keeps Solve's defaults).
 *
 * Returns false to bail (undecidable qqbar sign/equality, unclean Solve, a
 * resource-budget overrun, an identically-zero resultant, or FLINT off), leaving
 * every already-pushed *arr entry valid for the caller to free. */
bool rru_algebraic_fiber_roots(const Expr* factor, const Expr* var,
                               Expr** vv, Expr** vals, Expr** defs, int nlev,
                               Expr*** arr, int* n, int* cap,
                               int** prov, int factor_id,
                               const ReduceOpts* opts);

#endif /* REDUCE_ALGFIBER_H */
