/*
 * reduce_realfn.h
 *
 * Elementary-real-function support for `Reduce` over the Reals (Phase 9): the
 * preprocessing rewrites (Abs sign-splitting, Mod->Floor, integer-part
 * defining-inequality isolation) and the head->domain table that the general
 * univariate sign diagram (reduce_realdiag.c) consults.
 *
 * The polynomial engines (reduce_univar.c, reduce_cad.c) require every atom to
 * be a polynomial.  This module is what lets a univariate real statement built
 * from Abs, real radicals (u^(p/q)), Log, bounded-domain inverse-trig, and the
 * isolated integer-part / Mod forms reach the sign-diagram machinery: Abs is
 * eliminated by case analysis, Mod/Floor by their defining inequalities, and the
 * remaining partial-domain functions are described by the domain table (a real
 * domain constraint + the boundary whose roots are breakpoints).
 */
#ifndef REDUCE_REALFN_H
#define REDUCE_REALFN_H

#include "expr.h"
#include <stdbool.h>

/* One real-domain constraint:  poly >= 0  (strict == false)  or  poly > 0
 * (strict == true).  A partial-domain node contributes one or two of these; the
 * general engine roots each `poly` for breakpoints and tests each at every
 * sample point as the domain gate. */
typedef struct { Expr* poly; bool strict; } RDomCon;

/* True iff `e` contains, with `x` somewhere in the argument, a node the
 * polynomial engines cannot take directly: Abs, Floor/Ceiling/Round/IntegerPart,
 * Mod, Log, bounded-domain inverse-trig, or an even-order rational-power
 * radical.  Drives reduce.c's decision to route to the Reals and preprocess. */
bool reduce_stmt_has_realfn(const Expr* e, const Expr* x);

/* True iff `e` contains a selector piecewise head (Abs/Max/Min/Piecewise/Sign/
 * UnitStep/Ramp/Clip/HeavisideTheta/Boole/UnitBox) with one of the `nv` reduce
 * `vars` under it.  Drives the multivariate dispatch to case-split the piecewise
 * heads (reduce_piecewise_preprocess) into polynomial branches before FM/CAD. */
bool reduce_stmt_has_piecewise(const Expr* e, Expr** vars, int nv);

/* Rewrite the statement so the general real sign diagram can consume it:
 *   1. substitute every Mod[u,m] (m a positive constant) by u - m*Floor[u/m];
 *   2. expand a relational leaf that is linear in a single Floor/Ceiling/Round
 *      into its defining inequalities;
 *   3. eliminate every Abs[u] by sign-splitting into an Or of And-branches.
 * Radicals, Log and inverse-trig are left intact (the domain table handles
 * them).  Returns a freshly-owned rewritten Expr, or NULL when nothing fired
 * (the caller then keeps the original expression).  *changed reports whether a
 * rewrite fired. */
Expr* reduce_realfn_preprocess(const Expr* e, const Expr* x, bool* changed);

/* Multivariate (any variable count) piecewise preprocessing: the domain-agnostic
 * selector splits only (Abs, Min/Max, Piecewise/Sign/UnitStep/Ramp/Clip/
 * HeavisideTheta/Boole/UnitBox), iterated to a fixpoint.  Returns a freshly-owned
 * rewritten Expr, or NULL when nothing fired.  *changed reports whether a rewrite
 * fired.  The integer-part machinery is univariate and lives in
 * reduce_realfn_preprocess instead. */
Expr* reduce_piecewise_preprocess(const Expr* e, bool* changed);

/* Append the real-domain constraints of every partial-domain node in `e` whose
 * argument contains `x` to the growable array (*cons,*n,*cap).  Each RDomCon
 * owns its `poly`; the caller frees them.  Used per-atom by the general engine
 * to build both its domain gate and its domain-boundary breakpoints. */
void reduce_real_domain_collect(const Expr* e, const Expr* x,
                                RDomCon** cons, int* n, int* cap);

#endif /* REDUCE_REALFN_H */
