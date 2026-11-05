/*
 * hrep.h -- formula representations of hypergeometric functions.
 *
 * Implements Kelly Roach's algorithm ("Hypergeometric Function
 * Representations", ISSAC '96; "Meijer G Function Representations",
 * ISSAC '97). A pFq is written in closed form by walking from a known
 * representation at a nearby "origin" (a lookup certificate) with shift and
 * inverse-shift (contiguity) operators. Each step is an exact operation on a
 * coefficient vector over Q(parameters)(w), and the vector lives in the
 * C(w)-module spanned by the origin's companion basis
 * [f0, theta f0, ..., theta^(d-1) f0], theta = w d/dw.
 *
 * See hrep_core.c for the algebra and docs/spec/builtins/special-functions.md
 * ("Formula representations") for the user-level contract.
 */
#ifndef MATHILDA_HREP_H
#define MATHILDA_HREP_H

#include "expr.h"
#include <stddef.h>

/* Closed form of HypergeometricPFQ[{a1..ap}, {b1..bq}, z], or NULL when the
 * algorithm declines (no certificate in the residue class, an improper walk,
 * parameters outside Q(symbols), FLINT absent). Borrows every argument; the
 * result is owned and evaluated. Generic in symbolic parameters: an identity
 * is returned when it holds for all but isolated parameter values. */
Expr* hrep_pfq(Expr* const* a, size_t p, Expr* const* b, size_t q,
               const Expr* z);

/* The same for an evaluated HypergeometricPFQ[{..}, {..}, z] node; NULL when
 * `e` is not one or the algorithm declines. Borrows `e`. */
Expr* hrep_pfq_expr(const Expr* e);

/* Closed form of MeijerG[{{an}, {ap}}, {{bm}, {bq}}, z] by the G version of
 * the walk (Roach '97), or NULL. hrep_meijerg_expr takes the evaluated node
 * and runs as a probe region. */
Expr* hrep_meijerg(Expr* const* an, int n, Expr* const* ap, int np_,
                   Expr* const* bm, int m, Expr* const* bq, int nq_, const Expr* z);
Expr* hrep_meijerg_expr(const Expr* e);

/* Registers the diagnostic builtin HypergeometricPFQ`Representation. */
void hrep_init(void);

#endif /* MATHILDA_HREP_H */
