/* intrat.h — rational-function integrator (`Integrate`` package).
 *
 * This header declares the C entry points that back every public
 * symbol in the `Integrate`` context.  Each builtin is exposed under
 * its fully qualified name (e.g. `Integrate`HermiteReduce`) so that
 * package consumers can call them directly from the REPL for
 * end-to-end testing of the pipeline.
 *
 * Phase 1 of INTEGRATE_PLAN.md:
 *   * IntegrateRational    — top-level package entry (skeleton, expanded
 *                            in later phases).
 *   * IntegratePolynomial  — term-by-term integration of polynomials.
 *   * HermiteReduce        — Mack's linear Hermite reduction.
 *   * Helpers              — Content / Primitive / Monic /
 *                            LeadingCoefficient unit-test handles.
 *
 * Memory contract follows the standard picocas BuiltinFunc rule: the
 * caller (evaluator) owns `res`.  On success the builtin returns a new
 * Expr*; on failure it returns NULL and the evaluator preserves the
 * call unevaluated.
 */

#ifndef PICOCAS_INTRAT_H
#define PICOCAS_INTRAT_H

#include "expr.h"

/* Public package entry: corresponds to the Mathematica function
 * `IntegrateRational[f, x]`.  Phase 1 ships a thin implementation that
 * (a) tries the derivative-recognition fast path (`c*D'/D^k`),
 * (b) splits into HermiteReduce + IntegratePolynomial when possible,
 * and (c) returns the call unevaluated otherwise (the
 * Lazard-Rioboo-Trager log part lands in Phase 2). */
Expr* builtin_intrat_integraterational(Expr* res);

/* Polynomial-only term-by-term integration.  Expand[f,x] then maps
 * a*x^n -> a*x^(n+1)/(n+1) (with the n=-1 → Log[x] special case). */
Expr* builtin_intrat_integratepolynomial(Expr* res);

/* HermiteReduce[f, x]: returns {g, h} with f == D[g, x] + h and
 * Denominator[h] squarefree.  Direct port of the Mathematica
 * implementation at IntegrateRational.m:1303-1323. */
Expr* builtin_intrat_hermitereduce(Expr* res);

/* Helpers exposed for unit-testing of the lower-level building blocks. */
Expr* builtin_intrat_helpers_content(Expr* res);              /* gcd of CoefficientList */
Expr* builtin_intrat_helpers_primitive(Expr* res);            /* p / Content[p, x] */
Expr* builtin_intrat_helpers_monic(Expr* res);                /* p / lc[p, x]      */
Expr* builtin_intrat_helpers_leadingcoefficient(Expr* res);   /* coeff at top deg  */

/* Phase 2 — Lazard-Rioboo-Trager log part. */

/* IntRationalLogPart[A/D, x, t]: returns either a list of {Q_i(t),
 * S_i(t,x)} pairs (RootSum -> False, the default for tooling) or a
 * sum of RootSum heads (RootSum -> True). */
Expr* builtin_intrat_intrationallogpart(Expr* res);

/* SquareFree[p]: list of {factor, multiplicity} pairs indexed
 * densely by multiplicity 1..max.  IntegrateRational.m:1474-1487. */
Expr* builtin_intrat_helpers_squarefree(Expr* res);

/* ExtractConstants[f, x]: returns {const, simplified_f}. */
Expr* builtin_intrat_helpers_extractconstants(Expr* res);

/* ApartList[f, x, Extension -> alpha]: returns Apart's output as a
 * list of summands. */
Expr* builtin_intrat_helpers_apartlist(Expr* res);

/* Phase 3 — LogToAtan (Rioboo recursive). */
Expr* builtin_intrat_logtoatan(Expr* res);

/* Phase 4 — LogToReal (Rioboo: complex Log → real Log + ArcTan). */
Expr* builtin_intrat_logtoreal(Expr* res);

/* Phase 6 — Log → ArcTanh / combined Log post-processing. */
Expr* builtin_intrat_logtoarctanh(Expr* res);

/* Phase 8b — NaiveLogPart RootSum fallback.  Direct port of
 * IntegrateRational.m:1116-1124.  Returns the held-symbolic
 * RootSum form of the log part of Integrate[a/d, x]:
 *
 *   RootSum[Function[t, d(t)],
 *           Function[t, a(t) Log(x - t) / d'(t)]]
 *
 * Used as the universal fallback by the rational integrator when
 * LogToReal cannot close the log part to a real elementary
 * expression (Phase 8c). */
Expr* builtin_intrat_naive_log_part(Expr* res);

/* Register every Integrate` package symbol in the global symbol table.
 * Called from integrate_init() during core_init().  Idempotent (each
 * registration is a fresh symtab_add_builtin / attribute set). */
void intrat_init(void);

#endif /* PICOCAS_INTRAT_H */
