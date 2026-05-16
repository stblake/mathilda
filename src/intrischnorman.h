/* intrischnorman.h — parallel Risch / Risch-Norman heuristic
 * integrator (`Integrate`RischNorman[f, x]`).
 *
 * Implements Bronstein's "Poor Man's Integrator" (pmint, 2004),
 * which is the modern reformulation of Norman's 1976 heuristic.
 * The reference Maple source is `parallel_risch/pmint.maple`
 * (99 lines).  Our C port is a near line-for-line translation that
 * reuses Mathilda's existing CAS primitives (Together / Coefficient /
 * Factor / PolynomialGCD / RowReduce).
 *
 * The entry `Integrate`RischNorman[f, x]` is called by the
 * `Integrate` dispatcher (src/integrate.c) as a fall-through after
 * `Integrate`BronsteinRational` declines a non-rational integrand.
 * It returns either a fresh Expr* antiderivative on success or NULL
 * on failure (the caller bubbles the call symbolic).
 *
 * Memory contract follows the standard Mathilda BuiltinFunc rule:
 * the caller (evaluator) owns `res`.  On success the builtin returns
 * a freshly-allocated Expr*; on failure it returns NULL and the
 * evaluator preserves the call unevaluated.
 *
 * Implementation phases (see RISCH_NORMAN_PLAN.md):
 *   Phase 1: skeleton + dispatcher hook (always returns NULL).
 *   Phase 2: convert_to_tan + collect_indets + substitution maps.
 *   Phase 3: vector field + splitFactor + deflation + monomials.
 *   Phase 4: candidate ansatz + linear-system extract + RowReduce.
 *   Phase 5: log-candidate sum + getSpecial + K=I retry.
 *   Phase 6: dispatcher polish + post-hoc verification.
 */

#ifndef MATHILDA_INTRISCHNORMAN_H
#define MATHILDA_INTRISCHNORMAN_H

#include "expr.h"

/* Public package entry: corresponds to Bronstein's pmint(f, x).
 *
 * Returns a freshly-allocated antiderivative Expr* on success or
 * NULL on failure (the call bubbles back as unevaluated).
 */
Expr* builtin_rischnorman(Expr* res);

/* Register every Integrate`RischNorman* package symbol in the
 * global symbol table.  Called from integrate_init() during
 * core_init().  Idempotent.
 */
void intrischnorman_init(void);

#endif /* MATHILDA_INTRISCHNORMAN_H */
