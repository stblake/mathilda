/* risch_canonical.h — splitting factorization for the Bronstein Risch algorithm.
 *
 * Implements Bronstein, *Symbolic Integration I*, 2nd ed., §3.4–3.5:
 *
 *   - SplitFactor(p, D)            -> (p_n, p_s) with p = p_n p_s, p_s special,
 *                                     every squarefree factor of p_n normal
 *                                     (Theorem 3.5.1 / algorithm, p.100).
 *   - SplitSquarefreeFactor(p, D)  -> squarefree factorizations of p_n and p_s
 *                                     (algorithm, p.102).
 *
 * These are exposed as internal builtins (backtick-namespaced under `Risch``),
 * grounded in the differential-field primitives of risch_field.{c,h}.  The
 * derivation is passed as a List[Rule[var, Dvar], ...] naming the base variable
 * (Dvar = 1) and each monomial variable with its derivative.
 *
 *   Risch`SplitFactor[p, t, deriv]           -> {p_n, p_s}
 *   Risch`SplitSquarefreeFactor[p, t, deriv] -> {{N_1, ...}, {S_1, ...}}
 *   Risch`Derivation[p, deriv]               -> D[p]      (monomial derivation)
 *   Risch`NormalQ[p, t, deriv]               -> True | False
 *   Risch`SpecialQ[p, t, deriv]              -> True | False
 *
 * Memory contract: the standard Mathilda BuiltinFunc rule — the evaluator owns
 * the argument; a builtin returns a fresh Expr* on success or NULL on failure.
 */

#ifndef MATHILDA_RISCH_CANONICAL_H
#define MATHILDA_RISCH_CANONICAL_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"
#include "risch_field.h"

/* SplitFactor (Bronstein §3.5).  Writes owned p_n, p_s with p = p_n p_s, p_s
 * special, and every squarefree factor of p_n normal.  `t` is the monomial
 * variable and `d` its derivation. */
void risch_split_factor(const Expr* p, const Expr* t, const RischDeriv* d,
                        Expr** p_n, Expr** p_s);

/* CanonicalRepresentation (Bronstein §3.5, p.103).  Writes owned f_p, f_s, f_n
 * with f = f_p + f_s + f_n, f_p polynomial, f_s = b/d_s (special denominator,
 * proper), f_n = c/d_n (normal denominator, proper). */
void risch_canonical_representation(const Expr* f, const Expr* t, const RischDeriv* d,
                                    Expr** f_p, Expr** f_s, Expr** f_n);

/* Register the Risch` splitting-factorization builtins.  Called from
 * integrate_init() during core_init().  Idempotent. */
void integrate_risch_canonical_init(void);

#endif /* MATHILDA_RISCH_CANONICAL_H */
