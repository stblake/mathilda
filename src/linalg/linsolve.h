#ifndef MATSOL_H
#define MATSOL_H

#include <stdint.h>
#include "expr.h"

/* RowReduce[m]                                  -- default Method -> Automatic.
 * RowReduce[m, Method -> "<name>"]              -- explicit dispatcher.
 *
 * LinearSolve[m, b]                             -- default Method -> Automatic.
 * LinearSolve[m, b, Method -> "<name>"]         -- explicit dispatcher.
 *
 * Supported method names:
 *   "Automatic"                  -- aliased to "DivisionFreeRowReduction"
 *   "DivisionFreeRowReduction"   -- Bareiss-like fraction-free Gauss-Jordan
 *   "OneStepRowReduction"        -- classical Gauss-Jordan with division
 *   "CofactorExpansion"          -- Cramer's rule (LinearSolve); identity-if-
 *                                   invertible (RowReduce), with fallback to
 *                                   DivisionFreeRowReduction on
 *                                   singular / rectangular input.
 *                                   For Inverse: adjugate / det formula.
 *
 * Method -> Automatic (the symbol) is also accepted, matching Mathematica.
 *
 * The MatsolMethod enum and helpers below are shared with src/matinv.c
 * so that Inverse, RowReduce and LinearSolve all parse the same
 * Method-option grammar and reuse the same Together / exact-division
 * primitives used inside the OneStep workers.
 */
typedef enum {
    MATSOL_AUTOMATIC = 0,
    MATSOL_DIVFREE,
    MATSOL_ONESTEP,
    MATSOL_COFACTOR,
    MATSOL_INVALID
} MatsolMethod;

/* Parse a `Method -> <value>` rule.  Returns MATSOL_INVALID when the
 * argument isn't a Method rule with a recognised RHS. */
MatsolMethod matsol_parse_method_option(Expr* opt);

/* Resolve Method -> Automatic for the matrix argument `m`.
 *
 * Returns MATSOL_ONESTEP when `m` has any inexact (Real / MPFR) leaf and
 * MATSOL_DIVFREE otherwise. Fraction-free elimination is pointless over
 * floating point and its per-pivot polynomial GCD recursion overflowed the
 * stack on a 90x90 Real matrix; see the definition in linsolve.c. */
MatsolMethod matsol_resolve_automatic(const Expr* m);

/* True when `e` contains any inexact (Real / MPFR, including inside a Complex
 * or a nested List) leaf. Used both to resolve Method -> Automatic and to pick
 * the TYPE of the identity element in the OneStep eliminators: normalising a
 * pivot cell to the exact Integer 1 inside a Real computation puts an exact
 * value in an inexact result, and Inverse of a Real matrix came back with
 * {1, -2.0, 1} instead of {1.0, -2.0, 1.0}. */
bool matsol_is_inexact(const Expr* e);

/* Rate-limit a per-call warning so test loops don't spew.  `key` is
 * hashed to detect repeated invocations of the same call. */
void matsol_warn_once(uint64_t* last_hash, Expr* key, const char* msg);

/* Canonicalise a matrix entry via `Together` so subsequent
 * is_zero_poly checks behave correctly on symbolic rationals. */
Expr* matsol_canon_entry(Expr* e);

/* num / den, canonicalised.  Prefers exact polynomial division via
 * exact_div_wrapper when applicable, else builds Times[num, Power[den,
 * -1]] then Together-canonicalises. */
Expr* matsol_div_entry(Expr* num, Expr* den);

/* ZeroTest option (shared by RowReduce and NullSpace).
 *
 * RowReduce[m, ZeroTest -> f] / NullSpace[m, ZeroTest -> f]: f is a predicate
 * applied to each entry; the entry counts as zero iff TrueQ[f[entry]].  Both a
 * body (RootReduce[Together[#]] === 0 &) and a predicate head (PossibleZeroQ)
 * work. */

/* Extract a `ZeroTest -> f` rule: returns f (borrowed from `opt`) or NULL. */
Expr* matsol_parse_zerotest_option(Expr* opt);

/* Zero test for an entry: default (zt == NULL) is the structural is_zero_poly;
 * with a predicate, zero iff TrueQ[zt[e]].  `zt` and `e` are borrowed. */
int matsol_zt_is_zero(Expr* e, Expr* zt);

/* In-place exact reduced row echelon form of the row-major flat[rows*cols]
 * (owned Exprs), consulting `zt` for every zero decision and keeping each
 * updated entry canonical via Together.  A structural RowReduce would pivot on
 * an algebraic zero the predicate rejects, so the reduction itself uses zt. */
void matsol_rref_with_zerotest(Expr** flat, int rows, int cols, Expr* zt);

Expr* builtin_rowreduce(Expr* res);
Expr* builtin_linearsolve(Expr* res);
void  matsol_init(void);

#endif /* MATSOL_H */
