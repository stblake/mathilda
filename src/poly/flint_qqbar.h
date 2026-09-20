/* flint_qqbar.h — exact algebraic-number canonicalisation for RootReduce.
 *
 * Backs the WL-faithful behaviour of RootReduce on *constant* algebraic
 * numbers (no free symbols): an expression built from integers, rationals,
 * radicals Power[base, p/q], roots of unity (-1)^(p/q), the imaginary unit,
 * and Root[Function[t, poly], k] objects, combined by +, -, *, /, ^, is
 * canonicalised to a single representative — a rational number, a quadratic
 * radical expression, or a Root[Function[minpoly&], k] object.
 *
 * The engine is FLINT's `qqbar` (exact real/complex algebraic numbers via
 * minimal polynomial + isolating enclosure): rigorous, with no numeric zero
 * oracle. When FLINT is not compiled in (USE_FLINT off) every entry point is a
 * graceful no-op returning NULL / -1, matching the rest of the FLINT bridge.
 *
 * Parametric algebraic *functions* (radicals whose radicand carries a free
 * variable, e.g. the Goursat k^(1/3) towers) are NOT handled here — they stay
 * with flint_algebraic_field_canonical in flint_bridge.c. RootReduce dispatches
 * between the two by the presence of free symbols.
 */
#ifndef FLINT_QQBAR_H
#define FLINT_QQBAR_H

#include "expr.h"

/* RootReduce Method selector (see the WL Method -> option). */
typedef enum {
    QQBAR_METHOD_AUTOMATIC = 0,
    QQBAR_METHOD_RECURSIVE,
    QQBAR_METHOD_NUMBERFIELD
} QQBarMethod;

/* True if `e` has no free symbol and every head is one the qqbar converter can
 * handle (a necessary, cheap pre-filter; conversion may still fail on a degree
 * blow-up or an unsupported atom). Integers, rationals, radicals, roots of
 * unity, the imaginary unit and Root objects qualify; a bare symbol, Pi, E,
 * Log, Sin, ... do not. */
int flint_qqbar_is_constant_algebraic(const Expr* e);

/* Canonicalise the constant algebraic number `e` under the chosen Method.
 * Returns a fresh owned Expr (rational / quadratic radical / Root object) on
 * success, or NULL if `e` is not a constant algebraic number, exceeds the
 * degree cap, or FLINT is unavailable. Never mutates `e`. */
Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod method);

/* Exact algebraic equality of two constant algebraic numbers.
 * Returns 1 (equal), 0 (unequal), or -1 (undecided: either operand is not a
 * constant algebraic number, or FLINT is off). */
int flint_qqbar_equal(const Expr* a, const Expr* b);

/* Sign of (a - b) for real constant algebraic numbers: -1, 0, or 1;
 * returns -2 when undecided (non-real, non-algebraic, or FLINT off). */
int flint_qqbar_compare(const Expr* a, const Expr* b);

/* Realness of a constant algebraic number: 1 (a real algebraic number),
 * 0 (a genuinely non-real algebraic number), or -1 (undecided: `e` is not a
 * constant algebraic number, exceeds the degree cap, or FLINT is off).  Lets a
 * Reals-domain filter distinguish "drop this non-real branch" (0) from "cannot
 * decide -- leave the whole statement unevaluated" (-1). */
int flint_qqbar_is_real(const Expr* e);

/* ------------------------------------------------------------------ */
/*  AlgebraicNumber / ToNumberField Expr-level entry points            */
/*  (used by the builtins in algebraicnumber.c / tonumberfield.c).     */
/*  Every one returns a fresh owned Expr on success and NULL to        */
/*  decline — a non-constant-algebraic argument, malformed coeffs, a   */
/*  degree-cap overflow, or FLINT compiled out. None mutate an input.  */
/* ------------------------------------------------------------------ */

/* Canonicalise AlgebraicNumber[gen, coeffs] (coeffs a List of integers /
 * rationals): reduce gen to an algebraic integer and coeffs to the power basis
 * of its minimal polynomial. Returns the reduced AlgebraicNumber[g, {..}], or a
 * plain Integer/Rational when the value is rational. */
Expr* flint_qqbar_algebraic_number(const Expr* gen, const Expr* coeffs);

/* ToNumberField[a, theta]: express the constant algebraic number `a` in the
 * field Q(theta) as AlgebraicNumber[g, {..}] (g the algebraic-integer generator
 * of Q(theta)); returns the rational `a` when Q(theta) = Q. NULL if a is not in
 * Q(theta). */
Expr* flint_qqbar_to_number_field(const Expr* a, const Expr* theta);

/* ToNumberField[{a0..a_{n-1}}] (or `All`): build a common field Q(a0,..) via a
 * primitive element and return List[AlgebraicNumber[g, {..}], ..]. `smallest`
 * is advisory (Automatic and All both use the qqbar primitive element). */
Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n,
                                         int smallest);

/* ToNumberField[x]: express the single algebraic number x as an explicit
 * AlgebraicNumber in Q(x). */
Expr* flint_qqbar_to_number_field_self(const Expr* x);

/* Same-field arithmetic for the Plus/Times/Power combination pre-passes. `a`
 * (and `b`) are AlgebraicNumber objects; combining happens only when generators
 * are structurally equal (expr_eq). `r` is an Integer/Rational scalar. */
Expr* flint_qqbar_algnum_add(const Expr* a, const Expr* b);
Expr* flint_qqbar_algnum_mul(const Expr* a, const Expr* b);
Expr* flint_qqbar_algnum_pow(const Expr* a, long p);
Expr* flint_qqbar_algnum_add_rational(const Expr* a, const Expr* r);
Expr* flint_qqbar_algnum_scale_rational(const Expr* a, const Expr* r);

#endif /* FLINT_QQBAR_H */
