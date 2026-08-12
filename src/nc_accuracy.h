/* Mathilda — shared AccuracyGoal / PrecisionGoal handling for the numerical
 * calculus and root/solve functions.
 *
 * Every numerical routine (NLimit, NSum, NProduct, NIntegrate, NDSolve,
 * NResidue, ND, NSeries, FindRoot, NRoots, NSolve) seeks a result to a
 * requested number of digits. Two options govern that target:
 *
 *   AccuracyGoal  a  — absolute error allowed        (10^-a)
 *   PrecisionGoal p  — relative error allowed        (|x| 10^-p)
 *
 * The success criterion is the combined tolerance
 *
 *       estimated_error  <=  10^-a + |x| 10^-p
 *
 * exactly as documented for AccuracyGoal/PrecisionGoal. This header holds the
 * one shared parser and tolerance calculator so the eleven drivers do not each
 * re-derive them (they previously carried five near-identical copies), and so
 * the `MachinePrecision` setting is understood in one place.
 *
 * A driver's contract on a shortfall is: adaptively refine up to its
 * (user-overridable) resource cap, keep the best-of estimate, and if the goal
 * is still unmet at the cap, call nc_warn_goal() and RETURN THE BEST
 * APPROXIMATION — never NULL / $Failed. */
#ifndef NC_ACCURACY_H
#define NC_ACCURACY_H

#include <stdbool.h>
#include <math.h>
#include "expr.h"

/* Symbolic goal settings, encoded as digit counts so a plain `double` field
 * can carry any setting a user writes. */
#define NC_GOAL_AUTO    (-1.0)   /* Automatic -> resolve near full WorkingPrecision */
#define NC_GOAL_MACHINE (-2.0)   /* MachinePrecision -> track WorkingPrecision (default) */
#define NC_GOAL_OFF  INFINITY    /* Infinity  -> criterion disabled (term -> 0) */

/* Parse an AccuracyGoal / PrecisionGoal option value into a decimal-digit
 * count written to *digits_out:
 *   Automatic        -> NC_GOAL_AUTO
 *   Infinity         -> NC_GOAL_OFF
 *   MachinePrecision -> NC_GOAL_MACHINE (tracks WorkingPrecision, like Automatic)
 *   positive number  -> that value
 * Returns false (leaving *digits_out untouched) for any other value, so the
 * caller can report an invalid-option error. */
bool nc_parse_goal(const Expr* v, double* digits_out);

/* Resolve a goal (possibly NC_GOAL_AUTO) to a concrete digit count given the
 * working precision in decimal digits (machine precision ==
 * NUMERIC_MACHINE_PRECISION_DIGITS). Automatic -> wp_digits/2; every other
 * value — an explicit number or NC_GOAL_OFF (+Inf) — passes through. */
double nc_resolve_goal(double goal_digits, double wp_digits);

/* Combined absolute+relative tolerance for a result of magnitude `result_mag`:
 *     tol = 10^-a + |x| * 10^-p
 * with a = resolve(acc_digits), p = resolve(prec_digits). An Infinity goal
 * contributes 0 to its term (10^-Inf == 0), matching the documented meaning
 * that accuracy/precision is not used as a termination criterion. */
double nc_combined_tol(double acc_digits, double prec_digits,
                       double result_mag, double wp_digits);

/* Emit the standard "goal not met, returning best approximation" warning to
 * stderr in Mathematica's `Head::tag` style. Call once, only when the routine
 * has exhausted its refinement budget and still has achieved_err > tol. */
void nc_warn_goal(const char* head, double achieved_err, double tol);

#endif /* NC_ACCURACY_H */
