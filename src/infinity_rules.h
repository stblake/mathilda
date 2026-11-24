#ifndef MATHILDA_INFINITY_RULES_H
#define MATHILDA_INFINITY_RULES_H

/* Numeric special values: one mechanism for Indeterminate and the infinities
 * as arguments of a NumericFunction head (audit Wave 3, shared fix S11).
 *
 *   - An Indeterminate argument makes any NumericFunction call Indeterminate
 *     (the list-valued ReIm / QuotientRemainder give a list of them).
 *   - A per-head table answers directed infinities and ComplexInfinity where
 *     the value is systematic: Sin[ComplexInfinity] is Indeterminate,
 *     Log[-Infinity] is Infinity, BesselJ[n, Infinity] is 0, Abs/Re/Im/Sign/
 *     Arg/Conjugate read the direction, Power raises it.
 *
 * The evaluator consults the rule only when its argument classification saw
 * such an argument (eval.c, listable_arg_kinds), so the ordinary numeric call
 * pays nothing beyond that pass. */

#include "expr.h"
#include <stdbool.h>

/* Argument-kind bits, shared with eval.c's listable classification. */
#define INFR_ARG_INDET 32
#define INFR_ARG_INF   64

/* INFR_ARG_INDET / INFR_ARG_INF / 0 for one argument (the non-Listable
 * NumericFunction heads, whose arguments eval.c does not classify). */
int numeric_special_arg_kind(const Expr* arg);

/* The special value of the NumericFunction call `call` (borrowed), or NULL
 * when no rule applies. `has_indet` says an argument is Indeterminate. */
Expr* numeric_special_value_rule(const Expr* call, bool has_indet);

/* Direction class of an infinite quantity. */
typedef enum {
    INF_DIR_NONE = 0,   /* not an infinity */
    INF_DIR_POS,        /* Infinity */
    INF_DIR_NEG,        /* -Infinity */
    INF_DIR_I,          /* I Infinity */
    INF_DIR_NEG_I,      /* -I Infinity */
    INF_DIR_CINF,       /* ComplexInfinity */
    INF_DIR_OTHER       /* c Infinity, c a numeric quantity off the axes */
} InfDir;

/* Classify `e`; for INF_DIR_OTHER, *dir (if non-NULL) receives a new
 * expression for the (unnormalised) direction c. */
InfDir infinity_direction(const Expr* e, Expr** dir);

/* p_n at an infinity for the classical polynomial heads of explicit degree
 * (ChebyshevT/U, HermiteH, LegendreP, LaguerreL, BernoulliB), or NULL. */
Expr* infinity_poly_value(const Expr* call);

/* DirectedInfinity[z]: the canonical infinity in direction z. */
Expr* builtin_directedinfinity(Expr* res);

/* The infinity the exact call takes at the point an inexact call that hit a
 * pole names (Cot[0.] -> ComplexInfinity, ArcTanh[1.] -> Infinity), or NULL
 * when the arguments are not small dyadic points or the exact call is not
 * infinite there. Borrows `call`. */
Expr* numeric_pole_value(const Expr* call);

#endif
