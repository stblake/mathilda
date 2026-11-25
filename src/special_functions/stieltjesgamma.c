/* Mathilda -- StieltjesGamma, the Stieltjes constants gamma_n and their
 * generalisation gamma_n(a).
 *
 *   StieltjesGamma[n] = gamma_n, the coefficients of the Laurent expansion
 *   of the Riemann zeta function about s = 1:
 *
 *     zeta(s) = 1/(s-1) + Sum_{n>=0} ((-1)^n / n!) gamma_n (s-1)^n,
 *
 *   and StieltjesGamma[n, a] = gamma_n(a), the same coefficients for the
 *   Hurwitz zeta function zeta(s, a).
 *
 * gamma_0 is EulerGamma. The higher constants have no elementary closed form,
 * so an exact call stays symbolic, apart from StieltjesGamma[0] -> EulerGamma
 * and StieltjesGamma[n, 1] -> StieltjesGamma[n]. It is the natural output of
 * Series[Zeta[x], {x, 1, n}] (and the Taylor expansion of Zeta about 0).
 *
 * Numerically (an inexact a, or N[...]) both forms come from arb's
 * acb_dirichlet_stieltjes through flint_num_stieltjes; a at 0 or a negative
 * integer, where zeta(s, a) has no Laurent expansion about s = 1, gives
 * ComplexInfinity. The order n must be a non-negative machine integer
 * (StieltjesGamma::intnm otherwise, as in the reference).
 *
 * Attributes: Listable, NHoldFirst, Protected. NHoldFirst keeps N from turning
 * the order into a real; numeric_stieltjes_n is N's route to the value. */
#include "stieltjesgamma.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "arithmetic.h"
#include "attr.h"
#include "flint_num_bridge.h"
#include "message.h"
#include "numeric.h"
#include "print.h"
#include "symtab.h"

/* An exact Integer >= 0 (the order). */
static bool sg_order_ok(const Expr* n) {
    return n->type == EXPR_INTEGER && n->data.integer >= 0;
}

static Expr* sg_intnm(Expr* res) {
    char* s = expr_to_string(res);
    mth_message("StieltjesGamma", "intnm",
                "Non-negative machine-sized integer expected at position 1 in %s.",
                s ? s : "");
    free(s);
    return NULL;
}

/* A number whose value is 0 or a negative integer: a pole of zeta(s, a). */
static bool sg_is_pole(const Expr* a) {
    if (a->type == EXPR_INTEGER) return a->data.integer <= 0;
    if (a->type == EXPR_REAL) return a->data.real <= 0 && a->data.real == floor(a->data.real);
#ifdef USE_MPFR
    /* an arbitrary-precision 0 or negative integer (-2.`30, 0``10) is the same
     * pole (audit StieltjesGamma#2) */
    if (a->type == EXPR_MPFR) return mpfr_integer_p(a->data.mpfr) && mpfr_sgn(a->data.mpfr) <= 0;
#endif
    return false;
}

/* gamma_n(a) for an inexact (or exact, under N) a, with a's precision; the
 * machine result is lowered to a Real. NULL when arb declines. */
static Expr* sg_numeric(const Expr* n, const Expr* a) {
    if (sg_is_pole(a)) return expr_new_symbol(SYM_ComplexInfinity);
    Expr* r = flint_num_stieltjes(n, a);
    if (!r) return NULL;
#ifdef USE_MPFR
    /* An arbitrary-precision a below 53 bits keeps its own precision
     * (StieltjesGamma[1, 1.5`10] is a 10-digit number, as Sin[1.5`10] is --
     * StieltjesGamma#1 of the 2026-11 audit); a machine a lowers. */
    if (a->type == EXPR_MPFR && mpfr_get_prec(a->data.mpfr) < 53 && r->type == EXPR_MPFR) {
        Expr* t = expr_new_mpfr_bits(mpfr_get_prec(a->data.mpfr));
        mpfr_set(t->data.mpfr, r->data.mpfr, MPFR_RNDN);
        expr_free(r);
        return t;
    }
    if (numeric_min_inexact_bits(a) <= 53) r = numeric_lower_mpfr_to_machine(r);
#endif
    return r;
}

Expr* builtin_stieltjesgamma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    Expr* n = res->data.function.args[0];
    bool n_numeric = n->type == EXPR_INTEGER || n->type == EXPR_BIGINT ||
                     n->type == EXPR_REAL ||
#ifdef USE_MPFR
                     n->type == EXPR_MPFR ||
#endif
                     (n->type == EXPR_FUNCTION &&
                      n->data.function.head->type == EXPR_SYMBOL &&
                      (n->data.function.head->data.symbol.name == SYM_Rational ||
                       n->data.function.head->data.symbol.name == SYM_Complex));
    if (n_numeric && !sg_order_ok(n)) return sg_intnm(res);
    if (!n_numeric) return NULL;

    if (argc == 1) {
        if (n->data.integer == 0) return expr_new_symbol(SYM_EulerGamma);
        return NULL; /* inert; N[...] goes through numeric_stieltjes_n */
    }
    Expr* a = res->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == 1) {
        Expr* one[1] = { expr_copy(n) };
        return expr_new_function(expr_new_symbol(SYM_StieltjesGamma), one, 1);
    }
    if (sg_is_pole(a) && a->type == EXPR_INTEGER) return expr_new_symbol(SYM_ComplexInfinity);
    if (arg_is_inexact(a)) return sg_numeric(n, a);
    return NULL;
}

Expr* numeric_stieltjes_n(const Expr* call, NumericSpec spec) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    size_t argc = call->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* n = call->data.function.args[0];
    if (!sg_order_ok(n)) return NULL;
    /* The precision rides on `a`: N[..., p] numericalizes it (a = 1 for the
     * one-argument form) at p digits, so the result comes back at p. */
    Expr* one = expr_new_integer(1);
    const Expr* a_exact = argc == 2 ? call->data.function.args[1] : one;
    Expr* a = numericalize(a_exact, spec);
    expr_free(one);
    if (!a) return NULL;
    Expr* r = arg_is_inexact(a) ? sg_numeric(n, a) : NULL;
    expr_free(a);
    return r;
}

void stieltjesgamma_init(void) {
    symtab_add_builtin("StieltjesGamma", builtin_stieltjesgamma);
    symtab_get_def("StieltjesGamma")->attributes |=
        (ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
