/* Mathilda -- shared helpers for ChebyshevT and ChebyshevU (see the header for
 * the split rationale). Every routine here is family-agnostic: argument-shape
 * predicates, the explicit-polynomial assembler, the MPFR result builder, the
 * output-precision rule, and the value at x = 0. The family-specific pieces --
 * the recurrence seed, the negative-order reflection, and the non-integer-order
 * closed forms -- live in chebyshev_t.c / chebyshev_u.c. */
#include "chebyshev_common.h"
#include "sym_names.h"

#include <math.h>
#include <stdlib.h>

#include "arithmetic.h"      /* make_rational, make_complex, is_complex */
#include "eval.h"            /* eval_and_free */
#include "numeric.h"         /* numeric_min_inexact_bits, get_approx_mpfr */

#ifdef USE_MPFR
#define PRND MPFR_RNDN
#endif

/* ------------------------------------------------------------------ */
/* Argument-shape predicates / extractors                             */
/* ------------------------------------------------------------------ */

bool cheb_exact_int(const Expr* e, int64_t* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) {
        return expr_mpz_get_int64(e->data.bigint, out);
    }
    return false;
}

bool cheb_is_int(const Expr* e, int64_t v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

bool cheb_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

bool cheb_is_inexact(const Expr* e) {
    if (cheb_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return cheb_is_inexact_leaf(re) || cheb_is_inexact_leaf(im);
    return false;
}

bool cheb_is_number_leaf(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rational ||
         e->data.function.head->data.symbol.name == SYM_Complex)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!cheb_is_number_leaf(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

bool cheb_leaf_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_REAL)    { *out = e->data.real; return true; }
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    return false;
}

bool cheb_machine_xy(const Expr* x_e, double* xr, double* xi, bool* real) {
    if (x_e && x_e->type == EXPR_REAL) {
        *xr = x_e->data.real; *xi = 0.0; *real = true; return true;
    }
    Expr *re, *im;
    if (is_complex((Expr*)x_e, &re, &im) &&
        cheb_leaf_double(re, xr) && cheb_leaf_double(im, xi)) {
        *real = (*xi == 0.0);
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/* Explicit integer-coefficient polynomial assembly                   */
/* ------------------------------------------------------------------ */

mpz_t* cheb_mpz_array_new(size_t len) {
    mpz_t* a = (mpz_t*)malloc(len * sizeof(mpz_t));
    if (!a) return NULL;
    for (size_t i = 0; i < len; i++) mpz_init(a[i]);   /* value 0 */
    return a;
}

void cheb_mpz_array_free(mpz_t* a, size_t len) {
    if (!a) return;
    for (size_t i = 0; i < len; i++) mpz_clear(a[i]);
    free(a);
}

Expr* cheb_poly_from_coeffs(mpz_t* c, size_t deg, const Expr* x) {
    Expr** terms = (Expr**)malloc((deg + 1) * sizeof(Expr*));
    if (!terms) return NULL;
    size_t count = 0;

    for (size_t i = 0; i <= deg; i++) {
        if (mpz_sgn(c[i]) == 0) continue;
        Expr* coeff = expr_bigint_normalize(expr_new_bigint_from_mpz(c[i]));
        Expr* term;
        if (i == 0) {
            term = coeff;
        } else if (i == 1) {
            term = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ coeff, expr_copy((Expr*)x) }, 2);
        } else {
            Expr* xi = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy((Expr*)x),
                                     expr_new_integer((int64_t)i) }, 2);
            term = expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ coeff, xi }, 2);
        }
        terms[count++] = term;
    }

    Expr* poly;
    if (count == 0)      poly = expr_new_integer(0);
    else if (count == 1) poly = terms[0];
    else                 poly = expr_new_function(expr_new_symbol(SYM_Plus),
                                                  terms, count);
    free(terms);
    return eval_and_free(poly);
}

/* ------------------------------------------------------------------ */
/* Numeric result helpers (MPFR)                                       */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
Expr* cheb_make_result(const mpfr_t re, const mpfr_t im,
                       mpfr_prec_t out_prec, bool force_real) {
    bool real = force_real || mpfr_zero_p(im);
    if (out_prec <= 53) {
        double dr = mpfr_get_d(re, PRND);
        if (real) {
            if (isinf(dr) && !mpfr_inf_p(re)) return expr_new_mpfr_copy(re);
            return expr_new_real(dr);
        }
        double di = mpfr_get_d(im, PRND);
        bool overflow = (isinf(dr) && !mpfr_inf_p(re)) ||
                        (isinf(di) && !mpfr_inf_p(im));
        if (overflow)
            return make_complex(expr_new_mpfr_copy(re), expr_new_mpfr_copy(im));
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
    Expr* rr = expr_new_mpfr_bits(out_prec);
    mpfr_set(rr->data.mpfr, re, PRND);
    if (real) return rr;
    Expr* ii = expr_new_mpfr_bits(out_prec);
    mpfr_set(ii->data.mpfr, im, PRND);
    return make_complex(rr, ii);
}

mpfr_prec_t cheb_out_prec(const Expr* n_e, const Expr* x_e) {
    long bn = cheb_is_inexact(n_e) ? numeric_min_inexact_bits(n_e) : 0;
    long bx = cheb_is_inexact(x_e) ? numeric_min_inexact_bits(x_e) : 0;
    long bits;
    if (bn > 0 && bx > 0) bits = (bn < bx) ? bn : bx;
    else                  bits = (bn > bx) ? bn : bx;
    if (bits < 53) bits = 53;
    return (mpfr_prec_t)bits;
}

bool cheb_inexact_int_order(const Expr* n_e, int64_t* k) {
    double v;
    if (n_e->type == EXPR_REAL) v = n_e->data.real;
    else if (n_e->type == EXPR_MPFR) v = mpfr_get_d(n_e->data.mpfr, PRND);
    else return false;
    if (v != floor(v)) return false;
    if (v < -(double)CHEB_NUMERIC_DEG_CAP || v > (double)CHEB_NUMERIC_DEG_CAP)
        return false;
    *k = (int64_t)v;
    return true;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Shared special value                                                */
/* ------------------------------------------------------------------ */

Expr* cheb_value_at_zero(const Expr* n_e) {
    Expr* arg = expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ make_rational(1, 2), expr_copy((Expr*)n_e),
                              expr_new_symbol(SYM_Pi) }, 3);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Cos),
                   (Expr*[]){ arg }, 1));
}
