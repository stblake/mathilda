/* Mathilda -- ProductLog, the Lambert W function.
 *
 *   ProductLog[z]     principal branch W_0(z): the solution w of z = w e^w.
 *   ProductLog[k, z]  the k-th branch W_k(z), k any integer (k == 0 principal).
 *
 * Evaluation is layered so each kind of argument takes the cheapest accurate
 * route:
 *
 *   exact special values   ->  ProductLog[0] = 0, ProductLog[E] = 1,
 *                              ProductLog[-1/E] = -1, ProductLog[-Pi/2] = I Pi/2,
 *                              ProductLog[+-Infinity/ComplexInfinity] = Infinity,
 *                              ProductLog[k, 0] = -Infinity  (k != 0)
 *   numeric (real/complex)  ->  unified complex-MPFR Halley core; the result is
 *                              a Real / MPFR leaf when it is real-valued for the
 *                              chosen branch, otherwise Complex[..]
 *   everything else        ->  stays symbolic (return NULL)
 *
 * The numeric core (pl_core) builds on the shared `ncpx` complex-MPFR toolkit
 * (numeric_complex.h). It seeds an initial approximation by region --
 *   - branch-point series in p = sqrt(2(e z + 1)) near z = -1/e (branches 0,-1);
 *   - the Maclaurin seed z(1 - z + 3/2 z^2) for the principal branch near 0;
 *   - otherwise the asymptotic L1 - L2 + L2/L1 with L1 = log z + 2 pi i k,
 *     L2 = log L1
 * -- and refines it with Halley's cubically-convergent iteration
 *   w <- w - (w e^w - z) / (e^w (w+1) - (w+2)(w e^w - z)/(2w+2))   (Corless 1996).
 * A real seed keeps the whole iteration exactly real (every ncpx op preserves a
 * zero imaginary part), so real-valued branches return a real leaf with no
 * imaginary noise. Working precision carries guard bits above the requested
 * output precision.
 *
 * D[ProductLog[z], z] = ProductLog[z] / (z (1 + ProductLog[z]))  (calculus/deriv.c).
 * Series at 0, at the branch point -1/E, and at Infinity live in calculus/series.c.
 *
 * Attributes: Listable, NumericFunction, Protected, ReadProtected.
 */
#include "productlog.h"
#include "sym_names.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* is_rational, make_complex, is_complex, make_rational */
#include "numeric.h"      /* numeric_min_inexact_bits, get_approx_mpfr */
#include "numeric_complex.h" /* ncpx toolkit */
#include "attr.h"
#include "eval.h"         /* eval_and_free */
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI / M_E are POSIX, not C99 -- provide fallbacks (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool pl_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool pl_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries an inexact numeric leaf (also looking inside Complex). */
static bool pl_is_inexact(const Expr* e) {
    if (pl_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return pl_is_inexact_leaf(re) || pl_is_inexact_leaf(im);
    return false;
}

/* Construct the (already-evaluated) expression `Times[-1, Power[E, -1]]` = -1/E. */
static Expr* pl_make_neg_inv_e(void) {
    Expr* inv_e = expr_new_function(expr_new_symbol(SYM_Power),
                      (Expr*[]){ expr_new_symbol(SYM_E), expr_new_integer(-1) }, 2);
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ expr_new_integer(-1), inv_e }, 2);
    return eval_and_free(t);
}

/* Construct the (already-evaluated) expression -Pi/2 = Times[Rational[-1,2], Pi]. */
static Expr* pl_make_neg_half_pi(void) {
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ make_rational(-1, 2), expr_new_symbol(SYM_Pi) }, 2);
    return eval_and_free(t);
}

/* Construct the (already-evaluated) expression I Pi/2 = Times[Rational[1,2], I, Pi]. */
static Expr* pl_make_i_half_pi(void) {
    Expr* i = make_complex(expr_new_integer(0), expr_new_integer(1));
    Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                  (Expr*[]){ make_rational(1, 2), i, expr_new_symbol(SYM_Pi) }, 3);
    return eval_and_free(t);
}

#ifdef USE_MPFR
#define PRND MPFR_RNDN

/* Working precision: a comfortable guard above the requested output precision. */
static mpfr_prec_t pl_wp(mpfr_prec_t out_prec) {
    mpfr_prec_t wp = (out_prec < 64) ? 64 : out_prec;
    return wp + 24;
}

/* Build a numeric Expr from a real mpfr value at out_prec bits. */
static Expr* pl_real_result(const mpfr_t v, mpfr_prec_t out_prec) {
    if (out_prec <= 53) {
        double d = mpfr_get_d(v, PRND);
        if (isinf(d) && !mpfr_inf_p(v)) return expr_new_mpfr_copy(v);
        return expr_new_real(d);
    }
    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_set(out->data.mpfr, v, PRND);
    return out;
}

/* Build a result from an (re, im) mpfr pair: a real leaf if im is exactly zero,
 * otherwise Complex[..]. Promotes to MPFR on machine-precision overflow. */
static Expr* pl_make_result(const ncpx* w, mpfr_prec_t out_prec, bool force_real) {
    if (force_real || mpfr_zero_p(w->im)) return pl_real_result(w->re, out_prec);
    Expr *rr, *ii;
    if (out_prec <= 53) {
        double dr = mpfr_get_d(w->re, PRND), di = mpfr_get_d(w->im, PRND);
        bool overflow = (isinf(dr) && !mpfr_inf_p(w->re)) ||
                        (isinf(di) && !mpfr_inf_p(w->im));
        if (overflow) {
            rr = expr_new_mpfr_copy(w->re);
            ii = expr_new_mpfr_copy(w->im);
            return make_complex(rr, ii);
        }
        rr = expr_new_real(dr);
        ii = expr_new_real(di);
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, w->re, PRND);
        mpfr_set(ii->data.mpfr, w->im, PRND);
    }
    return make_complex(rr, ii);
}

/* ------------------------------------------------------------------ */
/* Numeric core                                                       */
/* ------------------------------------------------------------------ */

/* Set the ncpx `c` to the integer real value v. */
static void ncpx_set_si(ncpx* c, long v) {
    mpfr_set_si(c->re, v, PRND);
    mpfr_set_si(c->im, 0, PRND);
}

/* Initial approximation for W_k(z), written into `w` (init2'd at wp). */
static void pl_initial(ncpx* w, long k, const ncpx* z, mpfr_prec_t wp) {
    ncpx ep1, p, tmp, tmp2, L1, L2;
    ncpx_init(&ep1, wp); ncpx_init(&p, wp); ncpx_init(&tmp, wp);
    ncpx_init(&tmp2, wp); ncpx_init(&L1, wp); ncpx_init(&L2, wp);

    /* ep1 = e*z + 1  (zero exactly at the branch point z = -1/e). */
    mpfr_t e_const; mpfr_init2(e_const, wp);
    mpfr_set_si(e_const, 1, PRND);
    mpfr_exp(e_const, e_const, PRND);   /* e = exp(1) */
    ncpx_scale(&ep1, z, e_const);
    mpfr_add_ui(ep1.re, ep1.re, 1, PRND);

    mpfr_t mag_bp, mag_z, mag_L1; mpfr_init2(mag_bp, wp);
    mpfr_init2(mag_z, wp); mpfr_init2(mag_L1, wp);
    ncpx_abs(mag_bp, &ep1);
    ncpx_abs(mag_z, z);

    bool done = false;

    /* Branch-point region (branches 0 and -1 only): p = sqrt(2(e z + 1)),
     * w = -1 + p - p^2/3 + 11 p^3/72, with p negated for branch -1. */
    if ((k == 0 || k == -1) && mpfr_cmp_d(mag_bp, 0.3) < 0) {
        ncpx_set_si(&tmp, 2);
        ncpx_mul(&tmp, &tmp, &ep1, wp);      /* 2(e z + 1) */
        ncpx_sqrt(&p, &tmp, wp);             /* p = sqrt(...) */
        if (k == -1) ncpx_neg(&p, &p);
        /* w = -1 + p - p^2/3 + 11/72 p^3 */
        ncpx_set_si(w, -1);
        ncpx_add(w, w, &p);
        ncpx_mul(&tmp, &p, &p, wp);          /* p^2 */
        mpfr_set_d(e_const, -1.0 / 3.0, PRND);
        ncpx_scale(&tmp2, &tmp, e_const);    /* -p^2/3 */
        ncpx_add(w, w, &tmp2);
        ncpx_mul(&tmp, &tmp, &p, wp);        /* p^3 */
        mpfr_set_d(e_const, 11.0 / 72.0, PRND);
        ncpx_scale(&tmp2, &tmp, e_const);
        ncpx_add(w, w, &tmp2);
        done = true;
    }

    /* Principal branch near 0: Maclaurin seed w = z(1 - z + 3/2 z^2). */
    if (!done && k == 0 && mpfr_cmp_d(mag_z, 0.35) < 0) {
        ncpx_mul(&tmp, z, z, wp);            /* z^2 */
        mpfr_set_d(e_const, 1.5, PRND);
        ncpx_scale(w, &tmp, e_const);        /* 3/2 z^2 */
        ncpx_sub(w, w, z);                   /* 3/2 z^2 - z */
        mpfr_set_si(e_const, 1, PRND);       /* (1) */
        /* w = z*(1 - z + 3/2 z^2) = z + z*(3/2 z^2 - z) */
        ncpx_mul(w, w, z, wp);               /* z*(3/2 z^2 - z) */
        ncpx_add(w, w, z);                   /* + z */
        done = true;
    }

    /* Asymptotic seed: L1 = log z + 2 pi i k, L2 = log L1, w = L1 - L2 + L2/L1. */
    if (!done) {
        ncpx_log(&L1, z, wp);
        if (k != 0) {
            mpfr_t two_pi_k; mpfr_init2(two_pi_k, wp);
            mpfr_const_pi(two_pi_k, PRND);
            mpfr_mul_si(two_pi_k, two_pi_k, 2 * k, PRND);
            mpfr_add(L1.im, L1.im, two_pi_k, PRND);
            mpfr_clear(two_pi_k);
        }
        ncpx_abs(mag_L1, &L1);
        if (k == 0 && mpfr_cmp_d(mag_L1, 1.0) < 0) {
            /* z near 1: L1 ~ 0 makes log L1 blow up. Seed a benign real value. */
            mpfr_set_d(w->re, 0.5, PRND);
            mpfr_set_si(w->im, 0, PRND);
        } else {
            ncpx_log(&L2, &L1, wp);          /* L2 = log L1 */
            ncpx_div(&tmp, &L2, &L1, wp);    /* L2/L1 */
            ncpx_sub(w, &L1, &L2);           /* L1 - L2 */
            ncpx_add(w, w, &tmp);            /* + L2/L1 */
        }
    }

    mpfr_clear(e_const); mpfr_clear(mag_bp); mpfr_clear(mag_z); mpfr_clear(mag_L1);
    ncpx_clear(&ep1); ncpx_clear(&p); ncpx_clear(&tmp);
    ncpx_clear(&tmp2); ncpx_clear(&L1); ncpx_clear(&L2);
}

/* Halley refinement of W_k(z): converge `w` (already seeded, init2'd at wp). */
static void pl_core(ncpx* w, const ncpx* z, mpfr_prec_t wp) {
    ncpx ew, t, f, wp1, num, denom, dw, one;
    ncpx_init(&ew, wp); ncpx_init(&t, wp); ncpx_init(&f, wp);
    ncpx_init(&wp1, wp); ncpx_init(&num, wp); ncpx_init(&denom, wp);
    ncpx_init(&dw, wp); ncpx_init(&one, wp);
    ncpx_set_si(&one, 1);

    mpfr_t adw, aw, eps; mpfr_init2(adw, wp); mpfr_init2(aw, wp); mpfr_init2(eps, wp);
    /* eps = 2^-(wp-6): relative-step tolerance. */
    mpfr_set_ui(eps, 1, PRND);
    mpfr_div_2si(eps, eps, (long)wp - 6, PRND);

    for (int iter = 0; iter < 100; iter++) {
        ncpx_exp(&ew, w, wp);                /* e^w */
        ncpx_mul(&t, w, &ew, wp);            /* w e^w */
        ncpx_sub(&f, &t, z);                 /* f = w e^w - z */

        ncpx_add(&wp1, w, &one);             /* w + 1 */
        ncpx_mul(&denom, &ew, &wp1, wp);     /* e^w (w+1) */
        /* corr = (w+2) f / (2(w+1)) */
        ncpx_add(&num, &wp1, &one);          /* w + 2 */
        ncpx_mul(&num, &num, &f, wp);        /* (w+2) f */
        ncpx_add(&t, &wp1, &wp1);            /* 2(w+1) */
        ncpx_div(&num, &num, &t, wp);        /* corr */
        ncpx_sub(&denom, &denom, &num);      /* e^w(w+1) - corr */

        ncpx_div(&dw, &f, &denom, wp);       /* Halley step */
        ncpx_sub(w, w, &dw);                 /* w <- w - dw */

        ncpx_abs(adw, &dw);
        ncpx_abs(aw, w);
        mpfr_add_ui(aw, aw, 1, PRND);        /* 1 + |w| (mixed abs/rel tol) */
        mpfr_mul(aw, aw, eps, PRND);
        if (mpfr_cmp(adw, aw) <= 0) {
            /* One extra polishing step, then stop. */
            ncpx_exp(&ew, w, wp);
            ncpx_mul(&t, w, &ew, wp);
            ncpx_sub(&f, &t, z);
            ncpx_add(&wp1, w, &one);
            ncpx_mul(&denom, &ew, &wp1, wp);
            ncpx_add(&num, &wp1, &one);
            ncpx_mul(&num, &num, &f, wp);
            ncpx_add(&t, &wp1, &wp1);
            ncpx_div(&num, &num, &t, wp);
            ncpx_sub(&denom, &denom, &num);
            ncpx_div(&dw, &f, &denom, wp);
            ncpx_sub(w, w, &dw);
            break;
        }
    }

    mpfr_clear(adw); mpfr_clear(aw); mpfr_clear(eps);
    ncpx_clear(&ew); ncpx_clear(&t); ncpx_clear(&f); ncpx_clear(&wp1);
    ncpx_clear(&num); ncpx_clear(&denom); ncpx_clear(&dw); ncpx_clear(&one);
}

/* True iff W_k(x) is real for the real argument x (so a real leaf is returned).
 * Principal branch: x >= -1/e. Branch -1: -1/e <= x < 0. */
static bool pl_real_valued(long k, const mpfr_t x, mpfr_prec_t wp) {
    mpfr_t neg_inv_e; mpfr_init2(neg_inv_e, wp);
    mpfr_set_si(neg_inv_e, 1, PRND);
    mpfr_exp(neg_inv_e, neg_inv_e, PRND);            /* e */
    mpfr_si_div(neg_inv_e, -1, neg_inv_e, PRND);     /* -1/e */
    bool result = false;
    if (k == 0) {
        result = (mpfr_cmp(x, neg_inv_e) >= 0);
    } else if (k == -1) {
        result = (mpfr_cmp(x, neg_inv_e) >= 0) && (mpfr_sgn(x) < 0);
    }
    mpfr_clear(neg_inv_e);
    return result;
}

/* Numeric evaluation of W_k(z_arg) at out_prec bits. Returns NULL if z_arg is
 * not numeric (symbolic), letting the caller leave the call unevaluated. */
static Expr* pl_eval_numeric(long k, const Expr* z_arg, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = pl_wp(out_prec);
    ncpx z, w;
    ncpx_init(&z, wp); ncpx_init(&w, wp);
    bool inexact = false;
    if (!get_approx_mpfr(z_arg, z.re, z.im, &inexact)) {
        ncpx_clear(&z); ncpx_clear(&w);
        return NULL;
    }
    bool z_is_real = mpfr_zero_p(z.im);
    bool force_real = z_is_real && pl_real_valued(k, z.re, wp);

    pl_initial(&w, k, &z, wp);
    pl_core(&w, &z, wp);

    Expr* out = NULL;
    if (!mpfr_nan_p(w.re) && !mpfr_nan_p(w.im))
        out = pl_make_result(&w, out_prec, force_real);

    ncpx_clear(&z); ncpx_clear(&w);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Dispatch                                                           */
/* ------------------------------------------------------------------ */

/* Evaluate W_k(z): exact special values, then numeric, else NULL (symbolic). */
static Expr* productlog_eval(long k, Expr* z) {
    /* --- z == 0 (exact). --- */
    bool z_int_zero = (z->type == EXPR_INTEGER && z->data.integer == 0);
    bool z_real_zero = (z->type == EXPR_REAL && z->data.real == 0.0);
    if (z_int_zero || z_real_zero) {
        if (k == 0) return z_real_zero ? expr_new_real(0.0) : expr_new_integer(0);
        /* W_k(0) = -Infinity for k != 0. */
        return expr_new_function(expr_new_symbol(SYM_Times),
                   (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
    }

    /* --- Infinity / ComplexInfinity -> Infinity (any branch). --- */
    if (pl_is_symbol(z, "Infinity") || pl_is_symbol(z, "ComplexInfinity"))
        return expr_new_symbol(SYM_Infinity);

    /* --- ProductLog[E] = 1 (principal branch). --- */
    if (k == 0 && pl_is_symbol(z, "E"))
        return expr_new_integer(1);

    /* --- ProductLog[-1/E] = -1 (branches 0 and -1 meet at the branch point). --- */
    if (k == 0 || k == -1) {
        Expr* cand = pl_make_neg_inv_e();
        bool match = expr_eq(cand, z);
        expr_free(cand);
        if (match) return expr_new_integer(-1);
    }

    /* --- ProductLog[-Pi/2] = I Pi/2 (principal branch). --- */
    if (k == 0) {
        Expr* cand = pl_make_neg_half_pi();
        bool match = expr_eq(cand, z);
        expr_free(cand);
        if (match) return pl_make_i_half_pi();
    }

#ifdef USE_MPFR
    /* --- Numeric (real or complex, machine or arbitrary precision). --- */
    if (pl_is_inexact(z)) {
        long bits = numeric_min_inexact_bits(z);
        mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
        Expr* out = pl_eval_numeric(k, z, out_prec);
        if (out) return out;
    }
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-style diagnostic for a wrong argument count; returns NULL so the
 * evaluator leaves the call unevaluated. */
static Expr* productlog_emit_argt(size_t argc) {
    fprintf(stderr,
            "ProductLog::argt: ProductLog called with %zu arguments; "
            "1 or 2 arguments are expected.\n",
            argc);
    return NULL;
}

Expr* builtin_productlog(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 1)
        return productlog_eval(0, res->data.function.args[0]);

    if (argc == 2) {
        Expr* kexpr = res->data.function.args[0];
        Expr* z = res->data.function.args[1];
        /* The branch index must be an explicit (machine-range) integer. */
        if (kexpr->type != EXPR_INTEGER) return NULL;
        int64_t k = kexpr->data.integer;
        if (k < LONG_MIN || k > LONG_MAX) return NULL;
        return productlog_eval((long)k, z);
    }

    return productlog_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void productlog_init(void) {
    symtab_add_builtin("ProductLog", builtin_productlog);
    symtab_get_def("ProductLog")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED | ATTR_READPROTECTED);
    /* Docstring lives in info.c (info_init). */
}
