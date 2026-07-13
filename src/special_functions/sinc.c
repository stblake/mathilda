/* Mathilda -- the cardinal sine  Sinc[z] = Sin[z]/z  (Sinc[0] = 1).
 *
 * Sinc is entire and even, with a removable singularity at the origin. Each
 * kind of argument takes the cheapest route:
 *
 *   exact special values   ->  1 (at 0), 0 (at +-Infinity), Indeterminate
 *   machine real           ->  libm sin(x)/x
 *   arbitrary real (MPFR)  ->  mpfr_sin(x)/x at the input precision
 *   complex (any prec)     ->  sin(z)/z via the shared ncpx toolkit
 *   everything else        ->  stays symbolic (return NULL)
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "sinc.h"
#include "sym_names.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

#include "arithmetic.h"        /* is_rational, make_complex, is_complex */
#include "numeric.h"           /* numeric_min_inexact_bits */
#include "attr.h"
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#include "numeric_complex.h"   /* ncpx toolkit, numeric_mpfr_make_complex */
#define SRND MPFR_RNDN
#endif

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

static bool sinc_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool sinc_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

#ifndef USE_MPFR
static bool sinc_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}
#endif

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool sinc_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!sinc_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && sinc_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && sinc_is_symbol(a, "Infinity"))
        return true;
    return false;
}

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool sinc_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, SRND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          SRND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        SRND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          SRND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, SRND);
        mpfr_div_si(out, out, (long)d, SRND);
        return true;
    }
    return false;
}

/* Sinc(x) for a numeric real `arg` at out_prec bits (<= 53 yields a Real). */
static Expr* sinc_mpfr_real(const Expr* arg, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = out_prec + 16;
    mpfr_t x, s;
    mpfr_inits2(wp, x, s, (mpfr_ptr)0);
    if (!sinc_set_mpfr(x, arg)) { mpfr_clears(x, s, (mpfr_ptr)0); return NULL; }
    if (mpfr_zero_p(x)) {                       /* Sinc[0.] = 1 */
        mpfr_clears(x, s, (mpfr_ptr)0);
        return out_prec <= 53 ? expr_new_real(1.0)
                              : expr_new_mpfr_from_si(1, out_prec);
    }
    mpfr_sin(s, x, SRND);
    mpfr_div(s, s, x, SRND);                    /* sin(x)/x */
    Expr* out;
    if (out_prec <= 53) out = expr_new_real(mpfr_get_d(s, SRND));
    else { out = expr_new_mpfr_bits(out_prec); mpfr_set(out->data.mpfr, s, SRND); }
    mpfr_clears(x, s, (mpfr_ptr)0);
    return out;
}

/* Sinc(z) for a numeric complex `arg` at out_prec bits via ncpx: sin(z)/z. */
static Expr* sinc_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    mpfr_prec_t wp = out_prec + 24;
    ncpx z, sn;
    ncpx_init(&z, wp); ncpx_init(&sn, wp);
    Expr* out = NULL;
    if (sinc_set_mpfr(z.re, re) && sinc_set_mpfr(z.im, im)) {
        ncpx_sin(&sn, &z, wp);                  /* sin(z)      */
        ncpx_div(&sn, &sn, &z, wp);             /* sin(z)/z    */
        if (!mpfr_nan_p(sn.re) && !mpfr_nan_p(sn.im)) {
            if (out_prec <= 53) {
                Expr* rr = expr_new_real(mpfr_get_d(sn.re, SRND));
                Expr* ii = expr_new_real(mpfr_get_d(sn.im, SRND));
                out = mpfr_zero_p(sn.im) ? (expr_free(ii), rr) : make_complex(rr, ii);
            } else {
                out = numeric_mpfr_make_complex(sn.re, sn.im);
            }
        }
    }
    ncpx_clear(&z); ncpx_clear(&sn);
    return out;
}
#else
/* Machine-precision fallback (USE_MPFR=0 builds). */
static bool sinc_machine_complex(double complex z, double complex* out) {
    double complex v = csin(z) / z;
    if (!isfinite(creal(v)) || !isfinite(cimag(v))) return false;
    *out = v;
    return true;
}
#endif

/* ------------------------------------------------------------------ */
/* Sinc[z]                                                            */
/* ------------------------------------------------------------------ */

static Expr* sinc_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(1);                          /* Sinc[0] = 1 */
    if (sinc_is_symbol(arg, "Infinity"))       return expr_new_integer(0);
    if (sinc_is_neg_infinity(arg))             return expr_new_integer(0);
    if (sinc_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol(SYM_Indeterminate);
    if (sinc_is_symbol(arg, "Indeterminate"))   return expr_new_symbol(SYM_Indeterminate);

    /* 2. Numeric real (machine and arbitrary precision). */
    if (arg->type == EXPR_REAL) {
#ifdef USE_MPFR
        return sinc_mpfr_real(arg, 53);
#else
        double x = arg->data.real;
        return expr_new_real(x == 0.0 ? 1.0 : sin(x) / x);
#endif
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        return sinc_mpfr_real(arg, prec);
    }
#endif

    /* 3. Numeric complex (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (sinc_is_inexact(re) || sinc_is_inexact(im))) {
#ifdef USE_MPFR
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = sinc_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (sinc_to_double(re, &rr) && sinc_to_double(im, &ii)) {
                double complex v;
                if (sinc_machine_complex(rr + ii * I, &v)) {
                    if (cimag(v) == 0.0) return expr_new_real(creal(v));
                    return make_complex(expr_new_real(creal(v)),
                                        expr_new_real(cimag(v)));
                }
            }
#endif
        }
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

static Expr* sinc_emit_argx(size_t argc) {
    fprintf(stderr,
            "Sinc::argx: Sinc called with %zu arguments; 1 argument is expected.\n",
            argc);
    return NULL;
}

Expr* builtin_sinc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return sinc_one_arg(args[0]);
    return sinc_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void sinc_init(void) {
    symtab_add_builtin("Sinc", builtin_sinc);
    symtab_get_def("Sinc")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
