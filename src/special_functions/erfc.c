/* Mathilda -- the complementary error function.
 *
 *   Erfc[z]   complementary error function   erfc(z) = 1 - erf(z)
 *
 * erfc is an entire function (no branch cuts). It is the complement of erf;
 * unlike erf it has no symmetry that simplifies Erfc[-x] (erfc(-x) = 2 -
 * erfc(x), which Mathilda leaves unexpanded). Evaluation is layered so each
 * kind of argument takes the cheapest, most accurate route:
 *
 *   exact special values    ->  1, 0, 2, DirectedInfinity[-+I], ...
 *   machine real            ->  libm   erfc
 *   arbitrary real          ->  MPFR   mpfr_erfc   (cancellation-free even for
 *                               large positive z, where 1 - erf(z) would lose
 *                               all significance)
 *   complex (any precision) ->  1 - erf(z), with erf(z) from the
 *                               cancellation-aware Maclaurin series (DLMF
 *                               7.6.2) evaluated in MPFR with guard bits; the
 *                               complement is formed at working precision
 *                               before rounding, so even machine-precision
 *                               complex results carry full accuracy. A
 *                               double-complex series is the USE_MPFR=0
 *                               fallback.
 *   everything else         ->  stays symbolic (return NULL)
 *
 * The erf series  erf(z) = (2/sqrt(pi)) e^-z^2 Sum_{n>=0} t_n,
 *   t_0 = z,  t_n = t_{n-1} (2 z^2)/(2n+1),
 * is convergent for every z. For complex z the partial sums can reach
 * magnitude ~e^|z|^2 before the e^-z^2 prefactor brings them back, so the
 * MPFR path adds |z|^2/ln2 guard bits to absorb that cancellation exactly.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "erfc.h"
#include "sym_names.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* is_rational, make_complex, is_complex */
#include "numeric.h"      /* numeric_min_inexact_bits */
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI / M_LN2 are POSIX, not C99 -- provide fallbacks (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* 2/sqrt(pi), the leading coefficient of erf. */
#define ERFC_TWO_OVER_SQRT_PI 1.1283791670955126

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool erfc_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool erfc_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Coerce an exact-or-real leaf to a double. Succeeds for Integer, Real,
 * BigInt, Rational and (for sizing only) MPFR. */
static bool erfc_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool erfc_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!erfc_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && erfc_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && erfc_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* Match Times[Complex[0, s], Infinity] (s = +-1): +-I Infinity. Returns
 * +1 for I Infinity, -1 for -I Infinity, 0 if not of this form. */
static int erfc_directed_imag_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    if (!erfc_is_symbol(e->data.function.head, "Times")) return 0;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    Expr *cplx = NULL, *inf = NULL;
    if (erfc_is_symbol(b, "Infinity")) { cplx = a; inf = b; }
    else if (erfc_is_symbol(a, "Infinity")) { cplx = b; inf = a; }
    if (!inf) return 0;
    /* cplx must be Complex[0, +-1]. */
    Expr *re, *im;
    if (!is_complex(cplx, &re, &im)) return 0;
    if (!(re->type == EXPR_INTEGER && re->data.integer == 0)) return 0;
    if (im->type == EXPR_INTEGER && im->data.integer == 1)  return 1;
    if (im->type == EXPR_INTEGER && im->data.integer == -1) return -1;
    return 0;
}

/* Build DirectedInfinity[Complex[0, s]] for s = +-1. */
static Expr* erfc_make_directed_imag_infinity(int s) {
    Expr* dir = make_complex(expr_new_integer(0), expr_new_integer(s));
    return expr_new_function(expr_new_symbol(SYM_DirectedInfinity), &dir, 1);
}

#ifndef USE_MPFR
/* ------------------------------------------------------------------ */
/* Machine double-complex series (fallback path for USE_MPFR=0)       */
/* ------------------------------------------------------------------ */

/* erf(z) for a machine double-complex z via the convergent series.
 * Returns false if the result is not finite (overflow for very large
 * |z|), in which case the caller leaves the call symbolic. The
 * complementary value erfc(z) = 1 - erf(z) is formed by the caller. */
static bool erfc_erf_machine_complex(double complex z, double complex* out) {
    double complex z2 = z * z;
    double complex two_z2 = 2.0 * z2;
    double complex t = z;          /* t_0 */
    double complex s = z;          /* partial sum */
    double z2abs = cabs(z2);
    for (int n = 1; n < 100000; n++) {
        t *= two_z2 / (2.0 * n + 1.0);
        s += t;
        /* Only stop once past the term peak (n > |z^2|) and converged. */
        if ((double)n > z2abs && cabs(t) <= 1e-17 * cabs(s)) break;
    }
    double complex val = ERFC_TWO_OVER_SQRT_PI * cexp(-z2) * s;
    if (!isfinite(creal(val)) || !isfinite(cimag(val))) return false;
    *out = val;
    return true;
}
#endif /* !USE_MPFR */

#ifdef USE_MPFR
/* ------------------------------------------------------------------ */
/* Complex-MPFR toolkit (pairs of mpfr_t; no MPC library available).  */
/* Mirrors the `ecx` toolkit in erf.c -- each op runs at an explicit  */
/* working precision and is alias-safe (inputs read into temporaries  */
/* before any output component is written).                           */
/* ------------------------------------------------------------------ */

#define ERND MPFR_RNDN

typedef struct { mpfr_t re, im; } ecx;

static void ecx_init(ecx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
static void ecx_clear(ecx* z)               { mpfr_clear(z->re);    mpfr_clear(z->im);    }
static void ecx_set(ecx* d, const ecx* s)   { mpfr_set(d->re, s->re, ERND); mpfr_set(d->im, s->im, ERND); }

static void ecx_add(ecx* out, const ecx* a, const ecx* b) {
    mpfr_add(out->re, a->re, b->re, ERND);
    mpfr_add(out->im, a->im, b->im, ERND);
}

/* out = a * b. */
static void ecx_mul(ecx* out, const ecx* a, const ecx* b, mpfr_prec_t p) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(p, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, ERND);
    mpfr_mul(bd, a->im, b->im, ERND);
    mpfr_mul(ad, a->re, b->im, ERND);
    mpfr_mul(bc, a->im, b->re, ERND);
    mpfr_sub(out->re, ac, bd, ERND);
    mpfr_add(out->im, ad, bc, ERND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = z / r for a real scalar r (in place safe). */
static void ecx_div_d(ecx* out, const ecx* z, double r) {
    mpfr_div_d(out->re, z->re, r, ERND);
    mpfr_div_d(out->im, z->im, r, ERND);
}

/* |z| into mag. */
static void ecx_abs(mpfr_t mag, const ecx* z) { mpfr_hypot(mag, z->re, z->im, ERND); }

/* out = exp(z) = e^re (cos im + i sin im). */
static void ecx_exp(ecx* out, const ecx* z, mpfr_prec_t p) {
    mpfr_t ea, c, s;
    mpfr_inits2(p, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, ERND);
    mpfr_sin_cos(s, c, z->im, ERND);
    mpfr_mul(out->re, ea, c, ERND);
    mpfr_mul(out->im, ea, s, ERND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* erf(z) for complex z into already-init'd `out` (precision wp).
 * Returns false if the series fails to converge within the cap. */
static bool ecx_erf_series(ecx* out, const ecx* z, mpfr_prec_t wp, double z2abs) {
    ecx z2, two_z2, t, s, tmp;
    ecx_init(&z2, wp); ecx_init(&two_z2, wp);
    ecx_init(&t, wp); ecx_init(&s, wp); ecx_init(&tmp, wp);

    ecx_mul(&z2, z, z, wp);                  /* z^2 */
    mpfr_mul_2ui(two_z2.re, z2.re, 1, ERND); /* 2 z^2 */
    mpfr_mul_2ui(two_z2.im, z2.im, 1, ERND);
    ecx_set(&t, z);                          /* t_0 = z */
    ecx_set(&s, z);                          /* S    = z */

    mpfr_t mag, smag, eps;
    mpfr_inits2(wp, mag, smag, eps, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, ERND);
    mpfr_div_2ui(eps, eps, (unsigned long)(wp > 8 ? wp - 4 : 1), ERND); /* 2^-(wp-4) */

    /* Cap generously: term peak is near n ~ |z|^2, plus precision-driven tail. */
    unsigned long cap = (unsigned long)(2.0 * z2abs) + (unsigned long)wp + 1000;
    bool ok = false;
    for (unsigned long n = 1; n <= cap; n++) {
        ecx_mul(&tmp, &t, &two_z2, wp);
        ecx_div_d(&t, &tmp, (double)(2 * n + 1));   /* t_n = t_{n-1} (2 z^2)/(2n+1) */
        ecx_add(&s, &s, &t);
        if ((double)n > z2abs) {                    /* past the peak */
            ecx_abs(mag, &t);
            ecx_abs(smag, &s);
            mpfr_mul(smag, smag, eps, ERND);
            if (mpfr_cmp(mag, smag) < 0) { ok = true; break; }
        }
    }

    if (ok) {
        /* out = (2/sqrt(pi)) e^-z^2 S. */
        ecx nz2, ez;
        ecx_init(&nz2, wp); ecx_init(&ez, wp);
        mpfr_neg(nz2.re, z2.re, ERND);
        mpfr_neg(nz2.im, z2.im, ERND);
        ecx_exp(&ez, &nz2, wp);                     /* e^-z^2 */
        ecx_mul(out, &ez, &s, wp);

        mpfr_t coeff;
        mpfr_init2(coeff, wp);
        mpfr_const_pi(coeff, ERND);
        mpfr_sqrt(coeff, coeff, ERND);
        mpfr_ui_div(coeff, 2, coeff, ERND);         /* 2/sqrt(pi) */
        mpfr_mul(out->re, out->re, coeff, ERND);
        mpfr_mul(out->im, out->im, coeff, ERND);
        mpfr_clear(coeff);
        ecx_clear(&nz2); ecx_clear(&ez);
    }

    mpfr_clears(mag, smag, eps, (mpfr_ptr)0);
    ecx_clear(&z2); ecx_clear(&two_z2);
    ecx_clear(&t); ecx_clear(&s); ecx_clear(&tmp);
    return ok;
}

/* Set an already-init2'd mpfr from an exact-or-real leaf. */
static bool erfc_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, ERND); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          ERND); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        ERND); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          ERND); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, ERND);
        mpfr_div_si(out, out, (long)d, ERND);
        return true;
    }
    return false;
}

/* Build a complex result from (re, im) at out_prec: machine precision
 * (<= 53) yields Real parts, higher yields MPFR parts. make_complex drops
 * a zero imaginary part to a bare real. */
static Expr* erfc_complex_result(const mpfr_t re, const mpfr_t im, mpfr_prec_t out_prec) {
    Expr *rr, *ii;
    if (out_prec <= 53) {
        rr = expr_new_real(mpfr_get_d(re, ERND));
        ii = expr_new_real(mpfr_get_d(im, ERND));
    } else {
        rr = expr_new_mpfr_bits(out_prec);
        ii = expr_new_mpfr_bits(out_prec);
        mpfr_set(rr->data.mpfr, re, ERND);
        mpfr_set(ii->data.mpfr, im, ERND);
    }
    return make_complex(rr, ii);
}

/* Evaluate erfc(z) = 1 - erf(z) for a numeric complex `arg` (Complex[..]
 * with at least one inexact part) at out_prec bits. erf(z) comes from the
 * MPFR series; the complement 1 - erf(z) is formed at working precision so
 * the guard bits also absorb the 1 - erf cancellation that dominates for
 * large positive Re(z). Returns the result, or NULL if the parts are not
 * numeric / the series diverges. */
static Expr* erfc_mpfr_complex(Expr* re, Expr* im, mpfr_prec_t out_prec) {
    double red = 0.0, imd = 0.0;
    (void)erfc_to_double(re, &red);
    (void)erfc_to_double(im, &imd);
    /* Guard bits absorb the ~e^|z|^2 partial-sum cancellation exactly. */
    double z2abs = red * red + imd * imd;
    long guard = 64 + (long)(z2abs / M_LN2);
    mpfr_prec_t wp = out_prec + (mpfr_prec_t)guard;

    ecx z, g;
    ecx_init(&z, wp); ecx_init(&g, wp);
    Expr* out = NULL;
    if (erfc_set_mpfr(z.re, re) && erfc_set_mpfr(z.im, im) &&
        ecx_erf_series(&g, &z, wp, z2abs) &&
        !mpfr_nan_p(g.re) && !mpfr_nan_p(g.im) &&
        !mpfr_inf_p(g.re) && !mpfr_inf_p(g.im)) {
        /* erfc(z) = 1 - erf(z), formed at working precision. */
        mpfr_ui_sub(g.re, 1, g.re, ERND);    /* 1 - Re(erf) */
        mpfr_neg(g.im, g.im, ERND);          /*   - Im(erf) */
        out = erfc_complex_result(g.re, g.im, out_prec);
    }
    ecx_clear(&z); ecx_clear(&g);
    return out;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Erfc[z]                                                            */
/* ------------------------------------------------------------------ */

static Expr* erfc_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_integer(1);                       /* Erfc[0] = 1 */
    if (erfc_is_symbol(arg, "Infinity"))      return expr_new_integer(0);
    if (erfc_is_neg_infinity(arg))            return expr_new_integer(2);
    if (erfc_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol(SYM_ComplexInfinity);
    if (erfc_is_symbol(arg, "Indeterminate"))   return expr_new_symbol(SYM_Indeterminate);
    {
        /* erfc(z) = 1 - erf(z), so the directed imaginary infinities are
         * negated relative to erf: Erfc[I Infinity] = DirectedInfinity[-I]. */
        int s = erfc_directed_imag_infinity(arg);
        if (s != 0) return erfc_make_directed_imag_infinity(-s);
    }

    /* 2. Machine real. */
    if (arg->type == EXPR_REAL) {
        double r = erfc(arg->data.real);
        return expr_new_real(r);
    }

#ifdef USE_MPFR
    /* 3. Arbitrary-precision real (mpfr_erfc is cancellation-free). */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        mpfr_erfc(out->data.mpfr, arg->data.mpfr, ERND);
        return out;
    }
#endif

    /* 4. Complex argument (Complex[..] with an inexact part). */
    {
        Expr *re, *im;
        if (is_complex(arg, &re, &im) && (erfc_is_inexact(re) || erfc_is_inexact(im))) {
#ifdef USE_MPFR
            /* Unified MPFR series: machine precision (out_prec = 53) yields
             * Real parts with full accuracy; higher precision yields MPFR. */
            long bits = numeric_min_inexact_bits(arg);
            mpfr_prec_t out_prec = (bits && bits > 53) ? (mpfr_prec_t)bits : 53;
            Expr* out = erfc_mpfr_complex(re, im, out_prec);
            if (out) return out;
#else
            double rr, ii;
            if (erfc_to_double(re, &rr) && erfc_to_double(im, &ii)) {
                double complex v;
                if (erfc_erf_machine_complex(rr + ii * I, &v)) {
                    double complex w = 1.0 - v;   /* erfc = 1 - erf */
                    if (cimag(w) == 0.0) return expr_new_real(creal(w));
                    return make_complex(expr_new_real(creal(w)),
                                        expr_new_real(cimag(w)));
                }
            }
#endif
        }
    }

    /* 5. No symmetry rewrite (erfc is not odd); stay symbolic. */
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* erfc_emit_argx(size_t argc) {
    fprintf(stderr,
            "Erfc::argx: Erfc called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_erfc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return erfc_one_arg(args[0]);
    return erfc_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void erfc_init(void) {
    symtab_add_builtin("Erfc", builtin_erfc);
    symtab_get_def("Erfc")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
