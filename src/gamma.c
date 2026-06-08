/* Mathilda -- the Gamma function family.
 *
 *   Gamma[z]          Euler gamma function   Gamma(z) = Int_0^Inf t^(z-1) e^-t dt
 *   Gamma[a, z]       upper incomplete gamma Gamma(a,z) = Int_z^Inf t^(a-1) e^-t dt
 *   Gamma[a, z0, z1]  generalized incomplete = Gamma(a,z0) - Gamma(a,z1)
 *
 * Evaluation is layered so each kind of argument takes the cheapest exact
 * or fastest numeric route available:
 *
 *   exact integer / half-integer  ->  (z-1)! via the Factorial machinery
 *                                      (exact, BigInt, or rational*Sqrt[Pi])
 *   machine real        ->  libm   tgamma
 *   machine complex     ->  Lanczos approximation (double complex)
 *   arbitrary real      ->  MPFR   mpfr_gamma / mpfr_gamma_inc
 *   everything else     ->  stays symbolic (return NULL)
 *
 * Arbitrary-precision *complex* gamma is deliberately left symbolic: a
 * fixed-coefficient Lanczos series only carries ~15 correct digits, so
 * emitting it as an MPFR value would advertise a precision it does not
 * have. Reporting the input unevaluated is the honest behaviour.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "gamma.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "arithmetic.h"   /* is_rational, make_rational, is_complex, make_complex */
#include "attr.h"
#include "eval.h"          /* eval_and_free */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Small numeric-coercion helpers                                     */
/* ------------------------------------------------------------------ */

/* Coerce an exact-or-real leaf to a double. Succeeds for Integer, Real,
 * BigInt and Rational; fails (returns false) for symbols, complex values,
 * MPFR values and anything else. */
static bool gamma_to_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { return false; } /* handled by the MPFR path */
#endif
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    int64_t n, d;
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True if `e` is an inexact numeric leaf (Real or MPFR): its presence is
 * what turns an otherwise-symbolic call into a numeric one. */
static bool gamma_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` is exactly the symbol `name`. */
static bool gamma_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* True if `e` is -Infinity, represented as Times[-1, Infinity]. */
static bool gamma_is_neg_infinity(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return false;
    if (!gamma_is_symbol(e->data.function.head, "Times")) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (a->type == EXPR_INTEGER && a->data.integer == -1 && gamma_is_symbol(b, "Infinity"))
        return true;
    if (b->type == EXPR_INTEGER && b->data.integer == -1 && gamma_is_symbol(a, "Infinity"))
        return true;
    return false;
}

/* ------------------------------------------------------------------ */
/* Exact path: Gamma[z] = (z-1)! for integer / half-integer z         */
/* ------------------------------------------------------------------ */

/* Reuse the existing Factorial builtin (exact integers, BigInt, and the
 * Sqrt[Pi] half-integer rationals) by evaluating Factorial[z-1]. Only
 * called for is_rational(z) with denominator 1 or 2, the cases Factorial
 * actually closes. Returns NULL (and frees its scratch) if the Factorial
 * unexpectedly fails to reduce, leaving Gamma symbolic. */
static Expr* gamma_exact_via_factorial(int64_t n, int64_t d) {
    /* z = n/d, so z - 1 = (n - d)/d. */
    Expr* zm1 = make_rational(n - d, d);
    if (!zm1) return NULL;
    Expr* call = expr_new_function(expr_new_symbol("Factorial"), &zm1, 1);
    Expr* out = eval_and_free(call);
    /* If it came back still headed by Factorial, it did not reduce. */
    if (out && out->type == EXPR_FUNCTION &&
        gamma_is_symbol(out->data.function.head, "Factorial")) {
        expr_free(out);
        return NULL;
    }
    return out;
}

/* ------------------------------------------------------------------ */
/* Machine complex path: Lanczos approximation                        */
/* ------------------------------------------------------------------ */

/* Lanczos approximation (g = 7, n = 9), accurate to ~15 significant
 * digits across the complex plane. Uses the reflection formula for the
 * left half-plane where the series is ill-conditioned. */
static double complex gamma_lanczos(double complex z) {
    static const double g = 7.0;
    static const double c[9] = {
         0.99999999999980993,    676.5203681218851,    -1259.1392167224028,
       771.32342877765313,      -176.61502916214059,      12.507343278686905,
        -0.13857109526572012,      9.9843695780195716e-6,   1.5056327351493116e-7
    };
    if (creal(z) < 0.5) {
        /* Gamma(z) Gamma(1-z) = pi / sin(pi z). */
        return M_PI / (csin(M_PI * z) * gamma_lanczos(1.0 - z));
    }
    z -= 1.0;
    double complex x = c[0];
    for (int i = 1; i < 9; i++) x += c[i] / (z + (double)i);
    double complex t = z + g + 0.5;
    return csqrt(2.0 * M_PI) * cpow(t, z + 0.5) * cexp(-t) * x;
}

/* ------------------------------------------------------------------ */
/* MPFR helpers                                                       */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from an exact-or-real leaf. Succeeds for
 * Integer, Real, BigInt, MPFR and Rational; fails for complex / symbolic. */
static bool gamma_set_mpfr(mpfr_t out, const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          MPFR_RNDN); return true;
        default: break;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, MPFR_RNDN);
        mpfr_div_si(out, out, (long)d, MPFR_RNDN);
        return true;
    }
    return false;
}

/* Working precision for a numeric Gamma: the largest precision among any
 * MPFR operands, else 53 bits (machine). */
static mpfr_prec_t gamma_work_prec(const Expr* a, const Expr* b) {
    mpfr_prec_t p = 53;
    if (a && a->type == EXPR_MPFR && mpfr_get_prec(a->data.mpfr) > p) p = mpfr_get_prec(a->data.mpfr);
    if (b && b->type == EXPR_MPFR && mpfr_get_prec(b->data.mpfr) > p) p = mpfr_get_prec(b->data.mpfr);
    return p;
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* Gamma[z]  (the Euler gamma function)                               */
/* ------------------------------------------------------------------ */

static Expr* gamma_one_arg(Expr* arg) {
    /* 1. Exact integer / half-integer -> (z-1)!. */
    int64_t n, d;
    if (is_rational(arg, &n, &d) && (d == 1 || d == 2)) {
        Expr* exact = gamma_exact_via_factorial(n, d);
        if (exact) return exact;
    }

    /* 2. Symbolic infinities. */
    if (gamma_is_symbol(arg, "Infinity"))        return expr_new_symbol("Infinity");
    if (gamma_is_symbol(arg, "ComplexInfinity")) return expr_new_symbol("ComplexInfinity");
    if (gamma_is_symbol(arg, "Indeterminate"))   return expr_new_symbol("Indeterminate");
    if (gamma_is_neg_infinity(arg))              return expr_new_symbol("Indeterminate");

    /* 3. Machine real. */
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        double r = tgamma(v);
        if (isnan(r)) return expr_new_symbol("ComplexInfinity"); /* pole at <=0 integer */
        if (isinf(r)) {
            if (v <= 0.0) return expr_new_symbol("ComplexInfinity");
            return NULL; /* overflow for large positive z: stay symbolic */
        }
        return expr_new_real(r);
    }

#ifdef USE_MPFR
    /* 4. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        mpfr_gamma(out->data.mpfr, arg->data.mpfr, MPFR_RNDN);
        if (mpfr_inf_p(out->data.mpfr)) { expr_free(out); return expr_new_symbol("ComplexInfinity"); }
        if (mpfr_nan_p(out->data.mpfr)) { expr_free(out); return expr_new_symbol("ComplexInfinity"); }
        return out;
    }
#endif

    /* 5. Machine complex (Complex[..] with an inexact part, both parts
     *    representable as double). MPFR-complex stays symbolic. */
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        double rr, ii;
        bool inexact = gamma_is_inexact(re) || gamma_is_inexact(im);
        if (inexact && gamma_to_double(re, &rr) && gamma_to_double(im, &ii)) {
            double complex g = gamma_lanczos(rr + ii * I);
            if (cimag(g) == 0.0) return expr_new_real(creal(g));
            return make_complex(expr_new_real(creal(g)), expr_new_real(cimag(g)));
        }
    }

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Gamma[a, z]  (upper incomplete gamma)                              */
/* ------------------------------------------------------------------ */

static Expr* gamma_two_arg(Expr* a, Expr* z) {
    /* 1. Exact rewrites that hold for any z / a. */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        /* Gamma[a, 0] = Gamma[a]. */
        Expr* ga = expr_copy(a);
        return eval_and_free(expr_new_function(expr_new_symbol("Gamma"), &ga, 1));
    }
    if (gamma_is_symbol(z, "Infinity")) return expr_new_integer(0); /* Gamma[a, Inf] = 0 */
    if (a->type == EXPR_INTEGER && a->data.integer == 1) {
        /* Gamma[1, z] = E^-z = Exp[-z]. */
        Expr* nz   = eval_and_free(expr_new_function(
                        expr_new_symbol("Times"),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(z) }, 2));
        return eval_and_free(expr_new_function(expr_new_symbol("Exp"), &nz, 1));
    }

#ifdef USE_MPFR
    /* 2. Numeric real path -- needs at least one inexact operand so we do
     *    not silently turn exact Gamma[2, 3] into a float. */
    if (gamma_is_inexact(a) || gamma_is_inexact(z)) {
        mpfr_prec_t prec = gamma_work_prec(a, z);
        mpfr_t av, zv, rv;
        mpfr_init2(av, prec); mpfr_init2(zv, prec); mpfr_init2(rv, prec);
        bool ok = gamma_set_mpfr(av, a) && gamma_set_mpfr(zv, z);
        Expr* out = NULL;
        if (ok) {
            mpfr_gamma_inc(rv, av, zv, MPFR_RNDN);
            if (mpfr_nan_p(rv)) {
                out = NULL;
            } else if (a->type == EXPR_MPFR || z->type == EXPR_MPFR) {
                out = expr_new_mpfr_copy(rv);            /* arbitrary precision */
            } else {
                out = expr_new_real(mpfr_get_d(rv, MPFR_RNDN)); /* machine */
            }
        }
        mpfr_clear(av); mpfr_clear(zv); mpfr_clear(rv);
        if (out) return out;
    }
#else
    (void)gamma_is_inexact;
#endif

    return NULL; /* leave symbolic */
}

/* ------------------------------------------------------------------ */
/* Gamma[a, z0, z1] = Gamma[a, z0] - Gamma[a, z1]                      */
/* ------------------------------------------------------------------ */

static Expr* gamma_three_arg(Expr* a, Expr* z0, Expr* z1) {
    Expr* g0 = expr_new_function(expr_new_symbol("Gamma"),
                                 (Expr*[]){ expr_copy(a), expr_copy(z0) }, 2);
    Expr* g1 = expr_new_function(expr_new_symbol("Gamma"),
                                 (Expr*[]){ expr_copy(a), expr_copy(z1) }, 2);
    Expr* diff = expr_new_function(expr_new_symbol("Subtract"),
                                   (Expr*[]){ g0, g1 }, 2);
    return eval_and_free(diff);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

Expr* builtin_gamma(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return gamma_one_arg(args[0]);
    if (argc == 2) return gamma_two_arg(args[0], args[1]);
    if (argc == 3) return gamma_three_arg(args[0], args[1], args[2]);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void gamma_init(void) {
    symtab_add_builtin("Gamma", builtin_gamma);
    symtab_get_def("Gamma")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
