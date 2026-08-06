/* Mathilda -- BarnesG[z], the Barnes G-function.
 *
 *   G(1) = G(2) = 1,   G(z+1) = Gamma[z] G(z),
 *   integer:  G(n+1) = prod_{k=1}^{n-1} k!   (the superfactorial),
 *             G(m) = 0 for non-positive integer m (double zeros).
 *
 * Exact for integer orders (GMP); non-integer orders are left unevaluated (the
 * LogGamma/zeta'(-1) asymptotic continuation is not implemented).  N at an
 * integer order routes through the exact value and numericalize.  Used by
 * Product to recognise prod_{k=1}^{n-1} Gamma[k] = BarnesG[n].
 *
 * Memory: honours the builtin ownership contract.
 */

#include "barnesg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational, is_complex, make_rational */
#include "sym_names.h"
#include "expr.h"
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Runaway guard on the superfactorial product. */
#define BARNESG_MAX_N 2000

/* ---- numeric (real/complex/arbitrary-precision) continuation ---------------
 *
 * BarnesG has no elementary closed form off the integers, so it is evaluated
 * from the Barnes asymptotic expansion combined with the recurrence
 * G(w+1) = Gamma(w) G(w).  The whole thing is assembled as a Mathilda
 * expression and handed to N[], so it inherits Gamma's / Log's / Exp's complex
 * and MPFR kernels and the Glaisher-Kinkelin constant at full precision:
 *
 *   G(w) = exp( A(Z) ) / prod_{j=0}^{m-1} Gamma(w + j),   Z = w + m - 1,
 *
 * with m chosen so Re(Z) is large enough for the (truncated) asymptotic
 *
 *   log G(Z+1) = (Z^2/2 - 1/12) log Z - 3 Z^2/4 + (Z/2) log(2 pi)
 *              + 1/12 - log A + sum_{k>=1} B_{2k+2} / ((2k)(2k+2)) Z^{-2k}.
 *
 * Twelve series terms are carried, so the shift target scales with the
 * requested precision; accurate to ~40 digits (12 fixed Bernoulli terms). */

/* c_k = B_{2k+2} / ((2k)(2k+2)), k = 1..12, reduced (all fit int64).  Verified
 * numerically against exact Log[BarnesG[n]]; note the denominator is
 * (2k)(2k+2), NOT (2k)(2k+1)(2k+2). */
static const int64_t BARNES_C_NUM[12] = {
    -1, 1, -1, 1, -691, 1, -3617, 43867, -174611, 854513, -236364091, 657931 };
static const int64_t BARNES_C_DEN[12] = {
    240, 1008, 1440, 1056, 327600, 144, 114240, 229824, 118800, 60720,
    1441440, 288 };

static bool barnes_leaf_double(const Expr* e, double* out) {
    int64_t n, d;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;   return true;
        case EXPR_REAL:    *out = e->data.real;              return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (is_rational((Expr*)e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

/* Numeric real part (for choosing the shift). */
static double barnes_real_part(const Expr* e) {
    double v; Expr *re, *im;
    if (barnes_leaf_double(e, &v)) return v;
    if (is_complex((Expr*)e, &re, &im) && barnes_leaf_double(re, &v)) return v;
    return 0.0;
}

/* Working precision in bits (53 machine; the max over an MPFR complex's parts). */
static long barnes_target_bits(const Expr* e) {
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        long a = barnes_target_bits(re), b = barnes_target_bits(im);
        return a > b ? a : b;
    }
    return 53;
}

/* True for a machine/MPFR real, or a Complex[..] with an inexact part. */
static bool barnes_is_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return barnes_is_inexact(re) || barnes_is_inexact(im);
    return false;
}

static Expr* barnesg_numeric(Expr* w) {
    long bits = barnes_target_bits(w);
    double D = (double)bits * 0.30103;                 /* target decimal digits */
    long R = (long)ceil(pow(10.0, (D + 4.0) / 24.0)) + 1;  /* shift so 12 terms reach D */
    if (R < 12) R = 12;
    long m = (long)ceil((double)R - barnes_real_part(w) + 1.0);
    if (m < 0) m = 0;

    Expr* Z = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(w), expr_new_integer(m - 1) }, 2);      /* Z = w + m - 1 */

    /* A(Z): the asymptotic log G(Z+1), as a Plus of 5 + 12 terms. */
    Expr* Zsq = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy(Z), expr_new_integer(2) }, 2);
    Expr* polyA = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ make_rational(1, 2), Zsq }, 2),
                   make_rational(-1, 12) }, 2);
    Expr* termA = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ polyA, expr_new_function(expr_new_symbol(SYM_Log),
                       (Expr*[]){ expr_copy(Z) }, 1) }, 2);           /* (Z^2/2-1/12)Log[Z] */

    Expr* termB = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(-3, 4),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(Z), expr_new_integer(2) }, 2) }, 2);  /* -3/4 Z^2 */

    Expr* log2pi = expr_new_function(expr_new_symbol(SYM_Log),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(2), expr_new_symbol(SYM_Pi) }, 2) }, 1);
    Expr* termC = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(1, 2), expr_copy(Z), log2pi }, 3);  /* (Z/2)Log[2Pi] */

    Expr* termD = make_rational(1, 12);
    Expr* termE = expr_new_function(expr_new_symbol(SYM_Times),      /* -Log[Glaisher] */
        (Expr*[]){ expr_new_integer(-1),
                   expr_new_function(expr_new_symbol(SYM_Log),
                       (Expr*[]){ expr_new_symbol(SYM_Glaisher) }, 1) }, 2);

    size_t nterms = 5 + 12;
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * nterms);
    terms[0] = termA; terms[1] = termB; terms[2] = termC; terms[3] = termD; terms[4] = termE;
    for (int k = 1; k <= 12; k++) {
        Expr* Zpow = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(Z), expr_new_integer(-2 * k) }, 2);
        terms[4 + k] = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ make_rational(BARNES_C_NUM[k - 1], BARNES_C_DEN[k - 1]), Zpow }, 2);
    }
    Expr* asym = expr_new_function(expr_new_symbol(SYM_Plus), terms, nterms);
    free(terms);
    expr_free(Z);

    Expr* gexp = expr_new_function(expr_new_symbol(SYM_Exp), (Expr*[]){ asym }, 1);  /* G(w+m) */

    Expr* result_expr;
    if (m == 0) {
        result_expr = gexp;
    } else {
        Expr** gam = (Expr**)malloc(sizeof(Expr*) * (size_t)m);
        for (long j = 0; j < m; j++) {
            Expr* wj = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(w), expr_new_integer(j) }, 2);
            gam[j] = expr_new_function(expr_new_symbol(SYM_Gamma), (Expr*[]){ wj }, 1);
        }
        Expr* prod = (m == 1) ? gam[0]
            : expr_new_function(expr_new_symbol(SYM_Times), gam, (size_t)m);
        free(gam);
        Expr* invprod = expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ prod, expr_new_integer(-1) }, 2);
        result_expr = expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ gexp, invprod }, 2);                          /* / prod Gamma */
    }

    /* Force numericalisation of the symbolic constants (Log[2Pi], Log[Glaisher],
     * Pi) at the target precision: machine -> N[expr], MPFR -> N[expr, digits].
     * The asymptotic A(Z) is a sum of terms as large as ~(Z^2/2) log Z, so its
     * exponentiation loses ~log10(Z^2) digits; carry that many guard digits. */
    Expr* ncall;
    if (bits <= 53) {
        ncall = expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ result_expr }, 1);
    } else {
        long digits = (long)D + (long)ceil((D + 4.0) / 12.0) + 5;   /* + cancellation guard */
        ncall = expr_new_function(expr_new_symbol(SYM_N),
              (Expr*[]){ result_expr, expr_new_integer(digits) }, 2);
    }
    return eval_and_free(ncall);
}

Expr* builtin_barnesg(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t z, d;
    if (is_rational(arg, &z, &d) && d == 1) {
        if (z <= 0) return expr_new_integer(0);          /* G(0), G(-1), ... = 0 */
        if (z <= 2) return expr_new_integer(1);          /* G(1) = G(2) = 1 */
        if (z - 1 > BARNESG_MAX_N) return NULL;
        /* G(z) = prod_{k=1}^{z-2} k!  (z >= 3). */
        mpz_t result, fact;
        mpz_init_set_ui(result, 1);
        mpz_init_set_ui(fact, 1);
        for (int64_t k = 1; k <= z - 2; k++) {
            mpz_mul_ui(fact, fact, (unsigned long)k);    /* fact = k! */
            mpz_mul(result, result, fact);
        }
        mpz_clear(fact);
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(result));
        mpz_clear(result);
        return r;
    }

    /* Non-integer numeric argument (under N): asymptotic continuation. */
    if (barnes_is_inexact(arg)) {
        Expr* val = barnesg_numeric(arg);
        if (val && expr_is_numeric_like(val)) return val;
        expr_free(val);
    }

    return NULL;
}

void barnesg_init(void) {
    symtab_add_builtin("BarnesG", builtin_barnesg);
    symtab_get_def("BarnesG")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
