/* Mathilda -- the inverse complementary error function.
 *
 *   InverseErfc[s]   inverse of erfc: the z solving s = erfc(z)
 *
 * Since erfc maps the real line onto (0, 2) (decreasing from 2 at -Infinity to
 * 0 at +Infinity), explicit numerical values are produced only for *real*
 * s in [0, 2]; out-of-domain (s < 0 or s > 2) and complex inputs stay symbolic.
 * Evaluation mirrors InverseErf -- no complex machinery:
 *
 *   exact special values   ->  InverseErfc[0]=Infinity, InverseErfc[1]=0,
 *                              InverseErfc[2]=-Infinity
 *   machine real (0<s<2)   ->  Winitzki seed + Newton polish on libm erfc
 *   arbitrary real (0<s<2) ->  Newton with precision doubling on mpfr_erfc
 *   everything else        ->  stays symbolic (return NULL)
 *
 * Mathematically InverseErfc[s] = InverseErf[1 - s], but we do NOT route
 * through InverseErf: for small s (large z) forming 1 - s and inverting erf
 * near 1 loses all significance to cancellation. Instead Newton iterates on
 * erfc directly -- both libm and MPFR ship cancellation-free erfc -- with
 *   f(z) = erfc(z) - s,  f'(z) = -(2/sqrt(pi)) e^{-z^2}, i.e.
 *   z <- z + (erfc(z) - s) (sqrt(pi)/2) e^{z^2}.
 * Newton converges quadratically; the MPFR path doubles the working precision
 * each step so the final full-precision erfc dominates the cost.
 *
 * The derivative D[InverseErfc[z], z] = -(sqrt(pi)/2) e^{InverseErfc[z]^2}
 * lives in src/calculus/deriv.c; Series follows from it via the generic
 * Taylor-via-D fallback. (erfc is not odd, so unlike InverseErf there is no
 * auto-applied symmetry rewrite.)
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "inverfc.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* sqrt(pi), the reciprocal-derivative scale of erfc. */
#define INVERFC_SQRT_PI 1.77245385090551602730

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool ifc_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

/* Build -Infinity, represented (as elsewhere in Mathilda) as Times[-1, Infinity]. */
static Expr* ifc_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* ------------------------------------------------------------------ */
/* Numeric kernel: machine precision                                  */
/* ------------------------------------------------------------------ */

/* erfcinv(x) for x strictly inside (0, 2) at full double accuracy. The seed
 * is Winitzki's closed-form erfinv on w = 1 - x (since erfcinv(x)=erfinv(1-x));
 * Newton then polishes directly on the libm erfc, converging quadratically. */
static double inverfc_double(double x) {
    if (x == 1.0) return 0.0;
    const double a = 0.147;                 /* Winitzki's tuned constant */
    double w  = 1.0 - x;                    /* erf(z) = 1 - erfc(z) = w   */
    double ln = log(1.0 - w * w);           /* < 0 for |w| < 1            */
    double t  = 2.0 / (M_PI * a) + ln / 2.0;
    double sgn = (w < 0.0) ? -1.0 : 1.0;
    double y  = sgn * sqrt(sqrt(t * t - ln / a) - t);

    /* Newton polish: y <- y + (erfc(y) - x) (sqrt(pi)/2) e^{y^2}. */
    for (int i = 0; i < 4; i++) {
        double f = erfc(y) - x;
        if (f == 0.0) break;
        y += f * (INVERFC_SQRT_PI / 2.0) * exp(y * y);
    }
    return y;
}

/* ------------------------------------------------------------------ */
/* Numeric kernel: arbitrary precision (MPFR)                         */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* One Newton step at working precision `p`, refining y in place:
 *   y <- y + (erfc(y) - s) (sqrt(pi)/2) e^{y^2}.
 * Intermediates are held at precision p so the cost of the expensive
 * erfc/exp tracks p; the result rounds into y at y's own precision. */
static void inverfc_newton_step(mpfr_t y, const mpfr_t s, mpfr_prec_t p) {
    mpfr_t f, ex, c;
    mpfr_init2(f, p);
    mpfr_init2(ex, p);
    mpfr_init2(c, p);

    mpfr_erfc(f, y, MPFR_RNDN);             /* erfc(y)             */
    mpfr_sub(f, f, s, MPFR_RNDN);           /* erfc(y) - s         */
    mpfr_sqr(ex, y, MPFR_RNDN);             /* y^2                 */
    mpfr_exp(ex, ex, MPFR_RNDN);            /* e^{y^2}             */
    mpfr_const_pi(c, MPFR_RNDN);
    mpfr_sqrt(c, c, MPFR_RNDN);
    mpfr_div_ui(c, c, 2, MPFR_RNDN);        /* sqrt(pi)/2          */
    mpfr_mul(f, f, ex, MPFR_RNDN);
    mpfr_mul(f, f, c, MPFR_RNDN);           /* full correction     */
    mpfr_add(y, y, f, MPFR_RNDN);           /* y <- y + correction */

    mpfr_clears(f, ex, c, (mpfr_ptr)0);
}

/* erfcinv(s) for 0 < s < 2 into `out` (already init2'd at the target precision).
 * Newton with precision doubling: seed from the double kernel, then refine,
 * roughly doubling the correct-bit count each step until the target is met. */
static void inverfc_mpfr(mpfr_t out, const mpfr_t s, mpfr_prec_t target) {
    mpfr_prec_t wp = target + 32;           /* working precision incl. guard */

    mpfr_t y;
    mpfr_init2(y, wp);
    mpfr_set_d(y, inverfc_double(mpfr_get_d(s, MPFR_RNDN)), MPFR_RNDN);

    /* The double seed is good to ~50 bits; double the precision each step. */
    mpfr_prec_t cur = 50;
    while (cur < wp) {
        cur *= 2;
        mpfr_prec_t p = (cur < wp ? cur : wp) + 16;
        if (p > wp) p = wp;
        inverfc_newton_step(y, s, p);
    }
    inverfc_newton_step(y, s, wp);          /* final clean-up at full precision */

    mpfr_set(out, y, MPFR_RNDN);
    mpfr_clear(y);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* InverseErfc[s]                                                     */
/* ------------------------------------------------------------------ */

static Expr* inverfc_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer == 0) return expr_new_symbol(SYM_Infinity);   /* +Inf */
        if (arg->data.integer == 1) return expr_new_integer(0);           /* 0    */
        if (arg->data.integer == 2) return ifc_neg_infinity();            /* -Inf */
        return NULL;                          /* out of domain [0, 2] */
    }
    if (ifc_is_symbol(arg, "Indeterminate")) return expr_new_symbol(SYM_Indeterminate);

    /* 2. Machine real (only inside the real domain [0, 2]). */
    if (arg->type == EXPR_REAL) {
        double x = arg->data.real;
        if (x < 0.0 || x > 2.0) return NULL;             /* out of domain */
        if (x == 0.0) return expr_new_symbol(SYM_Infinity);
        if (x == 2.0) return ifc_neg_infinity();
        return expr_new_real(inverfc_double(x));
    }

#ifdef USE_MPFR
    /* 3. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) {
        int lo = mpfr_cmp_ui(arg->data.mpfr, 0);         /* s vs 0 */
        int hi = mpfr_cmp_ui(arg->data.mpfr, 2);         /* s vs 2 */
        if (lo < 0 || hi > 0) return NULL;               /* out of domain */
        if (lo == 0) return expr_new_symbol(SYM_Infinity);
        if (hi == 0) return ifc_neg_infinity();
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        inverfc_mpfr(out->data.mpfr, arg->data.mpfr, prec);
        return out;
    }
#endif

    return NULL; /* leave symbolic (exact rationals, s<0 or s>2, complex, symbols) */
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argx diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* inverfc_emit_argx(size_t argc) {
    fprintf(stderr,
            "InverseErfc::argx: InverseErfc called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_inverfc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return inverfc_one_arg(args[0]);
    return inverfc_emit_argx(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void inverfc_init(void) {
    symtab_add_builtin("InverseErfc", builtin_inverfc);
    symtab_get_def("InverseErfc")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
