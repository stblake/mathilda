/* Mathilda -- the inverse error function.
 *
 *   InverseErf[s]       inverse of erf: the z solving s = erf(z)
 *   InverseErf[z0, s]   inverse of the generalized error function:
 *                       the z solving s = Erf[z0, z] = erf(z) - erf(z0)
 *
 * Per Mathematica, explicit numerical values are produced only for *real*
 * s in [-1, 1]; complex and out-of-domain (|s| > 1) inputs stay symbolic.
 * Evaluation is therefore much simpler than Erf -- no complex machinery:
 *
 *   exact special values   ->  InverseErf[0]=0, InverseErf[1]=Infinity,
 *                              InverseErf[-1]=-Infinity
 *   symbolic odd argument  ->  InverseErf[-x] = -InverseErf[x]
 *   machine real (|s|<1)   ->  Winitzki seed + Newton polish on libm erf
 *   arbitrary real (|s|<1) ->  Newton with precision doubling on mpfr_erf
 *   everything else        ->  stays symbolic (return NULL)
 *
 * Neither C99 libm nor MPFR ships an inverse-erf, so the kernel is Newton's
 * method on f(z) = erf(z) - s, with f'(z) = (2/sqrt(pi)) e^{-z^2}, i.e.
 *   z <- z - (erf(z) - s) (sqrt(pi)/2) e^{z^2}.
 * Newton converges quadratically; the MPFR path doubles the working
 * precision each step so the final full-precision erf dominates the cost.
 *
 * The derivative D[InverseErf[z], z] = (sqrt(pi)/2) e^{InverseErf[z]^2}
 * lives in src/calculus/deriv.c; Series follows from it via the generic
 * Taylor-via-D fallback.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "inverf.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* is_rational */
#include "attr.h"
#include "eval.h"         /* eval_and_free */
#include "expr.h"
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* M_PI is POSIX, not C99 -- provide a fallback (see CLAUDE.md). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* sqrt(pi), the reciprocal-derivative scale of erf. */
#define INVERF_SQRT_PI 1.77245385090551602730

/* ------------------------------------------------------------------ */
/* Small leaf helpers                                                 */
/* ------------------------------------------------------------------ */

/* True if `e` is exactly the symbol `name`. */
static bool iee_is_symbol(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* True if `e` is a negative real number literal (Integer/Real/Rational/MPFR). */
static bool iee_is_neg_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL)    return e->data.real < 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_sgn(e->data.mpfr) < 0;
#endif
    int64_t n, d;
    if (is_rational(e, &n, &d)) return ((n < 0) ^ (d < 0)) != 0;
    return false;
}

/* True if `e` is a Times[c, ...] whose leading factor `c` is a negative
 * real number literal -- the trigger for the odd-symmetry rewrite. */
static bool iee_is_neg_leading_times(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count < 2) return false;
    if (!iee_is_symbol(e->data.function.head, "Times")) return false;
    return iee_is_neg_number(e->data.function.args[0]);
}

/* Recursively test whether `e` contains a subexpression headed by `head`. */
static bool iee_contains_head(const Expr* e, const char* head) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (iee_is_symbol(e->data.function.head, head)) return true;
    if (iee_contains_head(e->data.function.head, head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (iee_contains_head(e->data.function.args[i], head)) return true;
    return false;
}

/* Build -Infinity, represented (as elsewhere in Mathilda) as Times[-1, Infinity]. */
static Expr* iee_neg_infinity(void) {
    return expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
}

/* ------------------------------------------------------------------ */
/* Numeric kernel: machine precision                                  */
/* ------------------------------------------------------------------ */

/* erfinv(x) for x strictly inside (-1, 1) at full double accuracy.
 * Winitzki's closed-form approximation seeds Newton's iteration on the
 * libm erf, which converges quadratically to a fixed point. */
static double inverf_double(double x) {
    if (x == 0.0) return 0.0;
    const double a = 0.147;                 /* Winitzki's tuned constant */
    double ln = log(1.0 - x * x);           /* < 0 for |x| < 1 */
    double t  = 2.0 / (M_PI * a) + ln / 2.0;
    double sgn = (x < 0.0) ? -1.0 : 1.0;
    double y  = sgn * sqrt(sqrt(t * t - ln / a) - t);

    /* Newton polish: y <- y - (erf(y) - x) (sqrt(pi)/2) e^{y^2}. */
    for (int i = 0; i < 4; i++) {
        double f = erf(y) - x;
        if (f == 0.0) break;
        y -= f * (INVERF_SQRT_PI / 2.0) * exp(y * y);
    }
    return y;
}

/* ------------------------------------------------------------------ */
/* Numeric kernel: arbitrary precision (MPFR)                         */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* One Newton step at working precision `p`, refining y in place:
 *   y <- y - (erf(y) - s) (sqrt(pi)/2) e^{y^2}.
 * Intermediates are held at precision p so the cost of the expensive
 * erf/exp tracks p; the result rounds into y at y's own precision. */
static void inverf_newton_step(mpfr_t y, const mpfr_t s, mpfr_prec_t p) {
    mpfr_t f, ex, c;
    mpfr_init2(f, p);
    mpfr_init2(ex, p);
    mpfr_init2(c, p);

    mpfr_erf(f, y, MPFR_RNDN);              /* erf(y)              */
    mpfr_sub(f, f, s, MPFR_RNDN);           /* erf(y) - s          */
    mpfr_sqr(ex, y, MPFR_RNDN);             /* y^2                 */
    mpfr_exp(ex, ex, MPFR_RNDN);            /* e^{y^2}             */
    mpfr_const_pi(c, MPFR_RNDN);
    mpfr_sqrt(c, c, MPFR_RNDN);
    mpfr_div_ui(c, c, 2, MPFR_RNDN);        /* sqrt(pi)/2          */
    mpfr_mul(f, f, ex, MPFR_RNDN);
    mpfr_mul(f, f, c, MPFR_RNDN);           /* full correction     */
    mpfr_sub(y, y, f, MPFR_RNDN);           /* y <- y - correction */

    mpfr_clears(f, ex, c, (mpfr_ptr)0);
}

/* erfinv(s) for |s| < 1 into `out` (already init2'd at the target precision).
 * Newton with precision doubling: seed from the double kernel, then refine,
 * roughly doubling the correct-bit count each step until the target is met. */
static void inverf_mpfr(mpfr_t out, const mpfr_t s, mpfr_prec_t target) {
    mpfr_prec_t wp = target + 32;           /* working precision incl. guard */

    mpfr_t y;
    mpfr_init2(y, wp);
    mpfr_set_d(y, inverf_double(mpfr_get_d(s, MPFR_RNDN)), MPFR_RNDN);

    /* The double seed is good to ~50 bits; double the precision each step. */
    mpfr_prec_t cur = 50;
    while (cur < wp) {
        cur *= 2;
        mpfr_prec_t p = (cur < wp ? cur : wp) + 16;
        if (p > wp) p = wp;
        inverf_newton_step(y, s, p);
    }
    inverf_newton_step(y, s, wp);           /* final clean-up at full precision */

    mpfr_set(out, y, MPFR_RNDN);
    mpfr_clear(y);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------ */
/* InverseErf[s]                                                      */
/* ------------------------------------------------------------------ */

static Expr* inverf_one_arg(Expr* arg) {
    /* 1. Exact special values. */
    if (arg->type == EXPR_INTEGER) {
        if (arg->data.integer == 0)  return expr_new_integer(0);          /* 0 */
        if (arg->data.integer == 1)  return expr_new_symbol(SYM_Infinity);  /* +Inf */
        if (arg->data.integer == -1) return iee_neg_infinity();           /* -Inf */
        return NULL;                          /* |n| >= 2: out of domain */
    }
    if (iee_is_symbol(arg, "Indeterminate")) return expr_new_symbol(SYM_Indeterminate);

    /* 2. Machine real (only inside the real domain [-1, 1]). */
    if (arg->type == EXPR_REAL) {
        double x = arg->data.real;
        if (x > 1.0 || x < -1.0) return NULL;            /* out of domain */
        if (x == 1.0)  return expr_new_symbol(SYM_Infinity);
        if (x == -1.0) return iee_neg_infinity();
        return expr_new_real(inverf_double(x));
    }

#ifdef USE_MPFR
    /* 3. Arbitrary-precision real. */
    if (arg->type == EXPR_MPFR) {
        mpfr_t one;
        mpfr_init2(one, 2);
        mpfr_set_ui(one, 1, MPFR_RNDN);
        int c = mpfr_cmpabs(arg->data.mpfr, one);        /* |s| vs 1 */
        mpfr_clear(one);
        if (c > 0) return NULL;                          /* out of domain */
        if (c == 0)
            return mpfr_sgn(arg->data.mpfr) > 0 ? expr_new_symbol(SYM_Infinity)
                                                : iee_neg_infinity();
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* out = expr_new_mpfr_bits(prec);
        inverf_mpfr(out->data.mpfr, arg->data.mpfr, prec);
        return out;
    }
#endif

    /* 4. Odd symmetry: InverseErf[-x] = -InverseErf[x] for a symbolic
     *    negative-leading argument (numeric negatives handled above). */
    if (iee_is_neg_leading_times(arg)) {
        Expr* pos = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(arg) }, 2));
        Expr* inner = expr_new_function(expr_new_symbol(SYM_InverseErf), &pos, 1);
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), inner }, 2));
    }

    return NULL; /* leave symbolic (exact rationals, |s|>1, complex, symbols) */
}

/* ------------------------------------------------------------------ */
/* InverseErf[z0, s] = InverseErf[s + Erf[z0]]                        */
/* ------------------------------------------------------------------ */

static Expr* inverf_two_arg(Expr* z0, Expr* s) {
    /* s = erf(z) - erf(z0)  =>  z = InverseErf[s + Erf[z0]]. */
    Expr* erf_z0 = expr_new_function(expr_new_symbol(SYM_Erf),
                                     (Expr*[]){ expr_copy(z0) }, 1);
    Expr* sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                     (Expr*[]){ expr_copy(s), erf_z0 }, 2));
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol(SYM_InverseErf),
                                     &sum, 1));
    /* If Erf[z0] did not reduce (still an Erf head survives), keep the
     * two-argument form symbolic rather than committing to the rewrite. */
    if (iee_contains_head(res, "Erf")) {
        expr_free(res);
        return NULL;
    }
    return res;
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                */
/* ------------------------------------------------------------------ */

/* Mathematica-compatible argt diagnostic; returns NULL so the evaluator
 * leaves the call unevaluated. */
static Expr* inverf_emit_argt(size_t argc) {
    fprintf(stderr,
            "InverseErf::argt: InverseErf called with %zu argument%s; "
            "1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_inverf(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc == 1) return inverf_one_arg(args[0]);
    if (argc == 2) return inverf_two_arg(args[0], args[1]);
    return inverf_emit_argt(argc);
}

/* ------------------------------------------------------------------ */
/* Registration                                                       */
/* ------------------------------------------------------------------ */

void inverf_init(void) {
    symtab_add_builtin("InverseErf", builtin_inverf);
    symtab_get_def("InverseErf")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
