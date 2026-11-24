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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "sf_scalar.h"          /* sf_machine_scalar_args: the machine-scalar gate */
#include "sf_machine_error.h"   /* the shared raw-double kernels (ND / Compile too) */

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
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
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
/* Exposed (see sf_machine.h's rationale): this is already a pure
 * double algorithm, so the NDArray fast path and the Compile[] VM can use
 * it directly instead of sitting on a degrade sentinel. */
double inverfc_double(double x) {
    if (x == 1.0) return 0.0;
    {   /* the piecewise-polynomial kernel; the Winitzki + Newton fallback only
         * for what it declines */
        double v;
        if (sf_inverse_erfc_r(x, &v)) return v;
    }
    if (x > 0.0 && x < 1e-300) {
        /* Subnormal / tiny s: Winitzki's 1 - s rounds to 1 and the seed is
         * log(0). The leading asymptotics of erfc, e^{-z^2}/(z sqrt(pi)) = s,
         * give a seed the MPFR Newton refines (it used to start from NaN and
         * answer Indeterminate). */
        double L = -log(x);
        double z = sqrt(L);
        z = sqrt(L - log(INVERFC_SQRT_PI * z));
        return sqrt(L - log(INVERFC_SQRT_PI * z));
    }
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

/* erfcinv(t) for 0 < t <= 1 into `out` (already init2'd at the target
 * precision). Newton with precision doubling, from a seed that is good to ~50
 * bits whatever t is:
 *   - t in the normal double range: the double kernel on t's nearest double.
 *     Rounding t perturbs it by 2^-53 relative, which moves the root by
 *     ~2^-53 / (2 y) -- below the seed's own error.
 *   - t below it (an MPFR 10^-400 has no double): the leading asymptotics
 *     e^{-y^2} / (y sqrt(pi)) = t from log t (computed in MPFR, so it never
 *     underflows), then Newton at 64 bits until it settles. */
static void inverfc_mpfr_unit(mpfr_t out, const mpfr_t t, mpfr_prec_t target) {
    mpfr_prec_t wp = target + 32;           /* working precision incl. guard */

    mpfr_t y;
    mpfr_init2(y, wp);
    if (mpfr_get_exp(t) > -1000) {
        mpfr_set_d(y, inverfc_double(mpfr_get_d(t, MPFR_RNDN)), MPFR_RNDN);
    } else {
        mpfr_t lt;
        mpfr_init2(lt, 64);
        mpfr_log(lt, t, MPFR_RNDN);
        double L = -mpfr_get_d(lt, MPFR_RNDN);  /* > 693 here */
        mpfr_clear(lt);
        double z = sqrt(L);
        for (int i = 0; i < 4; i++) z = sqrt(L - log(INVERFC_SQRT_PI * z));
        mpfr_set_d(y, z, MPFR_RNDN);
        /* Relative error of the asymptotic seed is O(1/L^2): converge it to
         * double accuracy at 64 bits before the precision-doubling stage. */
        mpfr_t prev;
        mpfr_init2(prev, 64);
        for (int i = 0; i < 40; i++) {
            mpfr_set(prev, y, MPFR_RNDN);
            inverfc_newton_step(y, t, 64);
            mpfr_sub(prev, prev, y, MPFR_RNDN);
            if (mpfr_zero_p(prev) ||
                mpfr_get_exp(prev) < mpfr_get_exp(y) - 55) break;
        }
        mpfr_clear(prev);
    }

    /* The seed is good to ~50 bits; double the precision each step. */
    mpfr_prec_t cur = 50;
    while (cur < wp) {
        cur *= 2;
        mpfr_prec_t p = (cur < wp ? cur : wp) + 16;
        if (p > wp) p = wp;
        inverfc_newton_step(y, t, p);
    }
    inverfc_newton_step(y, t, wp);          /* final clean-up at full precision */

    mpfr_set(out, y, MPFR_RNDN);
    mpfr_clear(y);
}

/* erfcinv(s) for 0 < s < 2 into `out` (init2'd at the target precision). Above
 * 1 by the reflection erfcinv(s) = -erfcinv(2 - s): 2 - s is exact at s's own
 * precision (Sterbenz, 1 <= s < 2), so s = 2 - 10^-30 is solved on 10^-30
 * instead of on a seed from the double 2.0 (which was -Infinity, then NaN).
 * Shared with InverseErf near +-1 (erfinv(s) = erfcinv(1 - s)). */
void inverfc_mpfr(mpfr_t out, const mpfr_t s, mpfr_prec_t target) {
    if (mpfr_cmp_ui(s, 1) <= 0) { inverfc_mpfr_unit(out, s, target); return; }
    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(s));
    mpfr_ui_sub(t, 2, s, MPFR_RNDN);        /* exact */
    inverfc_mpfr_unit(out, t, target);
    mpfr_neg(out, out, MPFR_RNDN);
    mpfr_clear(t);
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
        {
            double v;
            if (sf_inverse_erfc_r(x, &v)) return expr_new_real(v);
        }
#ifdef USE_MPFR
        /* The kernel declines a subnormal s (its Halley residual cannot
         * resolve s's few bits; the Winitzki seed below breaks down there
         * too): solve at 53 bits in MPFR instead. */
        {
            mpfr_t s53, r53;
            mpfr_init2(s53, 53);
            mpfr_init2(r53, 53);
            mpfr_set_d(s53, x, MPFR_RNDN);
            inverfc_mpfr(r53, s53, 53);
            double v = mpfr_get_d(r53, MPFR_RNDN);
            mpfr_clears(s53, r53, (mpfr_ptr)0);
            if (isfinite(v)) return expr_new_real(v);
        }
#endif
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
    mth_message("InverseErfc", "argx", "InverseErfc called with %zu argument%s; 1 argument is expected.", argc, argc == 1 ? "" : "s");
    return NULL;
}

Expr* builtin_inverfc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Machine-scalar fast path (see builtin_inverf). */
    if (argc == 1) {
        SFMachineArgs ma;
        double v;
        if (sf_machine_scalar_args(res, &ma) && !ma.any_complex &&
            sf_inverse_erfc_r(ma.re[0], &v))
            return expr_new_real(v);
    }
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
