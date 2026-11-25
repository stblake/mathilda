/* Mathilda -- Hyperfactorial.
 *
 *   Hyperfactorial[n] = prod_{k=1}^{n} k^k   (H(0) = H(1) = 1).
 *
 * Exact for a non-negative integer order (GMP); non-positive-integer,
 * non-integer, or symbolic orders are left unevaluated (the analytic
 * K-function continuation is not implemented).  N at an integer order routes
 * through the exact value and numericalize, so machine and MPFR precision come
 * for free.  Used by Product to recognise prod k^k = Hyperfactorial[n].
 *
 * Memory: honours the builtin ownership contract (never frees res; returns a
 * fresh Expr* or NULL; clears every GMP temporary).
 */

#include "hyperfactorial.h"
#include "numeric.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational, is_complex */
#include "sym_names.h"
#include "expr.h"
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>

/* Runaway guard: k^k products past this are astronomical (Mathematica still
 * computes them, but we cap to protect memory, mirroring Factorial's stance). */
#define HYPERFACTORIAL_MAX_N 20000

/* True for a machine/MPFR real, or a Complex[..] with an inexact part. */
static bool hf_is_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return hf_is_inexact(re) || hf_is_inexact(im);
    return false;
}

/* Analytic continuation via the Barnes G-function, built as an expression and
 * evaluated (so it inherits BarnesG's / Gamma's numeric kernels):
 *
 *     Hyperfactorial[z] = Gamma[z+1]^z / BarnesG[z+1].
 *
 * Reproduces the integer values (H(4) = 27648) and extends to real/complex z. */
static Expr* hf_numeric(Expr* z) {
    Expr* gam = expr_new_function(expr_new_symbol(SYM_Gamma),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(z), expr_new_integer(1) }, 2) }, 1);   /* Gamma[z+1] */
    Expr* gpow = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ gam, expr_copy(z) }, 2);                             /* Gamma[z+1]^z */

    Expr* bg = expr_new_function(expr_new_symbol("BarnesG"),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_copy(z), expr_new_integer(1) }, 2) }, 1);   /* BarnesG[z+1] */
    Expr* inv_bg = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ bg, expr_new_integer(-1) }, 2);

    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ gpow, inv_bg }, 2);
    return eval_and_free(prod);
}

/* H(n) for an integer n, into `out`: the product prod_{k<=n} k^k for n >= 0;
 * for n = -m < 0 the recurrence H(z) = z^z H(z - 1) run downwards from
 * H(0) = 1 (0^0 = 1) gives H(-m) = (-1)^Floor[m/2] H(m - 1):
 * H(-1) = 1, H(-2) = -1, H(-3) = -4 (Hyperfactorial#1 of the 2026-11 audit).
 * False when |n| is past the size cap. */
static bool hf_integer(int64_t n, mpz_t out) {
    bool neg = n < 0;
    int64_t m = neg ? -(n + 1) : n;               /* H(-m') uses H(m' - 1) */
    if (m > HYPERFACTORIAL_MAX_N) return false;
    mpz_t term;
    mpz_set_ui(out, 1);
    mpz_init(term);
    for (int64_t k = 2; k <= m; k++) {
        mpz_ui_pow_ui(term, (unsigned long)k, (unsigned long)k);
        mpz_mul(out, out, term);
    }
    mpz_clear(term);
    if (neg && (((-n) / 2) & 1)) mpz_neg(out, out);
    return true;
}

Expr* builtin_hyperfactorial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t n, d;
    if (is_rational(arg, &n, &d) && d == 1) {
        mpz_t result;
        mpz_init(result);
        if (!hf_integer(n, result)) { mpz_clear(result); return NULL; }
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(result));
        mpz_clear(result);
        return r;
    }

    /* An inexact non-positive integer: the Barnes-G continuation is 0/0
     * there (Gamma[z+1] and BarnesG[z+1] both have poles), while H is
     * finite -- take the exact value (Hyperfactorial#2). */
    {
        double v = 0.0;
        bool real_int = false;
        if (arg->type == EXPR_REAL) {
            v = arg->data.real;
            real_int = (v <= 0.0 && v == floor(v) && v >= -1.0e6);
        }
#ifdef USE_MPFR
        else if (arg->type == EXPR_MPFR && mpfr_integer_p(arg->data.mpfr)
                 && mpfr_sgn(arg->data.mpfr) <= 0
                 && mpfr_cmp_si(arg->data.mpfr, -1000000) >= 0) {
            v = mpfr_get_d(arg->data.mpfr, MPFR_RNDN);
            real_int = true;
        }
#endif
        if (real_int) {
            mpz_t result;
            mpz_init(result);
            Expr* r = NULL;
            if (hf_integer((int64_t)v, result)) {
#ifdef USE_MPFR
                if (arg->type == EXPR_MPFR)
                    r = expr_new_mpfr_from_mpz(result, mpfr_get_prec(arg->data.mpfr));
                else
#endif
                if (mpz_sizeinbase(result, 2) < 1000)
                    r = expr_new_real(mpz_get_d(result));
#ifdef USE_MPFR
                else   /* past the double range: the extended-exponent real */
                    r = expr_new_mpfr_from_mpz(result, 53);
#endif
            }
            mpz_clear(result);
            if (r) return r;
        }
    }

    /* Non-integer numeric argument (under N): Barnes-G continuation. */
    if (hf_is_inexact(arg)) {
        Expr* val = hf_numeric(arg);
        if (val && expr_is_numeric_like(val)) return val;
        /* Hyperfactorial[50.] is ~10^1895: the machine continuation overflows a
         * double (and its intermediate can go NaN). Promote to an extended-
         * exponent MPFR real, or let the net make a past-range result Overflow[]. */
        val = numeric_promote_result_if_degenerate(val, res);
        if (val) return val;
    }

    return NULL;
}

void hyperfactorial_init(void) {
    symtab_add_builtin("Hyperfactorial", builtin_hyperfactorial);
    symtab_get_def("Hyperfactorial")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
