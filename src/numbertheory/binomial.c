/* binomial.c -- Binomial[] and binomial_polynomial helper.
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Binomial[n, k] expanded as a falling-factorial polynomial for a small
 * non-negative integer k and arbitrary n (symbolic, rational, complex,
 * real, integer, ...).  Yields
 *
 *     Binomial[n, k] = n (n - 1) (n - 2) ... (n - k + 1) / k!
 *
 * Returns a fully-evaluated expression.  Used for two cases:
 *   - concrete non-negative integer k <= 32, n anything non-integer:
 *     Binomial[n, 4]  -> n (n-1) (n-2) (n-3) / 24
 *     Binomial[1+I,5] -> -1/12 - I/12 (Times/Plus fold complex factors)
 *   - the symmetric-reduction case where Subtract[n, m] collapses to a
 *     small non-negative integer k: Binomial[9/2, 7/2] -> Binomial[9/2,1]
 *     -> 9/2; Binomial[n, n-1] -> Binomial[n, 1] -> n.
 *
 * Cap of 32 mirrors FactorialPower so a careless Binomial[n, 1000] does
 * not blow up the tree size.  Returns NULL for k outside [0, 32]. */
static Expr* binomial_polynomial(Expr* n, int64_t k) {
    if (k < 0 || k > 32) return NULL;
    if (k == 0) return expr_new_integer(1);
    if (k == 1) return expr_copy(n);

    /* numerator factors: n, n-1, ..., n-k+1, plus a trailing 1/k! coefficient. */
    Expr** factors = malloc(sizeof(Expr*) * (size_t)(k + 1));
    for (int64_t i = 0; i < k; i++) {
        if (i == 0) {
            factors[i] = expr_copy(n);
        } else {
            Expr* shift = expr_new_integer(-i);
            factors[i] = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(n), shift }, 2);
        }
    }

    /* 1 / k! coefficient.  k! is computed in GMP so values for k > 20
     * (where k! overflows int64) still round-trip exactly through a
     * Rational[1, BigInt] which the evaluator's is_rational_like folders
     * recognise. */
    mpz_t fact_z;
    mpz_init(fact_z);
    mpz_fac_ui(fact_z, (unsigned long)k);
    Expr* fact_expr = expr_bigint_normalize(expr_new_bigint_from_mpz(fact_z));
    mpz_clear(fact_z);

    Expr* inv;
    if (fact_expr->type == EXPR_INTEGER) {
        inv = make_rational(1, fact_expr->data.integer);
        expr_free(fact_expr);
    } else {
        Expr* one = expr_new_integer(1);
        inv = expr_new_function(expr_new_symbol(SYM_Rational),
            (Expr*[]){ one, fact_expr }, 2);
    }
    factors[k] = inv;

    Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
        factors, (size_t)(k + 1));
    free(factors);
    return eval_and_free(product);
}

/* Convert a real-numeric leaf (integer, bigint, rational, machine real, or
 * MPFR real) to a double.  Returns false for symbolic or complex operands, so
 * the machine Binomial path can decline rather than silently reading a
 * non-convertible operand (e.g. an exact Rational) as 0. */
static bool binomial_to_double(const Expr* e, double* out) {
    int64_t n, d;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;    return true;
        case EXPR_REAL:    *out = e->data.real;               return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);  return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (is_rational(e, &n, &d)) { *out = (double)n / (double)d; return true; }
    return false;
}

#ifdef USE_MPFR
/* Set an already-init2'd mpfr from a real-numeric leaf (integer, bigint,
 * rational, machine real, MPFR).  Returns false for symbolic / complex. */
static bool binomial_set_mpfr(mpfr_t out, const Expr* e) {
    int64_t n, d;
    switch (e->type) {
        case EXPR_INTEGER: mpfr_set_si(out, (long)e->data.integer, MPFR_RNDN); return true;
        case EXPR_BIGINT:  mpfr_set_z (out, e->data.bigint,        MPFR_RNDN); return true;
        case EXPR_REAL:    mpfr_set_d (out, e->data.real,          MPFR_RNDN); return true;
        case EXPR_MPFR:    mpfr_set   (out, e->data.mpfr,          MPFR_RNDN); return true;
        default: break;
    }
    if (is_rational(e, &n, &d)) {
        mpfr_set_si(out, (long)n, MPFR_RNDN);
        mpfr_div_si(out, out, (long)d, MPFR_RNDN);
        return true;
    }
    return false;
}

/* Arbitrary-precision generalised binomial via the Gamma quotient
 *
 *     Binomial[n, m] = Gamma[n+1] / (Gamma[m+1] Gamma[n-m+1]),
 *
 * evaluated in MPFR when at least one operand is an MPFR real (i.e. under
 * N[Binomial[..], prec]).  The working precision follows Mathematica's
 * contagion rule -- the minimum precision among the inexact operands, floored
 * at machine (53).  Returns NULL when an operand is not real-numeric, or when
 * the quotient overflows / is indeterminate, leaving the expression symbolic. */
static Expr* binomial_mpfr(const Expr* arg_n, const Expr* arg_m) {
    mpfr_prec_t out_prec = 0;
    if (arg_n->type == EXPR_MPFR) out_prec = mpfr_get_prec(arg_n->data.mpfr);
    if (arg_m->type == EXPR_MPFR) {
        mpfr_prec_t p = mpfr_get_prec(arg_m->data.mpfr);
        if (out_prec == 0 || p < out_prec) out_prec = p;
    }
    if (out_prec < 53) out_prec = 53;
    mpfr_prec_t wp = out_prec + 64;                 /* guard bits */

    mpfr_t nv, mv, g1, g2, g3, den;
    mpfr_inits2(wp, nv, mv, g1, g2, g3, den, (mpfr_ptr)0);
    if (!binomial_set_mpfr(nv, arg_n) || !binomial_set_mpfr(mv, arg_m)) {
        mpfr_clears(nv, mv, g1, g2, g3, den, (mpfr_ptr)0);
        return NULL;
    }
    mpfr_add_ui(g1, nv, 1, MPFR_RNDN);  mpfr_gamma(g1, g1, MPFR_RNDN);   /* Gamma[n+1]   */
    mpfr_add_ui(g2, mv, 1, MPFR_RNDN);  mpfr_gamma(g2, g2, MPFR_RNDN);   /* Gamma[m+1]   */
    mpfr_sub(g3, nv, mv, MPFR_RNDN);
    mpfr_add_ui(g3, g3, 1, MPFR_RNDN);  mpfr_gamma(g3, g3, MPFR_RNDN);   /* Gamma[n-m+1] */
    mpfr_mul(den, g2, g3, MPFR_RNDN);

    Expr* out = expr_new_mpfr_bits(out_prec);
    mpfr_div(out->data.mpfr, g1, den, MPFR_RNDN);
    mpfr_clears(nv, mv, g1, g2, g3, den, (mpfr_ptr)0);

    if (mpfr_nan_p(out->data.mpfr) || mpfr_inf_p(out->data.mpfr)) {
        expr_free(out);
        return NULL;
    }
    return out;
}
#endif /* USE_MPFR */

/* True for an inexact real leaf (machine or MPFR). */
static bool binomial_is_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True when e is a Complex[..] value carrying at least one inexact part --
 * i.e. a numeric complex produced by N[..], as opposed to an exact Gaussian
 * like 1 + I. */
static bool binomial_inexact_complex(Expr* e) {
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        return binomial_is_inexact(re) || binomial_is_inexact(im);
    }
    return false;
}

/* True for any Complex[..] operand, exact (1 + I) or inexact (2. + I). */
static bool binomial_is_complex(Expr* e) {
    Expr *re, *im;
    return is_complex(e, &re, &im);
}

/* True for an operand that makes the whole Binomial a numeric computation: an
 * inexact real (machine or MPFR) or a Complex[..] carrying an inexact part.
 * When either operand forces numeric AND a complex operand is present, the
 * generalised binomial is evaluated through the Gamma quotient (path 6) -- an
 * exact Gaussian sibling (1 + I) is then carried along by the Times/Plus
 * numeric contagion, so Binomial[1 + I, 5.] evaluates rather than falling
 * through to symbolic. */
static bool binomial_forces_numeric(Expr* e) {
    return binomial_is_inexact(e) || binomial_inexact_complex(e);
}

/* Build and evaluate  Gamma[n+1] / (Gamma[m+1] Gamma[n-m+1])  as an
 * expression.  This reuses Gamma's complex kernels (machine Lanczos and, under
 * USE_MPFR, the arbitrary-precision Spouge path) and the complex-arithmetic
 * folders rather than reimplementing complex Gamma here.  Returns the
 * evaluated expression; the caller checks it is actually numeric before
 * accepting it (a pole or a still-symbolic operand leaves an unevaluated
 * Gamma, which is rejected). */
static Expr* binomial_gamma_quotient(Expr* n, Expr* m) {
    Expr* np1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(n), expr_new_integer(1) }, 2);
    Expr* mp1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(m), expr_new_integer(1) }, 2);
    Expr* negm = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy(m) }, 2);
    Expr* nmm1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(n), negm, expr_new_integer(1) }, 3);

    Expr* gn  = expr_new_function(expr_new_symbol(SYM_Gamma), &np1,  1);
    Expr* gm  = expr_new_function(expr_new_symbol(SYM_Gamma), &mp1,  1);
    Expr* gnm = expr_new_function(expr_new_symbol(SYM_Gamma), &nmm1, 1);

    Expr* inv_gm  = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ gm,  expr_new_integer(-1) }, 2);
    Expr* inv_gnm = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ gnm, expr_new_integer(-1) }, 2);

    Expr* quotient = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ gn, inv_gm, inv_gnm }, 3);
    return eval_and_free(quotient);
}

/* Binomial[n, m] -- generalised binomial coefficient
 *
 *     Binomial[n, m] = Gamma[n+1] / (Gamma[m+1] * Gamma[n - m + 1])
 *
 * Evaluation strategy (first matching rule wins):
 *   1. Both args concrete integer-like: exact GMP via mpz_bin_ui, with
 *      Mathematica's Pascal extension Binomial[n, m] = (-1)^m Binomial[
 *      m - n - 1, m] for n < 0, m >= 0.
 *   2. Either arg machine-precision real: tgamma path -> EXPR_REAL.
 *   3. Symmetric identity: if Subtract[n, m] simplifies to a non-negative
 *      integer k <= 32, expand Binomial[n, k] as a falling-factorial
 *      polynomial.  This catches Binomial[n, n-1] -> n,
 *      Binomial[9/2, 7/2] -> 9/2, Binomial[n+1, n-1] -> n (n+1) / 2, ...
 *   4. Concrete non-negative integer m <= 32 with anything-n: expand
 *      Binomial[n, m] as a falling-factorial polynomial.  Handles
 *      Binomial[n, 4] -> n(n-1)(n-2)(n-3)/24 and lets the Times/Plus
 *      folders simplify complex / rational n (Binomial[1+I, 5]
 *      -> -1/12 - I/12).
 *   5. Arbitrary-precision real (MPFR) operand: the Gamma quotient in MPFR,
 *      so N[Binomial[7/3, 1/5], p] evaluates while the exact form stays
 *      symbolic.
 *   6. Complex numeric operand (a Complex[..] with an inexact part): the
 *      Gamma quotient evaluated as an expression, reusing Gamma's complex
 *      kernels -- N[Binomial[1/2 + I/3, 1/4]]. */
Expr* builtin_binomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* arg_n = res->data.function.args[0];
    Expr* arg_m = res->data.function.args[1];

    /* --- (1) exact integer / integer path. --- */
    if (expr_is_integer_like(arg_n) && expr_is_integer_like(arg_m)) {
        mpz_t n, m;
        expr_to_mpz(arg_n, n);
        expr_to_mpz(arg_m, m);

        if (mpz_sgn(m) < 0) {
            mpz_clears(n, m, NULL);
            return expr_new_integer(0);
        }
        if (mpz_sgn(n) >= 0 && mpz_cmp(m, n) > 0) {
            mpz_clears(n, m, NULL);
            return expr_new_integer(0);
        }
        /* mpz_bin_ui needs m to fit in unsigned long; pathologically
         * large m has no finite expansion we'd want anyway. */
        if (!mpz_fits_ulong_p(m)) {
            mpz_clears(n, m, NULL);
            return NULL;
        }
        unsigned long mk = mpz_get_ui(m);
        mpz_t r;
        mpz_init(r);
        if (mpz_sgn(n) < 0) {
            /* Pascal extension for negative n. */
            mpz_t shifted;
            mpz_init(shifted);
            mpz_neg(shifted, n);
            mpz_add(shifted, shifted, m);
            mpz_sub_ui(shifted, shifted, 1);
            mpz_bin_ui(r, shifted, mk);
            if (mk & 1UL) mpz_neg(r, r);
            mpz_clear(shifted);
        } else {
            mpz_bin_ui(r, n, mk);
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clears(n, m, r, NULL);
        return out;
    }

    /* --- (2) machine real path via tgamma.  Fires when either operand is a
     * machine real; the other may be any real-numeric leaf (integer, bigint,
     * rational, MPFR).  If an operand is symbolic/complex the path declines --
     * so Binomial[2.5, 1/5] evaluates numerically but Binomial[2.5, x] is left
     * symbolic, rather than both reading the sibling as 0. --- */
    if (arg_n->type == EXPR_REAL || arg_m->type == EXPR_REAL) {
        double nv, mv;
        if (binomial_to_double(arg_n, &nv) && binomial_to_double(arg_m, &mv)) {
            return expr_new_real(tgamma(nv + 1.0) / (tgamma(mv + 1.0) * tgamma(nv - mv + 1.0)));
        }
    }

    /* --- (3) symmetric identity: reduce when n - m is a small non-neg int.
     *
     * Warnings are muted around the Subtract evaluation: when n / m happen
     * to involve a 1/0 partial that lives only inside this exploratory
     * difference, the user is not subtracting -- we are -- so the
     * Power::infy diagnostic would be spurious. --- */
    {
        arith_warnings_mute_push();
        Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract),
            (Expr*[]){ expr_copy(arg_n), expr_copy(arg_m) }, 2);
        Expr* diff_eval = eval_and_free(diff);
        arith_warnings_mute_pop();

        int64_t k = -1;
        if (diff_eval && diff_eval->type == EXPR_INTEGER &&
            diff_eval->data.integer >= 0 && diff_eval->data.integer <= 32) {
            k = diff_eval->data.integer;
        }
        expr_free(diff_eval);

        if (k >= 0) {
            Expr* poly = binomial_polynomial(arg_n, k);
            if (poly) return poly;
        }
    }

    /* --- (4) concrete non-negative integer m: falling-factorial expand. --- */
    if (arg_m->type == EXPR_INTEGER && arg_m->data.integer >= 0 &&
        arg_m->data.integer <= 32) {
        return binomial_polynomial(arg_n, arg_m->data.integer);
    }

#ifdef USE_MPFR
    /* --- (5) arbitrary-precision real path via the Gamma quotient.  Placed
     * after (3)/(4) so an integer m keeps the exact falling-factorial form;
     * fires for N[Binomial[7/3, 1/5], prec] where both operands are MPFR. --- */
    if (arg_n->type == EXPR_MPFR || arg_m->type == EXPR_MPFR) {
        Expr* out = binomial_mpfr(arg_n, arg_m);
        if (out) return out;
    }
#endif

    /* --- (6) complex numeric path.  Fires when the computation is numeric --
     * an inexact operand is present (machine real 5., MPFR, or a Complex with
     * an inexact part) -- AND a complex operand is present, so the real-only
     * machine path (2) has already declined.  Evaluates the Gamma quotient,
     * reusing Gamma's machine-Lanczos / MPFR-Spouge complex kernels; any exact
     * complex operand (1 + I, 7 - 3 I) is numericalised by the Times/Plus
     * numeric contagion as the quotient folds.  So Binomial[1 + I, 5.] and
     * Binomial[2. + I, 7 - 3 I] both evaluate, while an exact Gaussian pair
     * (Binomial[1 + I, 2 + I]) has no inexact operand and stays symbolic.
     * Accepts the result only if it is numeric, so a pole or a symbolic
     * operand leaves Binomial symbolic. --- */
    if ((binomial_forces_numeric(arg_n) || binomial_forces_numeric(arg_m)) &&
        (binomial_is_complex(arg_n) || binomial_is_complex(arg_m))) {
        Expr* val = binomial_gamma_quotient(arg_n, arg_m);
        if (val && expr_is_numeric_like(val)) return val;
        expr_free(val);
    }

    return NULL;
}
