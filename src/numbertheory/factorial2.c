/* factorial2.c -- Factorial2[] (double factorial).
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
#include <gmp.h>

/* Factorial2[n] -- the double factorial.
 *
 * For non-negative integer n:
 *   n!! = n * (n-2) * (n-4) * ...    (terminating at 1 if n odd, 2 if n even).
 *   0!! = 1, (-1)!! = 1.
 * For negative odd integers n <= -3:
 *   n!! = ComplexInfinity (poles, like Factorial of negative integers).
 *
 * Symbolic / non-integer arguments return NULL (the call stays as-is)
 * so the surrounding evaluator can keep them in expression form.
 *
 * Implemented in two paths: int64 fast-path for small n (<= 30 -- fits in
 * int64 for all valid inputs since 31!! ~= 7.9e9 fits), and a GMP path
 * for larger n. */
/* True for an inexact numeric leaf: a machine/MPFR real, or a Complex[..] with
 * an inexact part.  Gates the analytic-continuation path so that under N a
 * non-integer argument gets a value, while exact rationals stay symbolic. */
static bool f2_is_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return f2_is_inexact(re) || f2_is_inexact(im);
    return false;
}

/* Analytic continuation of the double factorial, built as an expression and
 * evaluated (reusing Gamma / Power / Cos and their complex + MPFR kernels):
 *
 *     z!! = 2^(z/2) (Pi/2)^((Cos[z Pi] - 1)/4) Gamma[z/2 + 1].
 *
 * Reproduces the integer values (5!! = 15, 4!! = 8) and extends to real and
 * complex z. */
static Expr* f2_numeric(Expr* z) {
    Expr* half_z = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(1, 2), expr_copy(z) }, 2);          /* z/2 */

    Expr* term1 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_integer(2), expr_copy(half_z) }, 2);     /* 2^(z/2) */

    Expr* zpi = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_copy(z), expr_new_symbol(SYM_Pi) }, 2);      /* z Pi */
    Expr* cos_zpi = expr_new_function(expr_new_symbol(SYM_Cos), &zpi, 1);
    Expr* cos_minus1 = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ cos_zpi, expr_new_integer(-1) }, 2);
    Expr* expo2 = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ make_rational(1, 4), cos_minus1 }, 2);            /* (Cos[zPi]-1)/4 */

    Expr* pi_half = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_symbol(SYM_Pi), make_rational(1, 2) }, 2);
    Expr* term2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ pi_half, expo2 }, 2);                            /* (Pi/2)^... */

    Expr* garg = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ half_z, expr_new_integer(1) }, 2);               /* z/2 + 1 */
    Expr* term3 = expr_new_function(expr_new_symbol(SYM_Gamma), &garg, 1);

    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ term1, term2, term3 }, 3);
    return eval_and_free(prod);
}

Expr* builtin_factorial2(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Special-case the small literal -1 / 0 inputs. */
    if (arg->type == EXPR_INTEGER) {
        int64_t n = arg->data.integer;
        if (n == -1 || n == 0) return expr_new_integer(1);
        if (n < 0) {
            /* Even negative integers: ComplexInfinity (pole). Odd negative
             * integers below -1 also diverge. */
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        if (n <= 30) {
            int64_t f = 1;
            for (int64_t i = n; i >= 1; i -= 2) f *= i;
            return expr_new_integer(f);
        }
        /* GMP path. */
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (int64_t i = n; i >= 1; i -= 2) {
            mpz_set_si(k, (long)i);
            mpz_mul(r, r, k);
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    if (arg->type == EXPR_BIGINT) {
        /* For arguments large enough to be BigInt, use GMP. We only
         * handle non-negative inputs; negative BigInt is ComplexInfinity. */
        if (mpz_sgn(arg->data.bigint) < 0) {
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        if (mpz_sgn(arg->data.bigint) == 0) return expr_new_integer(1);
        /* Cap at a sane size: if mpz_fits_ulong_p fails the input is too
         * large for our explicit loop. Return NULL to leave it symbolic. */
        if (!mpz_fits_ulong_p(arg->data.bigint)) return NULL;
        unsigned long n = mpz_get_ui(arg->data.bigint);
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (unsigned long i = n; i >= 1; i -= 2) {
            mpz_set_ui(k, i);
            mpz_mul(r, r, k);
            if (i < 2) break; /* guard underflow when i==1 */
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    /* Non-integer numeric argument (under N): analytic continuation. */
    if (f2_is_inexact(arg)) {
        Expr* val = f2_numeric(arg);
        if (val && expr_is_numeric_like(val)) return val;
        expr_free(val);
    }

    return NULL; /* symbolic -- leave Factorial2[expr] alone */
}
