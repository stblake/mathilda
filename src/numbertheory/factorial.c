/* factorial.c -- Factorial[].
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
#include <gmp.h>

Expr* builtin_factorial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Machine Real: Factorial[x] = Gamma[x + 1] via libm tgamma. */
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        return expr_new_real(tgamma(v + 1.0));
    }
#ifdef USE_MPFR
    /* MPFR Real: same identity at full input precision. */
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        mpfr_t shifted;
        mpfr_init2(shifted, prec);
        mpfr_add_ui(shifted, arg->data.mpfr, 1, MPFR_RNDN);
        Expr* result = expr_new_mpfr_bits(prec);
        mpfr_gamma(result->data.mpfr, shifted, MPFR_RNDN);
        mpfr_clear(shifted);
        return result;
    }
#endif
    /* BigInt: factorial of a value that exceeds int64 is astronomical
     * (1e20! has ~10^21 digits) and would exhaust memory. Mathematica
     * leaves it symbolic for the same reason. */

    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        if (d == 1) {
            if (n < 0) return expr_new_symbol(SYM_ComplexInfinity);
            if (n <= 20) {
                int64_t f = 1;
                for (int64_t i = 2; i <= n; i++) f *= i;
                return expr_new_integer(f);
            } else {
                mpz_t result;
                mpz_init(result);
                mpz_fac_ui(result, (unsigned long)n);
                Expr* r = expr_new_bigint_from_mpz(result);
                mpz_clear(result);
                return r;
            }
        } else if (d == 2 || d == -2) {
            if (d == -2) { n = -n; }
            /* Half-integer Factorial: (n/2)! = coeff * Sqrt[Pi], where the
             * coefficient is built in GMP -- the odd double factorial and the
             * 2^k denominator both blow past int64 well before the argument is
             * large (e.g. Gamma[201/2] needs 199!! / 2^100 ~ 10^157). */
            mpz_t num, den;
            mpz_init_set_ui(num, 1);
            mpz_init_set_ui(den, 1);

            if (n > 0) {
                for (int64_t i = n; i >= 1; i -= 2) mpz_mul_si(num, num, (long)i);
                mpz_mul_2exp(den, den, (unsigned long)((n + 1) / 2)); /* 2^((n+1)/2) */
            } else {
                for (int64_t i = n + 2; i <= -1; i += 2) {
                    mpz_mul_ui(num, num, 2);
                    mpz_mul_si(den, den, (long)i);
                }
            }
            /* Keep the denominator positive (negative half-integers can flip it). */
            if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }

            Expr* coeff = mpz_pair_to_rational_expr(num, den);
            mpz_clears(num, den, NULL);
            if (!coeff) coeff = expr_new_integer(0);

            Expr* pi_sym = expr_new_symbol(SYM_Pi);
            Expr* half = make_rational(1, 2);
            Expr* sqrt_pi = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){pi_sym, half}, 2));
            
            if (coeff->type == EXPR_INTEGER && coeff->data.integer == 1) {
                expr_free(coeff);
                return sqrt_pi;
            } else {
                return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){coeff, sqrt_pi}, 2));
            }
        }
    }

    return NULL;
}
