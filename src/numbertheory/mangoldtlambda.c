/* mangoldtlambda.c -- MangoldtLambda[n], the von Mangoldt function.
 *
 *   Lambda(n) = Log[p]   if n = p^k for a prime p and k >= 1,
 *             = 0        otherwise (n = 1, n composite with two distinct
 *                        primes, and every n <= 0).
 *
 * The answer is EXACT (Log[p], not its value), as in the reference:
 * MangoldtLambda[8] = Log[2]; N[MangoldtLambda[9]] = 1.0986...  The prime-
 * power test is the shared mth_prime_power_base_{u64,mpz} (src/facint.c), the
 * same one PrimePowerQ uses, so a large n never needs a full factorisation:
 * trial division, a perfect-power reduction, then one primality test.
 *
 * Argument classes:
 *   integer n (Integer / BigInt)  -> Log[p] or 0
 *   inexact number (2., 2.`30)    -> MangoldtLambda::exact, unevaluated
 *   anything else (x, 9/2, 2 I)   -> unevaluated, silently
 *
 * No packed / NDArray / Compile surface, deliberately: the result is an exact
 * symbolic Log[p] (an integer buffer maps to {Log[2], Log[3], Log[2], ...}),
 * which no machine dtype holds, and an inexact argument is an error, not a
 * value -- so there is no machine function to give a kernel.  Recorded with a
 * reason in the audit tools' EXEMPT lists.
 *
 * Attributes: Listable, Protected (the reference's set, less ReadProtected). */

#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "facint.h"      /* mth_prime_power_base_u64 / _mpz */
#include "message.h"     /* mth_message: Quiet/Check funnel */
#include "print.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"      /* builtin_arg_error */
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>

/* Log[p] for an int64 prime p. */
static Expr* mangoldt_log_of(Expr* p) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Log), (Expr*[]){ p }, 1));
}

/* True for an inexact number (machine or MPFR real, or a Complex with an
 * inexact part): the reference's MangoldtLambda::exact class. */
static bool mangoldt_is_inexact_number(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return mangoldt_is_inexact_number(re) || mangoldt_is_inexact_number(im);
    return false;
}

Expr* builtin_mangoldtlambda(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("MangoldtLambda", argc, 1, 1);
    Expr* n = res->data.function.args[0];

    if (n->type == EXPR_INTEGER) {
        if (n->data.integer < 2) return expr_new_integer(0);
        uint64_t p;
        if (!mth_prime_power_base_u64((uint64_t)n->data.integer, &p))
            return expr_new_integer(0);
        return mangoldt_log_of(expr_new_integer((int64_t)p));
    }
    if (n->type == EXPR_BIGINT) {
        if (mpz_sgn(n->data.bigint) <= 0) return expr_new_integer(0);
        mpz_t p;
        mpz_init(p);
        Expr* out;
        if (mth_prime_power_base_mpz(n->data.bigint, p))
            out = mangoldt_log_of(expr_new_bigint_from_mpz(p));
        else
            out = expr_new_integer(0);
        mpz_clear(p);
        return out;
    }
    /* Gated on the arithmetic-probe mute: NSum samples its summand at
     * non-integer points to learn whether it is continuous, and those trial
     * evaluations must stay silent (and invisible to Check). */
    if (mangoldt_is_inexact_number(n) && !g_arith_warnings_muted) {
        char* s = expr_to_string(n);
        char* c = expr_to_string(res);
        mth_message_gated(g_arith_warnings_muted, "MangoldtLambda", "exact",
                          "Argument %s in %s is not an exact number.", s ? s : "", c ? c : "");
        free(s); free(c);
    }
    return NULL;
}

void mangoldtlambda_init(void) {
    symtab_add_builtin("MangoldtLambda", builtin_mangoldtlambda);
    symtab_get_def("MangoldtLambda")->attributes |= (ATTR_PROTECTED | ATTR_LISTABLE);
    symtab_set_docstring("MangoldtLambda",
        "MangoldtLambda[n]\n"
        "\tThe von Mangoldt function: Log[p] when n = p^k is a power of a prime p\n"
        "\t(k >= 1), and 0 otherwise (including n = 1 and every n <= 0).\n"
        "\tThe result is exact (MangoldtLambda[8] = Log[2]); an inexact argument\n"
        "\temits MangoldtLambda::exact and stays unevaluated. Listable.");
}
