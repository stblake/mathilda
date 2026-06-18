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
#include "arithmetic.h"   /* is_rational */
#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>

/* Runaway guard on the superfactorial product. */
#define BARNESG_MAX_N 2000

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
    return NULL;
}

void barnesg_init(void) {
    symtab_add_builtin("BarnesG", builtin_barnesg);
    symtab_get_def("BarnesG")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
