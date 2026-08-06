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
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"   /* is_rational, is_complex */
#include "sym_names.h"
#include "expr.h"
#include <gmp.h>
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

Expr* builtin_hyperfactorial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t n, d;
    if (is_rational(arg, &n, &d) && d == 1) {
        if (n < 0) return NULL;                         /* leave symbolic */
        if (n > HYPERFACTORIAL_MAX_N) return NULL;
        mpz_t result, term;
        mpz_init_set_ui(result, 1);
        mpz_init(term);
        for (int64_t k = 2; k <= n; k++) {
            mpz_ui_pow_ui(term, (unsigned long)k, (unsigned long)k);
            mpz_mul(result, result, term);
        }
        mpz_clear(term);
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(result));
        mpz_clear(result);
        return r;
    }

    /* Non-integer numeric argument (under N): Barnes-G continuation. */
    if (hf_is_inexact(arg)) {
        Expr* val = hf_numeric(arg);
        if (val && expr_is_numeric_like(val)) return val;
        expr_free(val);
    }

    return NULL;
}

void hyperfactorial_init(void) {
    symtab_add_builtin("Hyperfactorial", builtin_hyperfactorial);
    symtab_get_def("Hyperfactorial")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
}
