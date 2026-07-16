/* divisors.c -- Divisors[].
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

/* Emit the wrong-argument-count diagnostic for Divisors and return NULL so the
 * call is left unevaluated. */
static Expr* divisors_emit_argx(size_t npos) {
    if (!arith_warnings_muted()) {
        fprintf(stderr,
                "Divisors::argx: Divisors called with %zu argument%s; "
                "1 argument is expected.\n",
                npos, npos == 1 ? "" : "s");
    }
    return NULL;
}

/* Divisors[n] gives the ascending list of positive integers dividing n (sign
 * ignored).  Divisors[n, GaussianIntegers -> True] (or a non-real Gaussian
 * input) returns the divisors in Z[i], one first-quadrant representative per
 * associate class, sorted by (Re, Im).  Machine ints and BigInts are handled
 * uniformly through GMP.  Divisors is Listable, so list arguments are threaded
 * by the evaluator before this builtin runs.  Non-integer, zero, or 0-argument
 * calls are left unevaluated (the last with a Divisors::argx message). */
Expr* builtin_divisors(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    /* Separate the lone positional argument from a GaussianIntegers option. */
    bool gaussian = false, gaussian_set = false;
    Expr* n = NULL;
    size_t npos = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_FUNCTION &&
            a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_Rule &&
            a->data.function.arg_count == 2) {
            Expr* key = a->data.function.args[0];
            Expr* val = a->data.function.args[1];
            if (key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_GaussianIntegers &&
                val->type == EXPR_SYMBOL &&
                (val->data.symbol.name == SYM_True || val->data.symbol.name == SYM_False)) {
                gaussian = (val->data.symbol.name == SYM_True);
                gaussian_set = true;
                continue;
            }
            return NULL;  /* malformed / unknown option: leave unevaluated */
        }
        npos++;
        n = a;
    }

    if (npos != 1) return divisors_emit_argx(npos);

    /* Explicit option wins; otherwise auto-enable Z[i] for non-real input. */
    bool use_gaussian = gaussian_set ? gaussian : is_gaussian_integer(n);

    if (use_gaussian) {
        mpz_t a, b;
        if (!df_to_gaussian(n, a, b)) return NULL;
        Expr* out = divisors_gaussian(a, b);
        mpz_clears(a, b, NULL);
        return out;
    }

    /* Ordinary integer path. */
    if (!expr_is_integer_like(n)) return NULL;
    mpz_t v;
    expr_to_mpz(n, v);                 /* inits v */
    if (mpz_sgn(v) == 0) { mpz_clear(v); return NULL; }   /* Divisors[0] */
    mpz_abs(v, v);

    Expr* out;
    if (mpz_cmp_ui(v, 1) == 0) {
        Expr* items[1] = { expr_new_integer(1) };
        out = expr_new_function(expr_new_symbol(SYM_List), items, 1);
    } else {
        out = divisors_ordinary(v);
    }
    mpz_clear(v);
    return out;
}
