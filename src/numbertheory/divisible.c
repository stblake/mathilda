/* divisible.c -- Divisible[].
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

/* Registration, attributes, and docstrings for the number-theory builtins.
 * Called from core_init(); replaces the inline block that previously lived
 * in core.c. */
/* Emit the wrong-argument-count diagnostic for Divisible and return NULL so
 * the call is left unevaluated.  Too few arguments use Mathematica's `argm`
 * tag; too many use `argt`. */
static Expr* divisible_emit_argcount(size_t argc) {
    if (argc < 2) {
        fprintf(stderr,
                "Divisible::argm: Divisible called with %zu argument%s; "
                "2 or more arguments are expected.\n",
                argc, argc == 1 ? "" : "s");
    } else {
        fprintf(stderr,
                "Divisible::argt: Divisible called with %zu arguments; "
                "2 arguments are expected.\n",
                argc);
    }
    return NULL;
}


/* True when e is a concrete numeric quantity (exact number, recognised
 * constant such as Pi or E, or a NumericFunction applied to numeric
 * quantities).  Mirrors the file-static helper of the same name behind
 * NumericQ in core.c — inlined here, as in rationalize.c, because that helper
 * is not exported.  Used by Divisible to distinguish a manifestly
 * non-divisible numeric pair (returns False) from symbolic arguments (left
 * unevaluated). */
static bool divisible_is_numeric_quantity(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol.name;
        return name == SYM_Pi || name == SYM_E || name == SYM_I ||
               name == SYM_Infinity || name == SYM_ComplexInfinity ||
               name == SYM_EulerGamma || name == SYM_GoldenRatio ||
               name == SYM_Catalan || name == SYM_Degree ||
               name == SYM_GoldenAngle || name == SYM_Glaisher ||
               name == SYM_Khinchin;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head_name = e->data.function.head->data.symbol.name;
        if (head_name == SYM_Complex || head_name == SYM_Rational) return true;
        SymbolDef* def = symtab_get_def(head_name);
        if (def && (def->attributes & ATTR_NUMERICFUNCTION)) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!divisible_is_numeric_quantity(e->data.function.args[i])) return false;
            }
            return true;
        }
    }
    return false;
}

/* Divisible[n, m] yields True when n is an integer multiple of m and False
 * for numeric quantities that manifestly are not.  Over the integers this is
 * Mod[n, m] == 0, computed directly with GMP's mpz_divisible_p (handling
 * BigInts and the m == 0 edge: divisible by 0 iff n == 0).  For everything
 * else — Gaussian integers, rationals, and exact numeric quantities such as
 * 2 Pi or Sqrt[6] — it forms the quotient n/m, evaluates it exactly, and
 * reports True iff that quotient collapses to an integer or Gaussian integer.
 * Symbolic, non-numeric arguments leave the call unevaluated.  Divisible is
 * Listable, so list arguments are threaded by the evaluator before reaching
 * this builtin. */
Expr* builtin_divisible(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return divisible_emit_argcount(argc);

    Expr* n = res->data.function.args[0];
    Expr* m = res->data.function.args[1];

    /* Fast path: both integer-like.  mpz_divisible_p is exact at arbitrary
     * precision, so this also covers the large-integer cases. */
    if (expr_is_integer_like(n) && expr_is_integer_like(m)) {
        mpz_t a, b;
        mpz_init(a);
        mpz_init(b);
        expr_to_mpz(n, a);
        expr_to_mpz(m, b);
        bool divisible = mpz_divisible_p(a, b) != 0;
        mpz_clears(a, b, NULL);
        return expr_new_symbol(divisible ? SYM_True : SYM_False);
    }

    /* General path: build and evaluate the quotient n * m^-1. */
    Expr* inv_args[2] = { expr_copy(m), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
    Expr* q_args[2] = { expr_copy(n), inv };
    Expr* q = expr_new_function(expr_new_symbol(SYM_Times), q_args, 2);
    Expr* quotient = eval_and_free(q);

    bool divisible = expr_is_integer_like(quotient) || is_gaussian_integer(quotient);
    expr_free(quotient);
    if (divisible) return expr_new_symbol(SYM_True);

    /* Not an integer.  When both arguments are concrete numeric quantities
     * they are manifestly not divisible (False); otherwise (symbolic args)
     * leave Divisible unevaluated so user rules / pattern matching apply. */
    if (divisible_is_numeric_quantity(n) && divisible_is_numeric_quantity(m)) {
        return expr_new_symbol(SYM_False);
    }
    return NULL;
}
