/* divisible.c -- Divisible[].
 * Split from numbertheory.c; see numbertheory.h and
 * numbertheory_internal.h for the subsystem layout. */

#include "infinity_rules.h"   /* infinity_direction */
#include "numbertheory.h"
#include "numbertheory_internal.h"
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "symtab.h"
#include "attr.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "ndarray.h"
#include "ndinteger.h" /* ndint_binary_predicate */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

/* Registration, attributes, and docstrings for the number-theory builtins.
 * Called from core_init(); replaces the inline block that previously lived
 * in core.c. */
/* Divisible::argmu / ::argm -- fewer than two arguments. */
static Expr* divisible_emit_argcount(size_t argc) {
    if (argc == 1)
        mth_message_gated(g_arith_warnings_muted, "Divisible", "argmu",
            "Divisible called with 1 argument; 2 or more arguments are expected.");
    else
        mth_message_gated(g_arith_warnings_muted, "Divisible", "argm",
            "Divisible called with %zu arguments; 2 or more arguments are expected.", argc);
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

typedef enum { DIV_FALSE, DIV_TRUE, DIV_UNKNOWN } DivAnswer;

/* Is n an integer (or Gaussian-integer) multiple of m? DIV_UNKNOWN for a
 * symbolic pair, which leaves the call unevaluated. m must be nonzero. */
static DivAnswer divisible_pair(Expr* n, Expr* m) {
    /* Fast path: both integer-like. mpz_divisible_p is exact at arbitrary
     * precision, so this also covers the large-integer cases. */
    if (n->type == EXPR_INTEGER && m->type == EXPR_INTEGER &&
        m->data.integer != -1 && m->data.integer != 0)
        return (n->data.integer % m->data.integer == 0) ? DIV_TRUE : DIV_FALSE;
    if (expr_is_integer_like(n) && expr_is_integer_like(m)) {
        mpz_t a, b;
        expr_to_mpz(n, a);
        expr_to_mpz(m, b);
        bool divisible = mpz_divisible_p(a, b) != 0;
        mpz_clears(a, b, NULL);
        return divisible ? DIV_TRUE : DIV_FALSE;
    }

    /* General path: build and evaluate the quotient n * m^-1. */
    Expr* inv_args[2] = { expr_copy(m), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), inv_args, 2);
    Expr* q_args[2] = { expr_copy(n), inv };
    Expr* q = expr_new_function(expr_new_symbol(SYM_Times), q_args, 2);
    Expr* quotient = eval_and_free(q);
    bool divisible = expr_is_integer_like(quotient) || is_gaussian_integer(quotient);
    expr_free(quotient);
    if (divisible) return DIV_TRUE;

    /* Not an integer. Concrete numeric quantities are manifestly not
     * divisible (False); symbolic ones stay unevaluated so rules apply. */
    if (divisible_is_numeric_quantity(n) && divisible_is_numeric_quantity(m))
        return DIV_FALSE;
    return DIV_UNKNOWN;
}

static bool divisible_is_zero(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer == 0) ||
           (e->type == EXPR_BIGINT && mpz_sgn(e->data.bigint) == 0);
}

/* Divisible[n, m] yields True when n is an integer multiple of m and False for
 * numeric quantities that manifestly are not. Over the integers this is
 * mpz_divisible_p; for Gaussian integers, rationals and exact numeric
 * quantities such as 2 Pi or Sqrt[6] it evaluates the quotient n/m and asks
 * whether that is an integer or Gaussian integer.
 *
 * Divisible[n, m1, m2, ...] is a CHAIN: n divisible by m1,
 * m1 by m2, and so on (Divisible[12, 2, 3] is False, Divisible[8, 2, 2] True).
 *
 * Messages, each leaving the call unevaluated: a zero
 * divisor draws ::divz, an inexact argument ::exact, fewer than two arguments
 * ::argmu / ::argm. Symbolic arguments stay unevaluated silently. Listable, so
 * the evaluator threads list arguments first. */
Expr* builtin_divisible(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return divisible_emit_argcount(argc);
    Expr** args = res->data.function.args;

    /* A packed / visible int64 buffer: one pass, answering with a packed bool
     * List; any shape the kernel declines (a zero divisor included, so ::divz
     * still fires per element) degrades element by element. */
    for (size_t i = 0; i < argc; i++) {
        if (is_ndarray(args[i])) {
            Expr* nd = ndint_binary_predicate(res, NDBP_DIVISIBLE);
            return nd ? nd : ndarray_delist_and_reeval(res);
        }
    }

    /* An infinite argument is no integer to divide by, nor to divide: the
     * quotient n/Infinity = 0 is an integer only by accident of the limit, so
     * the call stays unevaluated (audit Divisible#1). */
    for (size_t i = 0; i < argc; i++)
        if (infinity_direction(args[i], NULL) != INF_DIR_NONE) return NULL;

    Expr* inexact = nt_first_inexact_arg(res);
    if (inexact) return nt_emit_exact("Divisible", inexact, res);

    for (size_t i = 1; i < argc; i++) {
        if (divisible_is_zero(args[i])) {
            /* ::divz only when every argument is a number; a symbolic call
             * (Divisible[x, 0]) stays unevaluated silently. */
            for (size_t j = 0; j < argc; j++)
                if (!divisible_is_numeric_quantity(args[j])) return NULL;
            char* call = expr_to_string(res);
            mth_message_gated(g_arith_warnings_muted, "Divisible", "divz",
                "The argument 0 in %s should be nonzero.", call ? call : "?");
            free(call);
            return NULL;
        }
    }

    bool unknown = false;
    for (size_t i = 0; i + 1 < argc; i++) {
        DivAnswer a = divisible_pair(args[i], args[i + 1]);
        if (a == DIV_FALSE) return expr_new_symbol(SYM_False);
        if (a == DIV_UNKNOWN) unknown = true;
    }
    return unknown ? NULL : expr_new_symbol(SYM_True);
}
