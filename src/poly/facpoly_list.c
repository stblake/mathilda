/*
 * facpoly_list.c -- FactorList[poly] / FactorList[poly, opts].
 *
 * A thin wrapper over Factor: factor `poly` via `Factor[poly, opts...]` (all
 * options are forwarded verbatim -- GaussianIntegers, Extension, ...), then
 * split the resulting product into {factor, exponent} pairs.
 *
 *   FactorList[x^2 - 1]              -> {{1, 1}, {-1 + x, 1}, {1 + x, 1}}
 *   FactorList[2 x^3 + 2 x^2 - ...]  -> {{2, 1}, {-1 + x, 1}, {1 + x, 2}}
 *   FactorList[x^4 - 2, Extension -> Sqrt[2]]
 *                                    -> {{1, 1}, {Sqrt[2] + x^2, 1}, {-Sqrt[2] + x^2, 1}}
 *
 * The first element is always the overall numerical factor {c, 1} (which is
 * {1, 1} when there is no numerical factor).  Denominator factors of a
 * rational function carry negative exponents.
 *
 * Parsing rules on the Factor output R (a Times, a bare factor, or a number):
 *   - a number literal (Integer / Rational / Real / Complex / ...) multiplies
 *     into the overall numerical factor `c`;
 *   - Power[base, e] with an *integer* e is the pair {base, e} (this is a
 *     factor raised to a multiplicity, positive or negative);
 *   - anything else -- including Power[base, 1/2] = Sqrt[base], which is an
 *     irreducible factor in its own right -- is the pair {factor, 1}.
 */

#include "facpoly_list.h"

#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===================================================================== */
/* Small helpers                                                         */
/* ===================================================================== */

static Expr* mk_int(long n) { return expr_new_integer(n); }

static Expr* mk_list2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ a, b }, 2);
}

static bool is_rule_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol.name == SYM_Rule ||
            e->data.function.head->data.symbol.name == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

static bool head_name_is(const Expr* e, const char* interned) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == interned;
}

/* A pure numeric literal: it belongs to the overall numerical factor. */
static bool is_number_literal(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default:
            break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        (head_name_is(e, SYM_Rational) || head_name_is(e, SYM_Complex))) {
        return true;
    }
    return false;
}

static bool is_integer_literal(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT);
}

/* ===================================================================== */
/* Diagnostics                                                           */
/* ===================================================================== */

static Expr* factorlist_emit_argx(size_t argc) {
    fprintf(stderr,
            "FactorList::argx: FactorList called with %zu arguments; "
            "1 argument is expected.\n", argc);
    return NULL;
}

static Expr* factorlist_emit_nonopt(Expr* bad, size_t pos, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "FactorList::nonopt: Options expected (instead of %s) beyond "
            "position %zu in %s. An option must be a rule or a list of rules.\n",
            bad_str ? bad_str : "?", pos, call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* ===================================================================== */
/* Factor-output -> {factor, exponent} pairs                             */
/* ===================================================================== */

/* Append the pairs for one multiplicative factor `f` to (bases,exps); or, when
 * f is a number literal, fold it into *numeric (consumed & replaced). */
static void absorb_factor(const Expr* f, Expr** numeric,
                          Expr*** bases, Expr*** exps, size_t* n, size_t* cap) {
    if (is_number_literal(f)) {
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                                       (Expr*[]){ *numeric, expr_copy((Expr*)f) }, 2);
        *numeric = eval_and_free(prod);
        return;
    }

    Expr* base; Expr* exp;
    if (head_name_is(f, SYM_Power) && f->data.function.arg_count == 2 &&
        is_integer_literal(f->data.function.args[1])) {
        base = expr_copy(f->data.function.args[0]);
        exp  = expr_copy(f->data.function.args[1]);
    } else {
        base = expr_copy((Expr*)f);
        exp  = mk_int(1);
    }

    if (*n == *cap) {
        *cap *= 2;
        *bases = realloc(*bases, *cap * sizeof(Expr*));
        *exps  = realloc(*exps,  *cap * sizeof(Expr*));
    }
    (*bases)[*n] = base;
    (*exps)[*n]  = exp;
    (*n)++;
}

/* ===================================================================== */
/* Entry point                                                           */
/* ===================================================================== */

Expr* builtin_factorlist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return factorlist_emit_argx(0);

    /* Positions 2+ must be options (rules).  Report the last offender. */
    Expr* last_bad = NULL;
    for (size_t i = 1; i < argc; i++) {
        if (!is_rule_head(res->data.function.args[i]))
            last_bad = res->data.function.args[i];
    }
    if (last_bad) return factorlist_emit_nonopt(last_bad, 1, res);

    /* Factor[poly, opts...] -- forward every argument verbatim. */
    Expr** fargs = malloc(argc * sizeof(Expr*));
    for (size_t i = 0; i < argc; i++)
        fargs[i] = expr_copy(res->data.function.args[i]);
    Expr* fac = eval_and_free(
        expr_new_function(expr_new_symbol(SYM_Factor), fargs, argc));
    free(fargs);
    if (!fac) return NULL;

    /* Split the factored form. */
    Expr* numeric = mk_int(1);
    size_t cap = 8, n = 0;
    Expr** bases = malloc(cap * sizeof(Expr*));
    Expr** exps  = malloc(cap * sizeof(Expr*));

    if (head_name_is(fac, SYM_Times)) {
        for (size_t i = 0; i < fac->data.function.arg_count; i++)
            absorb_factor(fac->data.function.args[i], &numeric,
                          &bases, &exps, &n, &cap);
    } else {
        absorb_factor(fac, &numeric, &bases, &exps, &n, &cap);
    }
    expr_free(fac);

    /* Assemble {{numeric, 1}, {base_i, exp_i}, ...}. */
    size_t total = n + 1;
    Expr** items = malloc(total * sizeof(Expr*));
    items[0] = mk_list2(numeric, mk_int(1));
    for (size_t i = 0; i < n; i++)
        items[i + 1] = mk_list2(bases[i], exps[i]);
    free(bases);
    free(exps);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, total);
    free(items);
    return out;
}

/* ===================================================================== */
/* Init                                                                  */
/* ===================================================================== */

void factorlist_init(void) {
    symtab_add_builtin("FactorList", builtin_factorlist);
    symtab_get_def("FactorList")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    /* Docstring lives in info.c. */
}
