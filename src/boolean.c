#include "boolean.h"
#include "symtab.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "ndarray.h"
#include "pack.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

Expr* builtin_not(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_SYMBOL) {
        if (arg->data.symbol.name == SYM_True) return expr_new_symbol(SYM_False);
        if (arg->data.symbol.name == SYM_False) return expr_new_symbol(SYM_True);
    }
    return NULL;
}

Expr* builtin_and(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 0) return expr_new_symbol(SYM_True);
    if (res->data.function.arg_count == 1) {
        /* And has HoldAll, so the single arg has not been evaluated yet. */
        return evaluate(res->data.function.args[0]);
    }

    bool changed = false;
    Expr** new_args = malloc(sizeof(Expr*) * res->data.function.arg_count);
    size_t new_count = 0;
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* evaluated_arg = evaluate(res->data.function.args[i]);
        if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol.name == SYM_False) {
            expr_free(evaluated_arg);
            for (size_t j = 0; j < new_count; j++) expr_free(new_args[j]);
            free(new_args);
            return expr_new_symbol(SYM_False);
        } else if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol.name == SYM_True) {
            expr_free(evaluated_arg);
            changed = true;
        } else {
            if (!expr_eq(evaluated_arg, res->data.function.args[i])) changed = true;
            new_args[new_count++] = evaluated_arg;
        }
    }
    if (new_count == 0) {
        free(new_args);
        return expr_new_symbol(SYM_True);
    }
    if (!changed) {
        for (size_t i = 0; i < new_count; i++) expr_free(new_args[i]);
        free(new_args);
        return NULL;
    }
    if (new_count == 1) {
        Expr* ret = new_args[0];
        free(new_args);
        return ret;
    }
    Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
    free(new_args);
    return ret;
}

Expr* builtin_or(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 0) return expr_new_symbol(SYM_False);
    if (res->data.function.arg_count == 1) {
        /* Or has HoldAll, so the single arg has not been evaluated yet. */
        return evaluate(res->data.function.args[0]);
    }

    bool changed = false;
    Expr** new_args = malloc(sizeof(Expr*) * res->data.function.arg_count);
    size_t new_count = 0;
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* evaluated_arg = evaluate(res->data.function.args[i]);
        if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol.name == SYM_True) {
            expr_free(evaluated_arg);
            for (size_t j = 0; j < new_count; j++) expr_free(new_args[j]);
            free(new_args);
            return expr_new_symbol(SYM_True);
        } else if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol.name == SYM_False) {
            expr_free(evaluated_arg);
            changed = true;
        } else {
            if (!expr_eq(evaluated_arg, res->data.function.args[i])) changed = true;
            new_args[new_count++] = evaluated_arg;
        }
    }
    if (new_count == 0) {
        free(new_args);
        return expr_new_symbol(SYM_False);
    }
    if (!changed) {
        for (size_t i = 0; i < new_count; i++) expr_free(new_args[i]);
        free(new_args);
        return NULL;
    }
    if (new_count == 1) {
        Expr* ret = new_args[0];
        free(new_args);
        return ret;
    }
    Expr* ret = expr_new_function(expr_copy(res->data.function.head), new_args, new_count);
    free(new_args);
    return ret;
}

/* Boole[expr]
 *   True  -> 1
 *   False -> 0
 *   anything else -> stays unevaluated
 *
 * With ATTR_LISTABLE assigned in attr.c, the evaluator threads Boole
 * over lists automatically, so this implementation only handles the
 * scalar 1-argument case. */
Expr* builtin_boole(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_SYMBOL) {
        if (arg->data.symbol.name == SYM_True)  return expr_new_integer(1);
        if (arg->data.symbol.name == SYM_False) return expr_new_integer(0);
    }
    /* A bool array is the natural input to Boole. A VISIBLE bool NDArray is an
     * atom, so the Listable threading above never fires for it and it would
     * otherwise stay unevaluated; map it to a packed int64 array (True->1,
     * False->0) directly off the buffer, inheriting the presentation. (A packed
     * bool List is threaded elementwise by the Listable step before we get here,
     * so this reads the visible surface; handling either is harmless.) */
    if (is_ndarray(arg) && arg->data.ndarray.dtype == NDT_BOOL) {
        size_t sz = ndarray_size(arg);
        if (sz == 0) return NULL;
        const uint8_t* in = (const uint8_t*)arg->data.ndarray.data;
        void* buf = NULL;
        Expr* out = ndbuild_open_like(arg, arg->data.ndarray.rank,
                                      arg->data.ndarray.dims, NDT_INT64, &buf);
        if (!out) return NULL;
        /* Present as a List (a packed int64 one), matching Boole of the plain
         * List of True/False -- same reasoning as the sign predicates. */
        out->data.ndarray.present_as = NDA_HEAD_LIST;
        pack_g_any_created = true;
        int64_t* ob = (int64_t*)buf;
        for (size_t i = 0; i < sz; i++) ob[i] = in[i] ? 1 : 0;
        return out;
    }
    return NULL;
}

/* ConditionalExpression[expr, cond]
 *   cond == True   -> expr
 *   cond == False  -> Undefined
 *   nested: ConditionalExpression[ConditionalExpression[e, c1], c2]
 *           -> ConditionalExpression[e, c1 && c2]
 *   otherwise stays unevaluated. */
Expr* builtin_conditional_expression(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* expr = res->data.function.args[0];
    Expr* cond = res->data.function.args[1];

    if (cond->type == EXPR_SYMBOL) {
        if (cond->data.symbol.name == SYM_True) {
            /* Steal the expr slot so the evaluator's free of res does not
             * double-free what we just returned. */
            res->data.function.args[0] = NULL;
            return expr;
        }
        if (cond->data.symbol.name == SYM_False) {
            return expr_new_symbol(SYM_Undefined);
        }
    }

    /* Nested flattening: ConditionalExpression[ConditionalExpression[e, c1], c2]
     * collapses to ConditionalExpression[e, And[c1, c2]]. */
    if (expr->type == EXPR_FUNCTION &&
        expr->data.function.head->type == EXPR_SYMBOL &&
        expr->data.function.head->data.symbol.name == SYM_ConditionalExpression &&
        expr->data.function.arg_count == 2) {
        Expr* inner_expr = expr_copy(expr->data.function.args[0]);
        Expr** and_args = malloc(sizeof(Expr*) * 2);
        and_args[0] = expr_copy(expr->data.function.args[1]);
        and_args[1] = expr_copy(cond);
        Expr* combined = expr_new_function(expr_new_symbol(SYM_And), and_args, 2);
        free(and_args);
        Expr* combined_eval = evaluate(combined);
        expr_free(combined);
        Expr** new_args = malloc(sizeof(Expr*) * 2);
        new_args[0] = inner_expr;
        new_args[1] = combined_eval;
        Expr* out = expr_new_function(expr_new_symbol(SYM_ConditionalExpression),
                                      new_args, 2);
        free(new_args);
        return out;
    }

    return NULL;
}

/* Xor[e1, e2, ...] -- exclusive or.  Flat/Orderless/OneIdentity, so the
 * evaluator has already flattened nested Xor and sorted the arguments; this
 * builtin folds the literal Booleans and cancels duplicate pairs (a Xor a =
 * False).  Xor[] -> False, Xor[a] -> a.  An odd number of consumed True
 * arguments negates the surviving core (Xor[True, a] -> Not[a]).  Returns NULL
 * (stays symbolic) when nothing simplifies. */
Expr* builtin_xor(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return expr_new_symbol(SYM_False);
    Expr** args = res->data.function.args;
    if (argc == 1) return expr_copy(args[0]);      /* OneIdentity collapse */

    int  parity  = 0;                              /* count of True, mod 2 */
    bool changed = false;
    Expr** keep = malloc(sizeof(Expr*) * argc);    /* borrowed pointers */
    int nk = 0;
    for (size_t i = 0; i < argc; i++) {
        Expr* a = args[i];
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_True)  { parity ^= 1; changed = true; continue; }
        if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_False) { changed = true; continue; }
        bool cancelled = false;                    /* a Xor a = False */
        for (int j = 0; j < nk; j++)
            if (keep[j] && expr_eq(keep[j], a)) { keep[j] = NULL; cancelled = true; changed = true; break; }
        if (!cancelled) keep[nk++] = a;
    }
    int m = 0;                                     /* compact out cancelled */
    for (int j = 0; j < nk; j++) if (keep[j]) keep[m++] = keep[j];

    if (m == 0) { free(keep); return expr_new_symbol(parity ? SYM_True : SYM_False); }
    if (!changed) { free(keep); return NULL; }     /* nothing simplified */

    Expr* core;
    if (m == 1) {
        core = expr_copy(keep[0]);
    } else {
        Expr** ca = malloc(sizeof(Expr*) * m);
        for (int j = 0; j < m; j++) ca[j] = expr_copy(keep[j]);
        core = expr_new_function(expr_new_symbol(SYM_Xor), ca, (size_t)m);
        free(ca);
    }
    free(keep);
    if (parity) {                                  /* odd Trues -> negate core */
        Expr* na[1] = { core };
        return expr_new_function(expr_new_symbol(SYM_Not), na, 1);
    }
    return core;
}

/* Implies[p, q] -- material implication (!p || q). */
Expr* builtin_implies(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* p = res->data.function.args[0];
    Expr* q = res->data.function.args[1];
    if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_False) return expr_new_symbol(SYM_True);  /* False => _ */
    if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_True)  return expr_copy(q);               /* True  => q */
    if (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True)  return expr_new_symbol(SYM_True);  /* _ => True  */
    if (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_False) {                                  /* p => False = !p */
        Expr* na[1] = { expr_copy(p) };
        return expr_new_function(expr_new_symbol(SYM_Not), na, 1);
    }
    if (expr_eq(p, q)) return expr_new_symbol(SYM_True);                                                /* p => p */
    return NULL;
}

void boolean_init(void) {
    symtab_add_builtin("And", builtin_and);
    symtab_add_builtin("Or", builtin_or);
    symtab_add_builtin("Not", builtin_not);
    symtab_add_builtin("Xor", builtin_xor);
    symtab_add_builtin("Implies", builtin_implies);
    symtab_add_builtin("Boole", builtin_boole);
    symtab_add_builtin("ConditionalExpression", builtin_conditional_expression);

    symtab_set_docstring("Xor",
        "Xor[e1, e2, ...]\n"
        "\tThe logical exclusive OR of the ei: True when an odd number of the\n"
        "\targuments are True.  Folds literal Booleans and cancels duplicate\n"
        "\targuments (a Xor a is False); Xor[] is False and Xor[e] is e.");
    symtab_set_docstring("Implies",
        "Implies[p, q]\n"
        "\tThe material implication p \\[Implies] q, equivalent to !p || q.\n"
        "\tImplies[False, q] and Implies[p, True] are True, Implies[True, q]\n"
        "\tis q, and Implies[p, False] is !p.");
}