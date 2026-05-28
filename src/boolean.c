#include "boolean.h"
#include "symtab.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>

Expr* builtin_not(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    if (arg->type == EXPR_SYMBOL) {
        if (arg->data.symbol == SYM_True) return expr_new_symbol("False");
        if (arg->data.symbol == SYM_False) return expr_new_symbol("True");
    }
    return NULL;
}

Expr* builtin_and(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count == 0) return expr_new_symbol("True");
    if (res->data.function.arg_count == 1) {
        /* And has HoldAll, so the single arg has not been evaluated yet. */
        return evaluate(res->data.function.args[0]);
    }

    bool changed = false;
    Expr** new_args = malloc(sizeof(Expr*) * res->data.function.arg_count);
    size_t new_count = 0;
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* evaluated_arg = evaluate(res->data.function.args[i]);
        if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol == SYM_False) {
            expr_free(evaluated_arg);
            for (size_t j = 0; j < new_count; j++) expr_free(new_args[j]);
            free(new_args);
            return expr_new_symbol("False");
        } else if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol == SYM_True) {
            expr_free(evaluated_arg);
            changed = true;
        } else {
            if (!expr_eq(evaluated_arg, res->data.function.args[i])) changed = true;
            new_args[new_count++] = evaluated_arg;
        }
    }
    if (new_count == 0) {
        free(new_args);
        return expr_new_symbol("True");
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
    if (res->data.function.arg_count == 0) return expr_new_symbol("False");
    if (res->data.function.arg_count == 1) {
        /* Or has HoldAll, so the single arg has not been evaluated yet. */
        return evaluate(res->data.function.args[0]);
    }

    bool changed = false;
    Expr** new_args = malloc(sizeof(Expr*) * res->data.function.arg_count);
    size_t new_count = 0;
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        Expr* evaluated_arg = evaluate(res->data.function.args[i]);
        if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol == SYM_True) {
            expr_free(evaluated_arg);
            for (size_t j = 0; j < new_count; j++) expr_free(new_args[j]);
            free(new_args);
            return expr_new_symbol("True");
        } else if (evaluated_arg->type == EXPR_SYMBOL && evaluated_arg->data.symbol == SYM_False) {
            expr_free(evaluated_arg);
            changed = true;
        } else {
            if (!expr_eq(evaluated_arg, res->data.function.args[i])) changed = true;
            new_args[new_count++] = evaluated_arg;
        }
    }
    if (new_count == 0) {
        free(new_args);
        return expr_new_symbol("False");
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
        if (arg->data.symbol == SYM_True)  return expr_new_integer(1);
        if (arg->data.symbol == SYM_False) return expr_new_integer(0);
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
        if (cond->data.symbol == SYM_True) {
            /* Steal the expr slot so the evaluator's free of res does not
             * double-free what we just returned. */
            res->data.function.args[0] = NULL;
            return expr;
        }
        if (cond->data.symbol == SYM_False) {
            return expr_new_symbol("Undefined");
        }
    }

    /* Nested flattening: ConditionalExpression[ConditionalExpression[e, c1], c2]
     * collapses to ConditionalExpression[e, And[c1, c2]]. */
    if (expr->type == EXPR_FUNCTION &&
        expr->data.function.head->type == EXPR_SYMBOL &&
        expr->data.function.head->data.symbol == SYM_ConditionalExpression &&
        expr->data.function.arg_count == 2) {
        Expr* inner_expr = expr_copy(expr->data.function.args[0]);
        Expr** and_args = malloc(sizeof(Expr*) * 2);
        and_args[0] = expr_copy(expr->data.function.args[1]);
        and_args[1] = expr_copy(cond);
        Expr* combined = expr_new_function(expr_new_symbol("And"), and_args, 2);
        free(and_args);
        Expr* combined_eval = evaluate(combined);
        expr_free(combined);
        Expr** new_args = malloc(sizeof(Expr*) * 2);
        new_args[0] = inner_expr;
        new_args[1] = combined_eval;
        Expr* out = expr_new_function(expr_new_symbol("ConditionalExpression"),
                                      new_args, 2);
        free(new_args);
        return out;
    }

    return NULL;
}

void boolean_init(void) {
    symtab_add_builtin("And", builtin_and);
    symtab_add_builtin("Or", builtin_or);
    symtab_add_builtin("Not", builtin_not);
    symtab_add_builtin("Boole", builtin_boole);
    symtab_add_builtin("ConditionalExpression", builtin_conditional_expression);
}