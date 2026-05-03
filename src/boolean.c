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

void boolean_init(void) {
    symtab_add_builtin("And", builtin_and);
    symtab_add_builtin("Or", builtin_or);
    symtab_add_builtin("Not", builtin_not);
}