#include "cond.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

Expr* builtin_if(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 4) return NULL;

    Expr* cond = res->data.function.args[0];

    if (cond->type == EXPR_SYMBOL) {
        if (cond->data.symbol == SYM_True) {
            return expr_copy(res->data.function.args[1]);
        } else if (cond->data.symbol == SYM_False) {
            if (argc >= 3) {
                return expr_copy(res->data.function.args[2]);
            } else {
                return expr_new_symbol("Null");
            }
        }
    }

    if (argc == 4) {
        return expr_copy(res->data.function.args[3]);
    }

    return NULL;
}

/*
 * builtin_which -- Which[t1, v1, t2, v2, ...]
 *
 * Which has the HoldAll attribute, so its arguments arrive unevaluated.
 * We evaluate the test expressions one pair at a time:
 *
 *   - if t_i evaluates to True   -> return v_i (held; the outer evaluator
 *                                   will continue evaluating it),
 *   - if t_i evaluates to False  -> drop this pair and continue,
 *   - otherwise                  -> return Which[t_i_eval, v_i, ...remaining...]
 *                                   with the inconclusive test in its
 *                                   evaluated form and the remaining
 *                                   arguments left unevaluated.
 *
 * If every test evaluates to False, Which returns Null. Which[] (no
 * arguments) likewise returns Null. An odd number of arguments is a
 * usage error: we leave the call unevaluated by returning NULL so the
 * user can see what they wrote.
 */
Expr* builtin_which(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) return expr_new_symbol("Null");
    if (argc % 2 != 0) return NULL;

    Expr** args = res->data.function.args;

    for (size_t i = 0; i < argc; i += 2) {
        Expr* test_eval = evaluate(args[i]);
        if (!test_eval) return NULL;

        bool is_true  = (test_eval->type == EXPR_SYMBOL &&
                         test_eval->data.symbol == SYM_True);
        bool is_false = (test_eval->type == EXPR_SYMBOL &&
                         test_eval->data.symbol == SYM_False);

        if (is_true) {
            expr_free(test_eval);
            return expr_copy(args[i + 1]);
        }

        if (is_false) {
            expr_free(test_eval);
            continue;
        }

        /* Inconclusive: return Which[t_i_eval, v_i, ...remaining...]. */
        size_t remaining = argc - i;
        Expr** new_args = malloc(sizeof(Expr*) * remaining);
        if (!new_args) {
            expr_free(test_eval);
            return NULL;
        }
        new_args[0] = test_eval;
        for (size_t j = 1; j < remaining; j++) {
            new_args[j] = expr_copy(args[i + j]);
        }
        Expr* head = expr_copy(res->data.function.head);
        Expr* out = expr_new_function(head, new_args, remaining);
        free(new_args);
        return out;
    }

    return expr_new_symbol("Null");
}

Expr* builtin_trueq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    
    Expr* cond = res->data.function.args[0];
    bool is_true = (cond->type == EXPR_SYMBOL && cond->data.symbol == SYM_True);
    
    return expr_new_symbol(is_true ? "True" : "False");
}

void cond_init(void) {
    symtab_add_builtin("If", builtin_if);
    symtab_get_def("If")->attributes |= ATTR_HOLDREST | ATTR_PROTECTED;

    symtab_add_builtin("Which", builtin_which);
    symtab_get_def("Which")->attributes |= ATTR_HOLDALL | ATTR_PROTECTED;

    symtab_add_builtin("TrueQ", builtin_trueq);
    symtab_get_def("TrueQ")->attributes |= ATTR_PROTECTED;
}