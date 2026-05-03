#include "cond.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "sym_names.h"
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

Expr* builtin_trueq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    
    Expr* cond = res->data.function.args[0];
    bool is_true = (cond->type == EXPR_SYMBOL && cond->data.symbol == SYM_True);
    
    return expr_new_symbol(is_true ? "True" : "False");
}

void cond_init(void) {
    symtab_add_builtin("If", builtin_if);
    symtab_get_def("If")->attributes |= ATTR_HOLDREST | ATTR_PROTECTED;
    
    symtab_add_builtin("TrueQ", builtin_trueq);
    symtab_get_def("TrueQ")->attributes |= ATTR_PROTECTED;
}