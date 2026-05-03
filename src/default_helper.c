#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>

Expr* get_default_value(Expr* pat_head, int pos, int total_pats) {
    if (!pat_head || pat_head->type != EXPR_SYMBOL) return NULL;
    
    // Construct Default[f, pos, total_pats]
    Expr* args3[3] = { expr_copy(pat_head), expr_new_integer(pos), expr_new_integer(total_pats) };
    Expr* def3 = expr_new_function(expr_new_symbol("Default"), args3, 3);
    Expr* res3 = evaluate(def3);
    if (!expr_eq(res3, def3)) {
        expr_free(def3);
        return res3;
    }
    expr_free(res3);
    expr_free(def3);
    
    // Construct Default[f, pos]
    Expr* args2[2] = { expr_copy(pat_head), expr_new_integer(pos) };
    Expr* def2 = expr_new_function(expr_new_symbol("Default"), args2, 2);
    Expr* res2 = evaluate(def2);
    if (!expr_eq(res2, def2)) {
        expr_free(def2);
        return res2;
    }
    expr_free(res2);
    expr_free(def2);
    
    // Construct Default[f]
    Expr* args1[1] = { expr_copy(pat_head) };
    Expr* def1 = expr_new_function(expr_new_symbol("Default"), args1, 1);
    Expr* res1 = evaluate(def1);
    if (!expr_eq(res1, def1)) {
        expr_free(def1);
        return res1;
    }
    expr_free(res1);
    expr_free(def1);
    
    // Built-in fallbacks
    const char* hn = pat_head->data.symbol;
    if (hn == SYM_Plus) return expr_new_integer(0);
    if (hn == SYM_Times || hn == SYM_Power) return expr_new_integer(1);
    
    return NULL;
}
