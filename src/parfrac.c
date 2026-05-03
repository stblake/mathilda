#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "parse.h"
#include "print.h"
#include "poly.h"
#include "attr.h"
#include "arithmetic.h"
#include "expand.h"
#include "core.h"
#include "rationalize.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

Expr* builtin_apart(Expr* res) {
    if (res->type != EXPR_FUNCTION || (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)) return NULL;

    /* Apart relies on PolynomialQuotient/Remainder and partial-fraction
     * decomposition over rationals — inexact coefficients break both. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_apart);
    }

    Expr* expr = res->data.function.args[0];
    
    // Check if it's an equation, inequality, list, or logical function to thread over
    if (expr->type == EXPR_FUNCTION) {
        const char* head = expr->data.function.head->type == EXPR_SYMBOL ? expr->data.function.head->data.symbol : "";
        if (strcmp(head, "List") == 0 || strcmp(head, "Equal") == 0 || strcmp(head, "Unequal") == 0 ||
            strcmp(head, "Less") == 0 || strcmp(head, "LessEqual") == 0 || strcmp(head, "Greater") == 0 || 
            strcmp(head, "GreaterEqual") == 0 || strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 || 
            strcmp(head, "Not") == 0) {
            Expr** args = malloc(sizeof(Expr*) * expr->data.function.arg_count);
            for (size_t i = 0; i < expr->data.function.arg_count; i++) {
                Expr** ap_args = res->data.function.arg_count == 1 ? (Expr*[]){expr_copy(expr->data.function.args[i])} : (Expr*[]){expr_copy(expr->data.function.args[i]), expr_copy(res->data.function.args[1])};
                args[i] = eval_and_free(expr_new_function(expr_new_symbol("Apart"), ap_args, res->data.function.arg_count));
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(expr->data.function.head), args, expr->data.function.arg_count));
            free(args);
            return ret;
        }
    }
    
    Expr* together = eval_and_free(expr_new_function(expr_new_symbol("Together"), (Expr*[]){expr_copy(expr)}, 1));
    
    Expr* var = NULL;
    if (res->data.function.arg_count == 2) {
        var = expr_copy(res->data.function.args[1]);
    } else {
        size_t v_count = 0, v_cap = 16;
        Expr** vars = malloc(sizeof(Expr*) * v_cap);
        collect_variables(together, &vars, &v_count, &v_cap);
        if (v_count > 0) {
            qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);
            var = expr_copy(vars[v_count - 1]);
        }
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
    }
    
    if (!var) {
        return eval_and_free(expr_new_function(expr_new_symbol("Expand"), (Expr*[]){together}, 1));
    }
    
    Expr* N = eval_and_free(expr_new_function(expr_new_symbol("Numerator"), (Expr*[]){expr_copy(together)}, 1));
    Expr* D = eval_and_free(expr_new_function(expr_new_symbol("Denominator"), (Expr*[]){expr_copy(together)}, 1));

    /* Apart's matrix-construction step assumes both N and D are polynomials
     * in `var`. When Together produces a numerator or denominator with
     * non-integer powers of var (e.g. Together[1/(y^(2/3)-1/y^(1/3))] =
     * y^(1/3)/(y-1)), get_coeff(y^(1/3), y, 0) returns 0 for every row,
     * yielding a zero solution and a wrong answer of 0. Detect this case
     * up front and return the Together'd form unchanged -- partial-fraction
     * decomposition is undefined when the numerator is not polynomial in
     * the chosen variable. */
    Expr* npq_args[2] = { expr_copy(N), expr_copy(var) };
    Expr* npq = eval_and_free(expr_new_function(expr_new_symbol("PolynomialQ"), npq_args, 2));
    Expr* dpq_args[2] = { expr_copy(D), expr_copy(var) };
    Expr* dpq = eval_and_free(expr_new_function(expr_new_symbol("PolynomialQ"), dpq_args, 2));
    bool n_poly = (npq->type == EXPR_SYMBOL && npq->data.symbol == SYM_True);
    bool d_poly = (dpq->type == EXPR_SYMBOL && dpq->data.symbol == SYM_True);
    expr_free(npq); expr_free(dpq);
    if (!n_poly || !d_poly) {
        expr_free(N); expr_free(D); expr_free(var);
        return together;
    }

    expr_free(together);

    if (get_degree_poly(D, var) == 0) {
        Expr* expanded = eval_and_free(expr_new_function(expr_new_symbol("Expand"), (Expr*[]){eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){N, eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){D, expr_new_integer(-1)}, 2))}, 2))}, 1));
        expr_free(var);
        return expanded;
    }
    
    Expr* Q = eval_and_free(expr_new_function(expr_new_symbol("PolynomialQuotient"), (Expr*[]){expr_copy(N), expr_copy(D), expr_copy(var)}, 3));
    Expr* R = eval_and_free(expr_new_function(expr_new_symbol("PolynomialRemainder"), (Expr*[]){expr_copy(N), expr_copy(D), expr_copy(var)}, 3));
    
    Expr* f_den = eval_and_free(expr_new_function(expr_new_symbol("Factor"), (Expr*[]){expr_copy(D)}, 1));
    expr_free(N); expr_free(D);
    
    size_t num_args = 1;
    Expr** args = &f_den;
    if (f_den->type == EXPR_FUNCTION && f_den->data.function.head->data.symbol == SYM_Times) {
        num_args = f_den->data.function.arg_count;
        args = f_den->data.function.args;
    }
    
    Expr* C = expr_new_integer(1);
    size_t m = 0;
    Expr** bases = malloc(sizeof(Expr*) * num_args);
    int64_t* ks = malloc(sizeof(int64_t) * num_args);
    
    for (size_t i = 0; i < num_args; i++) {
        Expr* factor = args[i];
        Expr* base = factor;
        int64_t k = 1;
        if (factor->type == EXPR_FUNCTION && factor->data.function.head->data.symbol == SYM_Power) {
            base = factor->data.function.args[0];
            if (factor->data.function.args[1]->type == EXPR_INTEGER) {
                k = factor->data.function.args[1]->data.integer;
            }
        }
        
        if (get_degree_poly(base, var) == 0) {
            Expr* term = (factor == base && k != 1) ? eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), expr_new_integer(k)}, 2)) : expr_copy(factor);
            C = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){C, term}, 2));
        } else {
            bases[m] = expr_copy(base);
            ks[m] = k;
            m++;
        }
    }
    
    if (m == 0) {
        Expr* C_inv = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){C, expr_new_integer(-1)}, 2));
        Expr* R_term = eval_and_free(expr_new_function(expr_new_symbol("Expand"), (Expr*[]){eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){R, C_inv}, 2))}, 1));
        Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){Q, R_term}, 2));
        expr_free(f_den); free(bases); free(ks); expr_free(var);
        return res;
    }
    
    Expr* D_main = expr_new_integer(1);
    for (size_t i = 0; i < m; i++) {
        Expr* term = ks[i] == 1 ? expr_copy(bases[i]) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(bases[i]), expr_new_integer(ks[i])}, 2));
        D_main = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){D_main, term}, 2));
    }
    
    Expr* exp_D_main = expr_expand(D_main);
    int S = get_degree_poly(exp_D_main, var);
    expr_free(exp_D_main);
    expr_free(D_main);
    
    Expr*** M = malloc(sizeof(Expr**) * S);
    for (int row = 0; row < S; row++) {
        M[row] = malloc(sizeof(Expr*) * (S + 1));
        for (int col = 0; col <= S; col++) M[row][col] = expr_new_integer(0);
    }
    
    int col = 0;
    for (size_t i = 0; i < m; i++) {
        Expr* exp_base = expr_expand(bases[i]);
        int d_B = get_degree_poly(exp_base, var);
        expr_free(exp_base);
        
        for (int j = 1; j <= ks[i]; j++) {
            for (int r = 0; r < d_B; r++) {
                Expr* M_ij = expr_new_integer(1);
                for (size_t p = 0; p < m; p++) {
                    int exp = (p == i) ? ks[p] - j : ks[p];
                    if (exp > 0) {
                        Expr* term = (exp == 1) ? expr_copy(bases[p]) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(bases[p]), expr_new_integer(exp)}, 2));
                        M_ij = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){M_ij, term}, 2));
                    }
                }
                
                Expr* var_pow = (r == 0) ? expr_new_integer(1) : ((r == 1) ? expr_copy(var) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(var), expr_new_integer(r)}, 2)));
                Expr* term_mult = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){M_ij, var_pow}, 2));
                Expr* expanded_term = expr_expand(term_mult);
                expr_free(term_mult);
                
                for (int row = 0; row < S; row++) {
                    expr_free(M[row][col]);
                    M[row][col] = get_coeff(expanded_term, var, row);
                }
                expr_free(expanded_term);
                col++;
            }
        }
    }
    
    Expr* C_inv = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){C, expr_new_integer(-1)}, 2));
    Expr* tmp_R = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){R, C_inv}, 2));
    Expr* R_rhs = expr_expand(tmp_R);
    expr_free(tmp_R);
    for (int row = 0; row < S; row++) {
        expr_free(M[row][S]);
        M[row][S] = get_coeff(R_rhs, var, row);
    }
    expr_free(R_rhs);
    
    Expr** row_args = malloc(sizeof(Expr*) * S);
    for (int row = 0; row < S; row++) {
        row_args[row] = expr_new_function(expr_new_symbol("List"), M[row], S + 1);
        free(M[row]);
    }
    free(M);
    
    Expr* matrix_expr = expr_new_function(expr_new_symbol("List"), row_args, S);
    free(row_args);
    
    Expr* reduced = eval_and_free(expr_new_function(expr_new_symbol("RowReduce"), (Expr*[]){matrix_expr}, 1));
    
    Expr* pfrac_sum = expr_new_integer(0);
    if (reduced->type == EXPR_FUNCTION && reduced->data.function.head->data.symbol == SYM_List && reduced->data.function.arg_count == (size_t)S) {
        col = 0;
        for (size_t i = 0; i < m; i++) {
            Expr* exp_base = expr_expand(bases[i]);
            int d_B = get_degree_poly(exp_base, var);
            expr_free(exp_base);
            
            for (int j = 1; j <= ks[i]; j++) {
                Expr* A_ij = expr_new_integer(0);
                for (int r = 0; r < d_B; r++) {
                    Expr* coeff = expr_copy(reduced->data.function.args[col]->data.function.args[S]);
                    
                    Expr* fact_coeff = eval_and_free(expr_new_function(expr_new_symbol("Factor"), (Expr*[]){coeff}, 1));
                    
                    Expr* var_pow = (r == 0) ? expr_new_integer(1) : ((r == 1) ? expr_copy(var) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(var), expr_new_integer(r)}, 2)));
                    Expr* term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){fact_coeff, var_pow}, 2));
                    
                    A_ij = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){A_ij, term}, 2));
                    col++;
                }
                
                Expr* denom = (j == 1) ? expr_copy(bases[i]) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(bases[i]), expr_new_integer(j)}, 2));
                Expr* fraction = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){A_ij, eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){denom, expr_new_integer(-1)}, 2))}, 2));
                
                pfrac_sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){pfrac_sum, fraction}, 2));
            }
        }
    } else {
        pfrac_sum = expr_new_integer(0); 
    }
    expr_free(reduced);
    
    Expr* final_res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){Q, pfrac_sum}, 2));
    
    for (size_t i = 0; i < m; i++) expr_free(bases[i]);
    free(bases); free(ks);
    expr_free(f_den); expr_free(var);
    
    return final_res;
}

void parfrac_init(void) {
    symtab_add_builtin("Apart", builtin_apart);
    symtab_get_def("Apart")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
}
