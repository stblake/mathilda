#include "rat.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "poly.h"
#include "expand.h"
#include "rationalize.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* Together's recursive walker.  Cancel and Together both share the
 * algebraic-generator (radical) pass exported from poly.c, but
 * together_recursive is internal because it embeds the Plus-merging
 * logic specific to combining over a common denominator. */
static Expr* together_recursive(Expr* e);

static Expr* negate_expr(Expr* e) {
    return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(e)}, 2));
}

static bool is_superficially_negative(Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer < 0;
    if (e->type == EXPR_REAL) return e->data.real < 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return n < 0;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol == SYM_Times) {
        if (e->data.function.arg_count > 0) {
            Expr* first = e->data.function.args[0];
            if (first->type == EXPR_INTEGER) return first->data.integer < 0;
            if (first->type == EXPR_REAL) return first->data.real < 0.0;
            if (is_rational(first, &n, &d)) return n < 0;
        }
    }
    return false;
}

static void extract_num_den(Expr* expr, Expr** num_out, Expr** den_out) {
    int64_t n, d;
    if (is_rational(expr, &n, &d)) {
        *num_out = expr_new_integer(n);
        *den_out = expr_new_integer(d);
        return;
    }

    Expr* re; Expr* im;
    if (is_complex(expr, &re, &im)) {
        int64_t rn = 0, rd = 1, in = 0, id = 1;
        is_rational(re, &rn, &rd);
        if (re->type == EXPR_INTEGER) { rn = re->data.integer; rd = 1; }
        is_rational(im, &in, &id);
        if (im->type == EXPR_INTEGER) { in = im->data.integer; id = 1; }

        int64_t common_den = lcm(rd, id);
        
        Expr* d_expr = expr_new_integer(common_den);
        *den_out = d_expr;
        *num_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(expr), expr_copy(d_expr)}, 2));
        return;
    }

    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && (expr->data.function.head->data.symbol == SYM_Power || expr->data.function.head->data.symbol == SYM_Exp)) {
        bool is_exp = expr->data.function.head->data.symbol == SYM_Exp;
        Expr* base = is_exp ? expr_new_symbol("E") : expr->data.function.args[0];
        Expr* exp = is_exp ? expr->data.function.args[0] : expr->data.function.args[1];

        if (exp->type == EXPR_FUNCTION && exp->data.function.head->type == EXPR_SYMBOL && exp->data.function.head->data.symbol == SYM_Plus) {
            size_t count = exp->data.function.arg_count;
            Expr** num_args = malloc(sizeof(Expr*) * count);
            Expr** den_args = malloc(sizeof(Expr*) * count);
            size_t n_c = 0, d_c = 0;
            for (size_t i = 0; i < count; i++) {
                Expr* arg = exp->data.function.args[i];
                if (is_superficially_negative(arg)) {
                    den_args[d_c++] = negate_expr(arg);
                } else {
                    num_args[n_c++] = expr_copy(arg);
                }
            }
            if (n_c == 0) *num_out = expr_new_integer(1);
            else if (n_c == 1) *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), num_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol("Plus"), num_args, n_c));
                *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), p}, 2));
            }
            if (d_c == 0) *den_out = expr_new_integer(1);
            else if (d_c == 1) *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), den_args[0]}, 2));
            else {
                Expr* p = eval_and_free(expr_new_function(expr_new_symbol("Plus"), den_args, d_c));
                *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), p}, 2));
            }
            free(num_args);
            free(den_args);
            if (is_exp) expr_free(base);
            return;
        } else {
            if (is_superficially_negative(exp)) {
                *num_out = expr_new_integer(1);
                *den_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), negate_expr(exp)}, 2));
            } else {
                if (is_exp) {
                    *num_out = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), expr_copy(exp)}, 2));
                } else {
                    *num_out = expr_copy(expr);
                }
                *den_out = expr_new_integer(1);
            }
            if (is_exp) expr_free(base);
            return;
        }
    }

    if (expr->type == EXPR_FUNCTION && expr->data.function.head->type == EXPR_SYMBOL && expr->data.function.head->data.symbol == SYM_Times) {
        size_t count = expr->data.function.arg_count;
        Expr** n_args = malloc(sizeof(Expr*) * count);
        Expr** d_args = malloc(sizeof(Expr*) * count);
        size_t n_c = 0, d_c = 0;
        for (size_t i = 0; i < count; i++) {
            Expr* n_out; Expr* d_out;
            extract_num_den(expr->data.function.args[i], &n_out, &d_out);
            if (!(n_out->type == EXPR_INTEGER && n_out->data.integer == 1)) n_args[n_c++] = n_out;
            else expr_free(n_out);
            if (!(d_out->type == EXPR_INTEGER && d_out->data.integer == 1)) d_args[d_c++] = d_out;
            else expr_free(d_out);
        }
        if (n_c == 0) *num_out = expr_new_integer(1);
        else if (n_c == 1) *num_out = n_args[0];
        else *num_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), n_args, n_c));

        if (d_c == 0) *den_out = expr_new_integer(1);
        else if (d_c == 1) *den_out = d_args[0];
        else *den_out = eval_and_free(expr_new_function(expr_new_symbol("Times"), d_args, d_c));

        free(n_args); free(d_args);
        return;
    }

    *num_out = expr_copy(expr);
    *den_out = expr_new_integer(1);
}

Expr* builtin_numerator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* n; Expr* d;
    extract_num_den(res->data.function.args[0], &n, &d);
    expr_free(d);
    return n;
}

Expr* builtin_denominator(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* n; Expr* d;
    extract_num_den(res->data.function.args[0], &n, &d);
    expr_free(n);
    return d;
}

static Expr* cancel_exact_div_wrapper(Expr* num, Expr* den) {
    if (is_zero_poly(num)) return expr_new_integer(0);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return expr_expand(num);

    Expr* exp_num = expr_expand(num);
    Expr* exp_den = expr_expand(den);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(exp_num, &vars, &v_count, &v_cap);
    collect_variables(exp_den, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    Expr* res = exact_poly_div(exp_num, exp_den, vars, v_count);

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);

    if (res) {
        expr_free(exp_num);
        expr_free(exp_den);
        return res;
    } else {
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){exp_num, t}, 2));
        return r;
    }
}

static Expr* cancel_recursive(Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        if (head == SYM_List || head == SYM_Plus ||
            head == SYM_Equal || head == SYM_Less ||
            head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And ||
            head == SYM_Or || head == SYM_Not) {
            
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = cancel_recursive(e->data.function.args[i]);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
            free(args);
            return ret;
        }
    }
    
    Expr* num; Expr* den;
    extract_num_den(e, &num, &den);
    
    if (den->type == EXPR_INTEGER && den->data.integer == 1) {
        expr_free(den);
        return num;
    }
    
    Expr* g = eval_and_free(expr_new_function(expr_new_symbol("PolynomialGCD"), (Expr*[]){expr_copy(num), expr_copy(den)}, 2));
    
    Expr* new_num = cancel_exact_div_wrapper(num, g);
    Expr* new_den = cancel_exact_div_wrapper(den, g);
    
    if (new_den && new_den->type == EXPR_INTEGER && new_den->data.integer < 0) {
        Expr* t1 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), new_num}, 2));
        Expr* t2 = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), new_den}, 2));
        new_num = expr_expand(t1);
        new_den = expr_expand(t2);
        expr_free(t1);
        expr_free(t2);
    }
    
    Expr* res;
    if (new_den && new_den->type == EXPR_INTEGER && new_den->data.integer == 1) {
        res = new_num;
        expr_free(new_den);
    } else if (new_num && new_den) {
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){new_den, expr_new_integer(-1)}, 2));
        res = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){new_num, inv_den}, 2));
    } else {
        if (new_num) expr_free(new_num);
        if (new_den) expr_free(new_den);
        res = expr_copy(e);
    }
    
    expr_free(g);
    expr_free(num);
    expr_free(den);
    
    return res;
}

Expr* builtin_cancel(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    /* Inexact coefficients block the rational-arithmetic cancellation
     * machinery (e.g. (x^2/9 - y^2/25.) — the 25. defeats Polynomial GCD).
     * Force-rationalise on the way in, numericalise on the way out. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_cancel);
    }

    Expr* arg = res->data.function.args[0];

    /* Algebraic-generator pass: if the input contains a sub-expression
     * u with fractional rational exponents (u^(p/q), q != 1), substitute
     * u -> g^m, run Together's recursive normalisation in g, then
     * substitute back.  together_recursive is preferred over
     * cancel_recursive directly because the substituted form may have
     * Plus terms with different g-exponents that need fraction
     * combination before GCD cancellation can fire. */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* combined = together_recursive(substituted);
            expr_free(substituted);
            Expr* final = poly_subst_radical_from_gen(combined, base, atom, m, gen);
            expr_free(combined);
            expr_free(base);
            expr_free(atom);
            free(gen);
            return final;
        }
    }

    return cancel_recursive(arg);
}

static Expr* together_recursive(Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* head = e->data.function.head->data.symbol;
        
        if (head == SYM_Plus) {
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = together_recursive(e->data.function.args[i]);
            }
            
            Expr** nums = malloc(sizeof(Expr*) * count);
            Expr** dens = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                extract_num_den(args[i], &nums[i], &dens[i]);
            }
            
            Expr* lcm_den = count > 0 ? expr_copy(dens[0]) : expr_new_integer(1);
            for (size_t i = 1; i < count; i++) {
                Expr* call_lcm = expr_new_function(expr_new_symbol("PolynomialLCM"), (Expr*[]){expr_copy(lcm_den), expr_copy(dens[i])}, 2);
                Expr* new_lcm = evaluate(call_lcm);
                expr_free(call_lcm);
                expr_free(lcm_den);
                lcm_den = new_lcm;
            }
            
            Expr** new_nums = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                Expr* q = cancel_exact_div_wrapper(lcm_den, dens[i]);
                new_nums[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){nums[i], q}, 2));
            }
            
            Expr* combined_num = eval_and_free(expr_new_function(expr_new_symbol("Plus"), new_nums, count));
            Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(lcm_den), expr_new_integer(-1)}, 2));
            Expr* combined_expr = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){combined_num, inv_den}, 2));
            
            Expr* res = cancel_recursive(combined_expr);
            expr_free(combined_expr);
            
            for (size_t i = 0; i < count; i++) {
                expr_free(args[i]);
                expr_free(dens[i]); 
            }
            free(args);
            free(nums);
            free(dens);
            free(new_nums);
            
            expr_free(lcm_den);
            
            return res;
        }
        
        if (head == SYM_List ||
            head == SYM_Equal || head == SYM_Less ||
            head == SYM_LessEqual || head == SYM_Greater ||
            head == SYM_GreaterEqual || head == SYM_And ||
            head == SYM_Or || head == SYM_Not) {
            
            size_t count = e->data.function.arg_count;
            Expr** args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                args[i] = together_recursive(e->data.function.args[i]);
            }
            Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
            free(args);
            return ret;
        }
        
        size_t count = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * count);
        for (size_t i = 0; i < count; i++) {
            args[i] = together_recursive(e->data.function.args[i]);
        }
        Expr* ret = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, count));
        free(args);
        
        Expr* cancelled = cancel_recursive(ret);
        expr_free(ret);
        return cancelled;
    }
    
    return cancel_recursive(e);
}

Expr* builtin_together(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    /* Inexact coefficients can't be combined by exact polynomial-LCM
     * machinery; rationalise inputs, run, then re-numericalise. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_together);
    }

    Expr* arg = res->data.function.args[0];

    /* Algebraic-generator pass (see builtin_cancel for the rationale). */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* combined = together_recursive(substituted);
            expr_free(substituted);
            Expr* final = poly_subst_radical_from_gen(combined, base, atom, m, gen);
            expr_free(combined);
            expr_free(base);
            expr_free(atom);
            free(gen);
            return final;
        }
    }

    return together_recursive(arg);
}

void rat_init(void) {
    symtab_add_builtin("Numerator", builtin_numerator);
    symtab_get_def("Numerator")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Denominator", builtin_denominator);
    symtab_get_def("Denominator")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Cancel", builtin_cancel);
    symtab_get_def("Cancel")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_add_builtin("Together", builtin_together);
    symtab_get_def("Together")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
}
