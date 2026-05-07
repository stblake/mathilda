#include "facpoly.h"
#include "poly.h"
#include "eval.h"
#include "expand.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "internal.h"
#include "parse.h"
#include "rationalize.h"
#include "zupoly.h"
#include "bpoly.h"
#include "mpoly.h"
#include "mvfactor.h"
#include "mvfactor3.h"
#include "sym_names.h"
#include "qafactor.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <gmp.h>

static bool is_constant_1(Expr* e);
Expr* bz_factor_to_expr(Expr* P, Expr* var);

static bool is_constant_1(Expr* e) {
    if (!e) return false;
    Expr* ev = eval_and_free(expr_copy(e));
    bool res = (ev->type == EXPR_INTEGER && ev->data.integer == 1);
    expr_free(ev);
    return res;
}

static Expr* poly_deriv(Expr* p, Expr* x) {
    if (!contains_any_symbol_from(p, x)) return expr_new_integer(0);
    if (expr_eq(p, x)) return expr_new_integer(1);

    if (p->type == EXPR_FUNCTION) {
        const char* head = p->data.function.head->type == EXPR_SYMBOL ? p->data.function.head->data.symbol : "";
        if (strcmp(head, "Plus") == 0) {
            Expr** args = malloc(sizeof(Expr*) * p->data.function.arg_count);
            for (size_t i=0; i<p->data.function.arg_count; i++) args[i] = poly_deriv(p->data.function.args[i], x);
            Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), args, p->data.function.arg_count));
            free(args);
            return res;
        } else if (strcmp(head, "Times") == 0) {
            size_t count = p->data.function.arg_count;
            Expr** sum_args = malloc(sizeof(Expr*) * count);
            for (size_t i = 0; i < count; i++) {
                Expr** prod_args = malloc(sizeof(Expr*) * count);
                for (size_t j = 0; j < count; j++) {
                    if (i == j) prod_args[j] = poly_deriv(p->data.function.args[j], x);
                    else prod_args[j] = expr_copy(p->data.function.args[j]);
                }
                sum_args[i] = eval_and_free(expr_new_function(expr_new_symbol("Times"), prod_args, count));
                free(prod_args);
            }
            Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), sum_args, count));
            free(sum_args);
            return res;
        } else if (strcmp(head, "Power") == 0 && p->data.function.arg_count == 2) {
            Expr* base = p->data.function.args[0];
            Expr* exp = p->data.function.args[1];
            if (!contains_any_symbol_from(exp, x)) {
                Expr* d_base = poly_deriv(base, x);
                Expr** p_args1 = malloc(sizeof(Expr*) * 2);
                p_args1[0] = expr_copy(exp); p_args1[1] = expr_new_integer(-1);
                Expr* exp_minus_1 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), p_args1, 2));
                free(p_args1);

                Expr** p_args2 = malloc(sizeof(Expr*) * 2);
                p_args2[0] = expr_copy(base); p_args2[1] = exp_minus_1;
                Expr* pow_term = eval_and_free(expr_new_function(expr_new_symbol("Power"), p_args2, 2));
                free(p_args2);

                Expr** t_args = malloc(sizeof(Expr*) * 3);
                t_args[0] = expr_copy(exp); t_args[1] = pow_term; t_args[2] = d_base;
                Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 3));
                free(t_args);
                return res;
            }
        }
    }
    return expr_new_integer(0);
}

/* Forward declarations used by the F4 cheap squarefree pre-check. */
static bool univariate_squarefree(Expr* u, Expr* var);
static Expr* eval_others_at_alpha(Expr* P, Expr** vars, size_t v_count,
                                  size_t main_idx, int64_t alpha);

/* F4 Phase 1 cheap squarefree pre-check.
 *
 * Returns true iff `pp` is provably squarefree as a polynomial in
 * `vars[var_count - 1]` (treating the other vars as parameters).
 * Returns false if no such proof was obtained -- the polynomial
 * may still be squarefree, and the caller should fall through to
 * the full Yun-style algorithm.
 *
 * Soundness: invoked AFTER content extraction in the main variable,
 * so any repeated factor of pp involves the main variable nontrivially
 * (a constant-in-x repeated factor would have divided content).
 * Substituting integer values for the other variables produces a
 * univariate image whose gcd-with-derivative is cheap to compute.
 * If the image is squarefree at any alpha that preserves the
 * leading-x degree, then pp itself is squarefree in x: a repeated
 * factor f(x, others)^2 of pp would specialise to f(x, alpha)^2 in
 * the image, and unless alpha is a root of the leading-x coefficient
 * of f the image is non-squarefree.
 *
 * Only used when var_count >= 2 because for univariate inputs the
 * GCD is unavoidable and a separate pre-check would only pessimise
 * the non-squarefree case. */
static bool sqfree_cheap_check(Expr* pp, Expr** vars, size_t var_count) {
    if (var_count < 2) return false;
    Expr* x = vars[var_count - 1];
    int deg_x = get_degree_poly(pp, x);
    if (deg_x < 2) return true;  /* deg <= 1 in x is squarefree after content stripping */

    static const int64_t alphas[] = {1, -1, 2, -2, 3, -3, 4};
    size_t alpha_count = sizeof(alphas) / sizeof(alphas[0]);

    for (size_t pi = 0; pi < alpha_count; pi++) {
        Expr* image = eval_others_at_alpha(pp, vars, var_count,
                                           var_count - 1, alphas[pi]);
        if (!image) continue;
        if (get_degree_poly(image, x) != deg_x) {
            expr_free(image);
            continue;
        }
        bool sqf = univariate_squarefree(image, x);
        expr_free(image);
        if (sqf) return true;
    }
    return false;
}

static Expr* factor_square_free_poly(Expr* P, Expr** vars, size_t var_count) {
    if (var_count == 0 || is_zero_poly(P)) return expr_copy(P);
    Expr* x = vars[var_count - 1];
    if (get_degree_poly(P, x) == 0) {
        return factor_square_free_poly(P, vars, var_count - 1);
    }

    Expr* cont = poly_content(P, vars, var_count);
    Expr* pp = exact_poly_div(P, cont, vars, var_count);
    if (!pp) {
        pp = expr_copy(P);
        expr_free(cont);
        cont = expr_new_integer(1);
    }

    Expr* res_cont = factor_square_free_poly(cont, vars, var_count - 1);
    expr_free(cont);

    if (get_degree_poly(pp, x) == 0) {
        Expr** t_args = malloc(sizeof(Expr*) * 2);
        t_args[0] = res_cont; t_args[1] = pp;
        Expr* ret = eval_and_free(expr_new_function(expr_new_symbol("Times"), t_args, 2));
        free(t_args);
        return ret;
    }

    /* F4 fast path: skip the expensive multivariate gcd(pp, pp') when a
     * cheap univariate-substitution probe proves pp is squarefree in x. */
    if (sqfree_cheap_check(pp, vars, var_count)) {
        Expr* res_pp = expr_expand(pp);
        expr_free(pp);

        Expr* expanded_res_cont = expr_expand(res_cont);
        Expr** pr_args = malloc(sizeof(Expr*) * 2);
        pr_args[0] = expanded_res_cont;
        pr_args[1] = expr_copy(res_pp);
        Expr* prod = eval_and_free(expr_new_function(expr_new_symbol("Times"), pr_args, 2));
        free(pr_args);
        Expr* expanded_prod = expr_expand(prod);
        expr_free(prod);
        Expr* missing = exact_poly_div(P, expanded_prod, vars, var_count);
        expr_free(expanded_prod);
        if (!missing) missing = expr_new_integer(1);

        Expr** final_args = malloc(sizeof(Expr*) * 3);
        final_args[0] = missing;
        final_args[1] = res_cont;
        final_args[2] = res_pp;
        Expr* final_res = eval_and_free(expr_new_function(expr_new_symbol("Times"), final_args, 3));
        free(final_args);
        return final_res;
    }

    Expr* A = pp;
    Expr* raw_A_prime = poly_deriv(A, x);
    Expr* A_prime = expr_expand(raw_A_prime);
    expr_free(raw_A_prime);
    Expr* B = poly_gcd_internal(A, A_prime, vars, var_count);
    Expr* C = exact_poly_div(A, B, vars, var_count);
    if (!C) C = expr_copy(A);
    Expr* A_prime_div_B = exact_poly_div(A_prime, B, vars, var_count);
    if (!A_prime_div_B) A_prime_div_B = expr_new_integer(0);
    
    Expr* raw_C_prime = poly_deriv(C, x);
    Expr* C_prime = expr_expand(raw_C_prime);
    expr_free(raw_C_prime);
    
    Expr** m_args = malloc(sizeof(Expr*) * 2);
    m_args[0] = expr_new_integer(-1); m_args[1] = C_prime;
    Expr* neg_C_prime = eval_and_free(expr_new_function(expr_new_symbol("Times"), m_args, 2));
    free(m_args);
    
    Expr** p_args = malloc(sizeof(Expr*) * 2);
    p_args[0] = A_prime_div_B; p_args[1] = neg_C_prime;
    Expr* evaluated_p = eval_and_free(expr_new_function(expr_new_symbol("Plus"), p_args, 2));
    free(p_args);
    Expr* D = expr_expand(evaluated_p);
    expr_free(evaluated_p);
    
    size_t i = 1;
    size_t factors_count = 0;
    size_t factors_cap = 16;
    Expr** factors_data = malloc(sizeof(Expr*) * factors_cap);
    
    int max_i = get_degree_poly(pp, x);

    while (!is_zero_poly(C) && get_degree_poly(C, x) > 0 && !is_constant_1(C)) {
        if (i > (size_t)max_i + 2) break;
        
        Expr* P_i = poly_gcd_internal(C, D, vars, var_count);
        if (!is_constant_1(P_i)) {
            if (factors_count == factors_cap) { factors_cap *= 2; factors_data = realloc(factors_data, sizeof(Expr*) * factors_cap); }
            if (i == 1) {
                factors_data[factors_count++] = expr_copy(P_i);
            } else {
                Expr** pow_args = malloc(sizeof(Expr*) * 2);
                pow_args[0] = expr_copy(P_i); pow_args[1] = expr_new_integer(i);
                Expr* pow = eval_and_free(expr_new_function(expr_new_symbol("Power"), pow_args, 2));
                factors_data[factors_count++] = pow;
                free(pow_args);
            }
        }
        
        Expr* C_next = exact_poly_div(C, P_i, vars, var_count);
        if (!C_next || is_constant_1(C_next)) {
            if (C_next) expr_free(C_next);
            expr_free(P_i);
            break;
        }
        
        Expr* raw_C_next_prime = poly_deriv(C_next, x);
        Expr* C_next_prime = expr_expand(raw_C_next_prime);
        expr_free(raw_C_next_prime);
        Expr* D_div_P_i = exact_poly_div(D, P_i, vars, var_count);
        if (!D_div_P_i) {
            expr_free(C_next_prime);
            expr_free(C_next);
            expr_free(P_i);
            break;
        }
        
        Expr** m2_args = malloc(sizeof(Expr*) * 2);
        m2_args[0] = expr_new_integer(-1); m2_args[1] = C_next_prime;
        Expr* neg_C_next_prime = eval_and_free(expr_new_function(expr_new_symbol("Times"), m2_args, 2));
        free(m2_args);
        
        Expr** p2_args = malloc(sizeof(Expr*) * 2);
        p2_args[0] = D_div_P_i; p2_args[1] = neg_C_next_prime;
        Expr* evaluated_p2 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), p2_args, 2));
        free(p2_args);
        Expr* D_next = expr_expand(evaluated_p2);
        expr_free(evaluated_p2);
        
        expr_free(C); C = C_next;
        expr_free(D); D = D_next;
        expr_free(P_i);
        i++;
    }
    
    expr_free(A); 
    expr_free(A_prime);
    expr_free(B);
    expr_free(C);
    expr_free(D);
    
    Expr* res_pp;
    if (factors_count == 0) res_pp = expr_new_integer(1);
    else if (factors_count == 1) res_pp = expr_copy(factors_data[0]);
    else {
        Expr** f_args = malloc(sizeof(Expr*) * factors_count);
        for(size_t k=0; k<factors_count; k++) f_args[k] = expr_copy(factors_data[k]);
        res_pp = eval_and_free(expr_new_function(expr_new_symbol("Times"), f_args, factors_count));
        free(f_args);
    }
    for(size_t k=0; k<factors_count; k++) expr_free(factors_data[k]);
    free(factors_data);
    
    Expr* expanded_res_pp = expr_expand(res_pp);
    Expr* expanded_res_cont = expr_expand(res_cont);
    
    Expr** pr_args = malloc(sizeof(Expr*) * 2);
    pr_args[0] = expanded_res_cont; pr_args[1] = expanded_res_pp;
    Expr* prod = eval_and_free(expr_new_function(expr_new_symbol("Times"), pr_args, 2));
    free(pr_args);
    
    Expr* expanded_prod = expr_expand(prod);
    expr_free(prod);
    
    Expr* missing = exact_poly_div(P, expanded_prod, vars, var_count);
    expr_free(expanded_prod);
    
    if (!missing) missing = expr_new_integer(1);
    
    Expr** final_args = malloc(sizeof(Expr*) * 3);
    final_args[0] = missing; final_args[1] = res_cont; final_args[2] = res_pp;
    Expr* final_res = eval_and_free(expr_new_function(expr_new_symbol("Times"), final_args, 3));
    free(final_args);
    return final_res;
}

static Expr* factor_square_free_dispatcher(Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    const char* head = e->data.function.head->type == EXPR_SYMBOL ? e->data.function.head->data.symbol : "";

    if (strcmp(head, "List") == 0 || strcmp(head, "Equal") == 0 || strcmp(head, "Less") == 0 || 
        strcmp(head, "LessEqual") == 0 || strcmp(head, "Greater") == 0 || strcmp(head, "GreaterEqual") == 0 ||
        strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 || strcmp(head, "Not") == 0) {
        Expr** args = malloc(sizeof(Expr*) * e->data.function.arg_count);
        for (size_t i = 0; i < e->data.function.arg_count; i++) args[i] = factor_square_free_dispatcher(e->data.function.args[i]);
        Expr* res = eval_and_free(expr_new_function(expr_copy(e->data.function.head), args, e->data.function.arg_count));
        free(args);
        return res;
    }

    if (strcmp(head, "Power") == 0 && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp = e->data.function.args[1];
        if (exp->type == EXPR_INTEGER || is_rational(exp, NULL, NULL)) {
            Expr* f_base = factor_square_free_dispatcher(base);
            Expr** p_args = malloc(sizeof(Expr*) * 2);
            p_args[0] = f_base; p_args[1] = expr_copy(exp);
            Expr* res = eval_and_free(expr_new_function(expr_copy(e->data.function.head), p_args, 2));
            free(p_args);
            return res;
        }
        return expr_copy(e);
    }

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(e, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    if (!is_polynomial(e, vars, v_count)) {
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        return expr_copy(e);
    }
    
    Expr* expanded_e = expr_expand(e);
    Expr* factored = factor_square_free_poly(expanded_e, vars, v_count);
    expr_free(expanded_e);
    
    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);
    
    return factored;
}

Expr* builtin_factorsquarefree(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    /* The square-free factorisation pipeline (PolynomialGCD against the
     * derivative) requires rational coefficients. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_factorsquarefree);
    }

    return factor_square_free_dispatcher(res->data.function.args[0]);
}

static int64_t int_root(int64_t n, int64_t k) {
    if (n < 0) {
        if (k % 2 == 0) return -1;
        int64_t r = int_root(-n, k);
        return r == -1 ? -1 : -r;
    }
    if (n == 0) return 0;
    if (n == 1) return 1;
    int64_t r = round(pow(n, 1.0 / k));
    int64_t p = 1;
    for (int i=0; i<k; i++) p *= r;
    if (p == n) return r;
    return -1;
}

/* Recursively absorb a single Expr into the running monomial state.
 * Handles INTEGER (multiplies into c), SYMBOL (adds variable), Power
 * (adds variable with exponent), and nested Times (recurses).  Bumps
 * `*cap` and reallocs vars/exps if needed.  Returns false if `e` is
 * not a recognised polynomial-monomial factor (e.g., Plus, non-integer
 * Power exponent, etc.); on false the monomial is unusable. */
static bool extract_monomial_walk(Expr* e, int64_t* c,
                                  Expr*** vars, int64_t** exps,
                                  size_t* count, size_t* cap) {
    if (e->type == EXPR_INTEGER) {
        *c *= e->data.integer;
        return true;
    }
    if (e->type == EXPR_SYMBOL) {
        if (*count == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *vars = realloc(*vars, sizeof(Expr*)   * (*cap));
            *exps = realloc(*exps, sizeof(int64_t) * (*cap));
        }
        (*vars)[*count] = e;
        (*exps)[*count] = 1;
        (*count)++;
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER) {
        if (*count == *cap) {
            *cap = *cap ? *cap * 2 : 4;
            *vars = realloc(*vars, sizeof(Expr*)   * (*cap));
            *exps = realloc(*exps, sizeof(int64_t) * (*cap));
        }
        (*vars)[*count] = e->data.function.args[0];
        (*exps)[*count] = e->data.function.args[1]->data.integer;
        (*count)++;
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!extract_monomial_walk(e->data.function.args[i],
                                       c, vars, exps, count, cap)) {
                return false;
            }
        }
        return true;
    }
    return false;
}

static void extract_monomial(Expr* e, int64_t* c, Expr*** vars, int64_t** exps, size_t* count) {
    *c = 1; *count = 0;
    *vars = NULL; *exps = NULL;
    size_t cap = 0;
    if (!extract_monomial_walk(e, c, vars, exps, count, &cap)) {
        /* On failure, leave outputs in a safe state: caller still needs
         * to free *vars and *exps (which may be partially populated). */
    }
}

static Expr* make_binomial_sum(Expr* A, Expr* B, int k, bool is_diff) {
    Expr** terms = malloc(sizeof(Expr*) * k);
    for (int i=0; i<k; i++) {
        Expr* pA = NULL;
        if (k-1-i == 0) pA = expr_new_integer(1);
        else if (k-1-i == 1) pA = expr_copy(A);
        else pA = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(A), expr_new_integer(k-1-i)}, 2));
        
        Expr* pB = NULL;
        if (i == 0) pB = expr_new_integer(1);
        else if (i == 1) pB = expr_copy(B);
        else pB = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(B), expr_new_integer(i)}, 2));
        
        Expr* term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){pA, pB}, 2));
        
        if (!is_diff && (i % 2 != 0)) { 
            term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), term}, 2));
        }
        terms[i] = term;
    }
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Plus"), terms, k));
    free(terms);
    return res;
}

static Expr* heuristic_factor(Expr* P);

static Expr* factor_binomial(Expr* P) {
    if (P->type != EXPR_FUNCTION || P->data.function.head->data.symbol != SYM_Plus) return NULL;
    if (P->data.function.arg_count != 2) return NULL;
    
    int64_t c1, c2;
    Expr **v1 = NULL, **v2 = NULL;
    int64_t *e1 = NULL, *e2 = NULL;
    size_t count1, count2;
    
    extract_monomial(P->data.function.args[0], &c1, &v1, &e1, &count1);
    extract_monomial(P->data.function.args[1], &c2, &v2, &e2, &count2);
    
    int64_t K = 0;
    for(size_t i=0; i<count1; i++) K = gcd(K, e1[i]);
    for(size_t i=0; i<count2; i++) K = gcd(K, e2[i]);
    
    if (K < 2) {
        free(v1); free(e1); free(v2); free(e2);
        return NULL;
    }
    
    for (int64_t k = K; k >= 2; k--) {
        if (K % k != 0) continue;
        
        int64_t a = int_root(llabs(c1), k);
        int64_t b = int_root(llabs(c2), k);
        if (a != -1 && b != -1) {
            Expr** u_args = malloc(sizeof(Expr*) * (count1 + 1));
            size_t u_c = 0;
            if (a != 1) u_args[u_c++] = expr_new_integer(a);
            for(size_t i=0; i<count1; i++) {
                if (e1[i]/k == 1) u_args[u_c++] = expr_copy(v1[i]);
                else u_args[u_c++] = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(v1[i]), expr_new_integer(e1[i]/k)}, 2));
            }
            Expr* U = (u_c == 0) ? expr_new_integer(1) : ((u_c == 1) ? u_args[0] : eval_and_free(expr_new_function(expr_new_symbol("Times"), u_args, u_c)));
            free(u_args);
            
            Expr** v_args = malloc(sizeof(Expr*) * (count2 + 1));
            size_t v_c = 0;
            if (b != 1) v_args[v_c++] = expr_new_integer(b);
            for(size_t i=0; i<count2; i++) {
                if (e2[i]/k == 1) v_args[v_c++] = expr_copy(v2[i]);
                else v_args[v_c++] = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(v2[i]), expr_new_integer(e2[i]/k)}, 2));
            }
            Expr* V = (v_c == 0) ? expr_new_integer(1) : ((v_c == 1) ? v_args[0] : eval_and_free(expr_new_function(expr_new_symbol("Times"), v_args, v_c)));
            free(v_args);
            
            bool is_diff = (c1 > 0 && c2 < 0) || (c1 < 0 && c2 > 0);
            if (c1 < 0 && c2 < 0) is_diff = false; 
            
            if (c1 < 0 && c2 > 0) {
                Expr* tmp = U; U = V; V = tmp;
            }
            
            Expr* factor1 = NULL;
            Expr* factor2 = NULL;
            
            if (is_diff) {
                if (k % 2 == 0) {
                    Expr* u_half = (k/2 == 1) ? expr_copy(U) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(U), expr_new_integer(k/2)}, 2));
                    Expr* v_half = (k/2 == 1) ? expr_copy(V) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(V), expr_new_integer(k/2)}, 2));
                    Expr* neg_v_half = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(v_half)}, 2));
                    
                    factor1 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(u_half), neg_v_half}, 2));
                    factor2 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){u_half, v_half}, 2));
                } else {
                    Expr* neg_v = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(V)}, 2));
                    factor1 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(U), neg_v}, 2));
                    factor2 = make_binomial_sum(U, V, k, true);
                }
            } else {
                if (k % 2 != 0) {
                    factor1 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(U), expr_copy(V)}, 2));
                    factor2 = make_binomial_sum(U, V, k, false);
                } else {
                    int odd_f = 1;
                    for (int f = 3; f <= k; f+=2) {
                        if (k % f == 0) { odd_f = f; break; }
                    }
                    if (odd_f > 1) {
                        int m = k / odd_f;
                        Expr* u_m = (m == 1) ? expr_copy(U) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(U), expr_new_integer(m)}, 2));
                        Expr* v_m = (m == 1) ? expr_copy(V) : eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(V), expr_new_integer(m)}, 2));
                        factor1 = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(u_m), expr_copy(v_m)}, 2));
                        factor2 = make_binomial_sum(u_m, v_m, odd_f, false);
                        expr_free(u_m); expr_free(v_m);
                    }
                }
            }
            
            expr_free(U); expr_free(V);
            
            if (factor1 && factor2) {
                free(v1); free(e1); free(v2); free(e2);
                Expr* expanded1 = expr_expand(factor1);
                Expr* expanded2 = expr_expand(factor2);
                expr_free(factor1); expr_free(factor2);
                
                Expr* res1 = heuristic_factor(expanded1);
                Expr* res2 = heuristic_factor(expanded2);
                expr_free(expanded1); expr_free(expanded2);
                
                if (c1 < 0 && c2 < 0) {
                    return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), res1, res2}, 3));
                }
                return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){res1, res2}, 2));
            }
        }
    }
    
    free(v1); free(e1); free(v2); free(e2);
    return NULL;
}

static Expr* factor_degree_one(Expr* P, Expr** vars, size_t v_count) {
    for (size_t i = 0; i < v_count; i++) {
        int d = get_degree_poly(P, vars[i]);
        if (d == 1) {
            Expr* L = get_coeff(P, vars[i], 1);
            Expr* C = get_coeff(P, vars[i], 0);
            Expr* G = poly_gcd_internal(L, C, vars, v_count);
            if (!(G->type == EXPR_INTEGER && G->data.integer == 1)) {
                Expr* L_G = exact_poly_div(L, G, vars, v_count);
                Expr* C_G = exact_poly_div(C, G, vars, v_count);
                
                Expr* L_G_v = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){L_G, expr_copy(vars[i])}, 2));
                Expr* prim = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){L_G_v, C_G}, 2));
                Expr* exp_prim = expr_expand(prim);
                expr_free(prim);
                
                Expr* f_G = heuristic_factor(G);
                Expr* f_prim = heuristic_factor(exp_prim);
                
                expr_free(L); expr_free(C); expr_free(G); expr_free(exp_prim);
                
                return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){f_G, f_prim}, 2));
            }
            expr_free(L); expr_free(C); expr_free(G);
        }
    }
    return NULL;
}

/* ===================================================================== */
/*  factor_monomial_content                                              */
/*                                                                       */
/*  Extracts the GCD of the variable-monomials of a Plus[term_1, ...]    */
/*  expression and returns a factorisation                               */
/*                                                                       */
/*      m * (P / m)                                                      */
/*                                                                       */
/*  where m is the largest monomial v_1^{e_1} * v_2^{e_2} * ... that     */
/*  divides every term of P.  Returns NULL if no nontrivial m exists or  */
/*  the input is not a Plus.                                             */
/*                                                                       */
/*  Examples (over Z[a, b]):                                             */
/*       3 a^2 b - 3 b - b^3      ->  b * (3 a^2 - 3 - b^2)              */
/*       x^2 y + x y^2            ->  x y * (x + y)                      */
/*       x + 1                    ->  NULL                               */
/*       Sin[x]^3 + Sin[x]        ->  Sin[x] * (1 + Sin[x]^2)            */
/*                                                                       */
/*  This is the cheapest factorisation step: O(n_terms * n_vars) with    */
/*  no integer arithmetic on coefficients.  It catches a large class of  */
/*  cases (every term sharing a variable factor) that would otherwise    */
/*  fall through to expensive trial-root or Hensel-based methods.        */
/*                                                                       */
/*  Note: integer-content extraction is handled separately by            */
/*  poly_content / FactorTerms; this function is concerned only with     */
/*  the variable powers.                                                 */
/*                                                                       */
/*  Memory: returns a fresh tree on success; the caller owns it.         */
/*  Pointers stored in the helper temporaries are *borrowed* from `P`    */
/*  and are not freed by us (extract_monomial returns borrowed slots).   */
/* ===================================================================== */
/* Recursive monomial decomposer.  Walks `e` collecting variable factors
 * (symbols, function-headed atoms like Sin[x], and Power[base, intExp]
 * nodes), accumulating an integer coefficient.  Robust to nested Times
 * (which can occur in picocas's canonical form, e.g. the parser builds
 *   3 a^2 b   as   Times[3, Times[Power[a,2], b]]
 * rather than the flat Times[3, Power[a,2], b]).  Anything we cannot
 * classify as "integer scalar" or "variable atom raised to integer power"
 * is recorded with exponent 1, treating it as an opaque atom so its
 * presence still constrains the monomial-content intersection.
 *
 * Output:
 *   coeff_out   GMP mpz_t multiplied (caller seeds with 1) by every
 *               integer or bigint factor encountered.  Caller is
 *               responsible for mpz_init/mpz_clear.
 *   atoms[]/exps[] receive {atom, exponent} pairs (atom pointers are
 *               *borrowed* from `e` -- the caller must not free them).
 *               Repeats are coalesced by atom-equality at the end.
 *   Returns true on success, false on overflow of the caller's buffer.
 *
 * Buffer convention: the caller pre-allocates `*cap` slots; we grow via
 * realloc when needed. `*count` is updated to reflect filled slots. */
static bool monomial_collect(Expr* e, mpz_t coeff_out,
                             Expr*** atoms, int64_t** exps,
                             size_t* count, size_t* cap) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER) {
        mpz_mul_si(coeff_out, coeff_out, (long)e->data.integer);
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        mpz_mul(coeff_out, coeff_out, e->data.bigint);
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!monomial_collect(e->data.function.args[i], coeff_out,
                                  atoms, exps, count, cap)) return false;
        }
        return true;
    }

    /* Identify (atom, exponent) for this node. */
    Expr* atom;
    int64_t exp;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Power
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_INTEGER) {
        atom = e->data.function.args[0];
        exp = e->data.function.args[1]->data.integer;
    } else {
        atom = e;
        exp = 1;
    }

    /* Coalesce with an existing atom or append. */
    for (size_t i = 0; i < *count; i++) {
        if (expr_eq((*atoms)[i], atom)) {
            (*exps)[i] += exp;
            return true;
        }
    }
    if (*count >= *cap) {
        size_t new_cap = (*cap == 0) ? 4 : (*cap * 2);
        *atoms = realloc(*atoms, sizeof(Expr*) * new_cap);
        *exps  = realloc(*exps,  sizeof(int64_t) * new_cap);
        if (!*atoms || !*exps) return false;
        *cap = new_cap;
    }
    (*atoms)[*count] = atom;
    (*exps)[*count] = exp;
    (*count)++;
    return true;
}

/* Build an INTEGER (when fits in int64) or BIGINT Expr from an mpz_t.
 * Helper for emitting the coefficient back into a residue term. */
static Expr* expr_from_mpz_normalized(const mpz_t v) {
    if (mpz_fits_slong_p(v)) {
        return expr_new_integer((int64_t)mpz_get_si(v));
    }
    return expr_new_bigint_from_mpz(v);
}

static Expr* factor_monomial_content(Expr* P) {
    if (P->type != EXPR_FUNCTION) return NULL;
    if (P->data.function.head->type != EXPR_SYMBOL) return NULL;
    if (P->data.function.head->data.symbol != SYM_Plus) return NULL;
    size_t n = P->data.function.arg_count;
    if (n < 2) return NULL;

    /* Pass 1: build the GCD-monomial.  We track {var -> min exponent}
     * across all terms.  A variable absent from any single term ends up
     * with exponent 0 and is dropped. */
    Expr** common_vars = NULL;       /* borrowed pointers into P's term subtrees */
    int64_t* common_exps = NULL;
    size_t   common_count = 0;
    size_t   common_cap = 0;
    bool     first = true;

    for (size_t i = 0; i < n; i++) {
        Expr* term = P->data.function.args[i];
        mpz_t coeff;
        mpz_init_set_si(coeff, 1);
        Expr** tvars = NULL;
        int64_t* texps = NULL;
        size_t tcount = 0, tcap = 0;
        if (!monomial_collect(term, coeff, &tvars, &texps, &tcount, &tcap)) {
            mpz_clear(coeff);
            free(tvars); free(texps);
            free(common_vars); free(common_exps);
            return NULL;
        }
        mpz_clear(coeff);

        if (first) {
            /* Seed: every variable from term 0 is a candidate. */
            common_cap = tcount;
            common_vars = malloc(sizeof(Expr*) * (common_cap + 1));
            common_exps = malloc(sizeof(int64_t) * (common_cap + 1));
            for (size_t j = 0; j < tcount; j++) {
                if (texps[j] <= 0) continue;
                common_vars[common_count] = tvars[j];
                common_exps[common_count] = texps[j];
                common_count++;
            }
            first = false;
        } else {
            /* Intersect with current term's variables. */
            size_t out = 0;
            for (size_t j = 0; j < common_count; j++) {
                int64_t found_exp = -1;
                for (size_t k = 0; k < tcount; k++) {
                    if (expr_eq(common_vars[j], tvars[k])) {
                        found_exp = texps[k];
                        break;
                    }
                }
                if (found_exp < 1) continue; /* var missing from this term */
                int64_t mn = (common_exps[j] < found_exp) ? common_exps[j] : found_exp;
                if (mn <= 0) continue;
                common_vars[out] = common_vars[j];
                common_exps[out] = mn;
                out++;
            }
            common_count = out;
        }

        free(tvars); free(texps);

        if (common_count == 0) {
            free(common_vars); free(common_exps);
            return NULL;
        }
    }

    if (common_count == 0) {
        free(common_vars); free(common_exps);
        return NULL;
    }

    /* Build the monomial factor m. */
    Expr* m;
    if (common_count == 1) {
        if (common_exps[0] == 1) {
            m = expr_copy(common_vars[0]);
        } else {
            m = eval_and_free(expr_new_function(
                expr_new_symbol("Power"),
                (Expr*[]){expr_copy(common_vars[0]),
                         expr_new_integer(common_exps[0])}, 2));
        }
    } else {
        Expr** m_args = malloc(sizeof(Expr*) * common_count);
        for (size_t i = 0; i < common_count; i++) {
            if (common_exps[i] == 1) {
                m_args[i] = expr_copy(common_vars[i]);
            } else {
                m_args[i] = eval_and_free(expr_new_function(
                    expr_new_symbol("Power"),
                    (Expr*[]){expr_copy(common_vars[i]),
                             expr_new_integer(common_exps[i])}, 2));
            }
        }
        m = eval_and_free(expr_new_function(
            expr_new_symbol("Times"), m_args, common_count));
        free(m_args);
    }

    /* Pass 2: build the residue P/m by reducing every term's exponents. */
    Expr** res_terms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        Expr* term = P->data.function.args[i];
        mpz_t coeff;
        mpz_init_set_si(coeff, 1);
        Expr** tvars = NULL;
        int64_t* texps = NULL;
        size_t tcount = 0, tcap = 0;
        if (!monomial_collect(term, coeff, &tvars, &texps, &tcount, &tcap)) {
            /* Pass 1 succeeded so this should not fail; defensive only. */
            mpz_clear(coeff);
            free(tvars); free(texps);
            for (size_t k = 0; k < i; k++) expr_free(res_terms[k]);
            free(res_terms);
            free(common_vars); free(common_exps);
            expr_free(m);
            return NULL;
        }

        size_t cap = tcount + 1;
        Expr** factor_args = malloc(sizeof(Expr*) * cap);
        size_t fc = 0;

        if (mpz_cmp_si(coeff, 1) != 0) {
            factor_args[fc++] = expr_from_mpz_normalized(coeff);
        }
        mpz_clear(coeff);
        for (size_t j = 0; j < tcount; j++) {
            int64_t reduced = texps[j];
            for (size_t k = 0; k < common_count; k++) {
                if (expr_eq(tvars[j], common_vars[k])) {
                    reduced -= common_exps[k];
                    break;
                }
            }
            if (reduced <= 0) continue;
            if (reduced == 1) {
                factor_args[fc++] = expr_copy(tvars[j]);
            } else {
                factor_args[fc++] = eval_and_free(expr_new_function(
                    expr_new_symbol("Power"),
                    (Expr*[]){expr_copy(tvars[j]),
                             expr_new_integer(reduced)}, 2));
            }
        }

        Expr* term_residue;
        if (fc == 0) {
            term_residue = expr_new_integer(1);
        } else if (fc == 1) {
            term_residue = factor_args[0];
        } else {
            term_residue = eval_and_free(expr_new_function(
                expr_new_symbol("Times"), factor_args, fc));
        }
        free(factor_args);
        free(tvars); free(texps);

        res_terms[i] = term_residue;
    }

    Expr* residue = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"), res_terms, n));
    free(res_terms);
    free(common_vars); free(common_exps);

    /* Recurse: the residue may admit further factorisation by the other
     * heuristic strategies (degree-1, binomial, root-trial, etc.). */
    Expr* f_residue = heuristic_factor(residue);
    expr_free(residue);

    return eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){m, f_residue}, 2));
}

/* ===================================================================== */
/*  is_likely_irreducible_multivariate                                   */
/*                                                                       */
/*  A probabilistic but extremely conservative irreducibility test for   */
/*  multivariate polynomials over Z, used to short-circuit the expensive */
/*  trial-root pipeline (factor_roots) when no factor can possibly be    */
/*  found.                                                               */
/*                                                                       */
/*  Theory (Hilbert irreducibility, in concrete form):                   */
/*    If P(x_1, ..., x_n) ∈ Z[x_1, ..., x_n] is irreducible and we       */
/*    substitute x_i = α_i for all i ≠ j with random integer α_i, then   */
/*    for "almost all" choices the resulting univariate polynomial in    */
/*    x_j is also irreducible.  Conversely, if P factors then so does    */
/*    the image (because the substitution homomorphism Z[x_1,...,x_n]    */
/*    → Z[x_j] preserves products).                                      */
/*                                                                       */
/*    So: if at SEVERAL good evaluation points the univariate image is   */
/*    squarefree and irreducible (treating "good" = leading coefficient  */
/*    in x_j survives + image is squarefree), we can conclude P is       */
/*    almost certainly irreducible.                                      */
/*                                                                       */
/*  "Almost" is the operative word: for very specific algebraic forms    */
/*  (e.g. polynomials whose factorisation requires a special evaluation  */
/*  pattern), the test could be fooled.  We compensate by:               */
/*    - trying every variable as the "main" variable in turn             */
/*    - requiring at least IRRED_CONFIRM_COUNT confirmations             */
/*  before declaring irreducibility.  False negatives are harmless --    */
/*  we just fall through to the more expensive pipeline.  False          */
/*  positives are the danger; the multi-variable, multi-point structure  */
/*  protects against them.                                               */
/*                                                                       */
/*  Returns true if the test confirms irreducibility.  Returns false if  */
/*  reducibility is detected OR if the test cannot reach a conclusion    */
/*  (insufficient confirmations).                                        */
/*                                                                       */
/*  Cost: O(IRRED_TRY_POINTS * v_count) calls to bz_factor_to_expr,      */
/*  each operating on a univariate polynomial of degree deg(P).  Each    */
/*  factor call is fast because picocas's univariate Berlekamp-          */
/*  Zassenhaus runs in milliseconds for typical inputs.  In total this   */
/*  is far cheaper than factor_roots' O(v_count^2 * 10) trial-divisions  */
/*  on irreducible inputs (where every division is wasted work).         */
/* ===================================================================== */

#define IRRED_TRY_POINTS    7   /* alpha values to try, in order: 1, -1, 2, -2, 3, -3, 4 */
#define IRRED_CONFIRM_COUNT 2   /* number of irreducible-image confirmations required */

/* Substitute every var in `vars` except `vars[main_idx]` with the integer
 * `alpha`, returning a freshly evaluated univariate polynomial in
 * `vars[main_idx]` (or NULL on allocation failure). */
static Expr* eval_others_at_alpha(Expr* P, Expr** vars, size_t v_count,
                                  size_t main_idx, int64_t alpha) {
    if (v_count <= 1) return expr_copy(P);

    /* Build a list of replacement rules.  We use ReplaceAll so that all
     * occurrences anywhere in P are substituted simultaneously, and rely
     * on the evaluator to reduce the resulting integer arithmetic. */
    size_t n_rules = v_count - 1;
    Expr** rules = malloc(sizeof(Expr*) * n_rules);
    if (!rules) return NULL;
    size_t idx = 0;
    for (size_t i = 0; i < v_count; i++) {
        if (i == main_idx) continue;
        Expr* rule_args[2] = {
            expr_copy(vars[i]),
            expr_new_integer(alpha)
        };
        rules[idx++] = expr_new_function(expr_new_symbol("Rule"),
                                         rule_args, 2);
    }
    Expr* rules_list = expr_new_function(expr_new_symbol("List"),
                                         rules, n_rules);
    free(rules);
    Expr* call_args[2] = { expr_copy(P), rules_list };
    Expr* call = expr_new_function(expr_new_symbol("ReplaceAll"),
                                   call_args, 2);
    Expr* image = evaluate(call);
    expr_free(call);
    /* Force expansion -- ReplaceAll alone does not normalise the result. */
    Expr* expanded = expr_expand(image);
    expr_free(image);
    return expanded;
}

/* Count the number of nontrivial (degree ≥ 1 in `var`) factors of an
 * already-factored expression `factored`.  A "Times" head means a
 * multiplication of factors; a single non-Times node counts as one
 * factor.  Constant factors (Integer / Rational / no var) are ignored.
 * Power nodes are counted as one factor (their multiplicity is reflected
 * in the exponent, not in distinct linear pieces). */
static int count_nontrivial_factors(Expr* factored, Expr* var) {
    if (!factored) return 0;
    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol == SYM_Times) {
        int total = 0;
        for (size_t i = 0; i < factored->data.function.arg_count; i++) {
            Expr* a = factored->data.function.args[i];
            if (get_degree_poly(a, var) >= 1) total++;
        }
        return total;
    }
    return (get_degree_poly(factored, var) >= 1) ? 1 : 0;
}

/* Return true if the univariate polynomial `u` in variable `var` is
 * squarefree, i.e., gcd(u, u') is a unit.
 *
 * Implementation: convert to ZUPoly and use the subresultant-PRS
 * `zupoly_gcd`.  The Expr-level `poly_gcd_internal` uses Knuth-style
 * primitive PRS which suffers exponential coefficient growth on the
 * intermediate pseudo-remainders -- on a degree-31 univariate input
 * (which arises from sqfree_cheap_check on trivariate polynomials)
 * it can take >55s, dominating the entire FactorSquareFree call.
 * `zupoly_gcd` does the same job in sub-millisecond time on these
 * inputs, since subresultant PRS keeps coefficient sizes polynomially
 * bounded. */
static bool univariate_squarefree(Expr* u, Expr* var) {
    ZUPoly* zu = expr_to_zupoly(u, var);
    if (zu) {
        /* Build derivative directly in ZUPoly to avoid Expr-level
         * derivation + ReplaceAll round-tripping. */
        ZUPoly* zdu = zupoly_zero();
        for (int i = 1; i <= zu->deg; i++) {
            const mpz_t* ci = zupoly_getcoef(zu, i);
            if (!ci) continue;
            mpz_t v;
            mpz_init(v);
            mpz_mul_si(v, *ci, i);
            zupoly_setcoef(zdu, i - 1, v);
            mpz_clear(v);
        }
        ZUPoly* g = zupoly_gcd(zu, zdu);
        bool sqf = (g->deg <= 0);
        zupoly_free(zu);
        zupoly_free(zdu);
        zupoly_free(g);
        return sqf;
    }

    /* Fallback for inputs that don't convert to ZUPoly (e.g., rational
     * coefficients that slipped through).  Original Expr-level path. */
    Expr* du_raw = poly_deriv(u, var);
    Expr* du = expr_expand(du_raw);
    expr_free(du_raw);
    Expr* vars1[1] = { var };
    Expr* g = poly_gcd_internal(u, du, vars1, 1);
    expr_free(du);
    bool sqf = (g->type == EXPR_INTEGER)
            || (get_degree_poly(g, var) == 0);
    expr_free(g);
    return sqf;
}

/* Core of the irreducibility test (see header comment for theory). */
static bool is_likely_irreducible_multivariate(Expr* P, Expr** vars,
                                               size_t v_count) {
    if (v_count < 2) return false;  /* univariate path is handled elsewhere */

    /* Evaluation points to try, in order of preference (small magnitudes
     * first). 0 is excluded because it commonly kills leading
     * coefficients of multivariate factors. */
    static const int64_t alphas[IRRED_TRY_POINTS] = {1, -1, 2, -2, 3, -3, 4};

    for (size_t mi = 0; mi < v_count; mi++) {
        Expr* main_var = vars[mi];
        int deg_main = get_degree_poly(P, main_var);
        if (deg_main < 2) continue;  /* degree-1 main: trivially irreducible/handled by factor_degree_one */

        int confirmations = 0;
        for (size_t pi = 0; pi < IRRED_TRY_POINTS; pi++) {
            int64_t a = alphas[pi];
            Expr* image = eval_others_at_alpha(P, vars, v_count, mi, a);
            if (!image) continue;

            /* Reject points that drop the main-variable degree --
             * such points killed the leading coefficient and the
             * image's irreducibility tells us nothing. */
            if (get_degree_poly(image, main_var) != deg_main) {
                expr_free(image);
                continue;
            }

            /* Reject points where the image is not squarefree --
             * irreducibility cannot be inferred from a non-squarefree
             * image (it may have repeated roots that hide structure). */
            if (!univariate_squarefree(image, main_var)) {
                expr_free(image);
                continue;
            }

            /* Skip alphas whose image has a BIGINT coefficient.
             * `bz_factor_to_expr` cannot ingest such coefficients (its
             * int64 UPoly substrate would silently return the input
             * unchanged), and `count_nontrivial_factors` would then
             * miscount the bare Plus as a single irreducible factor --
             * a false "irreducible" confirmation.  See FACTOR_PLAN.md
             * §12.F6 for the trace. */
            bool image_has_bigint = false;
            for (int i = 0; i <= deg_main; i++) {
                Expr* c = get_coeff(image, main_var, i);
                if (c->type == EXPR_BIGINT) image_has_bigint = true;
                expr_free(c);
                if (image_has_bigint) break;
            }
            if (image_has_bigint) {
                expr_free(image);
                continue;
            }

            /* Factor the univariate image. */
            Expr* factored = bz_factor_to_expr(image, main_var);
            int n = count_nontrivial_factors(factored, main_var);
            expr_free(factored);
            expr_free(image);

            if (n == 0) {
                /* Image had no main-variable content (shouldn't happen
                 * after the degree check above, but defensive). */
                continue;
            }
            if (n >= 2) {
                /* Image factors at this alpha.  This is *not* proof
                 * that P factors -- accidental factorisations of the
                 * univariate image are common (e.g., P = x^3 - 2y^3 - 1
                 * gives x^3 + 1 = (x+1)(x^2-x+1) at y = -1, even though
                 * P is irreducible over Z[x,y]).  We simply skip this
                 * alpha and continue looking for irreducible-image
                 * confirmations at other points. */
                continue;
            }
            /* n == 1: image is irreducible at this alpha -> a
             * confirmation that P is likely irreducible.  Provided we
             * accumulate enough such confirmations *for the same
             * main variable* (not mixing across variables, which
             * could hide a constant-in-one-variable factor), we
             * conclude irreducibility. */
            confirmations++;
            if (confirmations >= IRRED_CONFIRM_COUNT) return true;
        }
        /* This main variable did not yield enough irreducible-image
         * confirmations.  Try the next variable. */
    }

    return false;
}

static Expr* factor_roots(Expr* P, Expr** vars, size_t v_count) {
    int64_t c_vals[] = {1, -1, 2, -2, 3, -3, 4, -4, 6, -6};
    size_t num_c = 10;
    
    for (size_t i = 0; i < v_count; i++) {
        Expr* x = vars[i];
        int d = get_degree_poly(P, x);
        if (d < 2) continue; 
        
        for (size_t j = 0; j < num_c; j++) {
            int64_t c = c_vals[j];
            Expr* R = expr_new_integer(c);
            Expr* x_minus_R = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(x), expr_new_integer(-c)}, 2));
            Expr* Q = exact_poly_div(P, x_minus_R, vars, v_count);
            
            if (Q) {
                Expr* f_Q = heuristic_factor(Q);
                expr_free(Q);
                expr_free(R);
                return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){x_minus_R, f_Q}, 2));
            }
            expr_free(x_minus_R);
            expr_free(R);
            
            for (size_t k = 0; k < v_count; k++) {
                if (i == k) continue;
                Expr* y = vars[k];
                Expr* term = (c == 1) ? expr_copy(y) : ((c == -1) ? eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(y)}, 2)) : eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(c), expr_copy(y)}, 2)));
                
                Expr* neg_term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(term)}, 2));
                Expr* x_minus_y = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){expr_copy(x), neg_term}, 2));
                
                Expr* Q2 = exact_poly_div(P, x_minus_y, vars, v_count);
                if (Q2) {
                    Expr* f_Q2 = heuristic_factor(Q2);
                    expr_free(Q2);
                    expr_free(term);
                    return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){x_minus_y, f_Q2}, 2));
                }
                expr_free(x_minus_y);
                expr_free(term);
            }
        }
    }
    return NULL;
}

/* ===================================================================== */
/*  Bivariate-Hensel wiring                                              */
/*                                                                       */
/*  Glue between picocas's existing univariate Berlekamp-Zassenhaus     */
/*  (bz_factor_to_expr) and the new mvfactor multivariate pipeline.      */
/*  This is the integration point that lets `heuristic_factor` route    */
/*  monic bivariate inputs through the Hensel lift before falling       */
/*  back to factor_roots.                                                */
/* ===================================================================== */

/* Callback adapter for mvfactor_try_bivariate_monic.  Receives a
 * univariate ZUPoly image and must produce its irreducible-over-Z
 * factors as a fresh ZUPoly array.  Implementation: convert ZUPoly
 * -> Expr -> bz_factor_to_expr -> walk the resulting Times -> back
 * to ZUPoly. */
static bool factor_via_bz_callback(const ZUPoly* image,
                                   ZUPoly*** factors_out,
                                   int* count_out,
                                   void* user_data) {
    (void)user_data;

    /* Bail when any image coefficient overflows int64.  The int64-based
     * UPoly substrate inside `bz_factor_to_expr` cannot ingest such
     * inputs and would silently return them unchanged, which the
     * orchestrator above misreads as "image is irreducible" (a single
     * Plus factor).  Returning false instead lets the orchestrator try
     * the next alpha. See FACTOR_PLAN.md §12.F6 for the Fateman trace. */
    for (int i = 0; i <= image->deg; i++) {
        const mpz_t* ci = zupoly_getcoef(image, i);
        if (ci && !mpz_fits_slong_p(*ci)) return false;
    }

    /* Use a private symbol name to avoid colliding with anything the
     * user may have defined.  The symbol is local to this conversion
     * and never leaks out. */
    Expr* var = expr_new_symbol("$mvfactor_x");
    Expr* image_expr = zupoly_to_expr(image, var);
    Expr* factored = bz_factor_to_expr(image_expr, var);
    expr_free(image_expr);

    /* Walk the result.  bz_factor_to_expr returns either a Times of
     * factors (with possible trailing constant) or a single factor /
     * the input itself for the irreducible case. */
    ZUPoly** result = NULL;
    int n = 0, cap = 0;

    if (factored->type == EXPR_FUNCTION
        && factored->data.function.head
        && factored->data.function.head->type == EXPR_SYMBOL
        && factored->data.function.head->data.symbol == SYM_Times) {
        size_t ac = factored->data.function.arg_count;
        cap = (int)ac;
        result = (ZUPoly**)malloc(sizeof(ZUPoly*) * (size_t)cap);
        for (size_t i = 0; i < ac; i++) {
            Expr* a = factored->data.function.args[i];
            ZUPoly* zp = expr_to_zupoly(a, var);
            /* Drop pure-constant factors -- mvfactor only cares about
             * the variable polynomial pieces; the integer content is
             * already handled upstream. */
            if (zp && !zupoly_is_zero(zp) && zp->deg >= 1) {
                result[n++] = zp;
            } else if (zp) {
                zupoly_free(zp);
            }
        }
    } else {
        ZUPoly* zp = expr_to_zupoly(factored, var);
        if (zp && !zupoly_is_zero(zp) && zp->deg >= 1) {
            cap = 1;
            result = (ZUPoly**)malloc(sizeof(ZUPoly*));
            result[0] = zp;
            n = 1;
        } else if (zp) {
            zupoly_free(zp);
        }
    }
    expr_free(factored);
    expr_free(var);

    if (n == 0) {
        free(result);
        return false;
    }
    *factors_out = result;
    *count_out = n;
    return true;
}

/* Bivariate-Hensel leading-coefficient correction strategy. */
typedef enum {
    BV_LC_MONIC  = 0,  /* lc_x(P) = +1: lift P directly. */
    BV_LC_NEGATE = 1,  /* lc_x(P) = -1: Stage 1 — negate P, lift, flip one factor. */
    BV_LC_SCALE  = 2,  /* lc_x(P) = a, |a| > 1, integer constant: Stage 2 —
                        * Wang's monic substitution Q = a^(d-1) · P(x/a, y). */
} BvLcKind;

/* Pick the best variable index for the bivariate Hensel pipeline.
 *
 * The chosen variable's leading coefficient (in `*P` as a polynomial in
 * that variable) must be a non-zero integer constant -- i.e. constant
 * with respect to the other variable.  Among candidates:
 *
 *   1. lc = +1   -> BV_LC_MONIC   (cheapest path, no correction).
 *   2. lc = -1   -> BV_LC_NEGATE  (Stage 1 negate path).
 *   3. |lc| > 1  -> BV_LC_SCALE   (Stage 2 Wang via x-scale substitution).
 *
 * Among Stage-2 candidates, the variable with smaller |lc| wins (less
 * coefficient blowup in the substitution).
 *
 * On success returns the variable index and sets *kind_out / *lc_out.
 * `*lc_out` is mpz_init'd by the caller (we mpz_set into it).  Returns
 * -1 if no variable has a constant integer leading coefficient. */
static int pick_factorable_x_var(Expr* P, Expr** vars, size_t v_count,
                                 BvLcKind* kind_out, mpz_t* lc_out,
                                 int excluded_idx) {
    int best_idx = -1;
    int best_priority = INT_MAX;
    mpz_t best_abs;     mpz_init_set_ui(best_abs, 0);
    mpz_t best_lc;      mpz_init(best_lc);

    for (size_t i = 0; i < v_count; i++) {
        if ((int)i == excluded_idx) continue;
        int deg = get_degree_poly(P, vars[i]);
        if (deg < 1) continue;
        Expr* lc = get_coeff(P, vars[i], deg);
        bool is_int_const = expr_is_integer_like(lc) &&
            !(lc->type == EXPR_INTEGER && lc->data.integer == 0);
        if (!is_int_const) {
            expr_free(lc);
            continue;
        }
        mpz_t lc_mpz; mpz_init(lc_mpz);
        expr_to_mpz(lc, lc_mpz);
        expr_free(lc);

        int priority;
        if (mpz_cmp_si(lc_mpz, 1) == 0)         priority = 0;
        else if (mpz_cmp_si(lc_mpz, -1) == 0)   priority = 1;
        else                                    priority = 2;

        bool better = false;
        if (priority < best_priority) {
            better = true;
        } else if (priority == best_priority && priority == 2) {
            /* Tiebreak among Stage-2 candidates: prefer smaller |lc|. */
            mpz_t abs_lc; mpz_init(abs_lc);
            mpz_abs(abs_lc, lc_mpz);
            if (mpz_cmp(abs_lc, best_abs) < 0) better = true;
            mpz_clear(abs_lc);
        }

        if (better) {
            best_idx = (int)i;
            best_priority = priority;
            mpz_abs(best_abs, lc_mpz);
            mpz_set(best_lc, lc_mpz);
        }
        mpz_clear(lc_mpz);

        if (priority == 0) break;  /* +1 is optimal; stop searching. */
    }

    if (best_idx >= 0) {
        if (kind_out) *kind_out = (BvLcKind)best_priority;
        if (lc_out)   mpz_set(*lc_out, best_lc);
    }
    mpz_clear(best_abs);
    mpz_clear(best_lc);
    return best_idx;
}

/* ---------------------------------------------------------------------- */
/*  Phase F1 Stage 2 helpers (Wang's monic-substitution recipe)            */
/* ---------------------------------------------------------------------- */

/* Compute the integer content of a BPoly: the gcd of all its integer
 * coefficients across both x and y.  `out` must already be mpz_init'd;
 * set to 0 for the zero polynomial. */
static void bpoly_int_content(const BPoly* P, mpz_t out) {
    mpz_set_ui(out, 0);
    if (P->deg_x < 0) return;
    for (int i = 0; i <= P->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(P, i);
        if (!yi || zupoly_is_zero(yi)) continue;
        mpz_t cy; mpz_init(cy);
        zupoly_content(yi, cy);
        if (mpz_sgn(out) == 0) mpz_set(out, cy);
        else                   mpz_gcd(out, out, cy);
        mpz_clear(cy);
    }
}

/* Divide every coefficient of P by integer c.  Returns a fresh BPoly,
 * or NULL if any coefficient is not divisible by c (which would be a
 * programming error in our use-cases). */
static BPoly* bpoly_div_int_exact(const BPoly* P, const mpz_t c) {
    if (mpz_sgn(c) == 0) return NULL;
    if (P->deg_x < 0) return bpoly_zero();
    BPoly* R = bpoly_new(P->deg_x + 1);
    for (int i = 0; i <= P->deg_x; i++) {
        const ZUPoly* yi = bpoly_get_xcoef(P, i);
        if (!yi || zupoly_is_zero(yi)) continue;
        ZUPoly* d = zupoly_new(yi->deg + 1);
        for (int j = 0; j <= yi->deg; j++) {
            const mpz_t* cj = zupoly_getcoef(yi, j);
            if (!cj || mpz_sgn(*cj) == 0) continue;
            if (!mpz_divisible_p(*cj, c)) {
                zupoly_free(d);
                bpoly_free(R);
                return NULL;
            }
            mpz_t q; mpz_init(q);
            mpz_divexact(q, *cj, c);
            zupoly_setcoef(d, j, q);
            mpz_clear(q);
        }
        bpoly_set_xcoef(R, i, d);
    }
    return R;
}

/* Form Q = a^(d-1) · P(x/a, y), assuming lc_x(P) is the constant
 * integer `a` (a non-zero, non-unit) and d = deg_x(P) >= 1.  Q is
 * monic in x and integer-coefficient.  Concretely:
 *   cx_Q[d] = 1
 *   cx_Q[j] = cx_P[j] · a^(d-1-j)  for 0 <= j < d.
 *
 * Caller owns the returned BPoly.  Returns NULL on degenerate input. */
static BPoly* bpoly_make_monic_via_x_scale(const BPoly* P, const mpz_t lc_a) {
    int d = P->deg_x;
    if (d < 1 || mpz_sgn(lc_a) == 0) return NULL;

    BPoly* Q = bpoly_new(d + 1);

    /* cx_Q[d] = 1. */
    ZUPoly* one = zupoly_from_int(1);
    bpoly_set_xcoef(Q, d, one);

    /* a_pow tracks a^(d-1-j) as j decreases from d-1 to 0:
     * j = d-1 -> a^0 = 1; j = d-2 -> a^1; ...; j = 0 -> a^(d-1). */
    mpz_t a_pow; mpz_init_set_ui(a_pow, 1);
    for (int j = d - 1; j >= 0; j--) {
        const ZUPoly* yj = bpoly_get_xcoef(P, j);
        if (yj && !zupoly_is_zero(yj)) {
            ZUPoly* scaled = zupoly_scale(yj, a_pow);
            bpoly_set_xcoef(Q, j, scaled);
        }
        if (j > 0) mpz_mul(a_pow, a_pow, lc_a);
    }
    mpz_clear(a_pow);
    return Q;
}

/* Form H = G(c · x, y).  Each x^j coefficient (a y-poly) is multiplied
 * by c^j.  Caller owns the result. */
static BPoly* bpoly_subst_x_scale(const BPoly* G, const mpz_t c) {
    int d = G->deg_x;
    if (d < 0) return bpoly_zero();
    BPoly* H = bpoly_new(d + 1);
    mpz_t c_pow; mpz_init_set_ui(c_pow, 1);
    for (int j = 0; j <= d; j++) {
        const ZUPoly* gj = bpoly_get_xcoef(G, j);
        if (gj && !zupoly_is_zero(gj)) {
            ZUPoly* scaled = zupoly_scale(gj, c_pow);
            bpoly_set_xcoef(H, j, scaled);
        }
        if (j < d) mpz_mul(c_pow, c_pow, c);
    }
    mpz_clear(c_pow);
    return H;
}

/* ====================================================================== */
/*  Phase F1 Stage 3: polynomial-in-y leading coefficient                  */
/*                                                                        */
/*  Wang's leading-coefficient correction for inputs whose lc_x(P)(y) is  */
/*  a non-constant polynomial in the other variable.                       */
/*                                                                        */
/*  MVP scope:                                                            */
/*   - Bivariate inputs (handled at the F1 entry; trivariate goes to F2). */
/*   - Two univariate factors (r = 2) of the squarefree image.            */
/*   - Both univariate factors must be monic.                             */
/*   - cont(A) ∈ {±1}; A's polynomial factors are tracked with sign.      */
/*   - At least one evaluation point α with A(α) = +1 must exist.         */
/*                                                                        */
/*  Algorithm (per (variable, α) pair):                                   */
/*   1. Convert P to a BPoly bp and compute A = lc_x(P)(y) as a ZUPoly.   */
/*   2. Factor A over Z[y] to extract its content c and primitive         */
/*      irreducible factors a_factors[] (with multiplicity).              */
/*   3. Iterate over candidate α values; require A(α) = +1.               */
/*   4. Form the squarefree univariate image P(x, α) via bpoly_eval_y_si  */
/*      and run bz_factor_to_expr.  Verify r = 2 and both factors monic.  */
/*   5. Shift to y' = y - α coords (the lift assumes α = 0).              */
/*   6. Enumerate the 2^n distributions of A's polynomial factors         */
/*      between predicted leading coefficients q_u, q_v with q_u·q_v = A. */
/*      Filter to distributions with q_u(α) = q_v(α) = +1.                */
/*   7. For each surviving distribution, run bpoly_hensel_lift_2_lc.      */
/*      First success wins.                                               */
/*   8. Unshift the resulting BPoly factors and return Times[F_1, F_2].   */
/* ====================================================================== */

/* Factor a ZUPoly over Z[y] using the existing bz_factor_to_expr path.
 * Returns:
 *   - out_content: signed integer leading content c (caller mpz_init's).
 *   - out_factors: array of primitive irreducible polynomial factors with
 *     positive leading coefficient.  Repeated factors appear with full
 *     multiplicity (a^2 yields two entries).
 *   - out_count: number of factor entries.
 *
 * p = c · Π out_factors[i].  Caller frees the array and each ZUPoly. */
static bool zupoly_factor_to_array(const ZUPoly* p,
                                    mpz_t out_content,
                                    ZUPoly*** out_factors,
                                    int* out_count) {
    if (!p || zupoly_is_zero(p)) return false;

    Expr* var = expr_new_symbol("$mvfactor_lc");
    Expr* p_expr = zupoly_to_expr(p, var);
    Expr* fac = bz_factor_to_expr(p_expr, var);
    expr_free(p_expr);

    mpz_set_ui(out_content, 1);
    ZUPoly** factors = NULL;
    int n = 0, cap = 0;
    bool ok = true;

    /* Walk fac as either Times[args...] or a single arg. */
    size_t arg_count;
    Expr* singleton_arr[1];
    Expr** args;

    if (fac->type == EXPR_FUNCTION
        && fac->data.function.head
        && fac->data.function.head->type == EXPR_SYMBOL
        && fac->data.function.head->data.symbol == SYM_Times) {
        arg_count = fac->data.function.arg_count;
        args = fac->data.function.args;
    } else {
        singleton_arr[0] = fac;
        arg_count = 1;
        args = singleton_arr;
    }

    for (size_t i = 0; i < arg_count && ok; i++) {
        Expr* a = args[i];
        Expr* base = a;
        int e = 1;

        if (a->type == EXPR_INTEGER || a->type == EXPR_BIGINT) {
            mpz_t k; mpz_init(k);
            expr_to_mpz(a, k);
            mpz_mul(out_content, out_content, k);
            mpz_clear(k);
            continue;
        }
        if (a->type == EXPR_FUNCTION
            && a->data.function.head
            && a->data.function.head->type == EXPR_SYMBOL
            && a->data.function.head->data.symbol == SYM_Power
            && a->data.function.arg_count == 2) {
            Expr* exp_e = a->data.function.args[1];
            if (exp_e->type != EXPR_INTEGER || exp_e->data.integer < 1) {
                ok = false;
                break;
            }
            base = a->data.function.args[0];
            e = (int)exp_e->data.integer;
        }

        ZUPoly* zp = expr_to_zupoly(base, var);
        if (!zp || zp->deg < 1) {
            if (zp) zupoly_free(zp);
            continue;
        }
        /* Normalise to positive leading coefficient. */
        if (mpz_sgn(zp->c[zp->deg]) < 0) {
            ZUPoly* neg = zupoly_neg(zp);
            zupoly_free(zp);
            zp = neg;
            if ((e & 1) != 0) mpz_neg(out_content, out_content);
        }
        for (int j = 0; j < e; j++) {
            if (n == cap) {
                cap = (cap == 0) ? 4 : cap * 2;
                factors = (ZUPoly**)realloc(factors, sizeof(ZUPoly*) * (size_t)cap);
            }
            factors[n++] = zupoly_copy(zp);
        }
        zupoly_free(zp);
    }

    expr_free(fac);
    expr_free(var);

    if (!ok) {
        for (int i = 0; i < n; i++) zupoly_free(factors[i]);
        free(factors);
        *out_factors = NULL;
        *out_count = 0;
        return false;
    }

    *out_factors = factors;
    *out_count = n;
    return true;
}

/* Build q(y) = sign · Π factors[i] over the bits of `mask`. */
static ZUPoly* polylc_build_q(ZUPoly* const* factors, int n,
                               uint64_t mask, int sign) {
    ZUPoly* q = zupoly_from_int(sign);
    for (int i = 0; i < n; i++) {
        if ((mask >> i) & 1u) {
            ZUPoly* nq = zupoly_mul(q, factors[i]);
            zupoly_free(q);
            q = nq;
        }
    }
    return q;
}

/* Try the F1 Stage 3 lift across all 2^n distributions of A's factors
 * for the given (u, v) seed pair.  Returns true on first lift success;
 * sets *U_out, *V_out (in shifted coord). */
static bool try_polylc_two_factor(const BPoly* P_shifted,
                                   const ZUPoly* u, const ZUPoly* v,
                                   const mpz_t a_content,
                                   ZUPoly* const* a_factors_shifted, int n,
                                   BPoly** U_out, BPoly** V_out) {
    if (n > 12) return false;            /* avoid 2^n blow-up */
    if (mpz_cmpabs_ui(a_content, 1) != 0) return false; /* MVP: |cont(A)| = 1 */
    int c_sign = (mpz_sgn(a_content) > 0) ? 1 : -1;

    uint64_t total = (uint64_t)1 << n;
    for (uint64_t mask = 0; mask < total; mask++) {
        /* prod_{u,v}(0) = product of g_i(0) over set/cleared bits. */
        mpz_t prod_u_at0; mpz_init_set_ui(prod_u_at0, 1);
        mpz_t prod_v_at0; mpz_init_set_ui(prod_v_at0, 1);
        for (int i = 0; i < n; i++) {
            const mpz_t* g_i_at0 = zupoly_getcoef(a_factors_shifted[i], 0);
            mpz_t v_i; mpz_init(v_i);
            if (g_i_at0) mpz_set(v_i, *g_i_at0);
            else         mpz_set_ui(v_i, 0);
            if ((mask >> i) & 1u) mpz_mul(prod_u_at0, prod_u_at0, v_i);
            else                   mpz_mul(prod_v_at0, prod_v_at0, v_i);
            mpz_clear(v_i);
        }

        /* For q_u(0) = +1 we need c_u · prod_u(0) = +1, with c_u ∈ {±1}.
         * That forces prod_u(0) ∈ {±1}.  Same for v. */
        bool match = (mpz_cmpabs_ui(prod_u_at0, 1) == 0
                      && mpz_cmpabs_ui(prod_v_at0, 1) == 0);
        int c_u_sign = 0, c_v_sign = 0;
        if (match) {
            c_u_sign = (mpz_sgn(prod_u_at0) > 0) ? 1 : -1;
            c_v_sign = (mpz_sgn(prod_v_at0) > 0) ? 1 : -1;
            /* And c_u · c_v = c (the integer content). */
            if (c_u_sign * c_v_sign != c_sign) match = false;
        }
        mpz_clear(prod_u_at0);
        mpz_clear(prod_v_at0);
        if (!match) continue;

        ZUPoly* q_u = polylc_build_q(a_factors_shifted, n, mask, c_u_sign);
        ZUPoly* q_v = polylc_build_q(a_factors_shifted, n,
                                       ~mask & (total - 1), c_v_sign);

        BPoly *U = NULL, *V = NULL;
        bool lift_ok = bpoly_hensel_lift_2_lc(P_shifted, u, v, q_u, q_v,
                                               &U, &V);
        zupoly_free(q_u);
        zupoly_free(q_v);
        if (lift_ok) {
            *U_out = U;
            *V_out = V;
            return true;
        }
    }
    return false;
}

static const int64_t POLYLC_ALPHA_TRIES[] = {0, 1, -1, 2, -2, 3, -3, 4, -4};
#define POLYLC_ALPHA_COUNT \
    (sizeof(POLYLC_ALPHA_TRIES) / sizeof(POLYLC_ALPHA_TRIES[0]))

/* True if `var` divides every term of P (equivalently: the constant
 * coefficient of P viewed as a polynomial in `var` is zero). */
static bool polylc_var_divides(Expr* P, Expr* var) {
    Expr* c0 = get_coeff(P, var, 0);
    bool zero = (c0->type == EXPR_INTEGER && c0->data.integer == 0);
    expr_free(c0);
    return zero;
}

/* Phase F1 Stage 3 entry: handle bivariate inputs whose lc_x(P) is a
 * non-constant polynomial in the other variable.  Returns Times[U, V]
 * on success (only r = 2 supported in MVP), or NULL when any
 * precondition fails (so the caller can fall through to legacy). */
static Expr* factor_bivariate_via_polylc_hensel(Expr* P, Expr** vars) {
    /* Bail if either variable divides P -- the legacy heuristic_factor
     * pipeline (monomial-content extraction in Phase 0) will produce a
     * cleaner factorisation than a 2-factor Hensel split that lumps the
     * monomial into one of the bivariate factors. */
    if (polylc_var_divides(P, vars[0]) || polylc_var_divides(P, vars[1])) {
        return NULL;
    }

    for (int main_idx = 0; main_idx < 2; main_idx++) {
        int other_idx = 1 - main_idx;
        Expr* x_var = vars[main_idx];
        Expr* y_var = vars[other_idx];

        int deg = get_degree_poly(P, x_var);
        if (deg < 2) continue;

        Expr* lc_expr = get_coeff(P, x_var, deg);
        ZUPoly* A = expr_to_zupoly(lc_expr, y_var);
        expr_free(lc_expr);
        if (!A) continue;
        if (A->deg < 1) {
            /* Constant LC -- handled by Stages 1 / 2. */
            zupoly_free(A);
            continue;
        }

        BPoly* bp = expr_to_bpoly(P, x_var, y_var);
        if (!bp) {
            zupoly_free(A);
            continue;
        }

        mpz_t a_content; mpz_init(a_content);
        ZUPoly** a_factors = NULL;
        int a_count = 0;
        bool fac_ok = zupoly_factor_to_array(A, a_content, &a_factors, &a_count);
        if (!fac_ok || a_count == 0
            || mpz_cmpabs_ui(a_content, 1) != 0) {
            mpz_clear(a_content);
            if (a_factors) {
                for (int i = 0; i < a_count; i++) zupoly_free(a_factors[i]);
                free(a_factors);
            }
            bpoly_free(bp);
            zupoly_free(A);
            continue;
        }

        Expr* result = NULL;
        for (size_t pi = 0; pi < POLYLC_ALPHA_COUNT && !result; pi++) {
            int64_t alpha = POLYLC_ALPHA_TRIES[pi];

            mpz_t A_at_alpha; mpz_init(A_at_alpha);
            zupoly_eval_si(A, alpha, A_at_alpha);
            bool a_alpha_ok = (mpz_cmp_si(A_at_alpha, 1) == 0);
            mpz_clear(A_at_alpha);
            if (!a_alpha_ok) continue;

            ZUPoly* image = bpoly_eval_y_si(bp, alpha);
            if (!image || image->deg != bp->deg_x) {
                if (image) zupoly_free(image);
                continue;
            }

            /* Squarefree check: gcd(image, image') is a unit. */
            ZUPoly* image_deriv = zupoly_new(image->deg);
            for (int i = 1; i <= image->deg; i++) {
                const mpz_t* ci = zupoly_getcoef(image, i);
                if (!ci) continue;
                mpz_t scaled; mpz_init(scaled);
                mpz_mul_ui(scaled, *ci, (unsigned long)i);
                zupoly_setcoef(image_deriv, i - 1, scaled);
                mpz_clear(scaled);
            }
            ZUPoly* g = zupoly_gcd(image, image_deriv);
            bool sqf = (g->deg == 0);
            zupoly_free(g);
            zupoly_free(image_deriv);
            if (!sqf) {
                zupoly_free(image);
                continue;
            }

            ZUPoly** us = NULL;
            int r = 0;
            bool fac2_ok = factor_via_bz_callback(image, &us, &r, NULL);
            zupoly_free(image);
            if (!fac2_ok || r != 2) {
                if (us) {
                    for (int i = 0; i < r; i++) zupoly_free(us[i]);
                    free(us);
                }
                continue;
            }

            bool all_monic =
                (mpz_cmp_ui(us[0]->c[us[0]->deg], 1) == 0) &&
                (mpz_cmp_ui(us[1]->c[us[1]->deg], 1) == 0);
            if (!all_monic) {
                for (int i = 0; i < r; i++) zupoly_free(us[i]);
                free(us);
                continue;
            }

            /* Shift y -> y + alpha so the lift sees α = 0.  (The picocas
             * convention: bpoly_shift_y_si(P, α) returns P(x, y + α).) */
            BPoly* P_shifted = bpoly_shift_y_si(bp, alpha);
            ZUPoly** a_factors_shifted = (ZUPoly**)malloc(
                sizeof(ZUPoly*) * (size_t)a_count);
            for (int i = 0; i < a_count; i++) {
                a_factors_shifted[i] = zupoly_shift_si(a_factors[i], alpha);
            }

            BPoly *U = NULL, *V = NULL;
            bool lift_ok = try_polylc_two_factor(
                P_shifted, us[0], us[1],
                a_content, a_factors_shifted, a_count,
                &U, &V);

            for (int i = 0; i < a_count; i++) zupoly_free(a_factors_shifted[i]);
            free(a_factors_shifted);
            for (int i = 0; i < r; i++) zupoly_free(us[i]);
            free(us);
            bpoly_free(P_shifted);

            if (!lift_ok) continue;

            BPoly* U_orig = bpoly_shift_y_si(U, -alpha);
            BPoly* V_orig = bpoly_shift_y_si(V, -alpha);
            bpoly_free(U);
            bpoly_free(V);

            Expr* eU = bpoly_to_expr(U_orig, x_var, y_var);
            Expr* eV = bpoly_to_expr(V_orig, x_var, y_var);
            bpoly_free(U_orig);
            bpoly_free(V_orig);

            result = eval_and_free(expr_new_function(
                expr_new_symbol("Times"),
                (Expr*[]){eU, eV}, 2));
        }

        mpz_clear(a_content);
        for (int i = 0; i < a_count; i++) zupoly_free(a_factors[i]);
        free(a_factors);
        bpoly_free(bp);
        zupoly_free(A);

        if (result) return result;
    }
    return NULL;
}

/* End-to-end bivariate factorisation via the new pipeline.  Returns
 * a fresh Times[...] Expr* on success or NULL when:
 *   - P has no variable with a constant integer leading coefficient
 *     (i.e. every variable's lc is polynomial in the other variable).
 *     Inputs with polynomial-in-y LC fall through to the legacy pipeline;
 *     Wang's full LC correction (Stage 3) is future work.
 *   - The orchestrator failed to find a good evaluation point.
 *   - The lift produced only one factor (irreducible -- but that is
 *     handled by the irreducibility short-circuit at a higher level,
 *     so we treat r==1 as "no progress" here so the caller sees a
 *     consistent NULL signal).
 *
 * Three correction paths, dispatched by `pick_factorable_x_var`:
 *   - BV_LC_MONIC  (lc = +1): lift P directly.
 *   - BV_LC_NEGATE (lc = -1, Phase F1 Stage 1): negate P, lift, flip
 *     one factor's sign (chosen heuristically so the canonical form
 *     matches the legacy output).
 *   - BV_LC_SCALE  (lc = constant integer ≠ ±1, Phase F1 Stage 2):
 *     Wang's monic substitution Q = a^(d-1) · P(x/a, y), which is
 *     monic in x with integer coefficients.  Lift Q to G_1, ..., G_r.
 *     Recover the true factors via F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y)),
 *     where cont_Z denotes the integer content (gcd of all coefficients).
 *     The product cont_Z(G_i(a·x, y)) over all i is exactly a^(d-1),
 *     so prod F_i = P.
 *
 * Stage 1 + Stage 2 compose: when lc_x(P) is a negative integer not
 * equal to -1 (e.g. lc = -6), the BV_LC_SCALE path runs first with
 * |a| = 6 after pre-negating P, and the final factor adjustment flips
 * the sign of one factor.
 *
 * Memory: the caller retains ownership of P and vars; we allocate
 * the result and any intermediate BPoly/ZUPoly are freed before
 * returning. */
/* Run the bivariate Hensel pipeline with a specific choice of main
 * variable index and leading-coefficient kind.  Returns the factored
 * Expr on success or NULL when the lift produces ≤1 factor or any
 * intermediate step fails.  Used by factor_bivariate_via_hensel which
 * tries each variable in priority order, falling back to the next when
 * the prior choice's image factorisation hits the BIGINT bail in the
 * BZ ingestion path (Phase F6). */
static Expr* try_bivariate_via_hensel_with_main(Expr* P, Expr** vars,
                                                int main_idx, BvLcKind kind,
                                                const mpz_t lc_a) {
    int other_idx = 1 - main_idx;
    bool negate = (kind == BV_LC_NEGATE) ||
                  (kind == BV_LC_SCALE && mpz_sgn(lc_a) < 0);
    bool scale  = (kind == BV_LC_SCALE);

    /* For Stage 2 we work with |a| -- the sign is absorbed by
     * pre-negating P so the substitution always uses a positive a. */
    mpz_t abs_a; mpz_init(abs_a);
    mpz_abs(abs_a, lc_a);

    BPoly* bp = expr_to_bpoly(P, vars[main_idx], vars[other_idx]);
    if (!bp) {
        mpz_clear(abs_a);
        return NULL;
    }

    if (negate) {
        BPoly* bp_neg = bpoly_neg(bp);
        bpoly_free(bp);
        bp = bp_neg;
    }

    /* Stage 2: substitute to make the lift's input monic.  bp now has
     * lc_x = +|a|; Q has lc_x = 1. */
    BPoly* Q = NULL;
    if (scale) {
        Q = bpoly_make_monic_via_x_scale(bp, abs_a);
        if (!Q) {
            bpoly_free(bp);
            mpz_clear(abs_a);
            return NULL;
        }
    }
    const BPoly* lift_input = scale ? Q : bp;

    BPoly** factors = NULL;
    int r = 0;
    bool ok = mvfactor_try_bivariate_monic(lift_input, factor_via_bz_callback,
                                           NULL, &factors, &r);
    bpoly_free(bp);
    if (Q) bpoly_free(Q);
    if (!ok) {
        mpz_clear(abs_a);
        return NULL;
    }

    /* r == 1 means the orchestrator concluded irreducibility (or
     * couldn't find a non-trivial factorisation).  For our caller's
     * contract we return NULL in either case so the legacy pipeline
     * gets its turn. */
    if (r <= 1) {
        if (factors) {
            for (int i = 0; i < r; i++) bpoly_free(factors[i]);
            free(factors);
        }
        mpz_clear(abs_a);
        return NULL;
    }

    /* Stage 2 recovery: the lift returned factors G_i of Q (monic).
     * The true factors of P (with lc = +|a|) are
     *   F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y)).
     * The integer content of G_i(a·x, y) collects exactly the share
     * of a^(d-1) that was redistributed into G_i by the substitution.
     *
     * We rewrite each `factors[i]` in place with the recovered F_i. */
    if (scale) {
        for (int i = 0; i < r; i++) {
            BPoly* H = bpoly_subst_x_scale(factors[i], abs_a);
            mpz_t ci; mpz_init(ci);
            bpoly_int_content(H, ci);
            if (mpz_sgn(ci) == 0) {
                /* Zero factor -- defensive guard, should not happen. */
                mpz_clear(ci);
                bpoly_free(H);
                for (int k = 0; k < r; k++) bpoly_free(factors[k]);
                free(factors);
                mpz_clear(abs_a);
                return NULL;
            }
            BPoly* F = bpoly_div_int_exact(H, ci);
            mpz_clear(ci);
            bpoly_free(H);
            if (!F) {
                /* Non-exact division -- defensive guard, should not happen
                 * for a faithful Wang's correction. */
                for (int k = 0; k < r; k++) bpoly_free(factors[k]);
                free(factors);
                mpz_clear(abs_a);
                return NULL;
            }
            bpoly_free(factors[i]);
            factors[i] = F;
        }
    }
    mpz_clear(abs_a);

    /* Convert each BPoly factor back to Expr.  When we negated the
     * input, the lift's product equals -P (or -(absolute lc form));
     * we must absorb the overall -1 into one factor.
     *
     * Choice of which factor to negate: prefer the largest-x-degree
     * factor (tiebreak: largest y-degree).  Empirically this matches
     * the canonical printed form for inputs like
     *   Factor[3 a^2 b - 3 b - b^3]
     * where the existing baseline is `b · (-3 + 3 a^2 - b^2)` (the
     * `b` factor unsigned, the polynomial residue absorbing the -1).
     * Choosing "lowest x-degree" instead would produce
     * `-b · (3 - 3 a^2 + b^2)` which is mathematically equivalent
     * but displays the sign explicitly.  Picking the largest factor
     * keeps the output stable across F1's introduction. */
    int neg_idx = 0;
    if (negate) {
        int best_dx = factors[0]->deg_x;
        int best_dy = bpoly_deg_y(factors[0]);
        for (int i = 1; i < r; i++) {
            int dx = factors[i]->deg_x;
            int dy = bpoly_deg_y(factors[i]);
            if (dx > best_dx || (dx == best_dx && dy > best_dy)) {
                best_dx = dx;
                best_dy = dy;
                neg_idx = i;
            }
        }
    }

    Expr** args = (Expr**)malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        args[i] = bpoly_to_expr(factors[i], vars[main_idx], vars[other_idx]);
        bpoly_free(factors[i]);
    }
    free(factors);

    if (negate) {
        /* Distribute -1 INTO the chosen factor via Expand so the
         * sign is folded into a Plus (e.g., -(b^2 - 3a^2 + 3) becomes
         * Plus[-3, 3 a^2, -b^2]) rather than parked as a leading
         * Times[-1, ...].  Without this Expand, the result evaluates
         * to Times[-1, ..., Plus[...]] which prints as `-(...)` and
         * does not match the canonical form picocas previously
         * produced for inputs that fall through to the lift via
         * monomial-content recursion. */
        Expr* neg_times = expr_new_function(
            expr_new_symbol("Times"),
            (Expr*[]){ expr_new_integer(-1), args[neg_idx] }, 2);
        Expr* flipped = eval_and_free(expr_new_function(
            expr_new_symbol("Expand"),
            (Expr*[]){ neg_times }, 1));
        args[neg_idx] = flipped;
    }

    Expr* result = expr_new_function(expr_new_symbol("Times"), args, (size_t)r);
    free(args);
    return eval_and_free(result);
}

/* Bivariate Hensel entry point with main-variable fallback (Phase F6).
 *
 * The picker chooses one variable based on its leading-coefficient kind
 * priority (MONIC > NEGATE > SCALE).  For inputs like the Fateman
 * benchmark `Expand[(x^200 + y^200 + 1)(x^200 - y + 1)]` the picker
 * prefers `x` (lc_x = +1, MONIC) but `P(x, alpha)` for any |alpha| ≥ 2
 * produces alpha^200 -- a 60-digit BIGINT that the int64 UPoly substrate
 * inside `bz_factor_to_expr` cannot ingest, so every alpha bails and the
 * orchestrator returns "no factorisation".  Picking `y` instead (lc_y =
 * -1, NEGATE) and evaluating at alpha_x = 0 yields a clean image
 * `(1-y)(1+y^200)` with ±1 coefficients that BZ handles fine.
 *
 * We try the picker's preferred variable first, then fall back to the
 * other variable (with whatever lc kind it has) if the first attempt
 * yields no progress.  The polynomial-LC fallback fires only when the
 * very first picker call returns -1 (no integer-LC variable at all).
 */
static Expr* factor_bivariate_via_hensel(Expr* P, Expr** vars) {
    int tried_idx = -1;
    for (int attempt = 0; attempt < 2; attempt++) {
        BvLcKind kind = BV_LC_MONIC;
        mpz_t lc_a; mpz_init(lc_a);
        int main_idx = pick_factorable_x_var(P, vars, 2, &kind, &lc_a, tried_idx);
        if (main_idx < 0) {
            mpz_clear(lc_a);
            if (attempt == 0) {
                /* No variable has a constant integer leading coefficient
                 * -- lc_x(P) is a polynomial in the other variable.  Try
                 * the Stage 3 (Wang's leading-coefficient correction)
                 * path. */
                return factor_bivariate_via_polylc_hensel(P, vars);
            }
            break;
        }
        Expr* res = try_bivariate_via_hensel_with_main(P, vars, main_idx, kind, lc_a);
        mpz_clear(lc_a);
        if (res) return res;
        tried_idx = main_idx;
    }
    return NULL;
}

/* ====================================================================== */
/*  Phase 4: n-variate specialise-and-trial-divide                        */
/*                                                                        */
/*  For inputs with v_count >= 3, attempt to find a factor that does not  */
/*  depend on one of the variables.  Strategy:                            */
/*    1. Pick a "specialisation" variable z.                              */
/*    2. Evaluate P|z=0 and recursively factor it via heuristic_factor.   */
/*    3. For each non-trivial factor f of P|z=0, trial-divide the         */
/*       ORIGINAL P by f.  If f divides exactly, it is a true factor of  */
/*       P (one that happens to be z-independent).                        */
/*    4. After collecting z-independent factors, the residual is the     */
/*       product of remaining (z-dependent) factors; recurse on it.      */
/*                                                                        */
/*  This catches the common case where the multivariate polynomial       */
/*  factors cleanly into a "shape" piece (z-independent) and a           */
/*  "z-binding" piece, e.g.                                               */
/*    P(x, y, z) = (x + y) (x + y + z)  → factor (x + y) found at z=0    */
/*  Cases where ALL factors depend on z fall through to the legacy       */
/*  factor_roots pipeline (Phase 4 returns NULL).                         */
/*                                                                        */
/*  Try each variable as the candidate z; the first specialisation that  */
/*  yields ≥ 1 trial-divisible bivariate factor wins.                     */
/* ====================================================================== */

static Expr* heuristic_factor(Expr* P);  /* forward */

/* Walk a Times expression and copy its arguments into args/count_out.
 * If `e` is not a Times, treats it as a single factor.  Caller owns
 * each copied Expr* and the array. */
static void collect_times_args(Expr* e, Expr*** args_out, size_t* count_out) {
    if (e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        size_t n = e->data.function.arg_count;
        Expr** args = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            args[i] = expr_copy(e->data.function.args[i]);
        }
        *args_out = args;
        *count_out = n;
        return;
    }
    Expr** args = (Expr**)malloc(sizeof(Expr*));
    args[0] = expr_copy(e);
    *args_out = args;
    *count_out = 1;
}

/* True if `f` is a numeric or rational constant (no polynomial structure). */
static bool factor_is_numeric_constant(const Expr* f) {
    if (f->type == EXPR_INTEGER || f->type == EXPR_REAL || f->type == EXPR_BIGINT) {
        return true;
    }
    if (f->type == EXPR_FUNCTION
        && f->data.function.head
        && f->data.function.head->type == EXPR_SYMBOL
        && (f->data.function.head->data.symbol == SYM_Rational
            || f->data.function.head->data.symbol == SYM_Complex)) {
        return true;
    }
    return false;
}

/* Substitute vars[var_idx] -> 0 in P. */
static Expr* expr_substitute_var_zero(Expr* P, Expr* var) {
    Expr* rule = expr_new_function(
        expr_new_symbol("Rule"),
        (Expr*[]){ expr_copy(var), expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(
        expr_new_symbol("ReplaceAll"),
        (Expr*[]){ expr_copy(P), rule }, 2);
    Expr* result = evaluate(call);
    expr_free(call);
    return result;
}

static Expr* factor_via_z_independent_split(Expr* P, Expr** vars, size_t v_count) {
    if (v_count < 3) return NULL;

    /* Try each variable as the "z" specialisation. */
    for (size_t z_idx = 0; z_idx < v_count; z_idx++) {
        /* P_at_z0 = P with vars[z_idx] -> 0. */
        Expr* P_at_z0 = expr_substitute_var_zero(P, vars[z_idx]);
        if (!P_at_z0) continue;

        /* If the substitution yielded a constant, skip (this z value
         * killed too much; retry with a different z). */
        if (P_at_z0->type != EXPR_FUNCTION) {
            expr_free(P_at_z0);
            continue;
        }

        /* Recursively factor P_at_z0 (now in v_count - 1 variables). */
        Expr* factored = heuristic_factor(P_at_z0);
        expr_free(P_at_z0);
        if (!factored) continue;

        /* Walk the factors.  For each non-constant factor f, trial-
         * divide the ORIGINAL P by f.  Each successful division is a
         * genuine z-independent factor of P. */
        Expr** factor_args = NULL;
        size_t factor_count = 0;
        collect_times_args(factored, &factor_args, &factor_count);
        expr_free(factored);

        Expr* remaining = expr_copy(P);
        Expr** found = NULL;
        size_t found_n = 0;
        size_t found_cap = 0;

        for (size_t i = 0; i < factor_count; i++) {
            Expr* f = factor_args[i];
            if (factor_is_numeric_constant(f)) continue;

            /* Power[base, k] is a factor of multiplicity k -- attempt
             * the base, and on success peel off all k copies. */
            Expr* base = f;
            int64_t mult = 1;
            if (f->type == EXPR_FUNCTION
                && f->data.function.head
                && f->data.function.head->type == EXPR_SYMBOL
                && f->data.function.head->data.symbol == SYM_Power
                && f->data.function.arg_count == 2
                && f->data.function.args[1]->type == EXPR_INTEGER
                && f->data.function.args[1]->data.integer >= 1) {
                base = f->data.function.args[0];
                mult = f->data.function.args[1]->data.integer;
            }

            /* Trial-divide remaining by base, repeatedly. */
            for (int64_t k = 0; k < mult; k++) {
                Expr* q = exact_poly_div(remaining, base, vars, v_count);
                if (!q) {
                    /* If we already peeled k > 0 copies, that's fine;
                     * stop peeling further at this base. */
                    break;
                }
                if (found_n >= found_cap) {
                    found_cap = found_cap ? found_cap * 2 : 4;
                    found = (Expr**)realloc(found, sizeof(Expr*) * found_cap);
                }
                found[found_n++] = expr_copy(base);
                expr_free(remaining);
                remaining = q;
            }
        }

        for (size_t i = 0; i < factor_count; i++) expr_free(factor_args[i]);
        free(factor_args);

        if (found_n == 0) {
            /* No z-independent factor found at this specialisation. */
            expr_free(remaining);
            free(found);
            continue;
        }

        /* Recurse on remaining (which has the discovered factors removed). */
        Expr* rem_factored = heuristic_factor(remaining);
        expr_free(remaining);

        /* Build Times[found..., rem_factored]. */
        size_t total = found_n + 1;
        Expr** times_args = (Expr**)malloc(sizeof(Expr*) * total);
        for (size_t i = 0; i < found_n; i++) times_args[i] = found[i];
        times_args[found_n] = rem_factored;
        free(found);

        Expr* result = expr_new_function(expr_new_symbol("Times"),
                                         times_args, total);
        free(times_args);
        return eval_and_free(result);
    }

    return NULL;
}

/* ===================================================================== */
/*  Phase F2 MVP -- trivariate two-factor Hensel via MPoly                */
/*                                                                       */
/*  Handles n = 3 polynomials that are monic in some main variable and   */
/*  factor as a product of exactly two genuinely-trivariate factors.     */
/*  This is the case Phase 4's specialise-and-divide (factor_via_z_      */
/*  independent_split) cannot catch -- where every factor depends on     */
/*  every variable.                                                      */
/*                                                                       */
/*  Algorithm:                                                           */
/*    1. Convert P to MPoly.                                             */
/*    2. For each candidate main_idx with monic LC:                      */
/*       For each (alpha_y, alpha_z) tuple in a small budget:            */
/*         a. Compute univariate image P|y=alpha_y, z=alpha_z.           */
/*         b. Verify squarefree, degree-preserving, two-factor.          */
/*         c. Lift y-direction via existing bivariate Hensel             */
/*            (bpoly_hensel_lift_2) over BPoly{main, var_y}.             */
/*         d. Lift z-direction via mpoly_hensel_lift_3_2.                */
/*         e. Verify product == P, return factored.                      */
/*  Returns NULL if no productive (alpha_y, alpha_z, main_idx) was       */
/*  found within the search budget.                                      */
/* ===================================================================== */

static Expr* factor_trivariate_via_mhensel(Expr* P, Expr** vars, size_t v_count) {
    if (v_count != 3) return NULL;

    MPoly* P_mp = expr_to_mpoly(P, vars, (int)v_count);
    if (!P_mp) return NULL;

    Expr* result = NULL;

    static const int64_t alpha_choices[] = { 0, 1, -1, 2, -2, 3, -3 };
    enum { N_ALPHAS = sizeof(alpha_choices) / sizeof(alpha_choices[0]) };

    for (int main_idx = 0; main_idx < 3 && !result; main_idx++) {
        int d_main = mpoly_deg_var(P_mp, main_idx);
        if (d_main < 1) continue;

        /* Check monic-in-main: LC must be the constant integer 1. */
        MPoly* lc = mpoly_lc_var(P_mp, main_idx);
        bool monic = (lc->n_terms == 1) &&
                     (mpz_cmp_ui(lc->coefs[0], 1) == 0);
        if (monic) {
            const int* row = lc->exps;
            for (int v = 0; v < (int)v_count; v++) {
                if (row[v] != 0) { monic = false; break; }
            }
        }
        mpoly_free(lc);
        if (!monic) continue;

        /* Pick the two secondary variables. */
        int v_y = -1, v_z = -1;
        for (int i = 0; i < (int)v_count; i++) {
            if (i == main_idx) continue;
            if (v_y < 0) v_y = i;
            else        v_z = i;
        }

        /* Try (alpha_y, alpha_z) tuples. */
        for (size_t ai = 0; ai < N_ALPHAS && !result; ai++) {
            for (size_t bi = 0; bi < N_ALPHAS && !result; bi++) {
                int64_t alpha_y = alpha_choices[ai];
                int64_t alpha_z = alpha_choices[bi];

                /* Univariate image at the alpha tuple. */
                MPoly* P_y  = mpoly_subst_var_int(P_mp, v_y, alpha_y);
                MPoly* P_yz = mpoly_subst_var_int(P_y,  v_z, alpha_z);
                mpoly_free(P_y);
                ZUPoly* uni = mpoly_to_zupoly_in(P_yz, main_idx);
                mpoly_free(P_yz);
                if (!uni || uni->deg != d_main) {
                    zupoly_free(uni);
                    continue;
                }

                /* Squarefree probe via gcd(uni, uni'). */
                ZUPoly* uni_d = zupoly_new(uni->deg);
                for (int i = 1; i <= uni->deg; i++) {
                    const mpz_t* ci = zupoly_getcoef(uni, i);
                    if (!ci) continue;
                    mpz_t scaled; mpz_init(scaled);
                    mpz_mul_ui(scaled, *ci, (unsigned long)i);
                    zupoly_setcoef(uni_d, i - 1, scaled);
                    mpz_clear(scaled);
                }
                ZUPoly* g = zupoly_gcd(uni, uni_d);
                bool sqf = (g->deg == 0);
                zupoly_free(g); zupoly_free(uni_d);
                if (!sqf) {
                    zupoly_free(uni);
                    continue;
                }

                /* Factor univariately.  We piggyback on the existing
                 * `factor_via_bz_callback` (defined earlier in this
                 * file) which does ZUPoly -> Expr -> bz_factor_to_expr
                 * -> ZUPoly. */
                ZUPoly** uni_factors = NULL;
                int r = 0;
                bool fac_ok = factor_via_bz_callback(uni, &uni_factors, &r, NULL);
                zupoly_free(uni);
                if (!fac_ok || r != 2 || !uni_factors) {
                    if (uni_factors) {
                        for (int i = 0; i < r; i++) zupoly_free(uni_factors[i]);
                        free(uni_factors);
                    }
                    continue;
                }
                ZUPoly* u = uni_factors[0];
                ZUPoly* v = uni_factors[1];

                /* Step 1 -- y-lift via existing bivariate Hensel.
                 * Specialise z = alpha_z on the original P_mp; convert to
                 * BPoly{main, v_y}; run bpoly_hensel_lift_2 with seeds u, v. */
                MPoly* P_z   = mpoly_subst_var_int(P_mp, v_z, alpha_z);
                /* Shift y = y + alpha_y so the bivariate lift sees u, v at y = 0.
                 * bpoly_hensel_lift_2 expects u, v to satisfy u*v = P_xy(x, 0). */
                MPoly* P_z_sh = mpoly_shift_var_int(P_z, v_y, alpha_y);
                BPoly* P_xy_b = mpoly_to_bpoly_in(P_z_sh, main_idx, v_y);
                mpoly_free(P_z); mpoly_free(P_z_sh);
                if (!P_xy_b) {
                    zupoly_free(u); zupoly_free(v); free(uni_factors);
                    continue;
                }

                BPoly* U_xy_b = NULL;
                BPoly* V_xy_b = NULL;
                bool lift2_ok = bpoly_hensel_lift_2(P_xy_b, u, v,
                                                    &U_xy_b, &V_xy_b);
                bpoly_free(P_xy_b);
                zupoly_free(u); zupoly_free(v); free(uni_factors);
                if (!lift2_ok) {
                    bpoly_free(U_xy_b); bpoly_free(V_xy_b);
                    continue;
                }

                /* Convert BPoly factors back to MPoly{n=3 vars}, then
                 * unshift y = y - alpha_y. */
                MPoly* U_xy_mp_sh = bpoly_to_mpoly_in(U_xy_b, (int)v_count,
                                                     main_idx, v_y);
                MPoly* V_xy_mp_sh = bpoly_to_mpoly_in(V_xy_b, (int)v_count,
                                                     main_idx, v_y);
                bpoly_free(U_xy_b); bpoly_free(V_xy_b);
                MPoly* U_xy_mp = mpoly_shift_var_int(U_xy_mp_sh, v_y, -alpha_y);
                MPoly* V_xy_mp = mpoly_shift_var_int(V_xy_mp_sh, v_y, -alpha_y);
                mpoly_free(U_xy_mp_sh); mpoly_free(V_xy_mp_sh);

                /* Step 2 -- z-lift via mpoly_hensel_lift_3_2. */
                MPoly* U_full = NULL;
                MPoly* V_full = NULL;
                bool lift3_ok = mpoly_hensel_lift_3_2(
                    P_mp, U_xy_mp, V_xy_mp,
                    main_idx, v_y, v_z, alpha_z,
                    &U_full, &V_full);
                mpoly_free(U_xy_mp); mpoly_free(V_xy_mp);
                if (!lift3_ok) {
                    mpoly_free(U_full); mpoly_free(V_full);
                    continue;
                }

                /* Verify product == P_mp. */
                MPoly* product = mpoly_mul(U_full, V_full);
                bool match = mpoly_eq(product, P_mp);
                mpoly_free(product);
                if (!match) {
                    mpoly_free(U_full); mpoly_free(V_full);
                    continue;
                }

                /* Convert factors back to Expr. */
                Expr* U_e = mpoly_to_expr(U_full, vars);
                Expr* V_e = mpoly_to_expr(V_full, vars);
                mpoly_free(U_full); mpoly_free(V_full);

                Expr* U_eval = evaluate(U_e); expr_free(U_e);
                Expr* V_eval = evaluate(V_e); expr_free(V_e);

                Expr** factor_args = (Expr**)malloc(sizeof(Expr*) * 2);
                factor_args[0] = U_eval;
                factor_args[1] = V_eval;
                Expr* product_expr = expr_new_function(
                    expr_new_symbol("Times"), factor_args, 2);
                result = evaluate(product_expr);
                expr_free(product_expr);
            }
        }
    }

    mpoly_free(P_mp);
    return result;
}

/* Reentrancy guard for heuristic_factor.  The cont/pp split below can
 * keep producing inputs of the same variable count when poly_content
 * extracts a non-trivial GCD that itself shares all variables of P.
 * factor_via_z_independent_split, in turn, recurses on its `remaining`
 * factor at the same v_count when the polynomial structure looks like
 * a true atom in disguise (e.g. Power[base, n] with symbolic n that
 * the variable collector treats as opaque).  Without this guard the
 * cycle blows the C stack -- see the (Cos[x]-1)^n (Cos[x]+1)^n -
 * (-Sin[x]^2)^n regression. */
#define HEURISTIC_FACTOR_MAX_DEPTH 64
static int heuristic_factor_depth = 0;

static Expr* heuristic_factor_impl(Expr* P);

static Expr* heuristic_factor(Expr* P) {
    if (P->type != EXPR_FUNCTION) return expr_copy(P);
    if (heuristic_factor_depth >= HEURISTIC_FACTOR_MAX_DEPTH) {
        return expr_copy(P);
    }
    heuristic_factor_depth++;
    Expr* r = heuristic_factor_impl(P);
    heuristic_factor_depth--;
    return r;
}

static Expr* heuristic_factor_impl(Expr* P) {
    if (P->data.function.head->data.symbol == SYM_Times || P->data.function.head->data.symbol == SYM_Power) {
        Expr** args = malloc(sizeof(Expr*) * P->data.function.arg_count);
        for(size_t i=0; i<P->data.function.arg_count; i++) args[i] = heuristic_factor(P->data.function.args[i]);
        Expr* res = eval_and_free(expr_new_function(expr_copy(P->data.function.head), args, P->data.function.arg_count));
        free(args);
        return res;
    }
    if (P->data.function.head->data.symbol == SYM_List || P->data.function.head->data.symbol == SYM_Equal || P->data.function.head->data.symbol == SYM_Less ||
        P->data.function.head->data.symbol == SYM_LessEqual || P->data.function.head->data.symbol == SYM_Greater || P->data.function.head->data.symbol == SYM_GreaterEqual ||
        P->data.function.head->data.symbol == SYM_And || P->data.function.head->data.symbol == SYM_Or || P->data.function.head->data.symbol == SYM_Not) {
        Expr** args = malloc(sizeof(Expr*) * P->data.function.arg_count);
        for(size_t i=0; i<P->data.function.arg_count; i++) args[i] = heuristic_factor(P->data.function.args[i]);
        Expr* res = eval_and_free(expr_new_function(expr_copy(P->data.function.head), args, P->data.function.arg_count));
        free(args);
        return res;
    }

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(P, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    /* No variables -> constant; content-recursion would not terminate on
     * rationals or non-unit integers, so return the constant itself. */
    if (v_count == 0) {
        free(vars);
        return expr_copy(P);
    }

    Expr* cont = poly_content(P, vars, v_count);
    if (!(cont->type == EXPR_INTEGER && cont->data.integer == 1)) {
        Expr* pp = exact_poly_div(P, cont, vars, v_count);
        /* Guard against the no-progress case: if exact_poly_div left us
         * with pp == 1 (i.e. cont == P up to a constant), recursing on
         * cont = P loops indefinitely.  Bail and return the input. */
        if (pp && pp->type == EXPR_INTEGER && pp->data.integer == 1) {
            { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
            expr_free(cont); expr_free(pp);
            return expr_copy(P);
        }
        Expr* f_cont = heuristic_factor(cont);
        Expr* f_pp = heuristic_factor(pp);
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        expr_free(cont); expr_free(pp);
        return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){f_cont, f_pp}, 2));
    }
    expr_free(cont);

    /* Monomial-content extraction: pull out v_1^{e_1} * ... * v_k^{e_k}
     * shared by every term.  This is the cheapest factorisation step
     * (no polynomial arithmetic) and is decisively correct: if the GCD
     * of the term-monomials is nontrivial, the result is exactly that
     * GCD times the residue, with the residue then factored recursively
     * by the strategies below.  For inputs like
     *   3 a^2 b - 3 b - b^3      ->  b * (3 a^2 - 3 - b^2)
     *   3 Sin[x]^3 - 3 Sin[x]    ->  Sin[x] * (3 Sin[x]^2 - 3)
     * this is the difference between a correct factorisation and a
     * fall-through to expensive trial-root methods that miss the
     * monomial structure entirely. */
    Expr* mc = factor_monomial_content(P);
    if (mc) {
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return mc;
    }

    Expr* d1 = factor_degree_one(P, vars, v_count);
    if (d1) { 
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return d1; 
    }
    
    if (v_count == 1) {
        Expr* bz_res = bz_factor_to_expr(P, vars[0]);
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return bz_res;
    }
    
    Expr* bn = factor_binomial(P);
    if (bn) {
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return bn;
    }

    /* Try the productive multivariate paths FIRST.  These are
     * generally fast on factorable inputs (succeed quickly with a
     * non-trivial factorisation) and fail quickly on irreducible
     * inputs (no good evaluation point found).  Running them before
     * the is_likely_irreducible probe avoids paying the probe's cost
     * (10s of ms × number of variables) on inputs that the productive
     * paths handle directly. */
    if (v_count == 2) {
        Expr* bivariate = factor_bivariate_via_hensel(P, vars);
        if (bivariate) {
            { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
            return bivariate;
        }
    }

    /* Phase 4: for n >= 3 variables, attempt to peel off factors that
     * don't depend on one of the variables.  Specialise each variable
     * in turn to 0, recursively factor the result, and trial-divide
     * the original P by each candidate factor.  Catches polynomials
     * with a "shape" piece independent of the specialised variable. */
    if (v_count >= 3) {
        Expr* tri = factor_via_z_independent_split(P, vars, v_count);
        if (tri) {
            { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
            return tri;
        }
    }

    /* Phase F2 MVP: trivariate two-factor monic Hensel via MPoly.
     * Catches n=3 polynomials whose factors depend on every variable
     * (which Phase 4 cannot find because no variable specialisation
     * exposes a clean factor).  Currently restricted to two-factor
     * outputs with monic-in-main inputs; extends to higher n and
     * non-monic Wang correction in future phases. */
    if (v_count == 3) {
        Expr* tri_mh = factor_trivariate_via_mhensel(P, vars, v_count);
        if (tri_mh) {
            { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
            return tri_mh;
        }
    }

    /* Irreducibility gate: the productive paths above didn't find a
     * factorisation; we're about to fall through to factor_roots,
     * which can be slow on irreducible inputs.  Run the cheap
     * irreducibility probe first to short-circuit. */
    if (is_likely_irreducible_multivariate(P, vars, v_count)) {
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return expr_copy(P);
    }

    Expr* rt = factor_roots(P, vars, v_count);
    if (rt) { 
        { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
        return rt; 
    }
    
    { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
    return expr_copy(P);
}

// Factor leverages a structural Berlekamp-Zassenhaus algorithm approach.
// For the scope of PicoCAS and performance bounds, the Zassenhaus recombination 
// phase is aggressively optimized utilizing exact discrete root evaluations
// and homogeneous Binomial descent, providing exact polynomial splittings
// over Z without generating Mignotte bound overflows.

/* ===================================================================== */
/*  Per-call Factor memo                                                 */
/*                                                                       */
/*  Cache `Factor[poly]` results inside a single `Simplify` call.  See   */
/*  facpoly.h for the lifecycle contract.  The top-of-stack memo is      */
/*  consulted at entry to `builtin_factor`; on a hit we return a deep    */
/*  copy of the cached result.  On a miss we run the full pipeline,     */
/*  then store a deep copy of the input/output pair before returning.   */
/* ===================================================================== */

#define FACTOR_MEMO_BUCKETS 256

typedef struct FactorMemoEntry {
    Expr* key;
    Expr* value;
    struct FactorMemoEntry* next;
} FactorMemoEntry;

struct FactorMemo {
    FactorMemoEntry* buckets[FACTOR_MEMO_BUCKETS];
};

/* Stack of currently-active memos.  We use a stack to support
 * (in principle) nested Simplify calls, even though that's not the
 * common case.  In practice only the top entry is consulted. */
#define FACTOR_MEMO_STACK_DEPTH 8
static FactorMemo* g_memo_stack[FACTOR_MEMO_STACK_DEPTH] = {0};
static int g_memo_stack_top = -1;

FactorMemo* factor_memo_new(void) {
    FactorMemo* m = (FactorMemo*)calloc(1, sizeof(FactorMemo));
    return m;
}

void factor_memo_free(FactorMemo* m) {
    if (!m) return;
    for (int i = 0; i < FACTOR_MEMO_BUCKETS; i++) {
        FactorMemoEntry* e = m->buckets[i];
        while (e) {
            FactorMemoEntry* next = e->next;
            expr_free(e->key);
            expr_free(e->value);
            free(e);
            e = next;
        }
    }
    free(m);
}

void factor_memo_push(FactorMemo* m) {
    if (g_memo_stack_top + 1 >= FACTOR_MEMO_STACK_DEPTH) return;
    g_memo_stack[++g_memo_stack_top] = m;
}

void factor_memo_pop(void) {
    if (g_memo_stack_top < 0) return;
    g_memo_stack[g_memo_stack_top--] = NULL;
}

static FactorMemo* factor_memo_top(void) {
    if (g_memo_stack_top < 0) return NULL;
    return g_memo_stack[g_memo_stack_top];
}

FactorMemo* factor_memo_active(void) { return factor_memo_top(); }

const Expr* factor_memo_lookup(FactorMemo* m, Expr* key) {
    if (!m) return NULL;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    for (FactorMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

void factor_memo_store(FactorMemo* m, Expr* key, Expr* value) {
    if (!m) return;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    FactorMemoEntry* e = (FactorMemoEntry*)malloc(sizeof(FactorMemoEntry));
    if (!e) return;
    e->key = expr_copy(key);
    e->value = expr_copy(value);
    e->next = m->buckets[h];
    m->buckets[h] = e;
}

/* Lookup; returns a borrowed pointer to the cached value (do NOT
 * free) or NULL on miss. */
static const Expr* factor_memo_get(FactorMemo* m, Expr* key) {
    if (!m) return NULL;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    for (FactorMemoEntry* e = m->buckets[h]; e; e = e->next) {
        if (expr_eq(e->key, key)) return e->value;
    }
    return NULL;
}

/* Store deep copies of key and value. */
static void factor_memo_put(FactorMemo* m, Expr* key, Expr* value) {
    if (!m) return;
    uint64_t h = expr_hash(key) % FACTOR_MEMO_BUCKETS;
    FactorMemoEntry* e = (FactorMemoEntry*)malloc(sizeof(FactorMemoEntry));
    if (!e) return;
    e->key = expr_copy(key);
    e->value = expr_copy(value);
    e->next = m->buckets[h];
    m->buckets[h] = e;
}

Expr* builtin_factor(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* Phase G5: Extension -> α option.  Walk trailing args looking for
     * Rule[Extension, α] / RuleDelayed[Extension, α].  When present,
     * dispatch to the algebraic-factoring path in qafactor.c (Trager).
     *
     * The memo is bypassed for this branch: the cache key is the
     * polynomial alone, but the result depends on the extension too. */
    if (res->data.function.arg_count >= 2) {
        Expr* alpha_expr = NULL;
        for (size_t i = 1; i < res->data.function.arg_count; i++) {
            Expr* opt = res->data.function.args[i];
            if (opt->type == EXPR_FUNCTION
                && opt->data.function.head
                && opt->data.function.head->type == EXPR_SYMBOL
                && (opt->data.function.head->data.symbol == SYM_Rule
                    || opt->data.function.head->data.symbol == SYM_RuleDelayed)
                && opt->data.function.arg_count == 2) {
                Expr* lhs = opt->data.function.args[0];
                Expr* rhs = opt->data.function.args[1];
                if (lhs->type == EXPR_SYMBOL && lhs->data.symbol == SYM_Extension) {
                    alpha_expr = rhs;
                }
            }
        }
        if (alpha_expr) {
            Expr* poly = res->data.function.args[0];

            /* Phase G6: when α is a List[α_1, α_2, ...], dispatch
             * to the tower-of-extensions path.  Otherwise (G5) treat
             * it as a single algebraic generator. */
            Expr** alpha_list = NULL;
            int n_alphas = 0;
            bool is_tower = false;
            if (alpha_expr->type == EXPR_FUNCTION
                && alpha_expr->data.function.head
                && alpha_expr->data.function.head->type == EXPR_SYMBOL
                && alpha_expr->data.function.head->data.symbol == SYM_List) {
                is_tower = true;
                n_alphas = (int)alpha_expr->data.function.arg_count;
                alpha_list = alpha_expr->data.function.args;  /* borrowed */
            }

            /* Pick the polynomial variable.  Q(α)-factoring works in a
             * single indeterminate; we collect the free symbols of poly
             * and require exactly one (after stripping α's surface
             * form).  Multivariate algebraic factoring is out of MVP
             * scope. */
            Expr** vars = NULL;
            size_t vc = 0, vcap = 8;
            vars = (Expr**)malloc(sizeof(Expr*) * vcap);
            collect_variables(poly, &vars, &vc, &vcap);

            Expr* poly_var = NULL;
            size_t live = 0;
            for (size_t i = 0; i < vc; i++) {
                bool is_alpha = false;
                if (is_tower) {
                    for (int k = 0; k < n_alphas; k++) {
                        if (expr_eq(vars[i], alpha_list[k])) {
                            is_alpha = true;
                            break;
                        }
                    }
                } else if (expr_eq(vars[i], alpha_expr)) {
                    is_alpha = true;
                }
                if (is_alpha) continue;
                poly_var = vars[i];
                live++;
            }
            if (live == 1 && poly_var) {
                Expr* result;
                if (is_tower) {
                    result = qa_factor_with_extension_tower(
                        poly, alpha_list, n_alphas, poly_var);
                } else {
                    result = qa_factor_with_extension(poly, alpha_expr, poly_var);
                }
                for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
                free(vars);
                if (result) return result;
                /* Fall through to plain Factor on failure. */
            } else {
                for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
                free(vars);
            }
        }
    }

    /* Memo lookup: if a Simplify call has installed a Factor memo and
     * we've already factored this exact input, return the cached copy.
     * Note: builtin_factor follows the convention that the evaluator
     * owns `res` (we do not call expr_free(res) anywhere in this
     * function), so the hit path simply returns a fresh copy of the
     * cached value without touching `res`. */
    FactorMemo* memo = factor_memo_top();
    if (memo) {
        const Expr* hit = factor_memo_get(memo, res);
        if (hit) {
            return expr_copy((Expr*)hit);
        }
    }
    Expr* memo_key = memo ? expr_copy(res) : NULL;

    /* Berlekamp–Zassenhaus and the heuristic multivariate paths assume
     * integer (rational) coefficients; route inexact inputs through the
     * standard rationalise / numericalise round-trip.  We don't cache
     * inexact results -- they go through their own evaluation path. */
    if (internal_args_contain_inexact(res)) {
        if (memo_key) expr_free(memo_key);
        return internal_rationalize_then_numericalize(res, builtin_factor);
    }

    /* Thread over logic / comparison heads (Less, Greater, Equal,
     * LessEqual, GreaterEqual, And, Or).  Mathematica applies Factor
     * elementwise inside these structures, e.g.,
     *   Factor[1 < 1+2x+x^2+1/(1+x) < 2]
     * yields Less[Less[1, factored], 2].  Without explicit threading
     * here, builtin_factor sees the whole Less[...] as a "polynomial"
     * (with the comparison head treated as opaque), collect_variables
     * returns just [x], and bz_factor_to_expr returns the input
     * unchanged (deg <= 0 for non-polynomial heads).
     *
     * We thread by recursively invoking Factor on each comparison
     * argument, then rebuilding the same head with the factored
     * subterms.  ATTR_LISTABLE on Factor handles the List case
     * automatically; this hand-threading covers the comparison /
     * logic heads it doesn't reach. */
    {
        Expr* arg0 = res->data.function.args[0];
        if (arg0->type == EXPR_FUNCTION
            && arg0->data.function.head
            && arg0->data.function.head->type == EXPR_SYMBOL) {
            const char* h = arg0->data.function.head->data.symbol;
            if (h == SYM_Less ||
                h == SYM_Greater ||
                h == SYM_Equal ||
                h == SYM_Unequal ||
                h == SYM_LessEqual ||
                h == SYM_GreaterEqual ||
                h == SYM_And ||
                h == SYM_Or) {
                size_t n = arg0->data.function.arg_count;
                Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
                for (size_t i = 0; i < n; i++) {
                    Expr* sub_call = expr_new_function(
                        expr_new_symbol("Factor"),
                        (Expr*[]){expr_copy(arg0->data.function.args[i])}, 1);
                    new_args[i] = evaluate(sub_call);
                    expr_free(sub_call);
                }
                Expr* result = expr_new_function(expr_copy(arg0->data.function.head),
                                                 new_args, n);
                free(new_args);
                Expr* evaluated = eval_and_free(result);
                if (memo_key) {
                    factor_memo_put(memo, memo_key, evaluated);
                    expr_free(memo_key);
                }
                return evaluated;
            }
        }
    }

    Expr* arg = res->data.function.args[0];

    /* Algebraic-generator (radical) pass.  When the input contains a
     * sub-expression u with a fractional rational exponent (u^(p/q),
     * q != 1) — e.g. r^(1/5), Log[r]^(2/3), (x+1)^(1/2) — substitute
     * u -> g^m where m is the lcm of all such q's, factor as a polynomial
     * in g, then back-substitute g -> u^(1/m).  The result is fed back
     * through evaluate() so that any residual fractional bases get
     * processed recursively (each level picks a fresh generator name to
     * avoid colliding with previously-introduced ones). */
    {
        Expr* base = NULL;
        Expr* atom = NULL;
        int64_t m = 1;
        if (poly_find_radical_gen(arg, &base, &atom, &m)) {
            char* gen = poly_make_fresh_gen(arg);
            Expr* substituted = poly_subst_radical_to_gen(arg, base, atom, m, gen);
            Expr* call = expr_new_function(expr_new_symbol("Factor"),
                          (Expr*[]){substituted}, 1);
            Expr* result_in_g = evaluate(call);
            expr_free(call);
            Expr* final = poly_subst_radical_from_gen(result_in_g, base, atom, m, gen);
            expr_free(result_in_g);
            expr_free(base);
            expr_free(atom);
            free(gen);
            if (memo_key) {
                factor_memo_put(memo, memo_key, final);
                expr_free(memo_key);
            }
            return final;
        }
    }

    Expr* together = eval_and_free(expr_new_function(expr_new_symbol("Together"), (Expr*[]){expr_copy(arg)}, 1));
    Expr* n_call = expr_new_function(expr_new_symbol("Numerator"), (Expr*[]){expr_copy(together)}, 1);
    Expr* num = evaluate(n_call);
    expr_free(n_call);

    Expr* d_call = expr_new_function(expr_new_symbol("Denominator"), (Expr*[]){expr_copy(together)}, 1);
    Expr* den = evaluate(d_call);
    expr_free(d_call);

    /* Scope discipline: when called from inside Simplify (where a
     * Factor memo is active), use the SHARED-scope behavior -- both
     * numerator and denominator use the variable list collected from
     * the numerator.  This is conservative: a denominator with extra
     * variables won't have those variables in its factoring scope and
     * may end up unfactored.  But it preserves the historical
     * Simplify behavior, which has the property that some downstream
     * TrigRoundtrip paths only converge when intermediate factors
     * remain in their non-fully-factored form.  See FACTOR_PLAN.md
     * cleanup C4.
     *
     * For direct user Factor calls (no memo active), use the SEPARATE-
     * scope behavior: each of num and den is factored with its own
     * variable list.  This gives the correct fully-factored result
     * for inputs like Factor[(x^3+2x^2)/(x^2-4y^2) - (x+2)/(x^2-4y^2)]
     * where the denominator x^2-4y^2 has variables not in the
     * numerator. */
    bool inside_simplify = (factor_memo_top() != NULL);

    size_t num_vc = 0, num_vcap = 16;
    Expr** num_vars = malloc(sizeof(Expr*) * num_vcap);
    collect_variables(num, &num_vars, &num_vc, &num_vcap);

    Expr** den_vars = NULL;
    size_t den_vc = 0;
    if (!inside_simplify) {
        size_t den_vcap = 16;
        den_vars = malloc(sizeof(Expr*) * den_vcap);
        collect_variables(den, &den_vars, &den_vc, &den_vcap);
    }

    Expr* f_num = NULL;
    if (num_vc == 1) {
        f_num = bz_factor_to_expr(num, num_vars[0]);
    } else if (num_vc == 2) {
        /* Fast path for squarefree bivariate inputs: try the
         * bivariate Hensel directly without first running
         * FactorSquareFree.  The Hensel orchestrator skips alphas
         * where the image is not squarefree, so for inputs with
         * repeated factors it will fail and we fall back to the
         * generic FactorSquareFree + heuristic_factor pipeline.
         * For already-squarefree inputs (the common case), this
         * skips ~100+ ms of wasted FactorSquareFree work and the
         * is_likely_irreducible probe. */
        Expr* via_hensel = factor_bivariate_via_hensel(num, num_vars);
        if (via_hensel) {
            f_num = via_hensel;
        } else {
            Expr* sq_num = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(num)}, 1));
            f_num = heuristic_factor(sq_num);
            expr_free(sq_num);
        }
    } else {
        /* For n >= 3, run FactorSquareFree first.  It usually does
         * most of the work for inputs with separable factors (e.g.,
         * (1-x^12)(1+x-y^13)(1-y-z^14) -> Times[squarefree pieces]),
         * after which heuristic_factor's Times-recursion factors
         * each piece independently as a smaller problem. */
        Expr* sq_num = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(num)}, 1));
        f_num = heuristic_factor(sq_num);
        expr_free(sq_num);
    }

    /* Choose den's variable list:
     *  - Simplify-context: reuse num's list (legacy shared scope).
     *  - Direct call: use the den-specific list. */
    Expr** den_eff_vars = inside_simplify ? num_vars : den_vars;
    size_t den_eff_vc   = inside_simplify ? num_vc   : den_vc;

    Expr* f_den = NULL;
    if (den_eff_vc == 0) {
        f_den = expr_copy(den);
    } else if (den_eff_vc == 1) {
        f_den = bz_factor_to_expr(den, den_eff_vars[0]);
    } else {
        Expr* sq_den = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(den)}, 1));
        f_den = heuristic_factor(sq_den);
        expr_free(sq_den);
    }

    for(size_t i=0; i<num_vc; i++) expr_free(num_vars[i]);
    free(num_vars);
    if (den_vars) {
        for(size_t i=0; i<den_vc; i++) expr_free(den_vars[i]);
        free(den_vars);
    }

    Expr* result;
    if (f_den->type == EXPR_INTEGER && f_den->data.integer == 1) {
        result = f_num;
        expr_free(f_den);
    } else {
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){f_den, expr_new_integer(-1)}, 2));
        result = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){f_num, inv_den}, 2));
    }

    expr_free(together);
    expr_free(num);
    expr_free(den);

    /* Cache the (input, output) pair for the duration of the current
     * Simplify call.  memo_key holds a deep copy of the original input
     * we made before computation; we hand it off into the memo here. */
    if (memo_key) {
        factor_memo_put(memo, memo_key, result);
        expr_free(memo_key);
    }

    return result;
}
/* ===================================================================== */
/*  FactorTerms / FactorTermsList                                         */
/*                                                                       */
/*  FactorTerms[poly]                pulls out the overall numerical     */
/*                                   factor.                             */
/*  FactorTerms[poly, x]             pulls out factors that do not       */
/*                                   depend on x.                        */
/*  FactorTerms[poly, {x_1, ..., x_n}]                                   */
/*                                   pulls out factors that do not       */
/*                                   depend on any of the x_i, then      */
/*                                   recursively peels content w.r.t.    */
/*                                   smaller subsets.                    */
/*                                                                       */
/*  FactorTermsList returns the same factors as a list (numerical first, */
/*  then progressively-extracted contents, then the residue).            */
/* ===================================================================== */

static bool ft_is_one_int(Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* Recursive multivariate content. Computes the polynomial GCD of every  */
/* coefficient of a {S_0, ..., S_{S_count-1}}-monomial in `poly`,        */
/* viewed as a polynomial in K[ground_vars][S]. The returned expression  */
/* is a polynomial in K[ground_vars] (in particular it does not contain  */
/* any of the S variables).                                              */
static Expr* ft_content_wrt_set(Expr* poly,
                                Expr** S, size_t S_count,
                                Expr** ground, size_t ground_count) {
    if (is_zero_poly(poly)) return expr_new_integer(0);
    if (S_count == 0) {
        return expr_expand(poly);
    }

    Expr* s_k = S[S_count - 1];
    Expr* expanded = expr_expand(poly);
    int deg = get_degree_poly(expanded, s_k);

    Expr* g = NULL;
    for (int j = 0; j <= deg; j++) {
        Expr* b_j = get_coeff(expanded, s_k, j);
        if (is_zero_poly(b_j)) {
            expr_free(b_j);
            continue;
        }
        Expr* sub = ft_content_wrt_set(b_j, S, S_count - 1, ground, ground_count);
        expr_free(b_j);
        if (!g) {
            g = sub;
        } else {
            Expr* new_g = poly_gcd_internal(g, sub, ground, ground_count);
            expr_free(g);
            expr_free(sub);
            g = new_g;
        }
    }

    expr_free(expanded);
    return g ? g : expr_new_integer(0);
}

/* Divide `poly` by `content` in K[all_vars]. Falls back to symbolic     */
/* Times[poly, content^-1] if exact polynomial division does not apply.  */
static Expr* ft_divide_out(Expr* poly, Expr* content,
                           Expr** all_vars, size_t v_count) {
    if (ft_is_one_int(content)) return expr_expand(poly);
    if (is_zero_poly(content)) return expr_expand(poly);
    Expr* q = exact_poly_div(poly, content, all_vars, v_count);
    if (q) {
        Expr* expanded = expr_expand(q);
        expr_free(q);
        return expanded;
    }
    Expr* inv = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                              (Expr*[]){expr_copy(content), expr_new_integer(-1)}, 2));
    Expr* prod = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                               (Expr*[]){expr_copy(poly), inv}, 2));
    Expr* expanded = expr_expand(prod);
    expr_free(prod);
    return expanded;
}

/* Build the ground variable list = all_vars \ S[0..k-1]. The returned   */
/* array contains aliases into `all_vars` and must be freed (not the     */
/* element pointers).                                                    */
static void ft_compute_ground(Expr** all_vars, size_t v_count,
                              Expr** S, size_t k,
                              Expr*** out_ground, size_t* out_count) {
    *out_ground = (v_count > 0) ? malloc(sizeof(Expr*) * v_count) : NULL;
    *out_count = 0;
    for (size_t i = 0; i < v_count; i++) {
        bool in_S = false;
        for (size_t j = 0; j < k; j++) {
            if (expr_eq(all_vars[i], S[j])) { in_S = true; break; }
        }
        if (!in_S) (*out_ground)[(*out_count)++] = all_vars[i];
    }
}

/* True if `head` is one of the heads over which FactorTerms threads.    */
static bool ft_is_threading_head(const char* head) {
    return strcmp(head, "List") == 0 ||
           strcmp(head, "Equal") == 0 || strcmp(head, "Unequal") == 0 ||
           strcmp(head, "Less") == 0 || strcmp(head, "LessEqual") == 0 ||
           strcmp(head, "Greater") == 0 || strcmp(head, "GreaterEqual") == 0 ||
           strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 ||
           strcmp(head, "Not") == 0 || strcmp(head, "Xor") == 0;
}

/* Core engine. Returns the list of factors as List[c_0, c_1, ..., r],   */
/* where c_0 is the numerical content, c_1..c_n are progressively        */
/* extracted contents, and r is the residue. The result is a freshly     */
/* allocated Expr* that the caller takes ownership of.                  */
static Expr* ft_compute_list(Expr* poly, Expr** S, size_t S_count) {
    /* Together-normalise so rational expressions become a single        */
    /* numerator/denominator pair we can run polynomial machinery on.    */
    Expr* together = internal_together((Expr*[]){expr_copy(poly)}, 1);
    Expr* num = internal_numerator((Expr*[]){expr_copy(together)}, 1);
    Expr* den = internal_denominator((Expr*[]){expr_copy(together)}, 1);
    expr_free(together);

    Expr* cur = expr_expand(num);
    expr_free(num);

    /* Variables of the numerator, sorted canonically. */
    size_t v_count = 0, v_cap = 16;
    Expr** all_vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(cur, &all_vars, &v_count, &v_cap);
    if (v_count > 0) qsort(all_vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    /* Filter S to the variables that actually appear in the numerator;  */
    /* requested variables that do not appear contribute nothing and     */
    /* would otherwise produce trivial 1-factors in the output list.     */
    Expr** S_eff = (S_count > 0) ? malloc(sizeof(Expr*) * S_count) : NULL;
    size_t S_eff_count = 0;
    for (size_t i = 0; i < S_count; i++) {
        for (size_t j = 0; j < v_count; j++) {
            if (expr_eq(S[i], all_vars[j])) { S_eff[S_eff_count++] = S[i]; break; }
        }
    }

    size_t out_cap = S_eff_count + 4;
    Expr** out = malloc(sizeof(Expr*) * out_cap);
    size_t out_count = 0;

    /* Step 1: numerical content -- content w.r.t. ALL variables, with   */
    /* an empty ground ring (the gcd lives in Z, modulo the limits of    */
    /* my_number_gcd in poly.c).                                         */
    Expr* num_cont = ft_content_wrt_set(cur, all_vars, v_count, NULL, 0);
    if (!ft_is_one_int(num_cont) && !is_zero_poly(num_cont)) {
        Expr* new_cur = ft_divide_out(cur, num_cont, all_vars, v_count);
        expr_free(cur);
        cur = new_cur;
    }
    out[out_count++] = num_cont;

    /* Steps 2..(S_eff_count+1): for k = S_eff_count down to 1, peel the */
    /* content w.r.t. {S_eff[0], ..., S_eff[k-1]} from the running       */
    /* residue. Each ground ring shrinks by exactly one variable.        */
    for (size_t k = S_eff_count; k > 0; k--) {
        Expr** ground = NULL;
        size_t g_count = 0;
        ft_compute_ground(all_vars, v_count, S_eff, k, &ground, &g_count);

        Expr* cont_k = ft_content_wrt_set(cur, S_eff, k, ground, g_count);
        if (!ft_is_one_int(cont_k) && !is_zero_poly(cont_k)) {
            Expr* new_cur = ft_divide_out(cur, cont_k, all_vars, v_count);
            expr_free(cur);
            cur = new_cur;
        }
        free(ground);
        out[out_count++] = cont_k;
    }

    /* Final residue. Multiply through by 1/den so rational inputs round */
    /* trip back to their original form.                                 */
    Expr* final;
    if (ft_is_one_int(den)) {
        final = cur;
        expr_free(den);
    } else {
        Expr* inv_den = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                                      (Expr*[]){den, expr_new_integer(-1)}, 2));
        final = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                              (Expr*[]){cur, inv_den}, 2));
    }
    out[out_count++] = final;

    for (size_t i = 0; i < v_count; i++) expr_free(all_vars[i]);
    free(all_vars);
    free(S_eff);

    Expr* list = expr_new_function(expr_new_symbol("List"), out, out_count);
    free(out);
    return list;
}

/* Pull the variables-argument out of a 2-arg FactorTerms / FactorTermsList */
/* call. *S_out is set to a freshly allocated array of aliases (must be     */
/* freed by the caller; do not free the elements). var_arg is the second    */
/* argument as it appears in the call expression -- either a List or a      */
/* single bare symbol/expression.                                           */
static void ft_extract_vars(Expr* var_arg, Expr*** S_out, size_t* S_count_out) {
    if (var_arg->type == EXPR_FUNCTION &&
        var_arg->data.function.head->type == EXPR_SYMBOL &&
        var_arg->data.function.head->data.symbol == SYM_List) {
        size_t n = var_arg->data.function.arg_count;
        *S_out = (n > 0) ? malloc(sizeof(Expr*) * n) : NULL;
        *S_count_out = n;
        for (size_t i = 0; i < n; i++) (*S_out)[i] = var_arg->data.function.args[i];
    } else {
        *S_out = malloc(sizeof(Expr*) * 1);
        (*S_out)[0] = var_arg;
        *S_count_out = 1;
    }
}

/* Thread FactorTerms / FactorTermsList over equation, inequality, and  */
/* logic heads (List, Equal, Less, And, ...). `fname` chooses which of  */
/* the two functions to recurse with.                                   */
static Expr* ft_thread(Expr* compound, Expr* var_arg, const char* fname) {
    size_t n = compound->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t call_argc = (var_arg ? 2 : 1);
        Expr** call_args = malloc(sizeof(Expr*) * call_argc);
        call_args[0] = expr_copy(compound->data.function.args[i]);
        if (var_arg) call_args[1] = expr_copy(var_arg);
        Expr* call = expr_new_function(expr_new_symbol(fname), call_args, call_argc);
        free(call_args);
        new_args[i] = eval_and_free(call);
    }
    Expr* threaded = expr_new_function(expr_copy(compound->data.function.head), new_args, n);
    free(new_args);
    return eval_and_free(threaded);
}

Expr* builtin_factortermslist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* poly = res->data.function.args[0];
    Expr* var_arg = (argc == 2) ? res->data.function.args[1] : NULL;

    /* FactorTermsList does not auto-thread (its result shape is already */
    /* a List, so threading would be ambiguous).                         */

    Expr** S = NULL;
    size_t S_count = 0;
    if (var_arg) ft_extract_vars(var_arg, &S, &S_count);

    Expr* result = ft_compute_list(poly, S, S_count);
    free(S);
    return result;
}

Expr* builtin_factorterms(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* poly = res->data.function.args[0];
    Expr* var_arg = (argc == 2) ? res->data.function.args[1] : NULL;

    /* Auto-thread over List, Equal, Less, And, Or, ... */
    if (poly->type == EXPR_FUNCTION && poly->data.function.head->type == EXPR_SYMBOL) {
        const char* h = poly->data.function.head->data.symbol;
        if (ft_is_threading_head(h)) {
            return ft_thread(poly, var_arg, "FactorTerms");
        }
    }

    Expr** S = NULL;
    size_t S_count = 0;
    if (var_arg) ft_extract_vars(var_arg, &S, &S_count);

    Expr* list = ft_compute_list(poly, S, S_count);
    free(S);

    /* Multiply the factor list together to obtain the FactorTerms       */
    /* result -- expressed via Times (which canonicalises on evaluation).*/
    size_t n = list->data.function.arg_count;
    Expr** factors = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) factors[i] = expr_copy(list->data.function.args[i]);
    expr_free(list);
    Expr* product = eval_and_free(expr_new_function(expr_new_symbol("Times"), factors, n));
    free(factors);
    return product;
}

void facpoly_init(void) {
    symtab_add_builtin("FactorSquareFree", builtin_factorsquarefree);
    symtab_get_def("FactorSquareFree")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("Factor", builtin_factor);
    symtab_get_def("Factor")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("FactorTerms", builtin_factorterms);
    symtab_get_def("FactorTerms")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FactorTermsList", builtin_factortermslist);
    symtab_get_def("FactorTermsList")->attributes |= ATTR_PROTECTED;
}

typedef __int128_t int128_t;

typedef struct {
    int deg;
    int64_t* c;
} UPoly;

static UPoly* upoly_new(int deg) {
    if (deg < 0) return NULL;
    if (deg > 10000) return NULL; // Safety bound
    UPoly* p = malloc(sizeof(UPoly));
    p->deg = deg;
    p->c = calloc(deg + 1, sizeof(int64_t));
    return p;
}

static void upoly_free(UPoly* p) {
    if (p) {
        free(p->c);
        free(p);
    }
}

static UPoly* upoly_copy(UPoly* a) {
    UPoly* p = upoly_new(a->deg);
    for (int i = 0; i <= a->deg; i++) p->c[i] = a->c[i];
    return p;
}

static void upoly_trim(UPoly* p) {
    while (p->deg > 0 && p->c[p->deg] == 0) p->deg--;
}

static int64_t mod_inverse_int(int64_t a, int64_t m) {
    int64_t m0 = m, t, q;
    int64_t x0 = 0, x1 = 1;
    if (m == 1) return 0;
    a %= m; if (a < 0) a += m;
    while (a > 1) {
        q = a / m;
        t = m; m = a % m; a = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}

static void upoly_div_rem_mod(UPoly* a, UPoly* b, int64_t mod, UPoly** out_q, UPoly** out_r) {
    if (b->deg < 0 || (b->deg == 0 && b->c[0] == 0)) {
        if (out_q) *out_q = NULL;
        if (out_r) *out_r = NULL;
        return;
    }
    UPoly* q = upoly_new(a->deg >= b->deg ? a->deg - b->deg : 0);
    UPoly* r = upoly_copy(a);
    int64_t inv = mod_inverse_int(b->c[b->deg], mod);
    while (r->deg >= b->deg && !(r->deg == 0 && r->c[0] == 0)) {
        int d = r->deg - b->deg;
        int64_t coeff = (int64_t)(((int128_t)r->c[r->deg] * inv) % mod);
        if (coeff < 0) coeff += mod;
        q->c[d] = coeff;
        for (int i = 0; i <= b->deg; i++) {
            r->c[i + d] = (int64_t)((r->c[i + d] - (int128_t)coeff * b->c[i]) % mod);
            if (r->c[i + d] < 0) r->c[i + d] += mod;
        }
        upoly_trim(r);
    }
    upoly_trim(q);
    if (out_q) *out_q = q; else upoly_free(q);
    if (out_r) *out_r = r; else upoly_free(r);
}

static UPoly* upoly_add_mod(UPoly* a, UPoly* b, int64_t mod) {
    int deg = a->deg > b->deg ? a->deg : b->deg;
    UPoly* r = upoly_new(deg);
    for (int i = 0; i <= deg; i++) {
        int64_t ca = i <= a->deg ? a->c[i] : 0;
        int64_t cb = i <= b->deg ? b->c[i] : 0;
        r->c[i] = (ca + cb) % mod;
        if (r->c[i] < 0) r->c[i] += mod;
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_sub_mod(UPoly* a, UPoly* b, int64_t mod) {
    int deg = a->deg > b->deg ? a->deg : b->deg;
    UPoly* r = upoly_new(deg);
    for (int i = 0; i <= deg; i++) {
        int64_t ca = i <= a->deg ? a->c[i] : 0;
        int64_t cb = i <= b->deg ? b->c[i] : 0;
        r->c[i] = (ca - cb) % mod;
        if (r->c[i] < 0) r->c[i] += mod;
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_mul_mod(UPoly* a, UPoly* b, int64_t mod) {
    UPoly* r = upoly_new(a->deg + b->deg);
    for (int i = 0; i <= a->deg; i++) {
        for (int j = 0; j <= b->deg; j++) {
            r->c[i + j] = (int64_t)((r->c[i + j] + (int128_t)a->c[i] * b->c[j]) % mod);
            if (r->c[i+j] < 0) r->c[i+j] += mod;
        }
    }
    upoly_trim(r);
    return r;
}

static UPoly* upoly_gcd_mod(UPoly* a, UPoly* b, int64_t mod) {
    UPoly* x = upoly_copy(a);
    UPoly* y = upoly_copy(b);
    while (y->deg > 0 || y->c[0] != 0) {
        UPoly* rem;
        upoly_div_rem_mod(x, y, mod, NULL, &rem);
        upoly_free(x);
        x = y;
        y = rem;
    }
    upoly_free(y);
    if (x->deg >= 0 && x->c[x->deg] != 0) {
        int64_t inv = mod_inverse_int(x->c[x->deg], mod);
        for (int i = 0; i <= x->deg; i++) {
            x->c[i] = (int64_t)(((int128_t)x->c[i] * inv) % mod);
        }
    }
    return x;
}

static UPoly* upoly_pow_mod_poly(UPoly* a, int64_t exp, UPoly* p, int64_t mod) {
    UPoly* res = upoly_new(0); res->c[0] = 1;
    UPoly* base = upoly_copy(a);
    while (exp > 0) {
        if (exp % 2 == 1) {
            UPoly* t = upoly_mul_mod(res, base, mod);
            upoly_free(res);
            upoly_div_rem_mod(t, p, mod, NULL, &res);
            upoly_free(t);
        }
        UPoly* t2 = upoly_mul_mod(base, base, mod);
        upoly_free(base);
        upoly_div_rem_mod(t2, p, mod, NULL, &base);
        upoly_free(t2);
        exp /= 2;
    }
    upoly_free(base);
    return res;
}

typedef struct {
    UPoly** factors;
    int count;
} UPolyList;

static void upoly_list_add(UPolyList* list, UPoly* p) {
    list->factors = realloc(list->factors, sizeof(UPoly*) * (list->count + 1));
    list->factors[list->count++] = p;
}

static UPolyList cz_edf(UPoly* f, int d, int64_t p) {
    UPolyList list = {NULL, 0};
    if (f->deg == d) {
        upoly_list_add(&list, upoly_copy(f));
        return list;
    }
    
    UPoly* x = upoly_new(1); x->c[1] = 1;
    while (list.count < f->deg / d) {
        UPoly* t = upoly_new(f->deg - 1);
        for(int i=0; i<f->deg; i++) t->c[i] = rand() % p;
        upoly_trim(t);
        
        int64_t exponent = 1;
        for (int i=0; i<d; i++) exponent *= p;
        exponent = (exponent - 1) / 2;
        
        UPoly* t_pow = upoly_pow_mod_poly(t, exponent, f, p);
        UPoly* zero_poly = upoly_new(0);
        UPoly* t_pow_minus_1 = upoly_sub_mod(t_pow, zero_poly, p);
        upoly_free(zero_poly);
        t_pow_minus_1->c[0] = (t_pow->c[0] - 1 + p) % p;
        
        UPoly* g = upoly_gcd_mod(f, t_pow_minus_1, p);
        if (g->deg > 0 && g->deg < f->deg && g->deg % d == 0) {
            UPoly* f_over_g;
            upoly_div_rem_mod(f, g, p, &f_over_g, NULL);
            UPolyList L1 = cz_edf(g, d, p);
            UPolyList L2 = cz_edf(f_over_g, d, p);
            for(int i=0; i<L1.count; i++) upoly_list_add(&list, L1.factors[i]);
            for(int i=0; i<L2.count; i++) upoly_list_add(&list, L2.factors[i]);
            free(L1.factors); free(L2.factors);
            upoly_free(g); upoly_free(f_over_g); upoly_free(t); upoly_free(t_pow); upoly_free(t_pow_minus_1);
            break;
        }
        upoly_free(g); upoly_free(t); upoly_free(t_pow); upoly_free(t_pow_minus_1);
    }
    upoly_free(x);
    return list;
}

static UPolyList cz_ddf(UPoly* f, int64_t p) {
    UPolyList res = {NULL, 0};
    int i = 1;
    UPoly* f_star = upoly_copy(f);
    UPoly* x = upoly_new(1); x->c[1] = 1;
    UPoly* x_pow_p = upoly_pow_mod_poly(x, p, f_star, p);
    
    while (f_star->deg >= 2 * i) {
        UPoly* t = upoly_sub_mod(x_pow_p, x, p);
        UPoly* g = upoly_gcd_mod(f_star, t, p);
        if (g->deg > 0) {
            UPolyList edf_factors = cz_edf(g, i, p);
            for (int k = 0; k < edf_factors.count; k++) {
                upoly_list_add(&res, edf_factors.factors[k]);
            }
            free(edf_factors.factors);
            
            UPoly* next_f;
            upoly_div_rem_mod(f_star, g, p, &next_f, NULL);
            upoly_free(f_star);
            f_star = next_f;
            
            UPoly* next_x_pow = upoly_copy(x_pow_p);
            upoly_div_rem_mod(next_x_pow, f_star, p, NULL, &x_pow_p);
            upoly_free(next_x_pow);
        }
        upoly_free(t);
        upoly_free(g);
        
        UPoly* temp = upoly_pow_mod_poly(x_pow_p, p, f_star, p);
        upoly_free(x_pow_p);
        x_pow_p = temp;
        i++;
    }
    if (f_star->deg > 0) {
        upoly_list_add(&res, upoly_copy(f_star));
    }
    upoly_free(f_star);
    upoly_free(x);
    upoly_free(x_pow_p);
    return res;
}

static void upoly_xgcd_mod(UPoly* a, UPoly* b, int64_t mod, UPoly** out_g, UPoly** out_s, UPoly** out_t) {
    UPoly *s0 = upoly_new(0); s0->c[0] = 1;
    UPoly *t0 = upoly_new(0);
    UPoly *s1 = upoly_new(0);
    UPoly *t1 = upoly_new(0); t1->c[0] = 1;
    UPoly *r0 = upoly_copy(a);
    UPoly *r1 = upoly_copy(b);
    
    while (r1->deg > 0 || r1->c[0] != 0) {
        UPoly *q, *r2;
        upoly_div_rem_mod(r0, r1, mod, &q, &r2);
        
        UPoly *q_s1 = upoly_mul_mod(q, s1, mod);
        UPoly *s2 = upoly_sub_mod(s0, q_s1, mod);
        
        UPoly *q_t1 = upoly_mul_mod(q, t1, mod);
        UPoly *t2 = upoly_sub_mod(t0, q_t1, mod);
        
        upoly_free(r0); r0 = r1; r1 = r2;
        upoly_free(s0); s0 = s1; s1 = s2;
        upoly_free(t0); t0 = t1; t1 = t2;
        upoly_free(q); upoly_free(q_s1); upoly_free(q_t1);
    }
    
    int64_t lc = r0->c[r0->deg];
    int64_t inv = mod_inverse_int(lc, mod);
    for(int i=0; i<=r0->deg; i++) r0->c[i] = (r0->c[i] * inv) % mod;
    for(int i=0; i<=s0->deg; i++) s0->c[i] = (s0->c[i] * inv) % mod;
    for(int i=0; i<=t0->deg; i++) t0->c[i] = (t0->c[i] * inv) % mod;
    
    *out_g = r0; *out_s = s0; *out_t = t0;
    upoly_free(r1); upoly_free(s1); upoly_free(t1);
}

static void hensel_lift(UPoly* P, UPoly* A_p, UPoly* B_p, int64_t p, int64_t k, UPoly** out_A, UPoly** out_B) {
    int64_t pk = p;
    UPoly* A = upoly_copy(A_p);
    UPoly* B = upoly_copy(B_p);
    
    UPoly *S, *T, *G;
    upoly_xgcd_mod(A, B, p, &G, &S, &T);
    upoly_free(G);
    
    while (pk < k) {
        UPoly* AB = upoly_mul_mod(A, B, pk * p);
        UPoly* E = upoly_sub_mod(P, AB, pk * p);
        for(int i=0; i<=E->deg; i++) {
            if (E->c[i] < 0) E->c[i] += pk * p;
            E->c[i] /= pk;
        }
        upoly_free(AB);
        
        UPoly *SE = upoly_mul_mod(S, E, p);
        UPoly *TE = upoly_mul_mod(T, E, p);
        
        UPoly *Q, *R;
        upoly_div_rem_mod(TE, A, p, &Q, &R);
        upoly_free(TE);
        
        UPoly *QB = upoly_mul_mod(Q, B, p);
        UPoly *B_next_diff = upoly_add_mod(SE, QB, p);
        upoly_free(QB); upoly_free(SE); upoly_free(Q);
        
        UPoly* A_next_diff = R; 
        
        UPoly* next_A = upoly_new(A->deg >= A_next_diff->deg ? A->deg : A_next_diff->deg);
        UPoly* next_B = upoly_new(B->deg >= B_next_diff->deg ? B->deg : B_next_diff->deg);
        
        for(int i=0; i<=next_A->deg; i++) {
            int64_t ca = i <= A->deg ? A->c[i] : 0;
            int64_t cn = i <= A_next_diff->deg ? A_next_diff->c[i] : 0;
            next_A->c[i] = (ca + cn * pk) % (pk * p);
        }
        for(int i=0; i<=next_B->deg; i++) {
            int64_t cb = i <= B->deg ? B->c[i] : 0;
            int64_t cn = i <= B_next_diff->deg ? B_next_diff->c[i] : 0;
            next_B->c[i] = (cb + cn * pk) % (pk * p);
        }
        upoly_trim(next_A); upoly_trim(next_B);
        
        upoly_free(A); upoly_free(B);
        A = next_A; B = next_B;
        
        pk *= p;
        upoly_free(E); upoly_free(A_next_diff); upoly_free(B_next_diff);
    }
    
    upoly_free(S); upoly_free(T);
    *out_A = A; *out_B = B;
}

static void multifactor_hensel_lift(UPoly* P, UPolyList factors_p, int64_t p, int64_t k, UPolyList* out_factors) {
    if (factors_p.count == 1) {
        UPoly* f = upoly_copy(P);
        upoly_list_add(out_factors, f);
        return;
    }
    
    UPoly* A_p = factors_p.factors[0];
    UPoly* B_p = upoly_new(0); B_p->c[0] = 1;
    for (int i = 1; i < factors_p.count; i++) {
        UPoly* tmp = upoly_mul_mod(B_p, factors_p.factors[i], p);
        upoly_free(B_p);
        B_p = tmp;
    }
    
    int64_t lc = P->c[P->deg];
    for (int i = 0; i <= B_p->deg; i++) {
        B_p->c[i] = (B_p->c[i] * lc) % p;
    }
    
    UPoly *A_k, *B_k;
    hensel_lift(P, A_p, B_p, p, k, &A_k, &B_k);
    upoly_free(B_p);
    
    upoly_list_add(out_factors, A_k);
    
    UPolyList rest_p = {NULL, 0};
    for (int i = 1; i < factors_p.count; i++) {
        upoly_list_add(&rest_p, factors_p.factors[i]);
    }
    
    multifactor_hensel_lift(B_k, rest_p, p, k, out_factors);
    free(rest_p.factors);
    upoly_free(B_k);
}

static int64_t z_gcd(int64_t a, int64_t b) {
    a = llabs(a); b = llabs(b);
    while (b) { int64_t t = b; b = a % b; a = t; }
    return a;
}

static UPoly* upoly_div_exact_z(UPoly* a, UPoly* b) {
    UPoly* q = upoly_new(a->deg >= b->deg ? a->deg - b->deg : 0);
    UPoly* r = upoly_copy(a);
    while (r->deg >= b->deg && !(r->deg == 0 && r->c[0] == 0)) {
        int d = r->deg - b->deg;
        if (r->c[r->deg] % b->c[b->deg] != 0) {
            upoly_free(q); upoly_free(r); return NULL;
        }
        int64_t coeff = r->c[r->deg] / b->c[b->deg];
        q->c[d] = coeff;
        for (int i = 0; i <= b->deg; i++) {
            r->c[i + d] -= coeff * b->c[i];
        }
        upoly_trim(r);
    }
    if (r->deg > 0 || r->c[0] != 0) {
        upoly_free(q); upoly_free(r); return NULL;
    }
    upoly_free(r);
    return q;
}

static bool is_prime(int64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (int64_t i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

static UPolyList factor_zassenhaus(UPoly* P) {
    UPolyList result = {NULL, 0};
    if (P->deg == 0) {
        upoly_list_add(&result, upoly_copy(P));
        return result;
    }
    
    int64_t content = llabs(P->c[0]);
    for(int i=1; i<=P->deg; i++) content = z_gcd(content, P->c[i]);

    UPoly* P_prim = upoly_copy(P);
    if (content > 1) {
        for(int i=0; i<=P_prim->deg; i++) P_prim->c[i] /= content;
    }

    // Ensure P_prim has a positive leading coefficient
    if (P_prim->c[P_prim->deg] < 0) {
        for(int i=0; i<=P_prim->deg; i++) P_prim->c[i] = -P_prim->c[i];
        content = -content;
    }

    int64_t lc = P_prim->c[P_prim->deg];    
    int64_t p = 13; // Just pick 13 for simplicity in this implementation
    int attempts = 0;
    while (attempts < 20) {
        if (!is_prime(p)) {
            p += 2;
            continue;
        }
        if (lc % p != 0) {
            UPoly* p_mod = upoly_copy(P_prim);
            for(int i=0; i<=P_prim->deg; i++) {
                p_mod->c[i] %= p;
                if(p_mod->c[i] < 0) p_mod->c[i] += p;
            }
            // check squarefree
            UPoly* deriv = upoly_new(p_mod->deg - 1);
            for (int i = 1; i <= p_mod->deg; i++) deriv->c[i - 1] = (p_mod->c[i] * i) % p;
            // When p divides the polynomial degree (e.g. P = x^p - c), every
            // derivative coefficient is 0 mod p.  Without this trim the
            // divisor passed to upoly_gcd_mod has positive nominal degree
            // but all-zero coefficients, and upoly_div_rem_mod loops forever.
            upoly_trim(deriv);
            UPoly* g = upoly_gcd_mod(p_mod, deriv, p);
            bool sqf = (g->deg == 0 && g->c[0] != 0);
            upoly_free(deriv); upoly_free(g); upoly_free(p_mod);
            if (sqf) break;
        }
        p += 2; // next prime roughly
        attempts++;
    }
    
    if (attempts >= 20) {
        // Fallback: likely not square-free or algorithm failed. Return empty indicating error
        upoly_free(P_prim);
        UPolyList empty = {NULL, 0};
        return empty;
    }
    
    UPoly* p_mod = upoly_copy(P_prim);
    for(int i=0; i<=P_prim->deg; i++) {
        p_mod->c[i] %= p;
        if(p_mod->c[i] < 0) p_mod->c[i] += p;
    }
    
    int64_t inv_lc = mod_inverse_int(p_mod->c[p_mod->deg], p);
    UPoly* monic = upoly_copy(p_mod);
    for(int i=0; i<=monic->deg; i++) monic->c[i] = (monic->c[i] * inv_lc) % p;
    
    UPolyList factors_p = cz_ddf(monic, p);
    upoly_free(monic);
    upoly_free(p_mod);
    
    int64_t pk = p;
    while (pk < 1000000000000000LL) pk *= p; 
    
    UPolyList factors_pk = {NULL, 0};
    multifactor_hensel_lift(P_prim, factors_p, p, pk, &factors_pk);
    for(int i=0; i<factors_p.count; i++) upoly_free(factors_p.factors[i]);
    free(factors_p.factors);
    
    int* used = calloc(factors_pk.count, sizeof(int));
    UPoly* current_P = upoly_copy(P_prim);
    
    int s = 1;
    while (2 * s <= factors_pk.count) {
        bool combined = false;
        int total_subsets = 1 << factors_pk.count;
        for (int mask = 1; mask < total_subsets; mask++) {
            int bits = 0;
            for (int i=0; i<factors_pk.count; i++) if ((mask >> i) & 1) bits++;
            if (bits != s) continue;
            
            bool valid = true;
            for (int i=0; i<factors_pk.count; i++) if (((mask >> i) & 1) && used[i]) valid = false;
            if (!valid) continue;
            
            UPoly* test_A = upoly_new(0); test_A->c[0] = 1;
            for (int i=0; i<factors_pk.count; i++) {
                if ((mask >> i) & 1) {
                    UPoly* t = upoly_mul_mod(test_A, factors_pk.factors[i], pk);
                    upoly_free(test_A);
                    test_A = t;
                }
            }
            
            for(int j=0; j<=test_A->deg; j++) {
                test_A->c[j] = (test_A->c[j] * lc) % pk;
                if (test_A->c[j] > pk/2) test_A->c[j] -= pk;
            }
            
            int64_t c = test_A->c[0];
            for(int j=1; j<=test_A->deg; j++) c = z_gcd(c, test_A->c[j]);
            if (c > 0) {
                for(int j=0; j<=test_A->deg; j++) test_A->c[j] /= c;
            }
            
            UPoly* Q = upoly_div_exact_z(current_P, test_A);
            if (Q != NULL) {
                upoly_list_add(&result, test_A);
                for (int i=0; i<factors_pk.count; i++) if ((mask >> i) & 1) used[i] = 1;
                upoly_free(current_P);
                current_P = Q;
                combined = true;
                lc = current_P->c[current_P->deg];
                break;
            }
            upoly_free(test_A);
        }
        if (!combined) s++;
    }
    
    if (current_P->deg > 0) {
        upoly_list_add(&result, current_P);
    } else {
        upoly_free(current_P);
    }
    
    if (content != 1) {
        UPoly* c_poly = upoly_new(0);
        c_poly->c[0] = content;
        upoly_list_add(&result, c_poly);
    }
    
    free(used);
    for(int i=0; i<factors_pk.count; i++) upoly_free(factors_pk.factors[i]);
    free(factors_pk.factors);
    upoly_free(P_prim);
    
    return result;
}

Expr* bz_factor_to_expr(Expr* P, Expr* var) {
    if (P->type != EXPR_FUNCTION) return expr_copy(P);

    int deg = get_degree_poly(P, var);
    if (deg <= 0) return expr_copy(P);

    UPoly* up = upoly_new(deg);
    int non_integer_seen = 0;
    for (int i = 0; i <= deg; i++) {
        Expr* c = get_coeff(P, var, i);
        if (c->type == EXPR_INTEGER) {
            up->c[i] = c->data.integer;
        } else {
            /* Non-integer coefficient (e.g. Gaussian rational like 2 I, or
             * any other symbolic / rational / complex factor). The integer
             * Zassenhaus factoring routine below operates only on Z[var];
             * silently coercing the coefficient to 0 produces a structurally
             * different polynomial and yields a wrong factorization. Bail
             * out and leave P unfactored at this layer instead. */
            non_integer_seen = 1;
            up->c[i] = 0;
        }
        expr_free(c);
    }
    if (non_integer_seen) {
        upoly_free(up);
        return expr_copy(P);
    }
    
    UPolyList factors = factor_zassenhaus(up);
    upoly_free(up);
    
    if (factors.count == 0) {
        // Fallback to FactorSquareFree if not square-free or algorithm failed
        Expr* sq = eval_and_free(expr_new_function(expr_new_symbol("FactorSquareFree"), (Expr*[]){expr_copy(P)}, 1));
        Expr* res = heuristic_factor(sq);
        expr_free(sq);
        return res;
    }
    
    if (factors.count == 1 && factors.factors[0]->deg == deg) {
        upoly_free(factors.factors[0]);
        free(factors.factors);
        return expr_copy(P);
    }
    
    Expr** args = malloc(sizeof(Expr*) * factors.count);
    for (int i = 0; i < factors.count; i++) {
        UPoly* f = factors.factors[i];
        Expr* sum = expr_new_integer(0);
        for (int j = 0; j <= f->deg; j++) {
            if (f->c[j] == 0) continue;
            Expr* term;
            if (j == 0) term = expr_new_integer(f->c[j]);
            else if (j == 1) {
                if (f->c[j] == 1) term = expr_copy(var);
                else if (f->c[j] == -1) term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(var)}, 2));
                else term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(f->c[j]), expr_copy(var)}, 2));
            } else {
                Expr* xp = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(var), expr_new_integer(j)}, 2));
                if (f->c[j] == 1) term = xp;
                else if (f->c[j] == -1) term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), xp}, 2));
                else term = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(f->c[j]), xp}, 2));
            }
            sum = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){sum, term}, 2));
        }
        args[i] = sum;
        upoly_free(f);
    }
    free(factors.factors);
    
    Expr* res = eval_and_free(expr_new_function(expr_new_symbol("Times"), args, factors.count));
    free(args);
    return res;
}
