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
