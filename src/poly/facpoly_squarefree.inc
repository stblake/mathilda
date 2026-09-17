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
        const char* head = p->data.function.head->type == EXPR_SYMBOL ? p->data.function.head->data.symbol.name : "";
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

    const char* head = e->data.function.head->type == EXPR_SYMBOL ? e->data.function.head->data.symbol.name : "";

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

/* True for a pure numeric literal -- it belongs to the leading content of the
 * square-free product rather than to a factor. */
static bool sqf_is_number_literal(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:
        case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default:
            break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (strcmp(h, "Rational") == 0 || strcmp(h, "Complex") == 0) return true;
    }
    return false;
}

/* FactorSquareFree over an algebraic extension Q(alpha).
 *
 * Square-free factorisation is field-independent for rational-coefficient
 * inputs, but the user may supply coefficients in Q(alpha) or ask for the
 * decomposition to be performed there.  We obtain the irreducible-over-Q(alpha)
 * factorisation from Factor (the Trager engine in qafactor.c, reached through
 * builtin_factor's Extension option -- every argument, poly and options, is
 * forwarded verbatim), then reconstruct the square-free decomposition by
 * grouping the irreducible factors by multiplicity: for each distinct
 * multiplicity k, the square-free factor of multiplicity k is Expand[product of
 * all irreducibles that occur to the k-th power].  This needs no square-free-
 * over-Q(alpha) primitive of its own.  Returns NULL when Factor declines, so
 * the caller can fall back to the over-Q pipeline. */
static Expr* factor_square_free_over_extension(Expr* res) {
    size_t argc = res->data.function.arg_count;

    /* Factor[poly, opts...] -- forward every argument verbatim. */
    Expr** fargs = malloc(sizeof(Expr*) * argc);
    for (size_t i = 0; i < argc; i++) fargs[i] = expr_copy(res->data.function.args[i]);
    Expr* fac = eval_and_free(expr_new_function(expr_new_symbol("Factor"), fargs, argc));
    free(fargs);
    if (!fac) return NULL;

    /* Collect a numeric content plus (base, multiplicity) pairs. */
    Expr* content = expr_new_integer(1);
    size_t cap = 8, n = 0;
    Expr** bases = malloc(sizeof(Expr*) * cap);
    long*  mults = malloc(sizeof(long) * cap);

    bool fac_is_times = (fac->type == EXPR_FUNCTION && fac->data.function.head &&
                         fac->data.function.head->type == EXPR_SYMBOL &&
                         strcmp(fac->data.function.head->data.symbol.name, "Times") == 0);
    size_t nfac = fac_is_times ? fac->data.function.arg_count : 1;
    Expr** faclist = fac_is_times ? fac->data.function.args : &fac;

    for (size_t i = 0; i < nfac; i++) {
        Expr* f = faclist[i];
        if (sqf_is_number_literal(f)) {
            content = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ content, expr_copy(f) }, 2));
            continue;
        }
        Expr* base;
        long mult;
        if (f->type == EXPR_FUNCTION && f->data.function.head &&
            f->data.function.head->type == EXPR_SYMBOL &&
            strcmp(f->data.function.head->data.symbol.name, "Power") == 0 &&
            f->data.function.arg_count == 2 &&
            f->data.function.args[1]->type == EXPR_INTEGER) {
            base = expr_copy(f->data.function.args[0]);
            mult = (long)f->data.function.args[1]->data.integer;
        } else {
            base = expr_copy(f);
            mult = 1;
        }
        if (n == cap) {
            cap *= 2;
            bases = realloc(bases, sizeof(Expr*) * cap);
            mults = realloc(mults, sizeof(long) * cap);
        }
        bases[n] = base;
        mults[n] = mult;
        n++;
    }
    expr_free(fac);

    /* Group by multiplicity and assemble Times[content, Power[sqfree_k, k]...].
     * sqfree_k = Expand[product of bases sharing multiplicity k], matching the
     * expanded normal form of the over-Q FactorSquareFree path. */
    size_t term_cap = n + 1, term_n = 0;
    Expr** terms = malloc(sizeof(Expr*) * term_cap);
    terms[term_n++] = content;

    bool* used = calloc(n ? n : 1, sizeof(bool));
    for (size_t i = 0; i < n; i++) {
        if (used[i]) continue;
        long k = mults[i];
        size_t pc = 0, pcap = 4;
        Expr** pbases = malloc(sizeof(Expr*) * pcap);
        for (size_t j = i; j < n; j++) {
            if (used[j] || mults[j] != k) continue;
            used[j] = true;
            if (pc == pcap) { pcap *= 2; pbases = realloc(pbases, sizeof(Expr*) * pcap); }
            pbases[pc++] = expr_copy(bases[j]);
        }
        Expr* prod;
        if (pc == 1) {
            prod = pbases[0];
        } else {
            prod = eval_and_free(expr_new_function(expr_new_symbol("Times"), pbases, pc));
        }
        free(pbases);
        Expr* sqfree = expr_expand(prod);
        expr_free(prod);

        Expr* term;
        if (k == 1) {
            term = sqfree;
        } else {
            term = eval_and_free(expr_new_function(expr_new_symbol("Power"),
                        (Expr*[]){ sqfree, expr_new_integer(k) }, 2));
        }
        if (term_n == term_cap) { term_cap *= 2; terms = realloc(terms, sizeof(Expr*) * term_cap); }
        terms[term_n++] = term;
    }
    free(used);
    for (size_t i = 0; i < n; i++) expr_free(bases[i]);
    free(bases);
    free(mults);

    Expr* out = eval_and_free(expr_new_function(expr_new_symbol("Times"), terms, term_n));
    free(terms);
    return out;
}

Expr* builtin_factorsquarefree(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

    /* Parse a trailing Extension -> alpha option (mirrors builtin_factor).  Any
     * trailing argument that is not an Extension rule leaves the call
     * unevaluated (NULL), preserving the strict contract of the plain form. */
    size_t argc = res->data.function.arg_count;
    Expr* alpha_expr = NULL;  /* borrowed */
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (opt->type == EXPR_FUNCTION
            && opt->data.function.head
            && opt->data.function.head->type == EXPR_SYMBOL
            && (opt->data.function.head->data.symbol.name == SYM_Rule
                || opt->data.function.head->data.symbol.name == SYM_RuleDelayed)
            && opt->data.function.arg_count == 2
            && opt->data.function.args[0]->type == EXPR_SYMBOL
            && opt->data.function.args[0]->data.symbol.name == SYM_Extension) {
            alpha_expr = opt->data.function.args[1];
        } else {
            return NULL;  /* unknown / non-option trailing argument */
        }
    }

    /* A *real* extension is anything but Extension -> None (Automatic and an
     * explicit generator both route through Factor, which does its own
     * autodetect / over-Q fall-back when no generator applies). */
    bool real_ext = alpha_expr &&
        !(alpha_expr->type == EXPR_SYMBOL && alpha_expr->data.symbol.name
          && strcmp(alpha_expr->data.symbol.name, "None") == 0);

    if (real_ext) {
        Expr* out = factor_square_free_over_extension(res);
        if (out) return out;
        /* Fall through to the over-Q path on failure (graceful degradation). */
    }

    /* ---- Default over-Q path (unchanged; operates on args[0]) ---- */
#ifdef USE_FLINT
    /* Parametric radical field Q(t..)(sqrt k): squarefree-factor over the field
     * via FLINT (sqrt(k) -> S collapse). Only fires when a Sqrt[symbol] is
     * present; the classical pipeline blows up on a symbolic radicand. */
    {
        Expr* pf = flint_parametric_sqrt_factor_squarefree(res->data.function.args[0]);
        if (pf) return pf;
    }
#endif

    /* The square-free factorisation pipeline (PolynomialGCD against the
     * derivative) requires rational coefficients. */
    if (internal_args_contain_inexact(res)) {
        return internal_rationalize_then_numericalize(res, builtin_factorsquarefree);
    }

    return factor_square_free_dispatcher(res->data.function.args[0]);
}
