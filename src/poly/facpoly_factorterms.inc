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
    /* `poly` is nonzero here (the zero case returned at the top).  If no
     * monomial coefficient was found, the input is not a genuine polynomial in
     * the S variables -- e.g. it carries negative powers (a Laurent term such
     * as b^-1 left behind when an upstream Together failed to clear a
     * denominator), which the nonnegative-degree coefficient extraction cannot
     * see.  The content of a nonzero expression is a unit, NEVER 0; returning 0
     * here made FactorTerms multiply its factor list down to a spurious 0 (a
     * silent wrong answer, e.g. Simplify[(I - I Exp[-a b])/b^2] -> 0).  Return
     * the unit content 1 so nothing is factored out and the residue is intact. */
    return g ? g : expr_new_integer(1);
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
           strcmp(head, "Inequality") == 0 ||
           strcmp(head, "And") == 0 || strcmp(head, "Or") == 0 ||
           strcmp(head, "Not") == 0 || strcmp(head, "Xor") == 0;
}

/* Sign of the leading coefficient of the (expanded) polynomial `p` under
 * the lexicographic order of vars[0] > vars[1] > ...: +1, -1, or 0 when p is
 * zero or its leading coefficient is not a real number. */
static int ft_lex_lead_sign(Expr* p, Expr** vars, size_t nv) {
    if (nv == 0) {
        if (p->type == EXPR_INTEGER) return p->data.integer > 0 ? 1 : p->data.integer < 0 ? -1 : 0;
        if (p->type == EXPR_BIGINT) return mpz_sgn(p->data.bigint);
        if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL
            && p->data.function.head->data.symbol.name == SYM_Rational)
            return ft_lex_lead_sign(p->data.function.args[0], vars, 0);
        return 0;
    }
    int d = get_degree_poly(p, vars[0]);
    Expr* c = get_coeff(p, vars[0], d);
    int sg = ft_lex_lead_sign(c, vars + 1, nv - 1);
    expr_free(c);
    return sg;
}

/* The (positive) integer content of the polynomial `p`. */
static Expr* ft_number_content(Expr* p) {
    Expr* ex = expr_expand(p);
    size_t vc = 0, vcap = 8;
    Expr** vars = malloc(sizeof(Expr*) * vcap);
    collect_variables(ex, &vars, &vc, &vcap);
    Expr* c = ft_content_wrt_set(ex, vars, vc, NULL, 0);
    for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
    free(vars);
    expr_free(ex);
    if (!expr_is_integer_like(c)) { expr_free(c); return expr_new_integer(1); }
    return c;
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
    if (is_zero_poly(cur)) {                      /* FactorTermsList[0] = {1, 0} */
        expr_free(num_cont);
        num_cont = expr_new_integer(1);
    } else if (ft_lex_lead_sign(cur, all_vars, v_count) < 0) {
        /* The numerical factor carries the sign: the residue's leading term
         * (lexicographic, variables in canonical order) is positive, so
         * -6 x^2 - 12 x y is {-6, x^2 + 2 x y}. */
        Expr* neg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                  (Expr*[]){ expr_new_integer(-1), num_cont }, 2));
        num_cont = neg;
        Expr* nc = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                                 (Expr*[]){ expr_new_integer(-1), cur }, 2));
        cur = expr_expand(nc);
        expr_free(nc);
    }
    /* The integer content of the denominator belongs to the numerical factor
     * too: 3/4 x + 3/2 is {3/4, 2 + x}, (4 x + 4)/(6 y) is {2/3, ...}. */
    if (!ft_is_one_int(den)) {
        Expr* dc = expr_is_integer_like(den) ? expr_copy(den) : ft_number_content(den);
        if (!ft_is_one_int(dc)) {
            num_cont = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ num_cont, eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ expr_copy(dc), expr_new_integer(-1) }, 2)) }, 2));
            Expr* nd = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ den, eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ expr_copy(dc), expr_new_integer(-1) }, 2)) }, 2));
            den = nd;
        }
        expr_free(dc);
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
        Expr* q = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                                (Expr*[]){cur, inv_den}, 2));
        final = expr_expand(q);
        expr_free(q);
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
        var_arg->data.function.head->data.symbol.name == SYM_List) {
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
    bool is_ineq = (compound->data.function.head->type == EXPR_SYMBOL
                    && strcmp(compound->data.function.head->data.symbol.name, "Inequality") == 0);
    size_t n = compound->data.function.arg_count;
    Expr** new_args = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        if (is_ineq && (i & 1u) == 1) {
            new_args[i] = expr_copy(compound->data.function.args[i]);
            continue;
        }
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

/* Trailing options of FactorTerms / FactorTermsList: Modulus -> p (reduce
 * the coefficients mod p first) and Trig.  Returns the number of leading
 * non-option arguments, or 0 when the call stays unevaluated; *mod is the
 * modulus (borrowed) or NULL. */
static size_t ft_options(Expr* res, const char* head, Expr** mod) {
    *mod = NULL;
    size_t argc = res->data.function.arg_count;
    while (argc > 1 && fo_is_rule(res->data.function.args[argc - 1])) argc--;
    for (size_t i = argc; i < res->data.function.arg_count; i++) {
        Expr* r = res->data.function.args[i];
        const Expr* key = r->data.function.args[0];
        if (fo_sym_is(key, SYM_Modulus)) {
            if (!*mod) *mod = r->data.function.args[1];
        } else if (!fo_sym_is(key, SYM_Trig)) {
            char* kn = expr_to_string((Expr*)key);
            char* call = expr_to_string(res);
            mth_message(head, "optx", "Unknown option %s in %s.", kn ? kn : "?", call ? call : "?");
            free(kn); free(call);
            return 0;
        }
    }
    if (*mod && !(expr_is_integer_like(*mod)
                  && !(((*mod)->type == EXPR_INTEGER) && (*mod)->data.integer == 0)))
        *mod = NULL;
    return argc;
}

/* poly with its coefficients reduced mod m. */
static Expr* ft_reduce_mod(Expr* poly, Expr* m) {
    return eval_and_free(expr_new_function(expr_new_symbol("PolynomialMod"),
        (Expr*[]){ expr_copy(poly), expr_copy(m) }, 2));
}

Expr* builtin_factortermslist(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    Expr* mod = NULL;
    size_t argc = ft_options(res, "FactorTermsList", &mod);
    if (argc < 1 || argc > 2) return NULL;

    Expr* poly = res->data.function.args[0];
    Expr* var_arg = (argc == 2) ? res->data.function.args[1] : NULL;

    if (poly->type == EXPR_FUNCTION && poly->data.function.head->type == EXPR_SYMBOL
        && poly->data.function.head->data.symbol.name == SYM_List) {
        char* ps = expr_to_string(poly);
        mth_message("General", "poly", "%s is not a polynomial.", ps ? ps : "?");
        free(ps);
        return NULL;
    }

    /* FactorTermsList does not auto-thread (its result shape is already */
    /* a List, so threading would be ambiguous).                         */

    Expr** S = NULL;
    size_t S_count = 0;
    if (var_arg) ft_extract_vars(var_arg, &S, &S_count);

    Expr* reduced = mod ? ft_reduce_mod(poly, mod) : NULL;
    Expr* result = ft_compute_list(reduced ? reduced : poly, S, S_count);
    if (reduced) expr_free(reduced);
    free(S);
    return result;
}

Expr* builtin_factorterms(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    Expr* mod = NULL;
    size_t argc = ft_options(res, "FactorTerms", &mod);
    if (argc < 1 || argc > 2) return NULL;

    Expr* poly = res->data.function.args[0];
    Expr* var_arg = (argc == 2) ? res->data.function.args[1] : NULL;

    /* Auto-thread over List, Equal, Less, And, Or, ... */
    if (poly->type == EXPR_FUNCTION && poly->data.function.head->type == EXPR_SYMBOL) {
        const char* h = poly->data.function.head->data.symbol.name;
        if (ft_is_threading_head(h)) {
            return ft_thread(poly, var_arg, "FactorTerms");
        }
    }

    Expr** S = NULL;
    size_t S_count = 0;
    if (var_arg) ft_extract_vars(var_arg, &S, &S_count);

    Expr* reduced = mod ? ft_reduce_mod(poly, mod) : NULL;
    Expr* list = ft_compute_list(reduced ? reduced : poly, S, S_count);
    if (reduced) expr_free(reduced);
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
