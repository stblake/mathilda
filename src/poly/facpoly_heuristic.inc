#define HEURISTIC_FACTOR_MAX_DEPTH 64
static int heuristic_factor_depth = 0;

static Expr* heuristic_factor_impl(Expr* P);

static Expr* heuristic_factor(Expr* P) {
    /* Defensive: callers occasionally feed us NULL when an upstream
     * helper (e.g. exact_poly_div over a non-field coefficient ring)
     * has signalled "cannot proceed" via NULL.  Returning NULL here
     * keeps the cascade well-defined; callers must be prepared for
     * NULL just like the rest of the polynomial pipeline. */
    if (!P) return NULL;
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
        /* `exact_poly_div` returns NULL when the division is not
         * actually exact in the coefficient ring chosen by `vars`
         * (e.g. when P has radical coefficients like Sqrt[a] that
         * land outside Q or Q[i]).  Treat that as "no content can
         * be extracted at this level" and return P unchanged --
         * passing NULL down to heuristic_factor would deref NULL. */
        if (pp == NULL) {
            { for(size_t i=0; i<v_count; i++) expr_free(vars[i]); free(vars); }
            expr_free(cont);
            return expr_copy(P);
        }
        /* Guard against the no-progress case: if exact_poly_div left us
         * with pp == 1 (i.e. cont == P up to a constant), recursing on
         * cont = P loops indefinitely.  Bail and return the input. */
        if (pp->type == EXPR_INTEGER && pp->data.integer == 1) {
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
