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
/*  factor call is fast because Mathilda's univariate Berlekamp-          */
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
             * a false "irreducible" confirmation.  See plans/FACTOR_PLAN.md
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
