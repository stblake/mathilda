Expr* builtin_factor(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1) return NULL;

#ifdef USE_FLINT
    /* Parametric radical field Q(t..)(sqrt k), k a free symbol: factor over the
     * field via FLINT (sqrt(k) -> S collapse -> multivariate factoring over Q).
     * The classical algebraic-factoring path (qafactor/Trager) is built for
     * constant radicands and blows up / hangs on a symbolic radicand — this is
     * the Apart/Factor blocker in the parametric rational integrator. Only fires
     * when a Sqrt[symbol] is actually present, matching WL's default of factoring
     * over the field the input's radicals generate. */
    if (res->data.function.arg_count == 1) {
        Expr* pf = flint_parametric_sqrt_factor(res->data.function.args[0]);
        if (pf) return pf;
    }
#endif

    /* Phase G5: Extension -> α option.  Walk trailing args looking for
     * Rule[Extension, α] / RuleDelayed[Extension, α].  When present,
     * dispatch to the algebraic-factoring path in qafactor.c (Trager).
     *
     * Phase G9: when α is the symbol `Automatic`, run extension_autodetect
     * on the polynomial argument and use the resulting (tier-1) single
     * generator's surface form as α.
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
        /* Auto-detect when the user passed Extension -> Automatic. */
        Expr* alpha_auto = NULL;
        QATower* auto_tower = NULL;
        if (alpha_expr && alpha_expr->type == EXPR_SYMBOL
            && alpha_expr->data.symbol
            && strcmp(alpha_expr->data.symbol, "Automatic") == 0) {
            auto_tower = extension_autodetect(res->data.function.args[0]);
            if (auto_tower && auto_tower->n == 1) {
                alpha_auto = expr_copy(auto_tower->alpha_renders[0]);
                alpha_expr = alpha_auto;
            } else {
                /* No detectable generator, or multi-generator tower.
                 * Tier-1: fall through to plain Factor. */
                alpha_expr = NULL;
            }
        }
        /* Skip the qa-factor path when alpha_expr resolves to a no-op
         * sentinel (None or Automatic-with-no-detect). */
        if (alpha_expr && alpha_expr->type == EXPR_SYMBOL
            && alpha_expr->data.symbol
            && strcmp(alpha_expr->data.symbol, "None") == 0) {
            alpha_expr = NULL;
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
                if (result) {
                    if (alpha_auto) expr_free(alpha_auto);
                    if (auto_tower) qa_tower_free(auto_tower);
                    return result;
                }
                /* Fall through to plain Factor on failure. */
            } else {
                for (size_t i = 0; i < vc; i++) expr_free(vars[i]);
                free(vars);
            }
        }
        /* Cleanup auto-detect state before any fall-through to plain Factor. */
        if (alpha_auto) expr_free(alpha_auto);
        if (auto_tower) qa_tower_free(auto_tower);
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
     * LessEqual, GreaterEqual, Inequality, And, Or). Mathematica applies
     * Factor elementwise inside these structures, e.g.,
     *   Factor[1 < 1+2x+x^2+1/(1+x) < 2]
     * yields Inequality[1, Less, factored, Less, 2].  We thread by
     * recursively invoking Factor on each comparison argument, then
     * rebuilding the same head with the factored subterms.
     * ATTR_LISTABLE on Factor handles the List case automatically;
     * this hand-threading covers the comparison / logic heads it
     * doesn't reach.  For Inequality specifically, operator-symbol
     * slots (odd indices) are passed through verbatim. */
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
                h == SYM_Inequality ||
                h == SYM_And ||
                h == SYM_Or) {
                bool is_ineq = (h == SYM_Inequality);
                size_t n = arg0->data.function.arg_count;
                Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
                for (size_t i = 0; i < n; i++) {
                    if (is_ineq && (i & 1u) == 1) {
                        /* Operator-symbol slot: copy through unchanged. */
                        new_args[i] = expr_copy(arg0->data.function.args[i]);
                        continue;
                    }
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
     * remain in their non-fully-factored form.  See plans/FACTOR_PLAN.md
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
