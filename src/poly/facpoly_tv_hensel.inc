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
