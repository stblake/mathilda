/* ===================================================================== */
/*  Bivariate-Hensel wiring                                              */
/*                                                                       */
/*  Glue between Mathilda's existing univariate Berlekamp-Zassenhaus     */
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
     * the next alpha. See plans/FACTOR_PLAN.md §12.F6 for the Fateman trace. */
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

            /* Shift y -> y + alpha so the lift sees α = 0.  (The Mathilda
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
         * does not match the canonical form Mathilda previously
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
