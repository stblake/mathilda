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
        && e->data.function.head->data.symbol.name == SYM_Times) {
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
        && (f->data.function.head->data.symbol.name == SYM_Rational
            || f->data.function.head->data.symbol.name == SYM_Complex)) {
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

/* Nesting depth of the z-independent split.  The exploratory step (line:
 * `heuristic_factor(P_at_z0)` below) speculatively factors P specialised at
 * each of v_count variables -- and that speculative factorisation can itself
 * re-enter z_split, giving an O(v_count!) explosion: a single Simplify of a
 * derivative in Tanh[x/2]/Sech[x/2] reached v_count=6 and recursed four
 * levels deep, spending seconds re-specialising. The common, productive case
 * (peel one z-independent shape factor, then recurse on the smaller residual)
 * needs only shallow nesting, so cap the *speculative* nesting depth. Beyond
 * the cap z_split returns NULL and the caller falls through to the other
 * multivariate paths (Hensel, factor_roots) or returns P unfactored -- a
 * less-factored result, never an incorrect one. */
static int z_split_depth = 0;
#define Z_SPLIT_MAX_DEPTH 2

static Expr* factor_via_z_independent_split(Expr* P, Expr** vars, size_t v_count) {
    if (v_count < 3) return NULL;
    if (z_split_depth >= Z_SPLIT_MAX_DEPTH) return NULL;
    z_split_depth++;

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
                && f->data.function.head->data.symbol.name == SYM_Power
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
        z_split_depth--;
        return eval_and_free(result);
    }

    z_split_depth--;
    return NULL;
}
