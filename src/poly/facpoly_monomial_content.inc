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
 * (which can occur in Mathilda's canonical form, e.g. the parser builds
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
        && e->data.function.head->data.symbol.name == SYM_Times) {
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
        && e->data.function.head->data.symbol.name == SYM_Power
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
    if (P->data.function.head->data.symbol.name != SYM_Plus) return NULL;
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
