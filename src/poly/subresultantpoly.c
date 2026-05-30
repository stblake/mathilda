/* ====================================================================
 * subresultantpoly.c -- Subresultant polynomials.
 *
 *   SubresultantPolynomials[poly1, poly2, var]
 *     generates the list of subresultant polynomials {S_0, ..., S_m} of
 *     poly1 and poly2 with respect to var, where m = Exponent[poly2, var].
 *     The list has length m + 1, its first element is
 *     Resultant[poly1, poly2, var], and the coefficient of var^j in S_j is
 *     the j-th principal subresultant coefficient (i.e. Subresultants[
 *     poly1, poly2, var][[j+1]]).  Requires Exponent[poly1, var] >=
 *     Exponent[poly2, var] and exact coefficients.
 *
 * Each subresultant polynomial S_j is, by the fundamental theorem of
 * subresultants, either zero or a scalar multiple of a single member of
 * the subresultant polynomial remainder sequence (PRS).  We therefore
 * reuse the Bronstein gamma/beta/delta PRS (the same recurrence as
 * Subresultants in subresultants.c and Resultant in poly.c) and classify
 * every output index:
 *
 *   * Regular index  (j == deg(R_p) for a strict-drop chain step p):
 *       S_j = (psc_j / lc(R_p)) * R_p, the chain member rescaled so its
 *       leading coefficient is the j-th principal subresultant coefficient
 *       psc_j = clc_p^{delta_p}.  (delta_p == 1 leaves S_j = R_p.)
 *   * Defective index (j == deg(R_{p-1}) - 1 across a degree gap,
 *       delta_p > 1): a lower-degree polynomial computed directly from the
 *       determinant-polynomial definition.  Such indices sit high in the
 *       chain, so the associated Sylvester minor is small.
 *   * Otherwise: S_j = 0 (gap interior, or below the last chain degree).
 *
 * For algebraic-number coefficients (where the pseudo-remainder chain
 * bloats) we skip the PRS and build the whole list from the determinant-
 * polynomial definition, mirroring Resultant / Subresultants.
 *
 * Memory convention matches the rest of Mathilda: every helper returning
 * Expr* hands fresh ownership to the caller; builtin_subresultantpolynomials
 * leaves its input `res` alive for the evaluator to free.
 * ==================================================================== */

#include "internal.h"
#include "poly.h"
#include "expand.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "common.h"
#include "print.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

/* Descending coefficient array of an already-expanded polynomial:
 * out[k] = coeff of var^(n-k), so out[0] is the leading coefficient and
 * out[n] the constant term.  Each slot is owned by the caller. */
static Expr** desc_coeffs(Expr* expanded, Expr* var, int n) {
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)(n + 1));
    Expr** asc = NULL;
    if (get_all_coeffs_expanded(expanded, var, n, &asc)) {
        for (int k = 0; k <= n; k++) out[k] = asc[n - k];
        free(asc);
    } else {
        for (int k = 0; k <= n; k++) out[k] = get_coeff_expanded(expanded, var, n - k);
    }
    return out;
}

/* The degree-<=j truncation of var^shift * poly, returned as a polynomial
 * Sum_{d=0}^{j} coeff(var^(d-shift), poly) var^d.  pd is the descending
 * coefficient array (pd[k] = coeff of var^(deg-k)). */
static Expr* trunc_shift_poly(Expr** pd, int deg, int shift, int j, Expr* var) {
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(j + 1));
    int nt = 0;
    for (int d = 0; d <= j; d++) {
        int e = d - shift;                  /* exponent in the original poly */
        if (e < 0 || e > deg) continue;
        Expr* coeff = pd[deg - e];          /* coeff of var^e */
        if (d == 0) {
            terms[nt++] = expr_copy(coeff);
        } else {
            Expr* pw = internal_power((Expr*[]){expr_copy(var), expr_new_integer(d)}, 2);
            terms[nt++] = internal_times((Expr*[]){expr_copy(coeff), pw}, 2);
        }
    }
    Expr* poly = (nt == 0) ? expr_new_integer(0) : internal_plus(terms, (size_t)nt);
    free(terms);
    return poly;
}

/* The j-th subresultant polynomial S_j via the determinant-polynomial
 * definition.  p1d/p2d are descending coefficient arrays (length n+1 /
 * m+1).  dim = n + m - 2j; the first dim-1 columns hold the high-degree
 * coefficient columns (var^{n+m-j-1} .. var^{j+1}) as in subres_det, and
 * the last column holds each row's degree-<=j truncation as a polynomial,
 * so Det yields S_j directly.  An empty (0x0) matrix gives 1. */
static Expr* subres_poly_det(Expr** p1d, int n, Expr** p2d, int m, int j, Expr* var) {
    int dim = n + m - 2 * j;
    if (dim <= 0) return expr_new_integer(1);

    int rowsP = m - j;          /* shifts of poly1 */
    int rowsQ = n - j;          /* shifts of poly2 */

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
    for (int i = 0; i < rowsP; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
        for (int c = 0; c < dim - 1; c++) {
            int k = c - i;
            elems[c] = (k >= 0 && k <= n) ? expr_copy(p1d[k]) : expr_new_integer(0);
        }
        elems[dim - 1] = trunc_shift_poly(p1d, n, m - j - 1 - i, j, var);
        rows[i] = expr_new_function(expr_new_symbol("List"), elems, (size_t)dim);
        free(elems);
    }
    for (int i = 0; i < rowsQ; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
        for (int c = 0; c < dim - 1; c++) {
            int k = c - i;
            elems[c] = (k >= 0 && k <= m) ? expr_copy(p2d[k]) : expr_new_integer(0);
        }
        elems[dim - 1] = trunc_shift_poly(p2d, m, n - j - 1 - i, j, var);
        rows[rowsP + i] = expr_new_function(expr_new_symbol("List"), elems, (size_t)dim);
        free(elems);
    }

    Expr* matrix = expr_new_function(expr_new_symbol("List"), rows, (size_t)dim);
    free(rows);

    Expr* det_call = expr_new_function(expr_new_symbol("Det"), (Expr*[]){matrix}, 1);
    Expr* det = evaluate(det_call);
    expr_free(det_call);

    Expr* result = expr_expand(det);
    expr_free(det);
    return result;
}

/* Build the full subresultant-polynomial list {S_0, ..., S_m} by the
 * determinant-polynomial definition (the fallback / algebraic-coefficient
 * path). */
static Expr* subrespoly_determinant(Expr** p1d, int n, Expr** p2d, int m, Expr* var) {
    Expr** S = (Expr**)malloc(sizeof(Expr*) * (size_t)(m + 1));
    for (int j = 0; j <= m; j++) S[j] = subres_poly_det(p1d, n, p2d, m, j, var);
    Expr* list = expr_new_function(expr_new_symbol("List"), S, (size_t)(m + 1));
    free(S);
    return list;
}

/* Subresultant PRS -> subresultant polynomials.  e1/e2 are already
 * expanded; n = deg(e1, var) >= m = deg(e2, var) >= 0.  Returns the list
 * {S_0, ..., S_m}, or NULL to defer to the determinant path. */
static Expr* subrespoly_prs(Expr* e1, Expr* e2, int n, int m, Expr* var) {
    size_t cap = 8;
    Expr** R    = (Expr**)malloc(sizeof(Expr*) * cap);
    int*   deg  = (int*)  malloc(sizeof(int)   * cap);
    Expr** lc   = (Expr**)malloc(sizeof(Expr*) * cap);
    for (size_t t = 0; t < cap; t++) { R[t] = NULL; lc[t] = NULL; deg[t] = -1; }

    R[0] = expr_copy(e1);
    R[1] = expr_copy(e2);
    deg[0] = n; deg[1] = m;
    lc[0] = NULL;                                /* never read */
    lc[1] = get_coeff_expanded(R[1], var, m);

    Expr* gamma = expr_new_integer(-1);
    int prev_delta = n - m;
    Expr* beta_cur = ((prev_delta + 1) & 1) ? expr_new_integer(-1) : expr_new_integer(1);
    int i = 1;
    bool failed = false;

    while (!is_zero_poly(R[i])) {
        Expr* pr = pseudo_rem_standard(R[i - 1], R[i], var);
        Expr* Rn = pr ? poly_divide_by_scalar(pr, beta_cur, var) : NULL;
        if (pr) expr_free(pr);
        if (!Rn) { failed = true; break; }

        if ((size_t)(i + 2) >= cap) {
            size_t nc = cap * 2;
            R   = realloc(R,   sizeof(Expr*) * nc);
            deg = realloc(deg, sizeof(int)   * nc);
            lc  = realloc(lc,  sizeof(Expr*) * nc);
            for (size_t t = cap; t < nc; t++) { R[t] = NULL; lc[t] = NULL; deg[t] = -1; }
            cap = nc;
        }

        R[i + 1] = Rn;
        if (is_zero_poly(Rn)) { i = i + 1; break; }
        deg[i + 1] = get_degree_poly(Rn, var);
        lc[i + 1]  = get_coeff_expanded(Rn, var, deg[i + 1]);

        i = i + 1;

        /* gamma_i = (-lc_{i-1})^{delta_{i-1}} * gamma_{i-1}^{1-delta_{i-1}} */
        Expr* neg_r = internal_times((Expr*[]){expr_new_integer(-1), expr_copy(lc[i - 1])}, 2);
        Expr* term1 = internal_power((Expr*[]){neg_r, expr_new_integer(prev_delta)}, 2);
        Expr* term2 = internal_power((Expr*[]){expr_copy(gamma), expr_new_integer(1 - prev_delta)}, 2);
        Expr* gnew  = internal_times((Expr*[]){term1, term2}, 2);
        Expr* gsimp = internal_cancel((Expr*[]){gnew}, 1);
        expr_free(gamma);
        gamma = expr_expand(gsimp);
        expr_free(gsimp);

        int new_delta = deg[i - 1] - deg[i];     /* delta_i */
        prev_delta = new_delta;

        /* beta_i = -lc_{i-1} * gamma_i^{delta_i} (for the next division). */
        Expr* gpow  = internal_power((Expr*[]){expr_copy(gamma), expr_new_integer(new_delta)}, 2);
        Expr* bnew  = internal_times((Expr*[]){expr_new_integer(-1), expr_copy(lc[i - 1]), gpow}, 3);
        Expr* bsimp = internal_cancel((Expr*[]){bnew}, 1);
        expr_free(beta_cur);
        beta_cur = expr_expand(bsimp);
        expr_free(bsimp);
    }
    expr_free(beta_cur);
    expr_free(gamma);

    Expr* result = NULL;
    if (!failed) {
        /* Last nonzero chain index. */
        int t = i;
        while (t >= 0 && (R[t] == NULL || is_zero_poly(R[t]))) t--;

        Expr** S = (Expr**)malloc(sizeof(Expr*) * (size_t)(m + 1));
        for (int j = 0; j <= m; j++) S[j] = NULL;

        /* Walk the chain top-down.  Each member R_p (p >= 1) is the
         * defective subresultant of index deg(R_{p-1}) - 1; the regular
         * subresultant of degree deg(R_p) is R_p rescaled so its leading
         * coefficient is the principal subresultant coefficient at that
         * degree.  s_prev tracks that coefficient at the previous regular
         * degree, with the convention s at degree n = 1:
         *
         *   delta_p   = deg(R_{p-1}) - deg(R_p)
         *   S_{deg R_p}        = (lc(R_p) / s_prev)^{delta_p - 1} * R_p   (regular)
         *   S_{deg R_{p-1}-1}  = R_p                                       (defective)
         *   s at degree deg R_p = lc(R_p)^{delta_p} / s_prev^{delta_p - 1}
         *
         * Indices strictly inside a degree gap, and indices below the last
         * chain degree, are zero. */
        Expr* s_prev = expr_new_integer(1);
        for (int p = 1; p <= t; p++) {
            if (deg[p - 1] <= deg[p]) continue;          /* non-strict drop (p=1, n==m) */
            int delta = deg[p - 1] - deg[p];             /* delta_p >= 1 */
            int ec    = deg[p];                          /* regular index */

            /* Regular subresultant S_{deg R_p}. */
            if (ec >= 0 && ec <= m) {
                Expr* invsp = internal_power((Expr*[]){expr_copy(s_prev), expr_new_integer(-1)}, 2);
                Expr* ratio = internal_times((Expr*[]){expr_copy(lc[p]), invsp}, 2);
                Expr* scale = internal_power((Expr*[]){ratio, expr_new_integer(delta - 1)}, 2);
                Expr* prod  = internal_times((Expr*[]){scale, expr_copy(R[p])}, 2);
                Expr* sc    = internal_cancel((Expr*[]){prod}, 1);
                if (S[ec]) expr_free(S[ec]);
                S[ec] = expr_expand(sc);
                expr_free(sc);
            }

            /* Defective subresultant S_{deg R_{p-1} - 1} = R_p (degree gap). */
            if (delta > 1) {
                int dc = deg[p - 1] - 1;
                if (dc >= 0 && dc <= m) {
                    if (S[dc]) expr_free(S[dc]);
                    S[dc] = expr_expand(R[p]);
                }
            }

            /* s at degree deg(R_p) = lc(R_p)^delta / s_prev^{delta-1}. */
            Expr* num   = internal_power((Expr*[]){expr_copy(lc[p]), expr_new_integer(delta)}, 2);
            Expr* den   = internal_power((Expr*[]){expr_copy(s_prev), expr_new_integer(-(delta - 1))}, 2);
            Expr* snew  = internal_times((Expr*[]){num, den}, 2);
            Expr* ssimp = internal_cancel((Expr*[]){snew}, 1);
            expr_free(s_prev);
            s_prev = expr_expand(ssimp);
            expr_free(ssimp);
        }
        expr_free(s_prev);

        /* Equal input degrees: the top minor is empty -> S_m = 1. */
        if (n == m) {
            if (S[m]) expr_free(S[m]);
            S[m] = expr_new_integer(1);
        }

        for (int j = 0; j <= m; j++) if (!S[j]) S[j] = expr_new_integer(0);
        result = expr_new_function(expr_new_symbol("List"), S, (size_t)(m + 1));
        free(S);
    }

    for (size_t k = 0; k < cap; k++) {
        if (R[k])  expr_free(R[k]);
        if (lc[k]) expr_free(lc[k]);
    }
    free(R);
    free(deg);
    free(lc);
    return result;
}

/* SubresultantPolynomials::npolys -- inexact coefficients, a non-polynomial
 * argument, or deg(poly1) < deg(poly2).  Emits the diagnostic and leaves the
 * call unevaluated. */
static Expr* subrespoly_emit_npolys(Expr* p1, Expr* p2, Expr* var) {
    char* s1 = expr_to_string(p1);
    char* s2 = expr_to_string(p2);
    char* sv = expr_to_string(var);
    fprintf(stderr,
            "SubresultantPolynomials::npolys: %s and %s should be polynomials "
            "with exact coefficients and the degree of %s in %s should not be "
            "less than the degree of %s in %s.\n",
            s1 ? s1 : "?", s2 ? s2 : "?",
            s1 ? s1 : "?", sv ? sv : "?",
            s2 ? s2 : "?", sv ? sv : "?");
    free(s1);
    free(s2);
    free(sv);
    return NULL;
}

Expr* builtin_subresultantpolynomials(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* p1  = res->data.function.args[0];
    Expr* p2  = res->data.function.args[1];
    Expr* var = res->data.function.args[2];
    if (var->type != EXPR_SYMBOL) return NULL;

    /* Exact coefficients are required. */
    CommonInexactInfo ix1 = common_scan_inexact(p1);
    CommonInexactInfo ix2 = common_scan_inexact(p2);
    if (ix1.has_inexact || ix2.has_inexact) return subrespoly_emit_npolys(p1, p2, var);

    /* Both arguments must be polynomials in var. */
    Expr* pq1 = internal_polynomialq((Expr*[]){expr_copy(p1), expr_copy(var)}, 2);
    bool is_poly1 = (pq1->type == EXPR_SYMBOL && pq1->data.symbol == SYM_True);
    expr_free(pq1);
    Expr* pq2 = internal_polynomialq((Expr*[]){expr_copy(p2), expr_copy(var)}, 2);
    bool is_poly2 = (pq2->type == EXPR_SYMBOL && pq2->data.symbol == SYM_True);
    expr_free(pq2);
    if (!is_poly1 || !is_poly2) return subrespoly_emit_npolys(p1, p2, var);

    Expr* e1 = expr_expand(p1);
    Expr* e2 = expr_expand(p2);
    int n = get_degree_poly(e1, var);
    int m = get_degree_poly(e2, var);
    if (n < 0 || m < 0) {           /* identically-zero input: leave unevaluated */
        expr_free(e1);
        expr_free(e2);
        return NULL;
    }
    if (n < m) {                    /* requirement: deg(poly1) >= deg(poly2) */
        expr_free(e1);
        expr_free(e2);
        return subrespoly_emit_npolys(p1, p2, var);
    }

    Expr* result = NULL;

    /* Efficient path: subresultant PRS.  Skipped for algebraic-number
     * coefficients (where the chain bloats); the determinant definition
     * handles those, mirroring Subresultants. */
    if (!subres_has_algebraic(e1) && !subres_has_algebraic(e2)) {
        result = subrespoly_prs(e1, e2, n, m, var);
    }

    if (!result) {
        Expr** p1d = desc_coeffs(e1, var, n);
        Expr** p2d = desc_coeffs(e2, var, m);
        result = subrespoly_determinant(p1d, n, p2d, m, var);
        for (int k = 0; k <= n; k++) expr_free(p1d[k]);
        for (int k = 0; k <= m; k++) expr_free(p2d[k]);
        free(p1d);
        free(p2d);
    }

    expr_free(e1);
    expr_free(e2);
    return result;
}

void subresultantpolynomials_init(void) {
    symtab_add_builtin("SubresultantPolynomials", builtin_subresultantpolynomials);
    symtab_get_def("SubresultantPolynomials")->attributes |= ATTR_PROTECTED;
}
