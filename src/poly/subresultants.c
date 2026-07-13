/* ====================================================================
 * subresultants.c -- Principal subresultant coefficients.
 *
 *   Subresultants[poly1, poly2, var]
 *     generates the list of principal subresultant coefficients (PSCs)
 *     of poly1 and poly2 with respect to var.  The list has length
 *     Min[Exponent[poly1, var], Exponent[poly2, var]] + 1, its first
 *     element equals Resultant[poly1, poly2, var], and the first k PSCs
 *     vanish exactly when the polynomials share k roots (multiplicity
 *     counted).
 *
 * Two algorithms are provided:
 *
 *   1. subresultants_prs -- the efficient path.  Runs the subresultant
 *      polynomial remainder sequence (the same Bronstein gamma/beta/delta
 *      recurrence used by Resultant in poly.c), retaining each scaled
 *      remainder as a polynomial and placing it at its subresultant index
 *      via the fundamental theorem of subresultants.  PSC_j is the
 *      coefficient of var^j in the j-th subresultant polynomial S_j;
 *      defective (degree-gap) indices fall out as zero automatically.
 *
 *   2. subresultants_determinant -- the canonical definition and the
 *      fallback (used for algebraic-number coefficients, mirroring
 *      Resultant).  PSC_j = Det(M_j), where M_j is the Sylvester matrix
 *      restricted to its first (m-j) poly1-shift rows, first (n-j)
 *      poly2-shift rows, and first (n+m-2j) columns.  This generalises
 *      the Sylvester construction in resultant_internal.
 *
 * Memory convention matches the rest of Mathilda: every helper returning
 * Expr* hands fresh ownership to the caller; builtin_subresultants leaves
 * its input `res` alive for the evaluator to free.
 * ==================================================================== */

#include "internal.h"
#include "poly.h"
#include "expand.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <stdbool.h>

/* Descending coefficient array of an already-expanded polynomial:
 * returns a fresh array out[0..n] with out[k] = coeff of var^(n-k), so
 * out[0] is the leading coefficient and out[n] the constant term.  Each
 * slot is owned by the caller.  Uses the bulk extractor when possible. */
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

/* Det(M_j): the j-th principal subresultant coefficient via the Sylvester
 * minor.  p1d/p2d are descending coefficient arrays (length n+1 / m+1).
 * dim = n + m - 2j; an empty (0x0) matrix has determinant 1. */
static Expr* subres_det(Expr** p1d, int n, Expr** p2d, int m, int j) {
    int dim = n + m - 2 * j;
    if (dim <= 0) return expr_new_integer(1);

    int rowsP = m - j;          /* shifts of poly1 */
    int rowsQ = n - j;          /* shifts of poly2 */

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
    for (int i = 0; i < rowsP; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
        for (int c = 0; c < dim; c++) {
            int k = c - i;       /* index into descending coeff array */
            elems[c] = (k >= 0 && k <= n) ? expr_copy(p1d[k]) : expr_new_integer(0);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), elems, (size_t)dim);
        free(elems);
    }
    for (int i = 0; i < rowsQ; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)dim);
        for (int c = 0; c < dim; c++) {
            int k = c - i;
            elems[c] = (k >= 0 && k <= m) ? expr_copy(p2d[k]) : expr_new_integer(0);
        }
        rows[rowsP + i] = expr_new_function(expr_new_symbol(SYM_List), elems, (size_t)dim);
        free(elems);
    }

    Expr* matrix = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)dim);
    free(rows);

    Expr* det_call = expr_new_function(expr_new_symbol(SYM_Det), (Expr*[]){matrix}, 1);
    Expr* det = evaluate(det_call);
    expr_free(det_call);

    Expr* result = expr_expand(det);
    expr_free(det);
    return result;
}

/* Build the full PSC list {psc_0, ..., psc_L} (L = min(n,m)) by the
 * determinant definition. */
static Expr* subresultants_determinant(Expr** p1d, int n, Expr** p2d, int m) {
    int L = (n < m) ? n : m;
    Expr** psc = (Expr**)malloc(sizeof(Expr*) * (size_t)(L + 1));
    for (int j = 0; j <= L; j++) psc[j] = subres_det(p1d, n, p2d, m, j);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), psc, (size_t)(L + 1));
    free(psc);
    return list;
}

/* Subresultant polynomial remainder sequence -> principal subresultant
 * coefficients.  Runs the Bronstein gamma/beta/delta subresultant PRS
 * (identical recurrence to resultant_subresultant in poly.c), retaining
 * the degree and leading coefficient of each chain member.
 *
 * Fundamental theorem of subresultants: the nonzero principal subresultant
 * coefficients occur exactly at the degrees appearing in the PRS, and for a
 * chain member R_p (the (deg(R_{p-1})-1)-th subresultant polynomial),
 *
 *     psc[deg(R_p)] = lc(R_p) ^ (deg(R_{p-1}) - deg(R_p)).
 *
 * All other (defective / degree-gap) indices are zero.  The top index
 * L = min(n,m) is the regular subresultant of the input pair when the
 * degrees differ (handled by the p=1 step), and is 1 when the degrees are
 * equal (the empty 0x0 minor).
 *
 * e1/e2 are already expanded; n = deg(e1, var), m = deg(e2, var), both >= 0.
 * Returns the PSC list {psc_0, ..., psc_L}, or NULL to defer to the
 * determinant path. */
static Expr* subresultants_prs(Expr* e1, Expr* e2, int n, int m, Expr* var) {
    int L = (n < m) ? n : m;

    /* Orient so deg(A) >= deg(B); remember whether poly1/poly2 swapped. */
    Expr *A, *B;
    int dA, dB;
    bool swapped;
    if (n >= m) { A = e1; B = e2; dA = n; dB = m; swapped = false; }
    else        { A = e2; B = e1; dA = m; dB = n; swapped = true;  }

    size_t cap = 8;
    Expr** R    = (Expr**)malloc(sizeof(Expr*) * cap);
    int*   deg  = (int*)  malloc(sizeof(int)   * cap);
    Expr** lc   = (Expr**)malloc(sizeof(Expr*) * cap);
    for (size_t t = 0; t < cap; t++) { R[t] = NULL; lc[t] = NULL; deg[t] = -1; }

    R[0] = expr_copy(A);
    R[1] = expr_copy(B);
    deg[0] = dA; deg[1] = dB;
    lc[0] = NULL;                                /* never read */
    lc[1] = get_coeff_expanded(R[1], var, dB);

    /* Bronstein subresultant PRS (identical recurrence to
     * resultant_subresultant in poly.c, which is well-tested):
     *
     *   delta_i = deg(r_{i-1}) - deg(r_i)
     *   gamma_1 = -1,  beta_1 = (-1)^{delta_1+1}
     *   gamma_i = (-lc_{i-1})^{delta_{i-1}} * gamma_{i-1}^{1-delta_{i-1}}
     *   beta_i  = -lc_{i-1} * gamma_i^{delta_i}
     *   r_{i+1} = prem(r_{i-1}, r_i) / beta_i
     *
     * This recurrence yields chain members whose leading coefficients give
     * the Mathematica principal subresultant coefficients directly at normal
     * (delta=1) steps; degree-gap (delta>1) landings carry an extra
     * lc^{delta-1} factor which the extraction below deflates. */
    Expr* gamma = expr_new_integer(-1);
    int prev_delta = dA - dB;
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

        /* Principal subresultant coefficients via the subresultant chain.
         * Each chain member R_p is the defective subresultant of index
         * deg(R_{p-1}) - 1; the regular subresultant of degree deg(R_p) has
         * leading coefficient s_p, the p-th principal subresultant
         * coefficient.  With the convention s at degree deg(R_0) = 1,
         *
         *     delta_p = deg(R_{p-1}) - deg(R_p)
         *     s_p     = lc(R_p)^{delta_p} / s_{p-1}^{delta_p - 1}.
         *
         * Defective (degree-gap) indices have leading coefficient zero, so
         * only the regular degrees deg(R_p) receive a nonzero entry.  (The
         * earlier per-member deflation clc[p] = lc[p]/clc[p-1]^{delta-1} was
         * wrong across gaps whose preceding s_{p-1} is not a unit -- e.g.
         * 7 x^6 + 8 x - 9 vs 2 x^7 + ...; this recurrence carries the
         * cumulative s_{p-1} factor correctly.) */
        Expr** psc = (Expr**)malloc(sizeof(Expr*) * (size_t)(L + 1));
        for (int j = 0; j <= L; j++) psc[j] = NULL;

        Expr* s_prev = expr_new_integer(1);
        for (int p = 1; p <= t; p++) {
            if (deg[p - 1] <= deg[p]) continue;          /* non-strict drop (p=1, n==m) */
            int delta = deg[p - 1] - deg[p];             /* delta_p >= 1 */
            int ec    = deg[p];

            Expr* num   = internal_power((Expr*[]){expr_copy(lc[p]), expr_new_integer(delta)}, 2);
            Expr* den   = internal_power((Expr*[]){expr_copy(s_prev), expr_new_integer(-(delta - 1))}, 2);
            Expr* snew  = internal_times((Expr*[]){num, den}, 2);
            Expr* ssimp = internal_cancel((Expr*[]){snew}, 1);
            expr_free(s_prev);
            s_prev = expr_expand(ssimp);
            expr_free(ssimp);

            if (ec >= 0 && ec <= L) {
                if (psc[ec]) expr_free(psc[ec]);
                psc[ec] = expr_copy(s_prev);
            }
        }
        expr_free(s_prev);

        /* Equal input degrees: the top minor is empty -> psc_L = 1. */
        if (n == m) {
            if (psc[L]) expr_free(psc[L]);
            psc[L] = expr_new_integer(1);
        }

        /* Swap correction: S_j(e1,e2) = (-1)^{(n-j)(m-j)} S_j(e2,e1). */
        if (swapped) {
            for (int j = 0; j <= L; j++) {
                if (!psc[j]) continue;
                if ((((long)(n - j) * (long)(m - j)) & 1L)) {
                    Expr* neg = internal_times((Expr*[]){expr_new_integer(-1), psc[j]}, 2);
                    psc[j] = expr_expand(neg);
                    expr_free(neg);
                }
            }
        }

        for (int j = 0; j <= L; j++) if (!psc[j]) psc[j] = expr_new_integer(0);
        result = expr_new_function(expr_new_symbol(SYM_List), psc, (size_t)(L + 1));
        free(psc);
    }

    for (size_t t = 0; t < cap; t++) {
        if (R[t])  expr_free(R[t]);
        if (lc[t]) expr_free(lc[t]);
    }
    free(R);
    free(deg);
    free(lc);
    return result;
}

Expr* builtin_subresultants(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;

    Expr* p1  = res->data.function.args[0];
    Expr* p2  = res->data.function.args[1];
    Expr* var = res->data.function.args[2];
    if (var->type != EXPR_SYMBOL) return NULL;

    /* Both arguments must be polynomials in var (as Resultant requires). */
    Expr* pq1 = internal_polynomialq((Expr*[]){expr_copy(p1), expr_copy(var)}, 2);
    bool is_poly1 = (pq1->type == EXPR_SYMBOL && pq1->data.symbol.name == SYM_True);
    expr_free(pq1);
    Expr* pq2 = internal_polynomialq((Expr*[]){expr_copy(p2), expr_copy(var)}, 2);
    bool is_poly2 = (pq2->type == EXPR_SYMBOL && pq2->data.symbol.name == SYM_True);
    expr_free(pq2);
    if (!is_poly1 || !is_poly2) return NULL;

    Expr* e1 = expr_expand(p1);
    Expr* e2 = expr_expand(p2);
    int n = get_degree_poly(e1, var);
    int m = get_degree_poly(e2, var);
    if (n < 0 || m < 0) {   /* identically-zero input: leave unevaluated */
        expr_free(e1);
        expr_free(e2);
        return NULL;
    }

    Expr* result = NULL;

    /* Efficient path: subresultant PRS.  Skipped for algebraic-number
     * coefficients (Sqrt[N], cube roots, ...) where the pseudo-remainder
     * chain bloats -- the determinant definition handles those instead,
     * mirroring Resultant. */
    if (!subres_has_algebraic(e1) && !subres_has_algebraic(e2)) {
        result = subresultants_prs(e1, e2, n, m, var);
    }

    if (!result) {
        Expr** p1d = desc_coeffs(e1, var, n);
        Expr** p2d = desc_coeffs(e2, var, m);
        result = subresultants_determinant(p1d, n, p2d, m);
        for (int k = 0; k <= n; k++) expr_free(p1d[k]);
        for (int k = 0; k <= m; k++) expr_free(p2d[k]);
        free(p1d);
        free(p2d);
    }

    expr_free(e1);
    expr_free(e2);
    return result;
}

void subresultants_init(void) {
    symtab_add_builtin("Subresultants", builtin_subresultants);
    symtab_get_def("Subresultants")->attributes |= ATTR_PROTECTED;
}
