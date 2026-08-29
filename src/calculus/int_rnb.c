/* int_rnb.c — RischNormanBlake: parallel Risch-Norman integration over a
 * simple radical extension L = K(y), y^m = q(x), K = Q(x), Dx = 1.
 *
 * Reference: S. Blake, "Parallel Integration over Simple Radical Extensions"
 * (rn-radicals.tex); prototype rnrad.py / bounds.py.  Structure mirrors the
 * transcendental parallel-Risch engine intrischnorman.c: all polynomial /
 * rational-function algebra is carried on Expr* trees through the CAS
 * evaluator (Together / Cancel / Factor / Series / Solve / Det / ...), with
 * FLINT reached only for the exact-rational linear solve.
 *
 * Field layout (Proposition 3.5, rnrad.py:RadicalField):
 *   q = prod_j Q_j^j            (squarefree decomposition)
 *   E_i = prod_j Q_j^floor(ij/m),  w_i = y^i / E_i    (integral basis of O)
 *   w_i w_k = w_{(i+k) mod m} * prod_j Q_j^{floor((i+k)j/m)-floor(ij/m)-floor(kj/m)}
 *   D w_i = Lam_i w_i,  Lam_i = sum_j {ij/m} Q_j'/Q_j
 * An element of L is stored as an array a[0..m-1] of rational functions of x,
 * meaning sum_i a_i w_i.
 */

#include "int_rnb.h"

#include "expr.h"
#include "core.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "poly.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "risch_canonical.h"
#include "flint_bridge.h"
#include "print.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Wall-clock budget: a radical integrand with a complicated denominator can
 * drive the residue/nullspace/solve work arbitrarily deep.  Since this engine
 * runs inside the Integrate cascade, an unbounded run would hang the whole
 * dispatcher, so we cap it and decline (always a safe result) when exceeded. */
#define RNB_BUDGET_SEC 12.0
static clock_t g_rnb_deadline = 0;
static bool rnb_over_budget(void) {
    return g_rnb_deadline != 0 && clock() > g_rnb_deadline;
}

/* ------------------------------------------------------------------ */
/* Small expression constructors (structural, no evaluation).          */
/* ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_plus2(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
}
static Expr* mk_times2(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}
static Expr* mk_pow(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Power), args, 2);
}
static Expr* mk_div(Expr* a, Expr* b) {
    return mk_times2(a, mk_pow(b, mk_int(-1)));
}
static Expr* mk_rat(int64_t num, int64_t den) {
    /* num/den as an unevaluated Times so it survives to evaluation; simplest
     * is Rational-like via Times[num, Power[den,-1]] which evaluates to the
     * reduced rational. */
    return mk_div(mk_int(num), mk_int(den));
}
static Expr* mk_unary(const char* head, Expr* arg) {
    Expr* args[1] = { arg };
    return expr_new_function(expr_new_symbol(head), args, 1);
}
static Expr* mk_binary(const char* head, Expr* a1, Expr* a2) {
    Expr* args[2] = { a1, a2 };
    return expr_new_function(expr_new_symbol(head), args, 2);
}

/* ------------------------------------------------------------------ */
/* CAS-evaluator wrappers (each consumes its Expr* arguments).         */
/* ------------------------------------------------------------------ */

static Expr* eval_together(Expr* f) { return eval_and_free(mk_unary("Together", f)); }

static Expr* eval_expand(Expr* f)   { return eval_and_free(mk_unary("Expand", f)); }
static Expr* eval_cancel(Expr* f)   { return eval_and_free(mk_unary("Cancel", f)); }
static Expr* eval_numer(Expr* f)    { return eval_and_free(mk_unary("Numerator", f)); }
static Expr* eval_denom(Expr* f)    { return eval_and_free(mk_unary("Denominator", f)); }
MATHILDA_MAYBE_UNUSED static Expr* eval_factor(Expr* f)   { return eval_and_free(mk_unary("Factor", f)); }

/* D[f, x] via the evaluator.  Consumes f; x is copied. */
static Expr* eval_diff(Expr* f, Expr* x) {
    return eval_and_free(mk_binary("D", f, expr_copy(x)));
}
/* PolynomialGCD[a, b]; consumes both. */
MATHILDA_MAYBE_UNUSED static Expr* eval_poly_gcd(Expr* a, Expr* b) {
    return eval_and_free(mk_binary("PolynomialGCD", a, b));
}
/* PolynomialLCM[a, b]; consumes both. */
static Expr* eval_poly_lcm(Expr* a, Expr* b) {
    return eval_and_free(mk_binary("PolynomialLCM", a, b));
}
/* Degree of p in x, as Length[CoefficientList[p, x]] - 1 (0 if free of x).
 * Consumes p; x is copied. */
static int64_t eval_degree(Expr* p, Expr* x) {
    Expr* clist = eval_and_free(mk_binary("CoefficientList", p, expr_copy(x)));
    int64_t d = 0;
    if (clist && clist->type == EXPR_FUNCTION
        && clist->data.function.head
        && clist->data.function.head->type == EXPR_SYMBOL
        && clist->data.function.head->data.symbol.name == SYM_List
        && clist->data.function.arg_count > 0) {
        d = (int64_t)clist->data.function.arg_count - 1;
    }
    expr_free(clist);
    return d;
}

/* True iff `sym_name` (an interned symbol name) does not appear in expr. */
static bool expr_free_of_symbol(const Expr* expr, const char* sym_name) {
    if (!expr) return true;
    if (expr->type == EXPR_SYMBOL) return expr->data.symbol.name != sym_name;
    if (expr->type != EXPR_FUNCTION) return true;
    if (!expr_free_of_symbol(expr->data.function.head, sym_name)) return false;
    for (size_t i = 0; i < expr->data.function.arg_count; i++)
        if (!expr_free_of_symbol(expr->data.function.args[i], sym_name)) return false;
    return true;
}

/* Integer floor of a/b for b>0 (a may be negative). */
static int64_t ifloordiv(int64_t a, int64_t b) {
    int64_t qd = a / b, r = a % b;
    if (r != 0 && ((r < 0) != (b < 0))) qd -= 1;
    return qd;
}
static int64_t igcd(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* ------------------------------------------------------------------ */
/* Radical detection: find the single radical of x in f.               */
/* Sets *out_q (owned polynomial base) and *out_m (root order) so that  */
/* y = (*out_q)^(1/m).  Returns true on a clean single radical.         */
/* ------------------------------------------------------------------ */

/* Recursively search e for Power[base, exp] with base not free of x and exp a
 * non-integer rational p/den; collect the distinct base and the lcm of the
 * exponent denominators.  Rejects (returns false) on two distinct radical
 * bases (a tower — out of scope for n=1 single radical). */
static bool rnb_scan_radical(const Expr* e, const char* xname,
                             Expr** base_seen, int64_t* mden) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power && n == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (base && exp && !expr_free_of_symbol(base, xname)
            && expr_free_of_symbol(exp, xname)) {
            /* exponent p/den with den > 1 ?  Rationals are Rational[p, den]. */
            int64_t den = 0;
            if (exp->type == EXPR_FUNCTION
                && exp->data.function.head
                && exp->data.function.head->type == EXPR_SYMBOL
                && exp->data.function.head->data.symbol.name == SYM_Rational
                && exp->data.function.arg_count == 2
                && exp->data.function.args[1]->type == EXPR_INTEGER) {
                den = exp->data.function.args[1]->data.integer;
            }
            if (den > 1) {
                if (*base_seen == NULL) {
                    *base_seen = expr_copy((Expr*)base);
                    *mden = den;
                } else {
                    if (!expr_eq(*base_seen, base)) return false; /* two radicals */
                    *mden = (*mden / igcd(*mden, den)) * den;      /* lcm */
                }
                /* still recurse into base in case of nesting */
            }
        }
    }
    if (!rnb_scan_radical(e->data.function.head, xname, base_seen, mden)) return false;
    for (size_t i = 0; i < n; i++)
        if (!rnb_scan_radical(e->data.function.args[i], xname, base_seen, mden))
            return false;
    return true;
}

static bool rnb_find_radical(const Expr* f, const char* xname,
                             Expr** out_q, int64_t* out_m) {
    Expr* base = NULL; int64_t m = 0;
    if (!rnb_scan_radical(f, xname, &base, &m) || base == NULL || m < 2) {
        if (base) expr_free(base);
        return false;
    }
    *out_q = base;
    *out_m = m;
    return true;
}

/* ================================================================== */
/* RadicalField: integral basis, multiplication table, derivation.     */
/* ================================================================== */

typedef struct {
    int    m;
    Expr*  x;            /* borrowed variable symbol */
    const char* xname;
    Expr** Q;            /* size m: Q[j] (owned), Integer 1 if absent (j in 1..m-1) */
    bool*  present;      /* size m: Q[j] non-constant */
    Expr** E;            /* size m: E_i (owned) */
    Expr** Lam;          /* size m: Lam_i (owned rational function or Integer 0) */
    int*   mul_idx;      /* size m*m: (i+k) mod m */
    Expr** mul_coef;     /* size m*m: basis-product coefficient (owned) */
    Expr*  Qs;           /* squarefree part prod_j Q_j (owned) */
    Expr*  lead;         /* c = q / prod_j Q_j^j (leading constant); y^m = q */
} RadicalField;

MATHILDA_MAYBE_UNUSED static bool rf_is_one(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 1;
}

/* base^e evaluated, e >= 0.  base borrowed. */
static Expr* rf_powi(const Expr* base, int64_t e) {
    if (e == 0) return mk_int(1);
    if (e == 1) return expr_copy((Expr*)base);
    return eval_and_free(mk_pow(expr_copy((Expr*)base), mk_int(e)));
}

static void rf_free(RadicalField* F) {
    if (!F) return;
    if (F->Q)   { for (int j = 0; j < F->m; j++) if (F->Q[j]) expr_free(F->Q[j]); free(F->Q); }
    if (F->E)   { for (int i = 0; i < F->m; i++) if (F->E[i]) expr_free(F->E[i]); free(F->E); }
    if (F->Lam) { for (int i = 0; i < F->m; i++) if (F->Lam[i]) expr_free(F->Lam[i]); free(F->Lam); }
    if (F->mul_coef) { for (int t = 0; t < F->m * F->m; t++) if (F->mul_coef[t]) expr_free(F->mul_coef[t]); free(F->mul_coef); }
    if (F->mul_idx)  free(F->mul_idx);
    if (F->present)  free(F->present);
    if (F->Qs)       expr_free(F->Qs);
    if (F->lead)     expr_free(F->lead);
    memset(F, 0, sizeof(*F));
}

/* Build the field for y^m = q.  Returns false (and leaves F zeroed) when a
 * normalisation fails: N1 (q not m-th-power-free) or N2 (gcd(m,{j})!=1). */
static bool rf_setup(RadicalField* F, const Expr* q, int m, Expr* x) {
    memset(F, 0, sizeof(*F));
    if (m < 2) return false;
    F->m = m; F->x = x; F->xname = x->data.symbol.name;
    F->Q = (Expr**)calloc((size_t)m, sizeof(Expr*));
    F->present = (bool*)calloc((size_t)m, sizeof(bool));
    F->E = (Expr**)calloc((size_t)m, sizeof(Expr*));
    F->Lam = (Expr**)calloc((size_t)m, sizeof(Expr*));
    F->mul_idx = (int*)calloc((size_t)m * m, sizeof(int));
    F->mul_coef = (Expr**)calloc((size_t)m * m, sizeof(Expr*));
    for (int j = 0; j < m; j++) F->Q[j] = mk_int(1);

    /* squarefree decomposition q = prod_i sq[i]^{i+1} (sq[i] at multiplicity i+1) */
    size_t cnt = 0;
    Expr** sq = risch_squarefree_t(q, x, &cnt);
    bool n1_fail = false;
    for (size_t i = 0; i < cnt; i++) {
        int jmult = (int)i + 1;
        Expr* fac = sq[i];
        if (!fac) continue;
        int64_t degf = eval_degree(expr_copy(fac), x);
        if (degf >= 1) {
            if (jmult >= m) { n1_fail = true; }   /* q not m-th-power-free */
            else {
                expr_free(F->Q[jmult]);
                F->Q[jmult] = expr_copy(fac);
                F->present[jmult] = true;
            }
        }
    }
    for (size_t i = 0; i < cnt; i++) if (sq[i]) expr_free(sq[i]);
    free(sq);
    if (n1_fail) { rf_free(F); return false; }

    /* N2: gcd(m, {j : Q_j != 1}) == 1, and at least one Q_j present. */
    int64_t g = m; bool any = false;
    for (int j = 1; j < m; j++) if (F->present[j]) { g = igcd(g, j); any = true; }
    if (!any || g != 1) { rf_free(F); return false; }

    /* squarefree part Qs = prod_j Q_j */
    Expr* qs = mk_int(1);
    for (int j = 1; j < m; j++) if (F->present[j]) qs = mk_times2(qs, expr_copy(F->Q[j]));
    F->Qs = eval_expand(qs);

    /* leading constant c: risch_squarefree_t returns MONIC Q_j, so the actual
     * radicand is q = c * prod_j Q_j^j with c = q / prod_j Q_j^j a nonzero
     * constant.  y^m = q, so each m-th-power reduction y^m -> q contributes c;
     * this must ride the multiplication table (a sign for q = 1 - x^2). */
    Expr* prodQj = mk_int(1);
    for (int j = 1; j < m; j++)
        if (F->present[j]) prodQj = mk_times2(prodQj, rf_powi(F->Q[j], (int64_t)j));
    F->lead = eval_cancel(mk_div(expr_copy((Expr*)q), eval_expand(prodQj)));

    /* E_i = prod_j Q_j^floor(ij/m) */
    for (int i = 0; i < m; i++) {
        Expr* e = mk_int(1);
        for (int j = 1; j < m; j++) if (F->present[j]) {
            int64_t ex = ifloordiv((int64_t)i * j, m);
            if (ex > 0) e = mk_times2(e, rf_powi(F->Q[j], ex));
        }
        F->E[i] = eval_expand(e);
    }

    /* multiplication table */
    for (int i = 0; i < m; i++)
        for (int k = 0; k < m; k++) {
            F->mul_idx[i * m + k] = (i + k) % m;
            Expr* c = mk_int(1);
            for (int j = 1; j < m; j++) if (F->present[j]) {
                int64_t ex = ifloordiv((int64_t)(i + k) * j, m)
                           - ifloordiv((int64_t)i * j, m)
                           - ifloordiv((int64_t)k * j, m);
                if (ex > 0) c = mk_times2(c, rf_powi(F->Q[j], ex));
            }
            /* each wraparound i+k >= m applies one m-th-power reduction -> * lead */
            int64_t wraps = (i + k) / m;   /* 0 or 1 for i,k in [0,m) */
            for (int64_t w = 0; w < wraps; w++) c = mk_times2(c, expr_copy(F->lead));
            F->mul_coef[i * m + k] = eval_expand(c);
        }

    /* Lam_i = sum_j {ij/m} Q_j'/Q_j = sum_j ((ij mod m)/m) Q_j'/Q_j */
    for (int i = 0; i < m; i++) {
        Expr* lam = mk_int(0);
        for (int j = 1; j < m; j++) if (F->present[j]) {
            int r = (int)(((int64_t)i * j) % m);
            if (r == 0) continue;
            Expr* dq = eval_diff(expr_copy(F->Q[j]), x);
            Expr* term = mk_times2(mk_rat(r, m), mk_div(dq, expr_copy(F->Q[j])));
            lam = mk_plus2(lam, term);
        }
        F->Lam[i] = eval_cancel(eval_together(lam));
    }
    return true;
}

/* ---- elements: Expr* a[m], meaning sum_i a_i w_i --------------------- */

static Expr** elem_zero(int m) {
    size_t mm = (m < 1) ? 1u : (size_t)m;
    Expr** a = (Expr**)calloc(mm, sizeof(Expr*));
    for (int i = 0; i < m; i++) a[i] = mk_int(0);
    return a;
}
MATHILDA_MAYBE_UNUSED static Expr** elem_one(int m) {
    Expr** a = elem_zero(m);
    expr_free(a[0]); a[0] = mk_int(1);
    return a;
}
MATHILDA_MAYBE_UNUSED static Expr** elem_copy(Expr* const* a, int m) {
    Expr** b = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++) b[i] = expr_copy(a[i]);
    return b;
}
static void elem_free(Expr** a, int m) {
    if (!a) return;
    for (int i = 0; i < m; i++) if (a[i]) expr_free(a[i]);
    free(a);
}
MATHILDA_MAYBE_UNUSED static bool elem_is_zero(Expr* const* a, int m) {
    for (int i = 0; i < m; i++)
        if (!(a[i]->type == EXPR_INTEGER && a[i]->data.integer == 0)) return false;
    return true;
}

static Expr** rf_add(const RadicalField* F, Expr* const* a, Expr* const* b) {
    int m = F->m;
    Expr** r = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++)
        r[i] = eval_cancel(mk_plus2(expr_copy(a[i]), expr_copy(b[i])));
    return r;
}
static Expr** rf_scal(const RadicalField* F, const Expr* c, Expr* const* a) {
    int m = F->m;
    Expr** r = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++)
        r[i] = eval_cancel(mk_times2(expr_copy((Expr*)c), expr_copy(a[i])));
    return r;
}
static Expr** rf_mult(const RadicalField* F, Expr* const* a, Expr* const* b) {
    int m = F->m;
    Expr** r = elem_zero(m);
    for (int i = 0; i < m; i++) {
        if (a[i]->type == EXPR_INTEGER && a[i]->data.integer == 0) continue;
        for (int k = 0; k < m; k++) {
            if (b[k]->type == EXPR_INTEGER && b[k]->data.integer == 0) continue;
            int idx = F->mul_idx[i * m + k];
            Expr* prod = mk_times2(mk_times2(expr_copy(a[i]), expr_copy(b[k])),
                                   expr_copy(F->mul_coef[i * m + k]));
            Expr* acc = mk_plus2(r[idx], prod);   /* consumes old r[idx] */
            r[idx] = acc;
        }
    }
    for (int i = 0; i < m; i++) r[i] = eval_cancel(r[i]);
    return r;
}
static Expr** rf_D(const RadicalField* F, Expr* const* a) {
    int m = F->m;
    Expr** r = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++) {
        Expr* da = eval_diff(expr_copy(a[i]), F->x);
        Expr* lam_term = mk_times2(expr_copy(a[i]), expr_copy(F->Lam[i]));
        r[i] = eval_cancel(mk_plus2(da, lam_term));
    }
    return r;
}

/* matrix of multiplication-by-a in basis w (columns = a*w_k): List of rows. */
static Expr* rf_mult_matrix(const RadicalField* F, Expr* const* a) {
    int m = F->m;
    /* columns col_k = a * e_k */
    Expr*** cols = (Expr***)calloc((size_t)m, sizeof(Expr**));
    for (int k = 0; k < m; k++) {
        Expr** ek = elem_zero(m);
        expr_free(ek[k]); ek[k] = mk_int(1);
        cols[k] = rf_mult(F, a, ek);
        elem_free(ek, m);
    }
    Expr** rows = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++) {
        Expr** entries = (Expr**)calloc((size_t)m, sizeof(Expr*));
        for (int k = 0; k < m; k++) entries[k] = expr_copy(cols[k][i]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), entries, (size_t)m);
        free(entries);
    }
    for (int k = 0; k < m; k++) elem_free(cols[k], m);
    free(cols);
    Expr* M = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)m);
    free(rows);
    return M;
}
static Expr* rf_norm(const RadicalField* F, Expr* const* a) {
    Expr* M = rf_mult_matrix(F, a);
    return eval_cancel(mk_unary("Det", M));
}
/* a^{-1} as an element, via LinearSolve[mult_matrix(a), e0]. */
static Expr** rf_inv(const RadicalField* F, Expr* const* a) {
    int m = F->m;
    Expr* M = rf_mult_matrix(F, a);
    Expr** e0e = (Expr**)calloc((size_t)m, sizeof(Expr*));
    e0e[0] = mk_int(1);
    for (int i = 1; i < m; i++) e0e[i] = mk_int(0);
    Expr* e0 = expr_new_function(expr_new_symbol(SYM_List), e0e, (size_t)m);
    free(e0e);
    Expr* sol = eval_and_free(mk_binary("LinearSolve", M, e0));
    Expr** r = (Expr**)calloc((size_t)m, sizeof(Expr*));
    if (sol && sol->type == EXPR_FUNCTION
        && sol->data.function.head->type == EXPR_SYMBOL
        && sol->data.function.head->data.symbol.name == SYM_List
        && (int)sol->data.function.arg_count == m) {
        for (int i = 0; i < m; i++)
            r[i] = eval_cancel(expr_copy(sol->data.function.args[i]));
    } else {
        for (int i = 0; i < m; i++) r[i] = mk_int(0);  /* singular / failure */
    }
    if (sol) expr_free(sol);
    return r;
}

/* Render element a back into an expression in x and y = q^(1/m):
 * sum_i a_i * q^(i/m) / E_i.  q borrowed. */
static Expr* rf_to_y_expr(const RadicalField* F, Expr* const* a, const Expr* q) {
    int m = F->m;
    Expr* sum = mk_int(0);
    for (int i = 0; i < m; i++) {
        if (a[i]->type == EXPR_INTEGER && a[i]->data.integer == 0) continue;
        Expr* yi = (i == 0) ? mk_int(1)
                            : mk_pow(expr_copy((Expr*)q), mk_rat(i, m));
        Expr* term = mk_div(mk_times2(expr_copy(a[i]), yi), expr_copy(F->E[i]));
        sum = mk_plus2(sum, term);
    }
    return eval_and_free(sum);
}

/* ================================================================== */
/* Debug surface: Integrate`RNB`Info[f, x] — inspect the field.        */
/* ================================================================== */

static Expr* builtin_rnb_info(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* q = NULL; int64_t m = 0;
    if (!rnb_find_radical(f, x->data.symbol.name, &q, &m)) return NULL;

    RadicalField F;
    if (!rf_setup(&F, q, (int)m, x)) { expr_free(q); return NULL; }

    /* Build a descriptive List: {m, q, {Q_j...}, {E_i...}, {Lam_i...},
     * Norm[w_1]} — field data only; no integrand element needed here. */
    Expr** parts = (Expr**)calloc(6, sizeof(Expr*));
    parts[0] = mk_int((int64_t)m);
    parts[1] = expr_copy(q);
    /* Q_j list (present ones) */
    {
        Expr** qs = (Expr**)calloc((size_t)m, sizeof(Expr*)); size_t nq = 0;
        for (int j = 1; j < F.m; j++) if (F.present[j]) qs[nq++] = expr_copy(F.Q[j]);
        parts[2] = expr_new_function(expr_new_symbol(SYM_List), qs, nq);
        free(qs);
    }
    { Expr** es = (Expr**)calloc((size_t)m, sizeof(Expr*));
      for (int i = 0; i < F.m; i++) es[i] = expr_copy(F.E[i]);
      parts[3] = expr_new_function(expr_new_symbol(SYM_List), es, (size_t)m); free(es); }
    { Expr** ls = (Expr**)calloc((size_t)m, sizeof(Expr*));
      for (int i = 0; i < F.m; i++) ls[i] = expr_copy(F.Lam[i]);
      parts[4] = expr_new_function(expr_new_symbol(SYM_List), ls, (size_t)m); free(ls); }
    /* norm(w_1): element w_1 = (0,1,0,...) */
    { Expr** w1 = elem_zero(F.m); expr_free(w1[1 % F.m]); w1[1 % F.m] = mk_int(1);
      parts[5] = rf_norm(&F, w1); elem_free(w1, F.m); }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), parts, 6);
    free(parts);
    rf_free(&F);
    expr_free(q);
    return out;
}

/* ================================================================== */
/* Converting the integrand f into an element of O.                    */
/* ================================================================== */

/* True iff e's head is the interned symbol `name`. */
static bool head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Structurally rewrite every radical power Power[q, p/den] (den>1, den|m) of
 * the field base q into Power[Y, p*m/den], leaving all else intact.  Returns
 * a fresh tree; q, e borrowed. */
static Expr* rnb_sub_radical_to_Y(const Expr* e, const Expr* q, int m,
                                  const char* Yname) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    if (h && h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Power && n == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (expr_eq((Expr*)base, (Expr*)q)
            && exp->type == EXPR_FUNCTION
            && exp->data.function.head->type == EXPR_SYMBOL
            && exp->data.function.head->data.symbol.name == SYM_Rational
            && exp->data.function.arg_count == 2
            && exp->data.function.args[0]->type == EXPR_INTEGER
            && exp->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t p   = exp->data.function.args[0]->data.integer;
            int64_t den = exp->data.function.args[1]->data.integer;
            if (den > 1 && (m % den) == 0) {
                return mk_pow(expr_new_symbol(Yname), mk_int(p * (m / den)));
            }
        }
    }
    Expr* newhead = rnb_sub_radical_to_Y(h, q, m, Yname);
    Expr** newargs = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        newargs[i] = rnb_sub_radical_to_Y(e->data.function.args[i], q, m, Yname);
    Expr* r = expr_new_function(newhead, newargs, n);
    free(newargs);
    return r;
}

/* Convert a polynomial P in Y (rational-function coeffs in x) into an element:
 * Y^k = q^(k div m) E_{k mod m} w_{k mod m}.  P, q borrowed. */
static Expr** rf_from_y_poly(const RadicalField* F, const Expr* P,
                             const Expr* q, const char* Yname) {
    int m = F->m;
    Expr* Y = expr_new_symbol(Yname);
    Expr* clist = eval_and_free(mk_binary("CoefficientList", expr_copy((Expr*)P), Y));
    Expr** r = elem_zero(m);
    if (clist && head_is(clist, SYM_List)) {
        size_t nk = clist->data.function.arg_count;
        for (size_t k = 0; k < nk; k++) {
            Expr* ck = clist->data.function.args[k];
            if (ck->type == EXPR_INTEGER && ck->data.integer == 0) continue;
            int i = (int)(k % (size_t)m);
            int64_t qp = (int64_t)(k / (size_t)m);
            Expr* term = expr_copy(ck);
            if (qp > 0) term = mk_times2(term, rf_powi(q, qp));
            term = mk_times2(term, expr_copy(F->E[i]));
            r[i] = mk_plus2(r[i], term);
        }
        for (int i = 0; i < m; i++) r[i] = eval_cancel(r[i]);
    }
    if (clist) expr_free(clist);
    return r;
}

/* f (Expr, rational in x and the radical) -> element of O.  q borrowed. */
static Expr** rnb_f_to_element(const RadicalField* F, const Expr* f,
                               const Expr* q, const char* Yname) {
    int m = F->m;
    Expr* fsub = rnb_sub_radical_to_Y(f, q, m, Yname);
    Expr* fto  = eval_together(fsub);
    Expr* N    = eval_numer(expr_copy(fto));
    Expr* Dd   = eval_denom(fto);
    Expr** eN   = rf_from_y_poly(F, N, q, Yname);
    Expr** eD   = rf_from_y_poly(F, Dd, q, Yname);
    Expr** eDi  = rf_inv(F, eD);
    Expr** fe   = rf_mult(F, eN, eDi);
    expr_free(N); expr_free(Dd);
    elem_free(eN, m); elem_free(eD, m); elem_free(eDi, m);
    return fe;
}

/* Common denominator (a polynomial in x) of an element's coordinates. */
static Expr* rf_denominator_poly(const RadicalField* F, Expr* const* a) {
    int m = F->m;
    Expr* d = mk_int(1);
    for (int i = 0; i < m; i++) {
        Expr* di = eval_denom(eval_together(expr_copy(a[i])));
        d = eval_poly_lcm(d, di);
    }
    return d;
}

/* ================================================================== */
/* Heuristic parallel integration (given a logand list).               */
/* ================================================================== */

static uint64_t g_rnb_unk = 0;

/* --- algebraic-constant abstraction ---------------------------------------
 * Mathilda's Cancel/Together cancel a common SYMBOL factor but not a common
 * algebraic-constant factor (Sqrt[3/2]/(Sqrt[3/2](x^2-2)) stays uncancelled;
 * only the slow FullSimplify reduces it).  So before clearing denominators we
 * replace each distinct algebraic constant (Power[number, non-integer p/q],
 * e.g. Sqrt[2], (-1)^(1/3)) by a fresh symbol, which Cancel treats as an
 * ordinary variable, then substitute the constants back into the extracted
 * coefficients. */
static bool rnb_num_atom(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT
        || (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
            && e->data.function.head->data.symbol.name == SYM_Rational));
}
static bool rnb_is_alg_const(const Expr* e) {
    if (!head_is(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    return rnb_num_atom(e->data.function.args[0])
        && head_is(e->data.function.args[1], SYM_Rational);
}
static void rnb_collect_alg(const Expr* e, Expr*** ks, size_t* nk) {
    if (!e) return;
    if (rnb_is_alg_const(e)) {
        for (size_t i = 0; i < *nk; i++) if (expr_eq((*ks)[i], (Expr*)e)) return;
        *ks = (Expr**)realloc(*ks, (*nk + 1) * sizeof(Expr*));
        (*ks)[(*nk)++] = expr_copy((Expr*)e);
        return;   /* atomic: don't recurse into the radical */
    }
    if (e->type == EXPR_FUNCTION) {
        rnb_collect_alg(e->data.function.head, ks, nk);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rnb_collect_alg(e->data.function.args[i], ks, nk);
    }
}

/* Exact per-coordinate degree bounds for the numerator b_i (Cor. 4.7 /
 * bounds.py): B[i] = floor((k + e*Delta + nu_i)/e), computed from valuations
 * at the places at infinity.  Returns an owned int64_t[m]; d_in is den(f). */
static int64_t* rnb_exact_degree_bounds(const RadicalField* F, Expr* const* f_elem,
                                        const Expr* q, const Expr* d_in) {
    int m = F->m;
    int64_t N = eval_degree(expr_copy((Expr*)q), F->x);
    int64_t g = igcd(m, N); if (g < 1) g = 1;
    int64_t e = m / g, Ng = N / g;
    int64_t* degQ = (int64_t*)calloc((size_t)m, sizeof(int64_t));
    for (int j = 1; j < m; j++)
        if (F->present[j]) degQ[j] = eval_degree(expr_copy(F->Q[j]), F->x);
    int64_t* nu = (int64_t*)calloc((size_t)m, sizeof(int64_t));
    for (int i = 0; i < m; i++) {
        int64_t s = 0;
        for (int j = 1; j < m; j++)
            if (F->present[j]) s += ifloordiv((int64_t)i * j, m) * degQ[j];
        nu[i] = -(int64_t)i * Ng + e * s;
    }
    /* vf = min over nonzero coords of -e*deg(num_i)+e*deg(den_i)+nu_i */
    int64_t vf = 0; bool any = false;
    for (int i = 0; i < m; i++) {
        if (f_elem[i]->type == EXPR_INTEGER && f_elem[i]->data.integer == 0) continue;
        Expr* fi = eval_together(expr_copy(f_elem[i]));
        Expr* ni = eval_numer(expr_copy(fi));
        Expr* di = eval_denom(fi);
        int64_t val = -e * eval_degree(ni, F->x) + e * eval_degree(di, F->x) + nu[i];
        if (!any || val < vf) vf = val;
        any = true;
    }
    /* Delta = sum (j-1) deg(fac_j) over squarefree(d) */
    size_t sc = 0; Expr** sq = risch_squarefree_t(d_in, F->x, &sc);
    int64_t Delta = 0;
    for (size_t i = 0; i < sc; i++)
        if (sq[i]) Delta += (int64_t)i * eval_degree(expr_copy(sq[i]), F->x);
    for (size_t i = 0; i < sc; i++) if (sq[i]) expr_free(sq[i]);
    free(sq);
    int64_t k = e - vf; if (k < 0) k = 0;
    int64_t* B = (int64_t*)calloc((size_t)m, sizeof(int64_t));
    for (int i = 0; i < m; i++) B[i] = ifloordiv(k + e * Delta + nu[i], e);
    free(degQ); free(nu);
    return B;
}

/* Verification: substituting the solution (rules) and the zero-out (zlist)
 * into every extracted equation coefficient must give an algebraic zero.  The
 * coefficients are pure in the unknowns and algebraic CONSTANTS (no x), so
 * RootReduce is an exact, fast zero test — and unlike the residual it handles
 * products of the same radical (RootReduce knows Sqrt[3/2]^2 = 3/2). */
static bool rnb_eqs_satisfied(Expr* const* eqs, size_t neq, Expr* rules, Expr* zlist) {
    for (size_t t = 0; t < neq; t++) {
        Expr* lhs = eqs[t]->data.function.args[0];   /* the coefficient */
        Expr* s = eval_and_free(mk_binary("ReplaceAll", expr_copy(lhs), expr_copy(rules)));
        s = eval_and_free(mk_binary("ReplaceAll", s, expr_copy(zlist)));
        Expr* z = eval_and_free(mk_unary("RootReduce", s));
        bool zero = (z && z->type == EXPR_INTEGER && z->data.integer == 0);
        if (z) expr_free(z);
        if (!zero) return false;
    }
    return true;
}

/* Attempt v + sum_j gamma_j Log[u_j] with unknown numerator coeffs and gammas.
 * Returns the antiderivative expression (in x and y) on success, else NULL.
 * f_elem is the integrand as an element; logs[] are logand elements. */
static Expr* rnb_parallel_integrate(const RadicalField* F, Expr* const* f_elem,
                                    const Expr* q, const char* Yname,
                                    Expr*** logs, size_t nlog) {
    int m = F->m;
    (void)Yname;

    /* denominator of f and the Hermite ansatz denominator den_v */
    Expr* d = rf_denominator_poly(F, f_elem);
    size_t scnt = 0;
    Expr** sq = risch_squarefree_t(d, F->x, &scnt);
    Expr* den_v = mk_int(1);
    for (size_t i = 0; i < scnt; i++) {
        if (sq[i] && (int)i >= 1) den_v = mk_times2(den_v, rf_powi(sq[i], (int64_t)i));
    }
    den_v = eval_expand(den_v);
    for (size_t i = 0; i < scnt; i++) if (sq[i]) expr_free(sq[i]);
    free(sq);

    /* exact per-coordinate degree bounds B[i] (Cor. 4.7) */
    int64_t* B = rnb_exact_degree_bounds(F, f_elem, q, d);
    expr_free(d);

    /* unknown numerator coefficients mu[i][k] and log coeffs gamma[j] */
    Expr*** v = (Expr***)calloc(1, sizeof(Expr**));  /* v element */
    Expr** vel = (Expr**)calloc((size_t)m, sizeof(Expr*));
    /* collect all unknown symbols for zero-out */
    size_t max_unk = nlog + 4;
    for (int i = 0; i < m; i++) if (B[i] >= 0) max_unk += (size_t)(B[i] + 1);
    Expr** unknowns = (Expr**)calloc(max_unk, sizeof(Expr*));
    size_t nunk = 0;
    uint64_t base_id = g_rnb_unk;
    g_rnb_unk += max_unk;

    for (int i = 0; i < m; i++) {
        Expr* bi = mk_int(0);
        for (int64_t k = 0; k <= B[i]; k++) {
            char buf[48];
            snprintf(buf, sizeof(buf), "RNB$mu$%llu$%d$%lld",
                     (unsigned long long)base_id, i, (long long)k);
            Expr* mu = expr_new_symbol(intern_symbol(buf));
            unknowns[nunk++] = expr_copy(mu);
            Expr* term = mk_times2(mu, mk_pow(expr_copy(F->x), mk_int(k)));
            bi = mk_plus2(bi, term);
        }
        vel[i] = eval_cancel(mk_div(bi, expr_copy(den_v)));
    }
    v[0] = vel;
    free(B);

    Expr** gamma = (Expr**)calloc(nlog ? nlog : 1, sizeof(Expr*));
    for (size_t j = 0; j < nlog; j++) {
        char buf[48];
        snprintf(buf, sizeof(buf), "RNB$g$%llu$%zu", (unsigned long long)base_id, j);
        gamma[j] = expr_new_symbol(intern_symbol(buf));
        unknowns[nunk++] = expr_copy(gamma[j]);
    }

    /* Dg = D(v) + sum_j gamma_j D(u_j) u_j^{-1} */
    Expr** Dg = rf_D(F, vel);
    for (size_t j = 0; j < nlog; j++) {
        Expr** Du   = rf_D(F, logs[j]);
        Expr** uinv = rf_inv(F, logs[j]);
        Expr** dlog = rf_mult(F, Du, uinv);
        Expr** sc   = rf_scal(F, gamma[j], dlog);
        Expr** ndg  = rf_add(F, Dg, sc);
        elem_free(Dg, m); Dg = ndg;
        elem_free(Du, m); elem_free(uinv, m); elem_free(dlog, m); elem_free(sc, m);
    }

    /* residual f - Dg, then one equation per x-power per coordinate */
    Expr* negone = mk_int(-1);
    Expr** nDg = rf_scal(F, negone, Dg);
    expr_free(negone);
    Expr** resid = rf_add(F, f_elem, nDg);
    elem_free(nDg, m); elem_free(Dg, m);

    /* abstract algebraic constants -> fresh symbols so Cancel clears them */
    Expr** ks = NULL; size_t nk = 0;
    for (int i = 0; i < m; i++) rnb_collect_alg(resid[i], &ks, &nk);
    Expr** fwd = (Expr**)malloc((nk ? nk : 1) * sizeof(Expr*));
    Expr** bwd = (Expr**)malloc((nk ? nk : 1) * sizeof(Expr*));
    for (size_t t = 0; t < nk; t++) {
        char buf[32]; snprintf(buf, sizeof(buf), "RNB$a$%zu", t);
        Expr* s = expr_new_symbol(intern_symbol(buf));
        fwd[t] = mk_binary("Rule", expr_copy(ks[t]), expr_copy(s));
        bwd[t] = mk_binary("Rule", s, expr_copy(ks[t]));
    }
    Expr* fwd_rules = expr_new_function(expr_new_symbol(SYM_List), fwd, nk); free(fwd);
    Expr* bwd_rules = expr_new_function(expr_new_symbol(SYM_List), bwd, nk); free(bwd);

    size_t neq = 0, capeq = 0;
    Expr** eqs = NULL;
    for (int i = 0; i < m; i++) {
        Expr* ri = (nk > 0)
            ? eval_and_free(mk_binary("ReplaceAll", expr_copy(resid[i]), expr_copy(fwd_rules)))
            : expr_copy(resid[i]);
        Expr* num = eval_numer(eval_cancel(eval_together(ri)));
        Expr* clist = eval_and_free(mk_binary("CoefficientList",
                                    eval_expand(num), expr_copy(F->x)));
        if (clist && head_is(clist, SYM_List)) {
            size_t nc = clist->data.function.arg_count;
            for (size_t c = 0; c < nc; c++) {
                Expr* coef = clist->data.function.args[c];
                if (coef->type == EXPR_INTEGER && coef->data.integer == 0) continue;
                /* substitute the algebraic constants back */
                Expr* ec = (nk > 0)
                    ? eval_and_free(mk_binary("ReplaceAll", expr_copy(coef),
                                              expr_copy(bwd_rules)))
                    : expr_copy(coef);
                if (neq == capeq) {
                    capeq = capeq ? capeq * 2 : 16;
                    eqs = (Expr**)realloc(eqs, capeq * sizeof(Expr*));
                }
                eqs[neq++] = mk_binary("Equal", ec, mk_int(0));
            }
        }
        if (clist) expr_free(clist);
    }
    expr_free(fwd_rules); expr_free(bwd_rules);
    for (size_t t = 0; t < nk; t++) expr_free(ks[t]);
    free(ks);

    if (getenv("RNB_DEBUG")) {
        fprintf(stderr, "[rnb] nlog=%zu nunk=%zu neq=%zu\n", nlog, nunk, neq);
        for (size_t t = 0; t < neq && t < 40; t++) {
            char* s = expr_to_string(eqs[t]);
            fprintf(stderr, "  eq[%zu]: %s\n", t, s ? s : "?"); free(s);
        }
        for (size_t t = 0; t < nunk; t++) {
            char* s = expr_to_string(unknowns[t]);
            fprintf(stderr, "  unk[%zu]: %s\n", t, s ? s : "?"); free(s);
        }
    }

    /* candidate antiderivative G = to_y(v) + sum_j gamma_j Log[to_y(u_j)] */
    Expr* G = rf_to_y_expr(F, vel, q);
    for (size_t j = 0; j < nlog; j++) {
        Expr* uarg = rf_to_y_expr(F, logs[j], q);
        Expr* logu = mk_unary("Log", uarg);
        G = mk_plus2(G, mk_times2(expr_copy(gamma[j]), logu));
    }

    /* Solve the linear system */
    Expr* result = NULL;
    if (neq > 0 && !rnb_over_budget()) {
        Expr** eqcopy = (Expr**)malloc(neq * sizeof(Expr*));
        for (size_t t = 0; t < neq; t++) eqcopy[t] = expr_copy(eqs[t]);
        Expr* eqlist = expr_new_function(expr_new_symbol(SYM_List), eqcopy, neq);
        free(eqcopy);
        Expr** unkcopy = (Expr**)malloc((nunk ? nunk : 1) * sizeof(Expr*));
        for (size_t t = 0; t < nunk; t++) unkcopy[t] = expr_copy(unknowns[t]);
        Expr* unklist = expr_new_function(expr_new_symbol(SYM_List), unkcopy, nunk);
        free(unkcopy);
        /* An underdetermined ansatz legitimately leaves free unknowns
         * (Solve::svars, harmless here); we pin them to 0 below. */
        Expr* sol = eval_and_free(mk_binary("Solve", eqlist, unklist));
        if (getenv("RNB_DEBUG")) {
            char* s = expr_to_string(sol);
            fprintf(stderr, "[rnb] sol: %s\n", s ? s : "NULL"); free(s);
        }
        if (sol && head_is(sol, SYM_List) && sol->data.function.arg_count >= 1) {
            Expr* rules = sol->data.function.args[0];  /* first solution */
            /* zero-out rules for every unknown */
            Expr** zr = (Expr**)malloc((nunk ? nunk : 1) * sizeof(Expr*));
            for (size_t t = 0; t < nunk; t++)
                zr[t] = mk_binary("Rule", expr_copy(unknowns[t]), mk_int(0));
            Expr* zlist = expr_new_function(expr_new_symbol(SYM_List), zr, nunk);
            free(zr);
            /* verify the solution satisfies every equation before committing */
            if (rnb_eqs_satisfied(eqs, neq, rules, zlist)) {
                Expr* Gs = eval_and_free(mk_binary("ReplaceAll", expr_copy(G),
                                                   expr_copy(rules)));
                Expr* Gf = eval_and_free(mk_binary("ReplaceAll", Gs, expr_copy(zlist)));
                result = Gf;
            }
            expr_free(zlist);
        }
        if (sol) expr_free(sol);
    }
    elem_free(resid, m);
    for (size_t t = 0; t < neq; t++) expr_free(eqs[t]);
    free(eqs);

    for (size_t t = 0; t < nunk; t++) expr_free(unknowns[t]);
    free(unknowns);
    free(gamma);
    expr_free(G);
    expr_free(den_v);
    elem_free(vel, m);
    free(v);
    return result;
}

/* ================================================================== */
/* Exact tier (I): units at infinity for m=2 (polynomial Pell /        */
/* continued fraction of sqrt(q), Prop 4.9).                           */
/* ================================================================== */

#define RNB_CF_MAXSTEPS 64

/* Polynomial part at x=infinity of a rational-times-radical expression:
 * Normal[Series[e, {x, Infinity, 0}]] keeps exactly the x^k, k>=0 terms.
 * e is consumed. */
static Expr* rnb_polypart(const RadicalField* F, Expr* e) {
    Expr* spec_args[3] = { expr_copy(F->x), expr_new_symbol("Infinity"), mk_int(0) };
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 3);
    Expr* ser  = eval_and_free(mk_binary("Series", e, spec));
    Expr* norm = eval_and_free(mk_unary("Normal", ser));
    return eval_expand(norm);
}

/* Fundamental S-unit of O at the places at infinity for m=2 (deg q even):
 * the convergent {p, q} = p + q y from the first period of the continued
 * fraction of sqrt(q).  Returns an element {p,q} or NULL if none within the
 * step bound (odd deg q => O^x = F^*).  q borrowed. */
static Expr** rnb_cf_unit(const RadicalField* F, const Expr* q) {
    if (F->m != 2) return NULL;
    int64_t n = eval_degree(expr_copy((Expr*)q), F->x);
    if (n <= 0 || (n % 2) != 0) return NULL;

    Expr* P  = mk_int(0);
    Expr* Qd = mk_int(1);
    Expr* a0 = rnb_polypart(F, mk_pow(expr_copy((Expr*)q), mk_rat(1, 2)));
    Expr* p_prev = mk_int(1), *p_cur = expr_copy(a0);
    Expr* q_prev = mk_int(0), *q_cur = mk_int(1);
    Expr* alast  = a0;
    Expr** result = NULL;

    for (int step = 0; step < RNB_CF_MAXSTEPS; step++) {
        /* P <- alast*Qd - P */
        Expr* newP = eval_expand(mk_plus2(mk_times2(expr_copy(alast), expr_copy(Qd)),
                                          mk_times2(mk_int(-1), expr_copy(P))));
        expr_free(P); P = newP;
        /* Qd <- (q - P^2)/Qd  (exact) */
        Expr* num = mk_plus2(expr_copy((Expr*)q),
                             mk_times2(mk_int(-1), mk_pow(expr_copy(P), mk_int(2))));
        Expr* newQd = eval_cancel(mk_div(num, expr_copy(Qd)));
        expr_free(Qd); Qd = newQd;

        int64_t dQ = eval_degree(expr_copy(Qd), F->x);
        if (dQ == 0) {
            result = elem_zero(2);
            expr_free(result[0]); result[0] = expr_copy(p_cur);
            expr_free(result[1]); result[1] = expr_copy(q_cur);
            break;
        }
        /* an = polypart((P + sqrt q)/Qd) */
        Expr* pe = mk_div(mk_plus2(expr_copy(P),
                                   mk_pow(expr_copy((Expr*)q), mk_rat(1, 2))),
                          expr_copy(Qd));
        Expr* an = rnb_polypart(F, pe);
        /* convergents */
        Expr* np = eval_expand(mk_plus2(mk_times2(expr_copy(an), expr_copy(p_cur)),
                                        expr_copy(p_prev)));
        Expr* nq = eval_expand(mk_plus2(mk_times2(expr_copy(an), expr_copy(q_cur)),
                                        expr_copy(q_prev)));
        expr_free(p_prev); p_prev = p_cur; p_cur = np;
        expr_free(q_prev); q_prev = q_cur; q_cur = nq;
        expr_free(alast); alast = an;
    }

    expr_free(P); expr_free(Qd);
    expr_free(p_prev); expr_free(p_cur);
    expr_free(q_prev); expr_free(q_cur);
    expr_free(alast);
    return result;
}

/* ================================================================== */
/* Exact tier (II): places over the denominator, residues.             */
/* ================================================================== */

#define RNB_T "RNB$t"

/* ReplaceAll[e, var -> val]; e, var, val borrowed. */
static Expr* rnb_subst(const Expr* e, const Expr* var, const Expr* val) {
    Expr* rule = mk_binary("Rule", expr_copy((Expr*)var), expr_copy((Expr*)val));
    return eval_and_free(mk_binary("ReplaceAll", expr_copy((Expr*)e), rule));
}

/* Residue[expr, {var, at}].  expr consumed; var/at borrowed. */
static Expr* rnb_residue2(Expr* expr, const char* varname, const Expr* at) {
    Expr* spec_args[2] = { expr_new_symbol(varname), expr_copy((Expr*)at) };
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 2);
    return eval_and_free(mk_binary("Residue", expr, spec));
}

/* Roots of poly in `varname` (as algebraic exprs), via Solve.  Appends owned
 * root exprs to *out (realloc'd). */
static void rnb_poly_roots(const Expr* poly, const char* varname,
                           Expr*** out, size_t* nout) {
    Expr* var = expr_new_symbol(varname);
    Expr* eq  = mk_binary("Equal", expr_copy((Expr*)poly), mk_int(0));
    Expr* sol = eval_and_free(mk_binary("Solve", eq, var));
    if (sol && head_is(sol, SYM_List)) {
        for (size_t s = 0; s < sol->data.function.arg_count; s++) {
            Expr* soln = sol->data.function.args[s];   /* {var->r} */
            if (!head_is(soln, SYM_List)) continue;
            for (size_t rct = 0; rct < soln->data.function.arg_count; rct++) {
                Expr* rule = soln->data.function.args[rct];
                if (head_is(rule, SYM_Rule) && rule->data.function.arg_count == 2) {
                    *out = (Expr**)realloc(*out, (*nout + 1) * sizeof(Expr*));
                    (*out)[(*nout)++] = expr_copy(rule->data.function.args[1]);
                }
            }
        }
    }
    if (sol) expr_free(sol);
}

/* PolynomialRemainder[num, den, x] == 0 ? */
static bool rnb_divides(const Expr* den, const Expr* num, Expr* x) {
    Expr* args[3] = { expr_copy((Expr*)num), expr_copy((Expr*)den), expr_copy(x) };
    Expr* rem = eval_and_free(expr_new_function(
                    expr_new_symbol("PolynomialRemainder"), args, 3));
    bool z = (rem && rem->type == EXPR_INTEGER && rem->data.integer == 0);
    if (rem) expr_free(rem);
    return z;
}

/* A place over an affine prime. */
typedef struct { int ram; Expr* xi; Expr* eta; int j; } RnbPlace;

static void rnb_places_free(RnbPlace* p, size_t n) {
    if (!p) return;
    for (size_t i = 0; i < n; i++) { if (p[i].xi) expr_free(p[i].xi);
                                     if (p[i].eta) expr_free(p[i].eta); }
    free(p);
}

/* Numerical test: is e a real number?  Used to keep the m>=3 residue tier off
 * complex places, where the branch series (a root-of-a-complex-number Puiseux
 * expansion) is prohibitively slow — a documented scope limit, not a wrong
 * answer (the case simply declines). */
static bool rnb_is_real_num(const Expr* e) {
    Expr* im = eval_and_free(mk_unary("Im",
                   mk_binary("N", expr_copy((Expr*)e), mk_int(30))));
    Expr* av = eval_and_free(mk_unary("Abs", im));
    Expr* lt = mk_div(mk_int(1), mk_pow(mk_int(10), mk_int(15)));
    Expr* cmp = eval_and_free(mk_binary("Less", av, lt));
    bool r = (cmp->type == EXPR_SYMBOL && cmp->data.symbol.name == SYM_True);
    expr_free(cmp);
    return r;
}

/* Enumerate the affine places where f may have a pole (over the factors of
 * den(f)).  q borrowed; caller frees the array via rnb_places_free. */
static RnbPlace* rnb_collect_places(const RadicalField* F, Expr* const* f_elem,
                                    const Expr* q, size_t* nout) {
    *nout = 0;
    RnbPlace* places = NULL; size_t np = 0, cap = 0;
    Expr* d = rf_denominator_poly(F, f_elem);
    Expr* fl = eval_and_free(mk_unary("FactorList", d));  /* {{c,1},{fac,e},...} */
    if (!fl || !head_is(fl, SYM_List)) { if (fl) expr_free(fl); return NULL; }

    for (size_t fi = 0; fi < fl->data.function.arg_count; fi++) {
        Expr* pair = fl->data.function.args[fi];
        if (!head_is(pair, SYM_List) || pair->data.function.arg_count != 2) continue;
        Expr* fac = pair->data.function.args[0];
        int64_t dgf = eval_degree(expr_copy(fac), F->x);
        if (dgf < 1) continue;  /* constant / content */

        /* ramified? fac | Q_j for some present j (gcd(j,m)=1 by construction) */
        int ram_j = -1;
        for (int j = 1; j < F->m; j++)
            if (F->present[j] && rnb_divides(fac, F->Q[j], F->x)
                && igcd(j, F->m) == 1) { ram_j = j; break; }

        Expr** xis = NULL; size_t nxi = 0;
        rnb_poly_roots(fac, F->xname, &xis, &nxi);
        for (size_t r = 0; r < nxi; r++) {
            if (ram_j >= 0) {
                if (np == cap) { cap = cap ? cap*2 : 8;
                    places = (RnbPlace*)realloc(places, cap*sizeof(RnbPlace)); }
                places[np].ram = 1; places[np].xi = expr_copy(xis[r]);
                places[np].eta = NULL; places[np].j = ram_j; np++;
            } else {
                /* unramified: eta over Y^m - q(xi).  For m>=3 the branch series
                 * at a complex place is a root of a complex number and blows up;
                 * skip such places (the integral then declines, never wrong). */
                if (F->m >= 3 && !rnb_is_real_num(xis[r])) continue;
                Expr* q_at = rnb_subst(q, F->x, xis[r]);
                Expr* Ym   = mk_plus2(mk_pow(expr_new_symbol("RNB$Yr"), mk_int(F->m)),
                                      mk_times2(mk_int(-1), q_at));
                Expr** etas = NULL; size_t neta = 0;
                rnb_poly_roots(Ym, "RNB$Yr", &etas, &neta);
                expr_free(Ym);
                for (size_t e = 0; e < neta; e++) {
                    if (np == cap) { cap = cap ? cap*2 : 8;
                        places = (RnbPlace*)realloc(places, cap*sizeof(RnbPlace)); }
                    places[np].ram = 0; places[np].xi = expr_copy(xis[r]);
                    places[np].eta = expr_copy(etas[e]); places[np].j = -1; np++;
                }
                for (size_t e = 0; e < neta; e++) expr_free(etas[e]);
                free(etas);
            }
        }
        for (size_t r = 0; r < nxi; r++) expr_free(xis[r]);
        free(xis);
    }
    expr_free(fl);
    *nout = np;
    return places;
}

/* Branch series of y with y(xi)=eta at an unramified place, to given order in
 * t = x - xi: Normal[Series[eta*(q(xi+t)/q(xi))^(1/m), {t,0,order}]]. */
static Expr* rnb_branch_series(const RadicalField* F, const Expr* xi,
                               const Expr* eta, int64_t order, const Expr* q) {
    Expr* xit  = mk_plus2(expr_copy((Expr*)xi), expr_new_symbol(RNB_T));
    Expr* q_at = rnb_subst(q, F->x, xit); expr_free(xit);
    Expr* q0   = rnb_subst(q, F->x, xi);
    Expr* ratio = mk_div(q_at, q0);
    Expr* base = mk_times2(expr_copy((Expr*)eta), mk_pow(ratio, mk_rat(1, F->m)));
    Expr* spec_args[3] = { expr_new_symbol(RNB_T), mk_int(0), mk_int(order) };
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 3);
    Expr* ser = eval_and_free(mk_binary("Series", base, spec));
    return eval_and_free(mk_unary("Normal", ser));
}

/* Residue of f dx at an unramified place (xi, eta). */
static Expr* rnb_residue_unram(const RadicalField* F, Expr* const* f_elem,
                               const Expr* xi, const Expr* eta, const Expr* q,
                               int64_t ordr) {
    Expr* ys  = rnb_branch_series(F, xi, eta, ordr + 3, q);
    Expr* xit = mk_plus2(expr_copy((Expr*)xi), expr_new_symbol(RNB_T));
    Expr* tot = mk_int(0);
    for (int i = 0; i < F->m; i++) {
        if (f_elem[i]->type == EXPR_INTEGER && f_elem[i]->data.integer == 0) continue;
        Expr* fi  = rnb_subst(f_elem[i], F->x, xit);
        Expr* Ei  = rnb_subst(F->E[i], F->x, xit);
        Expr* ysi = (i == 0) ? mk_int(1) : mk_pow(expr_copy(ys), mk_int(i));
        Expr* term = mk_div(mk_times2(fi, ysi), Ei);
        tot = mk_plus2(tot, term);
    }
    expr_free(xit); expr_free(ys);
    Expr* zero = mk_int(0);
    Expr* res = rnb_residue2(tot, RNB_T, zero);
    expr_free(zero);
    return eval_and_free(mk_unary("Simplify", res));
}

/* Residue of f dx at the single ramified place over x=xi: residue of the
 * trace of the multiplication matrix. */
static Expr* rnb_residue_ram(const RadicalField* F, Expr* const* f_elem,
                             const Expr* xi) {
    Expr* M  = rf_mult_matrix(F, f_elem);
    Expr* tr = eval_cancel(mk_unary("Tr", M));
    /* Residue[0, ...] does not auto-simplify; short-circuit a zero trace. */
    if (tr && tr->type == EXPR_INTEGER && tr->data.integer == 0) return tr;
    Expr* res = rnb_residue2(tr, F->xname, xi);
    return eval_and_free(mk_unary("Simplify", res));
}

static Expr* rnb_residue_at(const RadicalField* F, Expr* const* f_elem,
                            const Expr* q, const RnbPlace* pl, int64_t ordr) {
    if (pl->ram) return rnb_residue_ram(F, f_elem, pl->xi);
    return rnb_residue_unram(F, f_elem, pl->xi, pl->eta, q, ordr);
}

/* ================================================================== */
/* Exact tier (III): elements with prescribed divisor -> logands.      */
/* ================================================================== */

#define RNB_NMAX    4
#define RNB_DEGMAX  3
#define RNB_MAXLOGS 16

static bool rnb_is_zero(const Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER && e->data.integer == 0) return true;
    Expr* z = eval_and_free(mk_unary("Simplify", expr_copy((Expr*)e)));
    bool zero = (z && z->type == EXPR_INTEGER && z->data.integer == 0);
    if (z) expr_free(z);
    return zero;
}

/* Conditions that sum_i alpha_i w_i vanishes to order >= r at unramified
 * (xi, eta): the t^0..t^{r-1} coefficients of the branch expansion. */
static Expr* rnb_valcond_unram(const RadicalField* F, Expr* const* alpha,
                               const Expr* xi, const Expr* eta, int r,
                               const Expr* q) {
    Expr* ys  = rnb_branch_series(F, xi, eta, r + 2, q);
    Expr* xit = mk_plus2(expr_copy((Expr*)xi), expr_new_symbol(RNB_T));
    Expr* tot = mk_int(0);
    for (int i = 0; i < F->m; i++) {
        Expr* ai  = rnb_subst(alpha[i], F->x, xit);
        Expr* Ei  = rnb_subst(F->E[i], F->x, xit);
        Expr* ysi = (i == 0) ? mk_int(1) : mk_pow(expr_copy(ys), mk_int(i));
        tot = mk_plus2(tot, mk_div(mk_times2(ai, ysi), Ei));
    }
    expr_free(xit); expr_free(ys);
    Expr* spec_args[3] = { expr_new_symbol(RNB_T), mk_int(0), mk_int(r) };
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List), spec_args, 3);
    Expr* ser  = eval_and_free(mk_unary("Normal",
                     eval_and_free(mk_binary("Series", tot, spec))));
    Expr** cells = (Expr**)calloc((size_t)(r > 0 ? r : 1), sizeof(Expr*));
    for (int k = 0; k < r; k++) {
        Expr* a3[3] = { expr_copy(ser), expr_new_symbol(RNB_T), mk_int(k) };
        cells[k] = eval_and_free(expr_new_function(
                       expr_new_symbol(SYM_Coefficient), a3, 3));
    }
    expr_free(ser);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)(r > 0 ? r : 0));
    free(cells);
    return out;
}

/* Conditions for vanishing to order >= r at the single ramified place over
 * x=xi (Q_j, e=m): e*v_p(alpha_i) + (ij mod m) >= r. */
static Expr* rnb_valcond_ram(const RadicalField* F, Expr* const* alpha,
                             const Expr* xi, int j, int r) {
    int m = F->m, e = m;
    Expr** cells = NULL; size_t nc = 0, cap = 0;
    for (int i = 0; i < m; i++) {
        int off = (int)(((int64_t)i * j) % m);
        int need = (r > off) ? ((r - off + e - 1) / e) : 0;
        for (int k = 0; k < need; k++) {
            Expr* dk;
            if (k == 0) dk = expr_copy(alpha[i]);
            else {
                Expr* dspec[2] = { expr_copy(F->x), mk_int(k) };
                Expr* dl = expr_new_function(expr_new_symbol(SYM_List), dspec, 2);
                dk = eval_and_free(mk_binary("D", expr_copy(alpha[i]), dl));
            }
            Expr* val = rnb_subst(dk, F->x, xi);
            expr_free(dk);
            if (nc == cap) { cap = cap ? cap*2 : 8;
                cells = (Expr**)realloc(cells, cap*sizeof(Expr*)); }
            cells[nc++] = val;
        }
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), cells, nc);
    free(cells);
    return out;
}

static uint64_t g_rnb_cid = 0;

/* Elements a in O with affine divisor >= N at each place of `pls`, of
 * coordinate degree <= degB: the nullspace of the vanishing conditions.
 * Returns owned array of elements (each Expr*[m]); count in *ne. */
static Expr*** rnb_find_element(const RadicalField* F, const RnbPlace* pls,
                                size_t npl, int N, int degB, const Expr* q,
                                size_t* ne) {
    *ne = 0;
    if (rnb_over_budget()) return NULL;
    int m = F->m;
    int ncols = m * (degB + 1);
    uint64_t g = g_rnb_cid++;
    Expr** cvar  = (Expr**)calloc((size_t)ncols, sizeof(Expr*));
    Expr** alpha = (Expr**)calloc((size_t)m, sizeof(Expr*));
    for (int i = 0; i < m; i++) {
        Expr* ai = mk_int(0);
        for (int k = 0; k <= degB; k++) {
            char buf[64];
            snprintf(buf, sizeof(buf), "RNB$c$%llu$%d$%d", (unsigned long long)g, i, k);
            Expr* c = expr_new_symbol(intern_symbol(buf));
            cvar[i * (degB + 1) + k] = expr_copy(c);
            ai = mk_plus2(ai, mk_times2(c, mk_pow(expr_copy(F->x), mk_int(k))));
        }
        alpha[i] = eval_expand(ai);
    }

    /* gather conditions */
    Expr** conds = NULL; size_t nc = 0, ccap = 0;
    for (size_t p = 0; p < npl; p++) {
        Expr* cl = pls[p].ram
            ? rnb_valcond_ram(F, alpha, pls[p].xi, pls[p].j, N)
            : rnb_valcond_unram(F, alpha, pls[p].xi, pls[p].eta, N, q);
        if (cl && head_is(cl, SYM_List))
            for (size_t t = 0; t < cl->data.function.arg_count; t++) {
                if (nc == ccap) { ccap = ccap ? ccap*2 : 16;
                    conds = (Expr**)realloc(conds, ccap*sizeof(Expr*)); }
                conds[nc++] = expr_copy(cl->data.function.args[t]);
            }
        if (cl) expr_free(cl);
    }

    /* matrix rows: coefficient of each cvar in each nonzero condition */
    Expr** rows = NULL; size_t nrow = 0;
    for (size_t r = 0; r < nc; r++) {
        Expr** entry = (Expr**)calloc((size_t)ncols, sizeof(Expr*));
        bool allzero = true;
        for (int col = 0; col < ncols; col++) {
            Expr* co = eval_and_free(mk_binary("Coefficient",
                           expr_copy(conds[r]), expr_copy(cvar[col])));
            entry[col] = co;
            if (!(co->type == EXPR_INTEGER && co->data.integer == 0)) allzero = false;
        }
        if (allzero) { for (int col = 0; col < ncols; col++) expr_free(entry[col]);
                       free(entry); continue; }
        rows = (Expr**)realloc(rows, (nrow + 1) * sizeof(Expr*));
        rows[nrow++] = expr_new_function(expr_new_symbol(SYM_List), entry, (size_t)ncols);
        free(entry);
    }

    Expr*** elems = NULL; size_t nel = 0;
    if (nrow > 0) {
        Expr** rc = (Expr**)malloc(nrow * sizeof(Expr*));
        for (size_t r = 0; r < nrow; r++) rc[r] = rows[r];
        Expr* M = expr_new_function(expr_new_symbol(SYM_List), rc, nrow);
        free(rc);
        Expr* ns = eval_and_free(mk_unary("NullSpace", M));
        if (ns && head_is(ns, SYM_List)) {
            for (size_t v = 0; v < ns->data.function.arg_count; v++) {
                Expr* vec = ns->data.function.args[v];
                if (!head_is(vec, SYM_List)
                    || (int)vec->data.function.arg_count != ncols) continue;
                Expr** a = elem_zero(m);
                for (int i = 0; i < m; i++) {
                    Expr* ai = mk_int(0);
                    for (int k = 0; k <= degB; k++) {
                        Expr* co = vec->data.function.args[i * (degB + 1) + k];
                        ai = mk_plus2(ai, mk_times2(expr_copy(co),
                                          mk_pow(expr_copy(F->x), mk_int(k))));
                    }
                    expr_free(a[i]); a[i] = eval_cancel(ai);
                }
                elems = (Expr***)realloc(elems, (nel + 1) * sizeof(Expr**));
                elems[nel++] = a;
            }
        }
        if (ns) expr_free(ns);
    }
    free(rows);

    for (size_t r = 0; r < nc; r++) expr_free(conds[r]);
    free(conds);
    for (size_t i = 0; i < (size_t)ncols; i++) expr_free(cvar[i]);
    free(cvar);
    for (int i = 0; i < m; i++) expr_free(alpha[i]);
    free(alpha);
    *ne = nel;
    return elems;
}

/* Certify that a's affine divisor is exactly N*(sum of the class places):
 * Norm(a) must equal a nonzero constant times prod_P (x - xi_P)^N.  (Used for
 * unramified classes, where Norm has a simple linear factor per place.) */
static bool rnb_affine_divisor_ok(const RadicalField* F, Expr* const* a,
                                  const RnbPlace* pls, size_t npl, int N) {
    Expr* Na = eval_numer(eval_together(rf_norm(F, a)));
    if (Na->type == EXPR_INTEGER && Na->data.integer == 0) { expr_free(Na); return false; }
    Expr* want = mk_int(1);
    for (size_t p = 0; p < npl; p++) {
        Expr* lin = mk_plus2(expr_copy(F->x),
                             mk_times2(mk_int(-1), expr_copy(pls[p].xi)));
        want = mk_times2(want, mk_pow(lin, mk_int(N)));
    }
    Expr* r  = eval_cancel(mk_div(Na, want));
    Expr* rn = eval_numer(eval_together(expr_copy(r)));
    Expr* rd = eval_denom(eval_together(r));
    bool ok = (eval_degree(rn, F->x) == 0) && (eval_degree(rd, F->x) == 0);
    return ok;
}

/* Append a logand element (ownership transferred) to a growing array. */
static void rnb_logs_push(Expr**** logs, size_t* nlog, size_t* cap, Expr** el, int m) {
    if (!el) return;
    if (*nlog >= RNB_MAXLOGS) { elem_free(el, m); return; }  /* pool cap */
    if (*nlog == *cap) { *cap = *cap ? *cap*2 : 8;
        *logs = (Expr***)realloc(*logs, *cap * sizeof(Expr**)); }
    (*logs)[(*nlog)++] = el;
}

/* Residue-driven S-unit logands (exact tier).  Appends elements to logs. */
static void rnb_exact_logands(const RadicalField* F, Expr* const* f_elem,
                              const Expr* q, Expr**** logs, size_t* nlog,
                              size_t* cap) {
    Expr* dd = rf_denominator_poly(F, f_elem);
    int64_t ordr = eval_degree(dd, F->x);
    size_t np = 0;
    RnbPlace* places = rnb_collect_places(F, f_elem, q, &np);
    if (!places) return;

    /* residue per place */
    Expr** resv = (Expr**)calloc(np ? np : 1, sizeof(Expr*));
    for (size_t i = 0; i < np; i++)
        resv[i] = rnb_over_budget() ? mk_int(0)
                                    : rnb_residue_at(F, f_elem, q, &places[i], ordr);

    /* group place indices by (nonzero) residue value */
    bool* used = (bool*)calloc(np ? np : 1, sizeof(bool));
    for (size_t i = 0; i < np; i++) {
        if (rnb_over_budget()) break;
        if (used[i] || rnb_is_zero(resv[i])) { used[i] = true; continue; }
        /* class = all j with residue == resv[i] */
        size_t* grp = (size_t*)calloc(np, sizeof(size_t)); size_t ng = 0;
        grp[ng++] = i; used[i] = true;
        for (size_t j = i + 1; j < np; j++) {
            if (used[j]) continue;
            Expr* diff = mk_plus2(expr_copy(resv[i]),
                                  mk_times2(mk_int(-1), expr_copy(resv[j])));
            if (rnb_is_zero(diff)) { grp[ng++] = j; used[j] = true; }
            expr_free(diff);
        }
        /* build the place subset and search for an element with divisor N*D_c */
        RnbPlace* sub = (RnbPlace*)calloc(ng, sizeof(RnbPlace));
        bool all_unram = true;
        for (size_t t = 0; t < ng; t++) { sub[t] = places[grp[t]];  /* shallow (borrowed) */
                                          if (sub[t].ram) all_unram = false; }
        bool found = false;
        for (int N = 1; N <= RNB_NMAX && !found; N++)
            for (int degB = 0; degB <= RNB_DEGMAX && !found; degB++) {
                size_t nel = 0;
                Expr*** els = rnb_find_element(F, sub, ng, N, degB, q, &nel);
                for (size_t e = 0; e < nel; e++) {
                    bool keep = !elem_is_zero((Expr* const*)els[e], F->m)
                        && (!all_unram
                            || rnb_affine_divisor_ok(F, (Expr* const*)els[e], sub, ng, N));
                    if (keep) { rnb_logs_push(logs, nlog, cap, els[e], F->m); found = true; }
                    else elem_free(els[e], F->m);
                }
                free(els);
            }
        free(sub); free(grp);
    }

    for (size_t i = 0; i < np; i++) if (resv[i]) expr_free(resv[i]);
    free(resv); free(used);
    rnb_places_free(places, np);
}

/* Debug: Integrate`RNB`Residues[f, x] -> {{ram?, xi, eta|j, residue}, ...}. */
static Expr* builtin_rnb_residues(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* q = NULL; int64_t m = 0;
    if (!rnb_find_radical(f, x->data.symbol.name, &q, &m)) return NULL;
    RadicalField F;
    if (!rf_setup(&F, q, (int)m, x)) { expr_free(q); return NULL; }
    const char* Yname = intern_symbol("RNB$Y");
    Expr** f_elem = rnb_f_to_element(&F, f, q, Yname);
    Expr* dd = rf_denominator_poly(&F, f_elem);
    int64_t ordr = eval_degree(dd, F.x);
    size_t np = 0;
    RnbPlace* places = rnb_collect_places(&F, f_elem, q, &np);

    Expr** rows = (Expr**)calloc(np ? np : 1, sizeof(Expr*));
    for (size_t i = 0; i < np; i++) {
        Expr* r = rnb_residue_at(&F, f_elem, q, &places[i], ordr);
        Expr* cell[4] = { mk_int(places[i].ram),
                          expr_copy(places[i].xi),
                          places[i].ram ? mk_int(places[i].j)
                                        : expr_copy(places[i].eta),
                          r };
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cell, 4);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, np);
    free(rows);
    rnb_places_free(places, np);
    elem_free(f_elem, F.m);
    rf_free(&F);
    expr_free(q);
    return out;
}

/* ================================================================== */
/* Diff-back verification.                                             */
/* ================================================================== */

/* True iff D[cand, x] - f is identically zero over Q(x)(radical).  Uses the
 * rigorous algebraic-field zero test (fast; the Groebner normal form behind
 * RootReduce) rather than FullSimplify, which is pathologically slow / hangs on
 * nested radicals.  Falls back to Cancel/Together for the purely rational case. */
MATHILDA_MAYBE_UNUSED static bool rnb_verify(const Expr* cand, const Expr* f, Expr* x) {
    Expr* dcand = eval_diff(expr_copy((Expr*)cand), x);
    Expr* diff  = mk_plus2(dcand, mk_times2(mk_int(-1), expr_copy((Expr*)f)));
    Expr* z = flint_algebraic_field_normalize(diff);
    if (z) {
        bool zero = (z->type == EXPR_INTEGER && z->data.integer == 0);
        expr_free(z); expr_free(diff);
        return zero;
    }
    Expr* zz = eval_and_free(mk_unary("Cancel", eval_together(diff)));
    bool zero = (zz && zz->type == EXPR_INTEGER && zz->data.integer == 0);
    if (zz) expr_free(zz);
    return zero;
}

/* Debug: Integrate`RNB`Logands[f, x] -> the candidate logands in y-form. */
static Expr* builtin_rnb_logands(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    Expr* q = NULL; int64_t m = 0;
    if (!rnb_find_radical(f, x->data.symbol.name, &q, &m)) return NULL;
    g_rnb_deadline = clock() + (clock_t)(RNB_BUDGET_SEC * CLOCKS_PER_SEC);
    RadicalField F;
    if (!rf_setup(&F, q, (int)m, x)) { expr_free(q); return NULL; }
    const char* Yname = intern_symbol("RNB$Y");
    Expr** f_elem = rnb_f_to_element(&F, f, q, Yname);
    Expr*** logs = NULL; size_t nlog = 0, caplog = 0;
    if (F.m == 2) rnb_logs_push(&logs, &nlog, &caplog, rnb_cf_unit(&F, q), F.m);
    rnb_exact_logands(&F, f_elem, q, &logs, &nlog, &caplog);
    Expr** cells = (Expr**)calloc(nlog ? nlog : 1, sizeof(Expr*));
    for (size_t j = 0; j < nlog; j++) cells[j] = rf_to_y_expr(&F, logs[j], q);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), cells, nlog);
    free(cells);
    for (size_t j = 0; j < nlog; j++) elem_free(logs[j], F.m);
    free(logs);
    elem_free(f_elem, F.m);
    rf_free(&F);
    expr_free(q);
    return out;
}

/* Final safety gate: a high-precision numerical diff-back at several rational
 * points.  The algebraic eqs-satisfied check only certifies the linear system
 * that was built; a numerical check independently catches a mis-built system
 * (e.g. a field-setup error).  Returns true iff some point verifies to < 1e-12
 * and none cleanly fails (> 1e-6); singular points are skipped. */
static bool rnb_verify_numeric(const Expr* cand, const Expr* f, Expr* x) {
    static const int pn[5] = { 53, 31, -23, 67, 141 };
    Expr* dcand = eval_diff(expr_copy((Expr*)cand), x);
    int passes = 0;
    for (int p = 0; p < 5; p++) {
        Expr* diff = mk_plus2(expr_copy(dcand),
                              mk_times2(mk_int(-1), expr_copy((Expr*)f)));
        Expr* rule = mk_binary("Rule", expr_copy(x), mk_div(mk_int(pn[p]), mk_int(100)));
        Expr* sub  = eval_and_free(mk_binary("ReplaceAll", diff, rule));
        Expr* val  = eval_and_free(mk_binary("N", sub, mk_int(30)));
        Expr* av   = eval_and_free(mk_unary("Abs", val));
        Expr* lo   = mk_div(mk_int(1), mk_pow(mk_int(10), mk_int(12)));
        Expr* hi   = mk_div(mk_int(1), mk_pow(mk_int(10), mk_int(6)));
        Expr* pass = eval_and_free(mk_binary("Less", expr_copy(av), lo));
        Expr* fail = eval_and_free(mk_binary("Greater", av, hi));
        bool is_pass = (pass->type == EXPR_SYMBOL && pass->data.symbol.name == SYM_True);
        bool is_fail = (fail->type == EXPR_SYMBOL && fail->data.symbol.name == SYM_True);
        expr_free(pass); expr_free(fail);
        if (is_pass) passes++;
        else if (is_fail) { expr_free(dcand); return false; }
    }
    expr_free(dcand);
    return passes >= 1;
}

/* ------------------------------------------------------------------ */
/* Engine driver.                                                      */
/* ------------------------------------------------------------------ */

static Expr* rnb_integrate(Expr* f, Expr* x) {
    g_rnb_deadline = clock() + (clock_t)(RNB_BUDGET_SEC * CLOCKS_PER_SEC);
    Expr* q = NULL; int64_t m = 0;
    if (!rnb_find_radical(f, x->data.symbol.name, &q, &m)) return NULL;

    RadicalField F;
    if (!rf_setup(&F, q, (int)m, x)) { expr_free(q); return NULL; }

    const char* Yname = intern_symbol("RNB$Y");
    Expr** f_elem = rnb_f_to_element(&F, f, q, Yname);

    /* Collect candidate logands (S-units of O). */
    Expr*** logs = NULL; size_t nlog = 0, caplog = 0;

    /* Units at infinity (m=2 continued fraction / Pell). */
    if (F.m == 2) rnb_logs_push(&logs, &nlog, &caplog, rnb_cf_unit(&F, q), F.m);

    /* Residue-driven S-units (branch-place / Jacobian divisor logands). */
    rnb_exact_logands(&F, f_elem, q, &logs, &nlog, &caplog);

    /* parallel_integrate returns element-verified candidates; a numerical
     * diff-back is the final guard against a mis-built system. */
    Expr* result = rnb_parallel_integrate(&F, f_elem, q, Yname, logs, nlog);
    if (result && !rnb_verify_numeric(result, f, x)) {
        expr_free(result); result = NULL;
    }

    for (size_t j = 0; j < nlog; j++) elem_free(logs[j], F.m);
    free(logs);
    elem_free(f_elem, F.m);
    rf_free(&F);
    expr_free(q);
    return result;
}

/* ------------------------------------------------------------------ */
/* Builtin entry + registration.                                       */
/* ------------------------------------------------------------------ */

Expr* builtin_rischnormanblake(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    /* This engine handles precisely the integrands the transcendental
     * pmint declines: those carrying a single radical of x.  If f has no
     * such radical, leave it to the other cascade stages. */
    Expr* q = NULL; int64_t m = 0;
    if (!rnb_find_radical(f, x->data.symbol.name, &q, &m)) return NULL;
    if (q) expr_free(q);

    return rnb_integrate(f, x);
}

static void install(const char* name, Expr* (*fn)(Expr*), const char* docstring) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (docstring) symtab_set_docstring(name, docstring);
}

void int_rnb_init(void) {
    install("Integrate`RischNormanBlake",
            builtin_rischnormanblake,
            "Integrate`RischNormanBlake[f, x] applies the parallel Risch-Norman\n"
            "method generalised to a simple radical extension L = K(y), y^m = q(x)\n"
            "(K = Q(x)), following S. Blake, \"Parallel Integration over Simple\n"
            "Radical Extensions\".  The antiderivative numerator ranges over the\n"
            "integral closure O with the Trager basis w_i = y^i / E_i, and the\n"
            "logands are the S-units of O (denominator factors, branch-place\n"
            "divisors, and units at infinity).  Returns an elementary\n"
            "antiderivative of f when one exists, else the call unevaluated.\n"
            "Handles exactly the single-radical-of-x integrands that\n"
            "Integrate`RischNorman declines.");

    install("Integrate`RNB`Info",
            builtin_rnb_info,
            "Integrate`RNB`Info[f, x] (debug) returns the RischNormanBlake field\n"
            "data for the radical y^m=q detected in f: {m, q, {Q_j}, {E_i},\n"
            "{Lam_i}, Norm[w_1]}.  For inspection and unit testing.");

    install("Integrate`RNB`Residues",
            builtin_rnb_residues,
            "Integrate`RNB`Residues[f, x] (debug) returns the residues of f dx at\n"
            "the affine places over den(f): {{ram?, xi, eta|j, residue}, ...}.");

    install("Integrate`RNB`Logands",
            builtin_rnb_logands,
            "Integrate`RNB`Logands[f, x] (debug) returns the candidate logands\n"
            "(S-units of O: units at infinity + residue-divisor elements) in\n"
            "y-form.  For inspection and unit testing.");
}
