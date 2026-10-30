/*
 * flint_polymat.c
 * ---------------
 * Polynomial matrices over K[x], K in {Q, GF(p), Q(theta)}, on FLINT's generic
 * rings: one implementation for the three coefficient fields (see
 * flint_polymat.h for the FLINT` builtins and their options).
 *
 * WHY.  The Trager integration package (src/internal/trager/) does its divisor
 * class arithmetic on ideals stored as Hermite forms of polynomial matrices
 * over K[x]: the reduced divisor classes of Jacobian.m, the torsion searches
 * over GF(p) of TorsionModP.m and the characteristic-0 arithmetic over the
 * residue field Q(theta).  Written in the package language those kernels pay an
 * expression-level Expand / PolynomialRemainder per coefficient operation; a
 * profile of the slow stress cases put 67% of every GF(p) class addition in the
 * Hermite form and ~40 ms per Q(theta) division.  Here the same algorithms run
 * on gr_poly with nmod / fmpq / nf coefficients.
 *
 * CONTRACT.  The results are exactly the canonical forms the package computes
 * (its routines remain the reference and the fallback):
 *   - Hermite form: upper triangular, monic pivots, every entry above a pivot
 *     of lower degree than the pivot, zero rows dropped; the module must have
 *     full column rank.  It is unique, so the native and the reference forms
 *     agree identically.  With "DeterminantMultiple" -> D (a multiple of the
 *     determinant of the module) the rows D e_k are adjoined and entries are
 *     reduced modulo D during the elimination, which bounds their degree.
 *   - Output entries are expanded polynomials in x (and t) with canonical
 *     coefficients: rationals over Q and Q(theta) (t-degree < deg mu),
 *     integers in [0, p) over GF(p).
 * Anything out of scope (a non-polynomial entry, another symbol, a rank-
 * deficient module, a non-prime or huge modulus, a reducible mu, a zero
 * divisor) leaves the call unevaluated, so the caller's fallback runs.
 */
#include "flint_polymat.h"

#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#ifdef USE_FLINT

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_mpoly.h>
#include <flint/ulong_extras.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/gr_mat.h>

#include "flint_bridge_internal.h"

/* gr_poly_xgcd arrived in FLINT 3.2 (see flint_bridge.c). */
#if defined(__FLINT_RELEASE) && __FLINT_RELEASE >= 30200
#define PM_GR_POLY_XGCD gr_poly_xgcd
#else
#define PM_GR_POLY_XGCD gr_poly_xgcd_euclidean
#endif

/* ------------------------------------------------------------------ */
/*  The coefficient field                                              */
/* ------------------------------------------------------------------ */

enum { PM_Q, PM_GF, PM_NF };

typedef struct {
    int kind;
    gr_ctx_t K;             /* the field                                     */
    int K_ready;
    ulong p;                /* GF(p)                                         */
    nf_t nf;                /* Q(theta): for element conversions only        */
    fmpq_poly_t mu;         /* Q(theta): the minimal polynomial              */
    slong n;                /* [K : prime field]                             */
    const char* tname;      /* Q(theta): the generator symbol                */
    const char* xname;      /* the polynomial variable                       */
    VarSet vs;              /* {x} or {x, t}: the fmpq_mpoly generators      */
    const fmpq_mpoly_ctx_struct* mctx;
    slong sz;               /* sizeof an element of K                        */
} PmField;

static const char* pm_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL) return NULL;
    return e->data.function.head->data.symbol.name;
}

static int pm_is_head(const Expr* e, const char* name) {
    const char* h = pm_head(e);
    return h && strcmp(h, name) == 0;
}

static int pm_varset_push(VarSet* vs, const char* name) {
    char** np = realloc(vs->names, (vs->count + 1) * sizeof(char*));
    if (!np) return 0;
    vs->names = np;
    vs->names[vs->count] = mathilda_strdup(name);
    if (!vs->names[vs->count]) return 0;
    vs->count++;
    vs->cap = vs->count;
    return 1;
}

/* mu irreducible over Q?  The last answer is cached: a caller works over one
 * residue field for many calls. */
static fmpq_poly_t g_pm_irred_mu;
static int g_pm_irred_ready = 0, g_pm_irred_ans = 0;

static int pm_irreducible(const fmpq_poly_t mu) {
    if (g_pm_irred_ready && fmpq_poly_equal(mu, g_pm_irred_mu)) return g_pm_irred_ans;
    fmpz_poly_t z; fmpz_poly_init(z);
    fmpq_poly_get_numerator(z, mu);
    fmpz_poly_factor_t f; fmpz_poly_factor_init(f);
    fmpz_poly_factor(f, z);
    int ans = (f->num == 1 && f->exp[0] == 1 &&
               fmpz_poly_degree(f->p + 0) == fmpz_poly_degree(z));
    fmpz_poly_factor_clear(f);
    fmpz_poly_clear(z);
    if (!g_pm_irred_ready) { fmpq_poly_init(g_pm_irred_mu); g_pm_irred_ready = 1; }
    fmpq_poly_set(g_pm_irred_mu, mu);
    g_pm_irred_ans = ans;
    return ans;
}

static void pm_field_clear(PmField* F) {
    if (F->K_ready) gr_ctx_clear(F->K);
    if (F->kind == PM_NF && F->K_ready) nf_clear(F->nf);
    if (F->kind == PM_NF) fmpq_poly_clear(F->mu);
    fb_varset_free(&F->vs);
}

/* The Hermite form's own options (string keys); NULL when absent. */
typedef struct {
    const Expr* D;          /* "DeterminantMultiple" -> D                    */
    const Expr* den;        /* "Denominator" -> den                          */
    const Expr* products;   /* "Products" -> {gens, table}                   */
} PmOpts;

/* Read the options (Rule[...] arguments) and set up the field.  String-keyed
 * options are accepted only when `xo` is given (the Hermite form).  Returns 1
 * on success; on failure everything is released. */
static int pm_field_init(PmField* F, const Expr* x, Expr* const* opts, size_t nopts,
                         PmOpts* xo) {
    memset(F, 0, sizeof *F);
    F->kind = PM_Q; F->n = 1;
    if (!x || x->type != EXPR_SYMBOL) return 0;
    F->xname = x->data.symbol.name;
    const Expr* ext = NULL;
    if (xo) memset(xo, 0, sizeof *xo);
    for (size_t i = 0; i < nopts; i++) {
        const Expr* r = opts[i];
        if (!pm_is_head(r, "Rule") || r->data.function.arg_count != 2) return 0;
        const Expr* k = r->data.function.args[0];
        const Expr* v = r->data.function.args[1];
        if (k->type == EXPR_SYMBOL && strcmp(k->data.symbol.name, "Modulus") == 0) {
            if (v->type != EXPR_INTEGER) return 0;
            if (v->data.integer == 0) continue;          /* Modulus -> 0: over Q */
            if (v->data.integer < 2 || !n_is_prime((ulong)v->data.integer)) return 0;
            F->kind = PM_GF; F->p = (ulong)v->data.integer;
        } else if (k->type == EXPR_SYMBOL && strcmp(k->data.symbol.name, "Extension") == 0) {
            if (!pm_is_head(v, "List") || v->data.function.arg_count != 2) return 0;
            ext = v;
        } else if (k->type == EXPR_STRING && xo &&
                   strcmp(k->data.string, "DeterminantMultiple") == 0) {
            xo->D = v;
        } else if (k->type == EXPR_STRING && xo &&
                   strcmp(k->data.string, "Denominator") == 0) {
            xo->den = v;
        } else if (k->type == EXPR_STRING && xo &&
                   strcmp(k->data.string, "Products") == 0) {
            xo->products = v;
        } else {
            return 0;
        }
    }
    if (ext && F->kind == PM_GF) return 0;
    if (!pm_varset_push(&F->vs, F->xname)) { fb_varset_free(&F->vs); return 0; }
    if (ext) {
        const Expr* t = ext->data.function.args[0];
        if (t->type != EXPR_SYMBOL || strcmp(t->data.symbol.name, F->xname) == 0 ||
            !pm_varset_push(&F->vs, t->data.symbol.name)) {
            fb_varset_free(&F->vs); return 0;
        }
        F->tname = t->data.symbol.name;
        F->kind = PM_NF;
        fmpq_poly_init(F->mu);
        /* mu as a univariate polynomial in t */
        const fmpq_mpoly_ctx_struct* c1 = fb_ctx_lex(1);
        VarSet tv = {0};
        int ok = pm_varset_push(&tv, F->tname);
        fmpq_mpoly_t m; fmpq_mpoly_init(m, c1);
        ok = ok && fb_to_mpoly(ext->data.function.args[1], m, c1, &tv);
        if (ok) ok = fmpq_mpoly_get_fmpq_poly(F->mu, m, 0, c1);
        fmpq_mpoly_clear(m, c1);
        fb_varset_free(&tv);
        F->n = fmpq_poly_degree(F->mu);
        if (!ok || F->n < 1 || !pm_irreducible(F->mu)) {
            fmpq_poly_clear(F->mu); fb_varset_free(&F->vs); return 0;
        }
        /* The gr context owns its own copy of Q[t]/(mu) (a context made "from
         * a reference" still frees it on clear); F->nf, built from the same mu
         * and so of the same element layout, serves the conversions
         * nf_elem <-> fmpq_poly only. */
        nf_init(F->nf, F->mu);
        gr_ctx_init_nf(F->K, F->mu);
    } else if (F->kind == PM_GF) {
        /* p is a checked prime >= 2.  gr_ctx_init_nmod returns void before FLINT
         * 3.3 and int (failing only for p == 0) after, so its result is unused. */
        gr_ctx_init_nmod(F->K, F->p);
    } else {
        gr_ctx_init_fmpq(F->K);
    }
    F->K_ready = 1;
    F->mctx = fb_ctx_lex((slong)F->vs.count);
    F->sz = gr_ctx_sizeof_elem(F->K);
    return 1;
}

#define PM_ENTRY(ptr, i, F) ((gr_ptr)((char*)(ptr) + (i) * (F)->sz))

/* ------------------------------------------------------------------ */
/*  Expr <-> K[x]                                                      */
/* ------------------------------------------------------------------ */

/* Set the K-element `c` from the t-polynomial `tp` (Q(theta)) or the rational
 * constant term of `tp` (Q, GF(p)). */
static int pm_set_elem(gr_ptr c, const fmpq_poly_t tp, PmField* F) {
    if (F->kind == PM_NF) {
        fmpq_poly_t r; fmpq_poly_init(r);
        fmpq_poly_rem(r, tp, F->mu);
        nf_elem_set_fmpq_poly((nf_elem_struct*)c, r, F->nf);
        fmpq_poly_clear(r);
        return 1;
    }
    fmpq_t q; fmpq_init(q);
    fmpq_poly_get_coeff_fmpq(q, tp, 0);
    int ok = gr_set_fmpq(c, q, F->K) == GR_SUCCESS;
    fmpq_clear(q);
    return ok;
}

/* Expr polynomial in x (and t) -> gr_poly over K.  0 when out of scope. */
static int pm_from_expr(gr_poly_t P, const Expr* e, PmField* F) {
    fmpq_mpoly_t A; fmpq_mpoly_init(A, F->mctx);
    int ok = fb_to_mpoly(e, A, F->mctx, &F->vs);
    if (ok) ok = gr_poly_zero(P, F->K) == GR_SUCCESS;
    if (ok && !fmpq_mpoly_is_zero(A, F->mctx)) {
        slong nv = (slong)F->vs.count;
        slong dx = fmpq_mpoly_degree_si(A, 0, F->mctx);
        if (F->kind != PM_NF && nv > 1) ok = 0;
        fmpq_poly_struct* tp = NULL;
        if (ok) {
            tp = malloc((size_t)(dx + 1) * sizeof(fmpq_poly_struct));
            if (!tp) ok = 0;
        }
        if (ok) {
            for (slong d = 0; d <= dx; d++) fmpq_poly_init(tp + d);
            slong exps[2] = {0, 0};
            fmpq_t c; fmpq_init(c);
            for (slong i = 0; i < fmpq_mpoly_length(A, F->mctx); i++) {
                fmpq_mpoly_get_term_exp_si(exps, A, i, F->mctx);
                fmpq_mpoly_get_term_coeff_fmpq(c, A, i, F->mctx);
                fmpq_t old; fmpq_init(old);
                slong tj = nv > 1 ? exps[1] : 0;
                fmpq_poly_get_coeff_fmpq(old, tp + exps[0], tj);
                fmpq_add(old, old, c);
                fmpq_poly_set_coeff_fmpq(tp + exps[0], tj, old);
                fmpq_clear(old);
            }
            fmpq_clear(c);
            gr_ptr cf = gr_heap_init(F->K);
            for (slong d = 0; d <= dx && ok; d++) {
                if (fmpq_poly_is_zero(tp + d)) continue;
                ok = pm_set_elem(cf, tp + d, F) &&
                     gr_poly_set_coeff_scalar(P, d, cf, F->K) == GR_SUCCESS;
            }
            gr_heap_clear(cf, F->K);
            for (slong d = 0; d <= dx; d++) fmpq_poly_clear(tp + d);
            free(tp);
        }
    }
    fmpq_mpoly_clear(A, F->mctx);
    return ok;
}

static Expr* pm_mono(const char* v, slong k) {
    if (k == 1) return expr_new_symbol(v);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_symbol(v), expr_new_integer((int64_t)k) }, 2);
}

/* Terms of the coefficient c of x^d, appended to terms[*nt]. */
static int pm_push_terms(Expr*** terms, size_t* nt, size_t* cap, gr_srcptr c, slong d,
                         PmField* F) {
    fmpq_poly_t tp; fmpq_poly_init(tp);
    if (F->kind == PM_NF) {
        nf_elem_get_fmpq_poly(tp, (const nf_elem_struct*)c, F->nf);
    } else if (F->kind == PM_GF) {
        fmpz_t z; fmpz_init_set_ui(z, *(const ulong*)c);
        fmpq_poly_set_fmpz(tp, z);
        fmpz_clear(z);
    } else {
        fmpq_poly_set_fmpq(tp, (const fmpq *)c);
    }
    fmpq_t q; fmpq_init(q);
    for (slong j = 0; j < fmpq_poly_length(tp); j++) {
        fmpq_poly_get_coeff_fmpq(q, tp, j);
        if (fmpq_is_zero(q)) continue;
        Expr* parts[3]; size_t np = 0;
        parts[np++] = fb_fmpq_to_expr(q);
        if (j > 0) parts[np++] = pm_mono(F->tname, j);
        if (d > 0) parts[np++] = pm_mono(F->xname, d);
        Expr* term = np == 1 ? parts[0]
                   : expr_new_function(expr_new_symbol(SYM_Times), parts, np);
        if (*nt == *cap) {
            size_t nc = *cap ? 2 * *cap : 16;
            Expr** g = realloc(*terms, nc * sizeof(Expr*));
            if (!g) { expr_free(term); fmpq_clear(q); fmpq_poly_clear(tp); return 0; }
            *terms = g; *cap = nc;
        }
        (*terms)[(*nt)++] = term;
    }
    fmpq_clear(q);
    fmpq_poly_clear(tp);
    return 1;
}

static Expr* pm_eval_free(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* gr_poly over K -> canonical expanded Expr. */
static Expr* pm_to_expr(const gr_poly_t P, PmField* F) {
    Expr** terms = NULL; size_t nt = 0, cap = 0;
    for (slong d = 0; d < P->length; d++) {
        gr_srcptr c = PM_ENTRY(P->coeffs, d, F);
        if (gr_is_zero(c, F->K) == T_TRUE) continue;
        if (!pm_push_terms(&terms, &nt, &cap, c, d, F)) {
            for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
            free(terms);
            return NULL;
        }
    }
    Expr* r;
    if (nt == 0) r = expr_new_integer(0);
    else if (nt == 1) r = terms[0];
    else r = expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return pm_eval_free(r);
}

static Expr* pm_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

/* ------------------------------------------------------------------ */
/*  Matrices of gr_poly                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    slong r, c;
    gr_poly_struct* e;      /* row-major, r*c entries */
} PmMat;

#define PM_AT(M, i, j) ((M)->e + (i) * (M)->c + (j))

static void pm_mat_init(PmMat* M, slong r, slong c, PmField* F) {
    M->r = r; M->c = c;
    M->e = malloc((size_t)(r * c > 0 ? r * c : 1) * sizeof(gr_poly_struct));
    for (slong i = 0; i < r * c; i++) gr_poly_init(M->e + i, F->K);
}

static void pm_mat_clear(PmMat* M, PmField* F) {
    for (slong i = 0; i < M->r * M->c; i++) gr_poly_clear(M->e + i, F->K);
    free(M->e);
}

/* A List of equal-length Lists of polynomials -> PmMat with `extra` zero rows
 * appended.  0 when out of scope. */
static int pm_mat_from_expr(PmMat* M, const Expr* m, slong extra, PmField* F) {
    if (!pm_is_head(m, "List") || m->data.function.arg_count == 0) return 0;
    slong r = (slong)m->data.function.arg_count;
    const Expr* r0 = m->data.function.args[0];
    if (!pm_is_head(r0, "List") || r0->data.function.arg_count == 0) return 0;
    slong c = (slong)r0->data.function.arg_count;
    pm_mat_init(M, r + extra, c, F);
    for (slong i = 0; i < r; i++) {
        const Expr* row = m->data.function.args[i];
        if (!pm_is_head(row, "List") || (slong)row->data.function.arg_count != c) {
            pm_mat_clear(M, F); return 0;
        }
        for (slong j = 0; j < c; j++)
            if (!pm_from_expr(PM_AT(M, i, j), row->data.function.args[j], F)) {
                pm_mat_clear(M, F); return 0;
            }
    }
    return 1;
}

static Expr* pm_mat_to_expr(const PmMat* M, slong rows, PmField* F) {
    Expr** R = malloc((size_t)(rows > 0 ? rows : 1) * sizeof(Expr*));
    for (slong i = 0; i < rows; i++) {
        Expr** row = malloc((size_t)M->c * sizeof(Expr*));
        for (slong j = 0; j < M->c; j++) {
            row[j] = pm_to_expr(PM_AT(M, i, j), F);
            if (!row[j]) {
                for (slong k = 0; k < j; k++) expr_free(row[k]);
                free(row);
                for (slong k = 0; k < i; k++) expr_free(R[k]);
                free(R);
                return NULL;
            }
        }
        R[i] = pm_list(row, (size_t)M->c);
        free(row);
    }
    Expr* out = pm_list(R, (size_t)rows);
    free(R);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Hermite normal form                                                */
/* ------------------------------------------------------------------ */

/* The algorithm of the package's hnf / hnfModDet (src/internal/trager/
 * Normalize.m), column by column: the nonzero entry of least degree becomes the
 * pivot; every row below is combined with the pivot row by the extended gcd
 * (a unimodular 2x2 step); the pivot is made monic.  With a modulus D (monic,
 * deg >= 1) the rows D e_k have been appended and every entry right of the
 * current column is reduced mod D after each step.  Finally the entries above
 * each pivot are reduced.  Returns 1 with the form in the first c rows, 0 on
 * rank deficiency or a FLINT failure. */
static int pm_hnf(PmMat* M, const gr_poly_t D, PmField* F) {
    gr_ctx_struct* K = F->K;
    slong nr = M->r, nc = M->c;
    int st = GR_SUCCESS;
    int modD = D != NULL;
    gr_poly_t g, s, t, a, b, q, r1, r2, tmp, tmp2;
    gr_poly_init(g, K); gr_poly_init(s, K); gr_poly_init(t, K); gr_poly_init(a, K);
    gr_poly_init(b, K); gr_poly_init(q, K); gr_poly_init(r1, K); gr_poly_init(r2, K);
    gr_poly_init(tmp, K); gr_poly_init(tmp2, K);
    gr_ptr lc = gr_heap_init(K);
    int ok = 1;
    for (slong j = 0; j < nc && ok; j++) {
        /* pivot: least degree in column j among rows j.. */
        slong piv = -1, best = 0;
        for (slong i = j; i < nr; i++) {
            slong len = PM_AT(M, i, j)->length;
            if (len > 0 && (piv < 0 || len < best)) { piv = i; best = len; }
        }
        if (piv < 0) { ok = 0; break; }
        if (piv != j)
            for (slong k = 0; k < nc; k++) gr_poly_swap(PM_AT(M, j, k), PM_AT(M, piv, k), K);
        for (slong i = j + 1; i < nr && ok; i++) {
            if (PM_AT(M, i, j)->length == 0) continue;
            st |= gr_poly_set(a, PM_AT(M, j, j), K);
            st |= gr_poly_set(b, PM_AT(M, i, j), K);
            st |= PM_GR_POLY_XGCD(g, s, t, a, b, K);      /* s a + t b = g */
            st |= gr_poly_divexact(a, a, g, K);          /* a/g */
            st |= gr_poly_divexact(b, b, g, K);          /* b/g */
            for (slong k = j; k < nc; k++) {
                gr_poly_struct* rj = PM_AT(M, j, k);
                gr_poly_struct* ri = PM_AT(M, i, k);
                st |= gr_poly_set(r1, rj, K);
                st |= gr_poly_set(r2, ri, K);
                st |= gr_poly_mul(tmp, s, r1, K);
                st |= gr_poly_mul(tmp2, t, r2, K);
                st |= gr_poly_add(rj, tmp, tmp2, K);
                st |= gr_poly_mul(tmp, a, r2, K);
                st |= gr_poly_mul(tmp2, b, r1, K);
                st |= gr_poly_sub(ri, tmp, tmp2, K);
                if (modD && k > j) {
                    st |= gr_poly_rem(rj, rj, D, K);
                    st |= gr_poly_rem(ri, ri, D, K);
                }
            }
            if (st != GR_SUCCESS) ok = 0;
        }
        if (!ok) break;
        /* monic pivot */
        gr_poly_struct* pv = PM_AT(M, j, j);
        if (pv->length == 0) { ok = 0; break; }
        st |= gr_inv(lc, PM_ENTRY(pv->coeffs, pv->length - 1, F), K);
        for (slong k = j; k < nc; k++)
            st |= gr_poly_mul_scalar(PM_AT(M, j, k), PM_AT(M, j, k), lc, K);
        if (st != GR_SUCCESS) ok = 0;
    }
    /* rows below the form must have vanished (always so modulo D) */
    for (slong i = nc; i < nr && ok; i++)
        for (slong k = 0; k < nc; k++)
            if (PM_AT(M, i, k)->length != 0) { ok = 0; break; }
    /* reduce above the pivots, column by column */
    for (slong j = 0; j < nc && ok; j++) {
        gr_poly_struct* pv = PM_AT(M, j, j);
        for (slong i = 0; i < j; i++) {
            gr_poly_struct* e = PM_AT(M, i, j);
            if (e->length < pv->length) continue;
            st |= gr_poly_divrem(q, tmp, e, pv, K);
            for (slong k = j; k < nc; k++) {
                st |= gr_poly_mul(tmp, q, PM_AT(M, j, k), K);
                st |= gr_poly_sub(PM_AT(M, i, k), PM_AT(M, i, k), tmp, K);
            }
        }
        if (st != GR_SUCCESS) ok = 0;
    }
    gr_heap_clear(lc, K);
    gr_poly_clear(g, K); gr_poly_clear(s, K); gr_poly_clear(t, K); gr_poly_clear(a, K);
    gr_poly_clear(b, K); gr_poly_clear(q, K); gr_poly_clear(r1, K); gr_poly_clear(r2, K);
    gr_poly_clear(tmp, K); gr_poly_clear(tmp2, K);
    return ok;
}

/* "Products" -> {gens, table}: the module spanned by the products r g (r a
 * row of A, g in gens) in the K[x]-algebra with basis w_0..w_{m-1} and
 * w_i w_j = f_ij w_{k_ij}, table[[i+1, j+1]] = {k_ij, f_ij}.  Replaces A by
 * those products.  0 when out of scope. */
static int pm_apply_products(PmMat* A, const Expr* spec, PmField* F) {
    if (!pm_is_head(spec, "List") || spec->data.function.arg_count != 2) return 0;
    const Expr* gl = spec->data.function.args[0];
    const Expr* tb = spec->data.function.args[1];
    slong m = A->c;
    if (!pm_is_head(tb, "List") || (slong)tb->data.function.arg_count != m) return 0;
    PmMat G;
    if (!pm_mat_from_expr(&G, gl, 0, F)) return 0;
    if (G.c != m) { pm_mat_clear(&G, F); return 0; }
    slong* tk = malloc((size_t)(m * m) * sizeof(slong));
    gr_poly_struct* tf = malloc((size_t)(m * m) * sizeof(gr_poly_struct));
    for (slong i = 0; i < m * m; i++) gr_poly_init(tf + i, F->K);
    int ok = 1;
    for (slong i = 0; i < m && ok; i++) {
        const Expr* row = tb->data.function.args[i];
        if (!pm_is_head(row, "List") || (slong)row->data.function.arg_count != m) { ok = 0; break; }
        for (slong j = 0; j < m && ok; j++) {
            const Expr* kf = row->data.function.args[j];
            if (!pm_is_head(kf, "List") || kf->data.function.arg_count != 2 ||
                kf->data.function.args[0]->type != EXPR_INTEGER) { ok = 0; break; }
            tk[i * m + j] = (slong)kf->data.function.args[0]->data.integer;
            if (tk[i * m + j] < 0 || tk[i * m + j] >= m) { ok = 0; break; }
            ok = pm_from_expr(tf + i * m + j, kf->data.function.args[1], F);
        }
    }
    if (ok) {
        PmMat P; pm_mat_init(&P, A->r * G.r, m, F);
        gr_poly_t u; gr_poly_init(u, F->K);
        int st = GR_SUCCESS;
        for (slong r = 0; r < A->r; r++)
            for (slong g = 0; g < G.r; g++) {
                slong row = r * G.r + g;
                for (slong i = 0; i < m; i++) {
                    if (PM_AT(A, r, i)->length == 0) continue;
                    for (slong j = 0; j < m; j++) {
                        if (PM_AT(&G, g, j)->length == 0) continue;
                        st |= gr_poly_mul(u, PM_AT(A, r, i), PM_AT(&G, g, j), F->K);
                        st |= gr_poly_mul(u, u, tf + i * m + j, F->K);
                        st |= gr_poly_add(PM_AT(&P, row, tk[i * m + j]),
                                          PM_AT(&P, row, tk[i * m + j]), u, F->K);
                    }
                }
            }
        gr_poly_clear(u, F->K);
        pm_mat_clear(A, F);
        *A = P;
        ok = st == GR_SUCCESS;
    }
    for (slong i = 0; i < m * m; i++) gr_poly_clear(tf + i, F->K);
    free(tf); free(tk);
    pm_mat_clear(&G, F);
    return ok;
}

/* The polynomial over the PRIME field (Q or GF(p)) of the K-polynomial P:
 * fails over Q(theta) when a coefficient is not rational. */
static int pm_prime_poly(fmpq_poly_t out, const gr_poly_t P, PmField* F) {
    fmpq_poly_zero(out);
    fmpq_t q; fmpq_init(q);
    int ok = 1;
    for (slong d = 0; d < P->length && ok; d++) {
        gr_srcptr c = PM_ENTRY(P->coeffs, d, F);
        if (F->kind == PM_GF) {
            fmpz_t z; fmpz_init_set_ui(z, *(const ulong*)c);
            fmpq_set_fmpz(q, z); fmpz_clear(z);
        } else if (F->kind == PM_NF) {
            /* rational iff the t-polynomial is constant */
            fmpq_poly_t tp; fmpq_poly_init(tp);
            nf_elem_get_fmpq_poly(tp, (const nf_elem_struct*)c, F->nf);
            ok = fmpq_poly_degree(tp) <= 0;
            if (ok) fmpq_poly_get_coeff_fmpq(q, tp, 0);
            fmpq_poly_clear(tp);
        } else {
            fmpq_set(q, (const fmpq *)c);
        }
        if (ok) fmpq_poly_set_coeff_fmpq(out, d, q);
    }
    fmpq_clear(q);
    return ok;
}

/* The largest factor g of den (a prime-field polynomial, as a K-polynomial)
 * dividing every entry of the first `rows` rows of M, monic, as a K-polynomial.
 * Over Q and GF(p) a gcd over K; over Q(theta) a gcd over Q of den and every
 * t-coefficient of every entry (a prime-field factor q divides an entry
 * Sum e_j(x) t^j exactly when it divides every e_j). */
static int pm_common_factor(gr_poly_t g, const gr_poly_t den, const PmMat* M, slong rows,
                            PmField* F) {
    int st = GR_SUCCESS;
    if (F->kind != PM_NF) {
        st |= gr_poly_set(g, den, F->K);
        for (slong i = 0; i < rows * M->c && g->length > 1; i++)
            if (M->e[i].length > 0) st |= gr_poly_gcd(g, g, M->e + i, F->K);
        if (g->length > 0) st |= gr_poly_make_monic(g, g, F->K);
        return st == GR_SUCCESS;
    }
    fmpq_poly_t G, e, tp; fmpq_poly_init(G); fmpq_poly_init(e); fmpq_poly_init(tp);
    int ok = pm_prime_poly(G, den, F);
    fmpq_t q; fmpq_init(q);
    for (slong i = 0; ok && i < rows * M->c && fmpq_poly_degree(G) > 0; i++) {
        const gr_poly_struct* E = M->e + i;
        if (E->length == 0) continue;
        /* t-coefficients of the entry: e_j(x) = Sum_d coeff_t^j(c_d) x^d */
        for (slong j = 0; j < F->n && fmpq_poly_degree(G) > 0; j++) {
            fmpq_poly_zero(e);
            for (slong d = 0; d < E->length; d++) {
                nf_elem_get_fmpq_poly(tp, (const nf_elem_struct*)PM_ENTRY(E->coeffs, d, F), F->nf);
                fmpq_poly_get_coeff_fmpq(q, tp, j);
                if (!fmpq_is_zero(q)) fmpq_poly_set_coeff_fmpq(e, d, q);
            }
            if (!fmpq_poly_is_zero(e)) fmpq_poly_gcd(G, G, e);
        }
    }
    if (ok) {
        if (!fmpq_poly_is_zero(G)) fmpq_poly_make_monic(G, G);
        st |= gr_poly_zero(g, F->K);
        gr_ptr c = gr_heap_init(F->K);
        for (slong d = 0; d <= fmpq_poly_degree(G); d++) {
            fmpq_poly_get_coeff_fmpq(q, G, d);
            st |= gr_set_fmpq(c, q, F->K);
            st |= gr_poly_set_coeff_scalar(g, d, c, F->K);
        }
        gr_heap_clear(c, F->K);
    }
    fmpq_clear(q);
    fmpq_poly_clear(G); fmpq_poly_clear(e); fmpq_poly_clear(tp);
    return ok && st == GR_SUCCESS;
}

/* FLINT`PolynomialHermiteForm[M, x, opts]
 *   "Products" -> {gens, table}  the module of the products (pm_apply_products)
 *   "DeterminantMultiple" -> D  entries reduced modulo D during the elimination
 *   "Denominator" -> den        the module M/den: common factors of den and
 *                               every entry are divided out and {H, den'} is
 *                               returned (den' monic over the prime field).
 *                               A Hermite form divided by a monic common factor
 *                               is again one, so no second pass is needed. */
static Expr* builtin_flint_polynomialhermiteform(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr** a = res->data.function.args;
    PmField F;
    PmOpts xo;
    if (!pm_field_init(&F, a[1], a + 2, res->data.function.arg_count - 2, &xo))
        return NULL;
    Expr* out = NULL;
    gr_poly_t D, den, g; gr_poly_init(D, F.K); gr_poly_init(den, F.K); gr_poly_init(g, F.K);
    int haveD = 0, ok = 1;
    if (xo.D) {
        ok = pm_from_expr(D, xo.D, &F);
        if (ok && D->length >= 2) {          /* deg D >= 1: a usable modulus */
            haveD = gr_poly_make_monic(D, D, F.K) == GR_SUCCESS;
            ok = haveD;
        }
    }
    if (ok && xo.den) {
        ok = pm_from_expr(den, xo.den, &F) && den->length > 0;
        if (ok) ok = gr_poly_make_monic(den, den, F.K) == GR_SUCCESS;
        if (ok && F.kind == PM_NF) {         /* den must lie over Q */
            fmpq_poly_t chk; fmpq_poly_init(chk);
            ok = pm_prime_poly(chk, den, &F);
            fmpq_poly_clear(chk);
        }
    }
    PmMat M;
    if (ok && pm_mat_from_expr(&M, a[0], 0, &F)) {
        if (xo.products) ok = pm_apply_products(&M, xo.products, &F);
        slong c = M.c;
        if (ok && haveD) {
            /* append D e_k, reducing the given rows mod D first */
            slong r0 = M.r;
            PmMat W; pm_mat_init(&W, r0 + c, c, &F);
            int st = GR_SUCCESS;
            for (slong i = 0; i < r0; i++)
                for (slong k = 0; k < c; k++)
                    st |= gr_poly_rem(PM_AT(&W, i, k), PM_AT(&M, i, k), D, F.K);
            for (slong k = 0; k < c; k++) st |= gr_poly_set(PM_AT(&W, r0 + k, k), D, F.K);
            pm_mat_clear(&M, &F);
            M = W;
            ok = st == GR_SUCCESS;
        }
        ok = ok && M.r >= c && pm_hnf(&M, haveD ? D : NULL, &F);
        if (ok && xo.den) {
            /* divide out the common factors (repeatedly: a square in den) */
            int st = GR_SUCCESS;
            while (ok && den->length > 1) {
                if (!pm_common_factor(g, den, &M, c, &F)) { ok = 0; break; }
                if (g->length <= 1) break;
                for (slong i = 0; i < c * c; i++)
                    if (M.e[i].length > 0) st |= gr_poly_divexact(M.e + i, M.e + i, g, F.K);
                st |= gr_poly_divexact(den, den, g, F.K);
                if (st != GR_SUCCESS) ok = 0;
            }
        }
        if (ok) {
            Expr* H = pm_mat_to_expr(&M, c, &F);
            if (H && xo.den) {
                Expr* de = pm_to_expr(den, &F);
                if (de) out = pm_list((Expr*[]){ H, de }, 2);
                else expr_free(H);
            } else {
                out = H;
            }
        }
        pm_mat_clear(&M, &F);
    }
    gr_poly_clear(D, F.K); gr_poly_clear(den, F.K); gr_poly_clear(g, F.K);
    pm_field_clear(&F);
    return out;
}

/* FLINT`PolynomialDetAdjugate[M, x, opts] -> {det, adj} (adj . M = det Id). */
static Expr* builtin_flint_polynomialdetadjugate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 2) return NULL;
    Expr** a = res->data.function.args;
    PmField F;
    PmOpts xo;
    if (!pm_field_init(&F, a[1], a + 2, res->data.function.arg_count - 2, &xo)) return NULL;
    Expr* out = NULL;
    PmMat M;
    if (xo.D || xo.den) { pm_field_clear(&F); return NULL; }
    if (pm_mat_from_expr(&M, a[0], 0, &F)) {
        int prod_ok = !xo.products || pm_apply_products(&M, xo.products, &F);
        if (prod_ok && M.r == M.c) {
            slong n = M.r;
            gr_ctx_t Rx; gr_ctx_init_gr_poly(Rx, F.K);
            gr_mat_t A, adj;
            gr_mat_init(A, n, n, Rx); gr_mat_init(adj, n, n, Rx);
            gr_poly_t det; gr_poly_init(det, F.K);
            int st = GR_SUCCESS;
            for (slong i = 0; i < n; i++)
                for (slong j = 0; j < n; j++)
                    st |= gr_poly_set((gr_poly_struct*)gr_mat_entry_ptr(A, i, j, Rx),
                                      PM_AT(&M, i, j), F.K);
            if (n == 1) {
                st |= gr_poly_set(det, PM_AT(&M, 0, 0), F.K);
                st |= gr_poly_one((gr_poly_struct*)gr_mat_entry_ptr(adj, 0, 0, Rx), F.K);
            } else {
                st |= gr_mat_adjugate_cofactor(adj, det, A, Rx);
            }
            if (st == GR_SUCCESS) {
                PmMat Adj; pm_mat_init(&Adj, n, n, &F);
                for (slong i = 0; i < n; i++)
                    for (slong j = 0; j < n; j++)
                        st |= gr_poly_set(PM_AT(&Adj, i, j),
                              (const gr_poly_struct*)gr_mat_entry_ptr(adj, i, j, Rx), F.K);
                Expr* de = st == GR_SUCCESS ? pm_to_expr(det, &F) : NULL;
                Expr* ae = de ? pm_mat_to_expr(&Adj, n, &F) : NULL;
                if (de && ae) out = pm_list((Expr*[]){ de, ae }, 2);
                else { if (de) expr_free(de); if (ae) expr_free(ae); }
                pm_mat_clear(&Adj, &F);
            }
            gr_poly_clear(det, F.K);
            gr_mat_clear(A, Rx); gr_mat_clear(adj, Rx);
            gr_ctx_clear(Rx);
        }
        pm_mat_clear(&M, &F);
    }
    pm_field_clear(&F);
    return out;
}

/* FLINT`PolynomialQuotientRemainder[a, b, x, opts] -> {q, r}. */
static Expr* builtin_flint_polynomialquotientremainder(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** a = res->data.function.args;
    PmField F;
    if (!pm_field_init(&F, a[2], a + 3, res->data.function.arg_count - 3, NULL)) return NULL;
    Expr* out = NULL;
    gr_poly_t A, B, Q, R;
    gr_poly_init(A, F.K); gr_poly_init(B, F.K); gr_poly_init(Q, F.K); gr_poly_init(R, F.K);
    if (pm_from_expr(A, a[0], &F) && pm_from_expr(B, a[1], &F) && B->length > 0 &&
        gr_poly_divrem(Q, R, A, B, F.K) == GR_SUCCESS) {
        Expr* qe = pm_to_expr(Q, &F);
        Expr* re = qe ? pm_to_expr(R, &F) : NULL;
        if (qe && re) out = pm_list((Expr*[]){ qe, re }, 2);
        else if (qe) expr_free(qe);
    }
    gr_poly_clear(A, F.K); gr_poly_clear(B, F.K); gr_poly_clear(Q, F.K); gr_poly_clear(R, F.K);
    pm_field_clear(&F);
    return out;
}

/* FLINT`PolynomialExtendedGCD[a, b, x, opts] -> {g, {s, t}}, g monic (0 when
 * a = b = 0), s a + t b = g. */
static Expr* builtin_flint_polynomialextendedgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** a = res->data.function.args;
    PmField F;
    if (!pm_field_init(&F, a[2], a + 3, res->data.function.arg_count - 3, NULL)) return NULL;
    Expr* out = NULL;
    gr_poly_t A, B, G, S, T;
    gr_poly_init(A, F.K); gr_poly_init(B, F.K); gr_poly_init(G, F.K);
    gr_poly_init(S, F.K); gr_poly_init(T, F.K);
    if (pm_from_expr(A, a[0], &F) && pm_from_expr(B, a[1], &F)) {
        int st = GR_SUCCESS;
        if (A->length == 0 && B->length == 0) {
            st |= gr_poly_zero(G, F.K); st |= gr_poly_zero(S, F.K); st |= gr_poly_zero(T, F.K);
        } else {
            st |= PM_GR_POLY_XGCD(G, S, T, A, B, F.K);
            /* normalise: g monic (FLINT already does; make it a contract) */
            gr_ptr lc = gr_heap_init(F.K);
            st |= gr_inv(lc, PM_ENTRY(G->coeffs, G->length - 1, &F), F.K);
            st |= gr_poly_mul_scalar(G, G, lc, F.K);
            st |= gr_poly_mul_scalar(S, S, lc, F.K);
            st |= gr_poly_mul_scalar(T, T, lc, F.K);
            gr_heap_clear(lc, F.K);
        }
        if (st == GR_SUCCESS) {
            Expr* ge = pm_to_expr(G, &F);
            Expr* se = ge ? pm_to_expr(S, &F) : NULL;
            Expr* te = se ? pm_to_expr(T, &F) : NULL;
            if (ge && se && te)
                out = pm_list((Expr*[]){ ge, pm_list((Expr*[]){ se, te }, 2) }, 2);
            else { if (ge) expr_free(ge); if (se) expr_free(se); }
        }
    }
    gr_poly_clear(A, F.K); gr_poly_clear(B, F.K); gr_poly_clear(G, F.K);
    gr_poly_clear(S, F.K); gr_poly_clear(T, F.K);
    pm_field_clear(&F);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Shifted weak Popov form                                            */
/* ------------------------------------------------------------------ */

/* The key of entry j of a row: scale deg(e_j) + shift_j (-infinity for 0).
 * The leading column of a row is the column of greatest key (the first one on
 * a tie); -1 for the zero row. */
static slong pm_lead(const PmMat* M, slong i, slong scale, const slong* shift, slong* key) {
    slong col = -1, best = 0;
    for (slong j = 0; j < M->c; j++) {
        slong len = PM_AT(M, i, j)->length;
        if (len == 0) continue;
        slong k = scale * (len - 1) + shift[j];
        if (col < 0 || k > best) { col = j; best = k; }
    }
    *key = best;
    return col;
}

/* The package's reducedRows (src/internal/trager/Jacobian.m), step for step:
 * while two rows share a leading column, take the FIRST such pair (i < j in
 * lexicographic order, as Subsets[Range[n], {2}] lists them), and cancel the
 * leading term of the row of greater key (lo) by f x^s times the other (hi),
 * s = (key_lo - key_hi)/scale, f the ratio of the leading coefficients.  The
 * form is not unique; following the same sequence makes the result identical
 * to the package's.  Returns 0 past `cap` steps or on a FLINT failure. */
static int pm_weak_popov(PmMat* M, slong scale, const slong* shift, slong cap, PmField* F) {
    gr_ctx_struct* K = F->K;
    slong n = M->r;
    slong* col = malloc((size_t)(n > 0 ? n : 1) * sizeof(slong));
    slong* key = malloc((size_t)(n > 0 ? n : 1) * sizeof(slong));
    gr_ptr f = gr_heap_init(K);
    gr_poly_t tmp; gr_poly_init(tmp, K);
    int ok = 1, st = GR_SUCCESS;
    for (slong i = 0; i < n; i++) col[i] = pm_lead(M, i, scale, shift, key + i);
    for (slong it = 0; ok; it++) {
        if (it > cap) { ok = 0; break; }
        slong ci = -1, cj = -1;
        for (slong i = 0; i < n && ci < 0; i++)
            for (slong j = i + 1; j < n; j++)
                if (col[i] >= 0 && col[i] == col[j]) { ci = i; cj = j; break; }
        if (ci < 0) break;
        slong lo = key[ci] >= key[cj] ? ci : cj, hi = lo == ci ? cj : ci;
        slong c = col[ci];
        slong diff = key[lo] - key[hi];
        if (diff % scale != 0) { ok = 0; break; }
        slong sh = diff / scale;
        gr_poly_struct* el = PM_AT(M, lo, c);
        gr_poly_struct* eh = PM_AT(M, hi, c);
        st |= gr_div(f, PM_ENTRY(el->coeffs, el->length - 1, F),
                     PM_ENTRY(eh->coeffs, eh->length - 1, F), K);
        for (slong k = 0; k < M->c; k++) {
            gr_poly_struct* rh = PM_AT(M, hi, k);
            if (rh->length == 0) continue;
            st |= gr_poly_shift_left(tmp, rh, sh, K);
            st |= gr_poly_mul_scalar(tmp, tmp, f, K);
            st |= gr_poly_sub(PM_AT(M, lo, k), PM_AT(M, lo, k), tmp, K);
        }
        if (st != GR_SUCCESS) { ok = 0; break; }
        col[lo] = pm_lead(M, lo, scale, shift, key + lo);
    }
    gr_poly_clear(tmp, K);
    gr_heap_clear(f, K);
    free(col); free(key);
    return ok;
}

/* FLINT`PolynomialWeakPopovForm[M, x, shifts, opts] */
static Expr* builtin_flint_polynomialweakpopovform(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** a = res->data.function.args;
    const Expr* sh = a[2];
    if (!pm_is_head(sh, "List")) return NULL;
    /* split our own option ("DegreeScale" -> k) from the field options */
    size_t no = res->data.function.arg_count - 3;
    Expr** fo = malloc((no ? no : 1) * sizeof(Expr*));
    size_t nfo = 0;
    slong scale = 1;
    for (size_t i = 0; i < no; i++) {
        Expr* r = a[3 + i];
        if (pm_is_head(r, "Rule") && r->data.function.arg_count == 2 &&
            r->data.function.args[0]->type == EXPR_STRING &&
            strcmp(r->data.function.args[0]->data.string, "DegreeScale") == 0) {
            const Expr* v = r->data.function.args[1];
            if (v->type != EXPR_INTEGER || v->data.integer < 1) { free(fo); return NULL; }
            scale = (slong)v->data.integer;
        } else {
            fo[nfo++] = r;
        }
    }
    PmField F;
    int okf = pm_field_init(&F, a[1], fo, nfo, NULL);
    free(fo);
    if (!okf) return NULL;
    Expr* out = NULL;
    PmMat M;
    if (pm_mat_from_expr(&M, a[0], 0, &F)) {
        slong c = M.c;
        slong* shift = malloc((size_t)c * sizeof(slong));
        int ok = (slong)sh->data.function.arg_count == c;
        slong degsum = 0;
        for (slong j = 0; ok && j < c; j++) {
            const Expr* v = sh->data.function.args[j];
            if (v->type != EXPR_INTEGER) ok = 0;
            else shift[j] = (slong)v->data.integer;
        }
        for (slong i = 0; ok && i < M.r * c; i++)
            if (M.e[i].length > 0) degsum += M.e[i].length - 1;
        /* the package's cap: 4 (total degree + rows) + 200 */
        if (ok && pm_weak_popov(&M, scale, shift, 4 * (degsum + M.r) + 200, &F))
            out = pm_mat_to_expr(&M, M.r, &F);
        free(shift);
        pm_mat_clear(&M, &F);
    }
    pm_field_clear(&F);
    return out;
}

/* ------------------------------------------------------------------ */
/*  gcd over K[x]: modular over Q(theta)                               */
/* ------------------------------------------------------------------ */

/* The Euclidean gcd over a number field is fraction-based: its remainders'
 * coefficients explode with the field degree (a gcd of a sextic and a degree-36
 * polynomial over a degree-36 field took seconds).  The bridge's certified
 * modular engine (Encarnacion; fb_field_mpoly_gcd) does not, but wants an
 * algebraic INTEGER generator: theta' = d theta with d the lcm of the
 * denominators of mu, whose minimal polynomial d^n mu(T/d) is monic integral.
 * A, B are mapped t -> T/d, the gcd taken in Q[x, T]/(mu_d), and mapped back.
 * Returns 1 with G monic, 0 on decline (the caller falls back). */
static int pm_nf_modular_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, PmField* F) {
    slong n = F->n;
    fmpz_t d; fmpz_init(d);
    fmpz_one(d);
    for (slong i = 0; i <= n; i++) {
        fmpq_t c; fmpq_init(c);
        fmpq_poly_get_coeff_fmpq(c, F->mu, i);
        fmpz_lcm(d, d, fmpq_denref(c));
        fmpq_clear(c);
    }
    /* mu_d(T) = d^n mu(T/d): coefficient of T^i is c_i d^(n-i) (c_n = 1) */
    Expr** cs = malloc((size_t)(n + 1) * sizeof(Expr*));
    fmpz_t z, dp; fmpz_init(z); fmpz_init(dp);
    fmpq_t c, w; fmpq_init(c); fmpq_init(w);
    int ok = 1;
    for (slong i = 0; i <= n; i++) {
        fmpq_poly_get_coeff_fmpq(c, F->mu, i);
        fmpz_pow_ui(dp, d, (ulong)(n - i));
        fmpq_mul_fmpz(w, c, dp);
        if (!fmpz_is_one(fmpq_denref(w))) ok = 0;
        fmpq_t one; fmpq_init(one); fmpq_set_fmpz(one, fmpq_numref(w));
        cs[i] = fb_fmpq_to_expr(one);
        fmpq_clear(one);
    }
    Expr* mp = expr_new_function(expr_new_symbol(SYM_List), cs, (size_t)(n + 1));
    free(cs);
    const fmpq_mpoly_ctx_struct* ctx2 = fb_ctx_lex(2);    /* gen0 = x, gen1 = T */
    fmpq_mpoly_t M, Am, Bm, Gm;
    fmpq_mpoly_init(M, ctx2); fmpq_mpoly_init(Am, ctx2);
    fmpq_mpoly_init(Bm, ctx2); fmpq_mpoly_init(Gm, ctx2);
    ok = ok && fb_field_build_minpoly(M, mp, ctx2, 2, 1);
    /* K-poly -> Q[x, T]: t^j = T^j / d^j */
    fmpq_poly_t tp; fmpq_poly_init(tp);
    for (int which = 0; ok && which < 2; which++) {
        const gr_poly_struct* P = which ? B : A;
        fmpq_mpoly_struct* Q = which ? Bm : Am;
        ulong ex[2];
        for (slong k = 0; k < P->length; k++) {
            nf_elem_get_fmpq_poly(tp, (const nf_elem_struct*)PM_ENTRY(P->coeffs, k, F), F->nf);
            for (slong j = 0; j < fmpq_poly_length(tp); j++) {
                fmpq_poly_get_coeff_fmpq(c, tp, j);
                if (fmpq_is_zero(c)) continue;
                fmpz_pow_ui(dp, d, (ulong)j);
                fmpq_div_fmpz(w, c, dp);
                ex[0] = (ulong)k; ex[1] = (ulong)j;
                fmpq_mpoly_set_coeff_fmpq_ui(Q, w, ex, ctx2);
            }
        }
        if (fmpq_mpoly_is_zero(Q, ctx2)) ok = 0;
    }
    ok = ok && fb_field_mpoly_gcd(Gm, Am, Bm, M, ctx2, mp, 1, 1, n);
    if (ok) {
        /* back: T^j = d^j t^j, coefficient of x^k collected per k */
        slong dx = fmpq_mpoly_degree_si(Gm, 0, ctx2);
        fmpq_poly_struct* cf = malloc((size_t)(dx + 1) * sizeof(fmpq_poly_struct));
        for (slong k = 0; k <= dx; k++) fmpq_poly_init(cf + k);
        ulong ex[2];
        for (slong i = 0; i < fmpq_mpoly_length(Gm, ctx2); i++) {
            fmpq_mpoly_get_term_exp_ui(ex, Gm, i, ctx2);
            fmpq_mpoly_get_term_coeff_fmpq(c, Gm, i, ctx2);
            fmpz_pow_ui(dp, d, ex[1]);
            fmpq_mul_fmpz(w, c, dp);
            fmpq_poly_set_coeff_fmpq(cf + ex[0], (slong)ex[1], w);
        }
        int st = gr_poly_zero(G, F->K);
        gr_ptr e = gr_heap_init(F->K);
        for (slong k = 0; k <= dx; k++) {
            if (fmpq_poly_is_zero(cf + k)) continue;
            pm_set_elem(e, cf + k, F);
            st |= gr_poly_set_coeff_scalar(G, k, e, F->K);
        }
        gr_heap_clear(e, F->K);
        for (slong k = 0; k <= dx; k++) fmpq_poly_clear(cf + k);
        free(cf);
        if (G->length > 0) st |= gr_poly_make_monic(G, G, F->K);
        ok = st == GR_SUCCESS && G->length > 0;
    }
    fmpq_poly_clear(tp);
    fmpq_mpoly_clear(M, ctx2); fmpq_mpoly_clear(Am, ctx2);
    fmpq_mpoly_clear(Bm, ctx2); fmpq_mpoly_clear(Gm, ctx2);
    expr_free(mp);
    fmpq_clear(c); fmpq_clear(w); fmpz_clear(z); fmpz_clear(dp); fmpz_clear(d);
    return ok;
}

/* FLINT`PolynomialFieldGCD[a, b, x, opts] -> the monic gcd over K[x] (0 when
 * a = b = 0).  Over Q(theta) the certified modular engine when both operands
 * have positive degree; the Euclidean gcd otherwise or on its decline. */
static Expr* builtin_flint_polynomialfieldgcd(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 3) return NULL;
    Expr** a = res->data.function.args;
    PmField F;
    if (!pm_field_init(&F, a[2], a + 3, res->data.function.arg_count - 3, NULL)) return NULL;
    Expr* out = NULL;
    gr_poly_t A, B, G;
    gr_poly_init(A, F.K); gr_poly_init(B, F.K); gr_poly_init(G, F.K);
    if (pm_from_expr(A, a[0], &F) && pm_from_expr(B, a[1], &F)) {
        int ok = 0;
        if (F.kind == PM_NF && A->length > 1 && B->length > 1)
            ok = pm_nf_modular_gcd(G, A, B, &F);
        if (!ok) {
            ok = gr_poly_gcd(G, A, B, F.K) == GR_SUCCESS;
            if (ok && G->length > 0) ok = gr_poly_make_monic(G, G, F.K) == GR_SUCCESS;
        }
        if (ok) out = pm_to_expr(G, &F);
    }
    gr_poly_clear(A, F.K); gr_poly_clear(B, F.K); gr_poly_clear(G, F.K);
    pm_field_clear(&F);
    return out;
}

static void pm_cleanup(void) {
    if (g_pm_irred_ready) { fmpq_poly_clear(g_pm_irred_mu); g_pm_irred_ready = 0; }
}

void flint_polymat_init(void) {
    atexit(pm_cleanup);

    symtab_add_builtin(SYM_FLINT_PolynomialHermiteForm, builtin_flint_polynomialhermiteform);
    symtab_get_def(SYM_FLINT_PolynomialHermiteForm)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialHermiteForm,
        "FLINT`PolynomialHermiteForm[m, x, opts] gives the Hermite normal form of "
        "the module over K[x] spanned by the rows of the polynomial matrix m: "
        "upper triangular, monic pivots, entries above each pivot of lower "
        "degree, zero rows dropped (the module must have full column rank). "
        "K is Q, GF(p) with Modulus -> p, or Q[t]/(mu) with Extension -> {t, mu}; "
        "\"DeterminantMultiple\" -> D (a multiple of the determinant of the "
        "module) reduces entries modulo D during the elimination; \"Products\" -> "
        "{gens, table} takes the module spanned by the products of the rows with "
        "gens in the algebra with basis w_0..w_(m-1) and w_i w_j = f w_k "
        "(table[[i+1, j+1]] = {k, f}); \"Denominator\" -> den gives {H, den'} "
        "for the module m/den in lowest terms (factors of den, over the prime "
        "field, dividing every entry removed). Computed via FLINT gr_poly. "
        "Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_PolynomialWeakPopovForm, builtin_flint_polynomialweakpopovform);
    symtab_get_def(SYM_FLINT_PolynomialWeakPopovForm)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialWeakPopovForm,
        "FLINT`PolynomialWeakPopovForm[m, x, shifts, opts] reduces the rows of "
        "the polynomial matrix m over K[x] (same module) until their leading "
        "columns are distinct: the leading column of a row maximises "
        "k deg(e_j) + shifts[[j]], k given by \"DegreeScale\" -> k (default 1), "
        "and two rows with the same leading column are reduced by cancelling the "
        "leading term of the one of greater key against x^s times the other. "
        "K and the field opts as for FLINT`PolynomialHermiteForm. Computed via "
        "FLINT gr_poly. Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_PolynomialDetAdjugate, builtin_flint_polynomialdetadjugate);
    symtab_get_def(SYM_FLINT_PolynomialDetAdjugate)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialDetAdjugate,
        "FLINT`PolynomialDetAdjugate[m, x, opts] gives {det, adj} for the square "
        "polynomial matrix m over K[x]: its determinant and its adjugate, with "
        "adj . m == det IdentityMatrix. K and the field opts as for "
        "FLINT`PolynomialHermiteForm, and \"Products\" -> {gens, table} as there "
        "(m replaced by the products of its rows with gens). Computed via FLINT "
        "gr_mat. Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_PolynomialQuotientRemainder,
                       builtin_flint_polynomialquotientremainder);
    symtab_get_def(SYM_FLINT_PolynomialQuotientRemainder)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialQuotientRemainder,
        "FLINT`PolynomialQuotientRemainder[a, b, x, opts] gives {q, r} with "
        "a == q b + r and deg r < deg b, over K[x]: K is Q, GF(p) with "
        "Modulus -> p, or Q[t]/(mu) with Extension -> {t, mu} (coefficients "
        "polynomials in t, reduced modulo mu). Computed via FLINT gr_poly. "
        "Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_PolynomialFieldGCD, builtin_flint_polynomialfieldgcd);
    symtab_get_def(SYM_FLINT_PolynomialFieldGCD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialFieldGCD,
        "FLINT`PolynomialFieldGCD[a, b, x, opts] gives the monic greatest common "
        "divisor of a and b over K[x]. K and opts as for "
        "FLINT`PolynomialQuotientRemainder; over Q[t]/(mu) the certified modular "
        "algorithm is used (no coefficient growth with the field degree). "
        "Computed via FLINT. Returns unevaluated if out of scope.");

    symtab_add_builtin(SYM_FLINT_PolynomialExtendedGCD, builtin_flint_polynomialextendedgcd);
    symtab_get_def(SYM_FLINT_PolynomialExtendedGCD)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolynomialExtendedGCD,
        "FLINT`PolynomialExtendedGCD[a, b, x, opts] gives {g, {s, t}} with g the "
        "monic greatest common divisor of a and b over K[x] and s a + t b == g. "
        "K and opts as for FLINT`PolynomialQuotientRemainder. Computed via FLINT "
        "gr_poly. Returns unevaluated if out of scope.");
}

#else /* !USE_FLINT */

void flint_polymat_init(void) {}

#endif /* USE_FLINT */
