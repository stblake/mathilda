/*
 * reduce_univar.c
 *
 * Univariate real sign-diagram engine for `Reduce` (REDUCE_PLAN.md, Phase 2).
 * See reduce_univar.h for the algorithm.  The exact sign oracle, rational sample
 * selection, polynomial test, and Solve-based real-root isolation live in the
 * shared reduce_real_util.{c,h} (also used by the multivariate CAD); this file
 * is the 1-D cell construction, truth evaluation, and segment emission on top.
 */
#include "reduce_univar.h"
#include "reduce_real_util.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ *
 *  Truth of the formula at a sample point                             *
 * ------------------------------------------------------------------ */

/* Truth of one atom at the sample: 1 true, 0 false, -1 undecided. */
static int atom_truth_at(const RAtom* a, const Expr* x, const Expr* sample) {
    if (a->rel == R_ELEM) return -1;
    int s = rru_poly_sign_at(a->poly, x, sample);
    if (s == -2) return -1;
    switch (a->rel) {
        case R_EQ: return s == 0;
        case R_NE: return s != 0;
        case R_LT: return s < 0;
        case R_LE: return s <= 0;
        default:   return -1;
    }
}

/* Truth of the whole DNF formula at the sample: 1, 0, or -1 (undecided). */
static int form_truth_at(const RForm* F, const Expr* x, const Expr* sample) {
    if (F->is_true) return 1;
    for (int i = 0; i < F->n; i++) {
        RConj* c = F->c[i];
        int conj = 1;
        for (int k = 0; k < c->n; k++) {
            int t = atom_truth_at(&c->a[k], x, sample);
            if (t == -1) return -1;
            if (t == 0) { conj = 0; break; }
        }
        if (conj) return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ *
 *  Segment emission                                                   *
 * ------------------------------------------------------------------ */

static Expr* seg_to_expr(const Expr* lo, bool lo_open,
                         const Expr* hi, bool hi_open, const Expr* x) {
    bool has_lo = (lo != NULL), has_hi = (hi != NULL);
    if (!has_lo && !has_hi) return expr_new_symbol(SYM_True);
    if (!has_lo) {
        const char* op = hi_open ? SYM_Less : SYM_LessEqual;
        return expr_new_function(expr_new_symbol(op),
            (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)hi) }, 2);
    }
    if (!has_hi) {
        const char* op = lo_open ? SYM_Greater : SYM_GreaterEqual;
        return expr_new_function(expr_new_symbol(op),
            (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)lo) }, 2);
    }
    if (expr_eq(lo, hi)) {                  /* isolated point: x == b */
        return expr_new_function(expr_new_symbol(SYM_Equal),
            (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)lo) }, 2);
    }
    const char* op1 = lo_open ? SYM_Less : SYM_LessEqual;
    const char* op2 = hi_open ? SYM_Less : SYM_LessEqual;
    return expr_new_function(expr_new_symbol(SYM_Inequality),
        (Expr*[]){ expr_copy((Expr*)lo), expr_new_symbol(op1),
                   expr_copy((Expr*)x),  expr_new_symbol(op2),
                   expr_copy((Expr*)hi) }, 5);
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

/* Gather the distinct real breakpoints (roots of every polynomial in F) in
 * ascending order, deduped.  Returns an owned array of length *m_out (caller
 * frees the entries and the array); sets *ok=false to decline (a non-polynomial
 * atom, an ELEM, a cleared variable denominator, or an undecidable ordering). */
static Expr** collect_breakpoints(const RForm* F, const Expr* x, int* m_out, bool* ok) {
    *ok = true; *m_out = 0;
    Expr** polys = NULL; int np = 0, pcap = 0;
    for (int i = 0; *ok && i < F->n; i++) {
        RConj* c = F->c[i];
        for (int k = 0; *ok && k < c->n; k++) {
            RAtom* a = &c->a[k];
            if (a->rel == R_ELEM || a->nonconst_denom) { *ok = false; break; }
            bool dup = false;
            for (int j = 0; j < np; j++) if (expr_eq(polys[j], a->poly)) { dup = true; break; }
            if (dup) continue;
            if (!rru_is_polynomial(a->poly, x)) { *ok = false; break; }
            if (np == pcap) { pcap = pcap ? pcap * 2 : 8; polys = realloc(polys, (size_t)pcap * sizeof(Expr*)); }
            polys[np++] = a->poly;          /* borrowed */
        }
    }
    Expr** roots = NULL; int nr = 0, rcap = 0;
    for (int j = 0; *ok && j < np; j++)
        if (!rru_collect_roots(polys[j], x, &roots, &nr, &rcap, NULL, 0)) *ok = false;
    free(polys);
    if (!*ok) { for (int i = 0; i < nr; i++) expr_free(roots[i]); free(roots); return NULL; }
    for (int i = 0; *ok && i < nr; i++) {
        int mn = i;
        for (int j = i + 1; j < nr; j++) {
            int c = rru_sign_compare(roots[j], roots[mn]);
            if (c == -2) { *ok = false; break; }
            if (c < 0) mn = j;
        }
        if (*ok && mn != i) { Expr* t = roots[i]; roots[i] = roots[mn]; roots[mn] = t; }
    }
    int m = 0;
    Expr** bp = malloc((size_t)(nr + 1) * sizeof(Expr*));
    for (int i = 0; *ok && i < nr; i++) {
        if (m == 0) { bp[m++] = roots[i]; continue; }
        int c = rru_sign_compare(bp[m - 1], roots[i]);
        if (c == -2) {
            *ok = false;
            expr_free(roots[i]);
            for (int j = i + 1; j < nr; j++) expr_free(roots[j]);
            break;
        }
        if (c == 0) { expr_free(roots[i]); continue; }
        bp[m++] = roots[i];
    }
    free(roots);
    if (!*ok) { for (int i = 0; i < m; i++) expr_free(bp[i]); free(bp); return NULL; }
    *m_out = m;
    return bp;
}

Expr* reduce_univar(const RForm* F, const Expr* x, Expr** vars, int nv) {
    (void)vars; (void)nv;
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    bool ok = true; int m = 0;
    Expr** bp = collect_breakpoints(F, x, &m, &ok);
    if (!ok) return NULL;

    /* 4. No breakpoints: the whole line is one cell. */
    if (m == 0) {
        Expr* sample = expr_new_integer(0);
        int t = form_truth_at(F, x, sample);
        expr_free(sample); free(bp);
        if (t == -1) return NULL;
        return expr_new_symbol(t ? SYM_True : SYM_False);
    }

    /* 5. Truth at one sample per cell (2m+1 cells, alternating interval/point). */
    int ncells = 2 * m + 1;
    int* truth = malloc(sizeof(int) * (size_t)ncells);
    for (int idx = 0; ok && idx < ncells; idx++) {
        Expr* sample;
        if (idx % 2 == 1) {                 /* point P_j */
            sample = expr_copy(bp[(idx + 1) / 2 - 1]);
        } else {                            /* interval I_j */
            int j = idx / 2;
            const Expr* lo = (j == 0) ? NULL : bp[j - 1];
            const Expr* hi = (j == m) ? NULL : bp[j];
            sample = rru_rational_between(lo, hi);
            if (!sample) { ok = false; break; }
        }
        int t = form_truth_at(F, x, sample);
        expr_free(sample);
        if (t == -1) { ok = false; break; }
        truth[idx] = t;
    }
    if (!ok) { for (int i = 0; i < m; i++) expr_free(bp[i]); free(bp); free(truth); return NULL; }

    /* 6. Emit the union of satisfying cells. */
    bool all_true = true, all_false = true;
    for (int i = 0; i < ncells; i++) { if (truth[i]) all_false = false; else all_true = false; }

    Expr* result;
    if (all_true) {
        result = expr_new_symbol(SYM_True);
    } else if (all_false) {
        result = expr_new_symbol(SYM_False);
    } else {
        bool all_int_true = true;
        for (int j = 0; j <= m; j++) if (!truth[2 * j]) { all_int_true = false; break; }

        if (all_int_true) {
            /* cofinite: everything except finitely many excluded points. */
            Expr** parts = NULL; int npar = 0;
            for (int j = 1; j <= m; j++) {
                if (truth[2 * j - 1]) continue;
                Expr* ne = expr_new_function(expr_new_symbol(SYM_Unequal),
                    (Expr*[]){ expr_copy((Expr*)x), expr_copy(bp[j - 1]) }, 2);
                parts = realloc(parts, (size_t)(npar + 1) * sizeof(Expr*));
                parts[npar++] = ne;
            }
            result = (npar == 1) ? parts[0]
                : expr_new_function(expr_new_symbol(SYM_And), parts, (size_t)npar);
            free(parts);
        } else {
            /* merge maximal runs of true cells into interval/point segments. */
            Expr** segs = NULL; int nseg = 0;
            bool active = false; const Expr* lo = NULL; bool lo_open = false;
            for (int idx = 0; idx < ncells; idx++) {
                bool is_int = (idx % 2 == 0);
                if (truth[idx] && !active) {
                    active = true;
                    if (is_int) {
                        int j = idx / 2;
                        if (j == 0) { lo = NULL; lo_open = false; }
                        else { lo = bp[j - 1]; lo_open = true; }
                    } else { lo = bp[(idx + 1) / 2 - 1]; lo_open = false; }
                } else if (!truth[idx] && active) {
                    active = false;
                    const Expr* hi; bool hi_open;
                    if (is_int) { hi = bp[idx / 2 - 1]; hi_open = false; }
                    else { hi = bp[(idx + 1) / 2 - 1]; hi_open = true; }
                    Expr* seg = seg_to_expr(lo, lo_open, hi, hi_open, x);
                    segs = realloc(segs, (size_t)(nseg + 1) * sizeof(Expr*));
                    segs[nseg++] = seg;
                }
            }
            if (active) {
                Expr* seg = seg_to_expr(lo, lo_open, NULL, false, x);
                segs = realloc(segs, (size_t)(nseg + 1) * sizeof(Expr*));
                segs[nseg++] = seg;
            }
            result = (nseg == 1) ? segs[0]
                : expr_new_function(expr_new_symbol(SYM_Or), segs, (size_t)nseg);
            free(segs);
        }
    }

    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp); free(truth);
    return eval_and_free(result);   /* flatten And/Or, fuse Inequality chains */
}

/* ------------------------------------------------------------------ *
 *  Integer domain (bounded enumeration)                               *
 * ------------------------------------------------------------------ */

Expr* reduce_univar_integers(const RForm* F, const Expr* x, Expr** vars, int nv) {
    (void)vars; (void)nv;
    if (F->is_true) return expr_new_symbol(SYM_True);   /* every integer */
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    bool ok = true; int m = 0;
    Expr** bp = collect_breakpoints(F, x, &m, &ok);
    if (!ok) return NULL;

    /* No sign changes: the formula is constant over the whole line. */
    if (m == 0) {
        Expr* z = expr_new_integer(0);
        int t = form_truth_at(F, x, z);
        expr_free(z); free(bp);
        if (t == -1) return NULL;
        return expr_new_symbol(t ? SYM_True : SYM_False);
    }

    /* The integer solution set is bounded only if BOTH unbounded end-cells are
     * unsatisfied; otherwise it is infinite and we decline (a later phase can
     * emit an `x >= k`-style form). */
    bool ap, bp_ok;
    double lo = rru_approx_double(bp[0], &ap);
    double hi = rru_approx_double(bp[m - 1], &bp_ok);
    if (!ap || !bp_ok) { for (int i = 0; i < m; i++) expr_free(bp[i]); free(bp); return NULL; }
    long a = (long)floor(lo) - 1;
    long b = (long)ceil(hi) + 1;

    Expr* sl = expr_new_integer(a);
    Expr* sr = expr_new_integer(b);
    int tl = form_truth_at(F, x, sl);
    int tr = form_truth_at(F, x, sr);
    expr_free(sl); expr_free(sr);
    if (tl == -1 || tr == -1 || tl == 1 || tr == 1) {
        for (int i = 0; i < m; i++) expr_free(bp[i]);
        free(bp);
        return NULL;                    /* undecidable or unbounded -> decline */
    }
    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp);

    /* Enumerate the integers in [a, b] and keep the satisfying ones. */
    Expr** sols = NULL; int ns = 0;
    for (long n = a; n <= b; n++) {
        Expr* s = expr_new_integer(n);
        int t = form_truth_at(F, x, s);
        if (t == -1) {
            expr_free(s);
            for (int i = 0; i < ns; i++) expr_free(sols[i]);
            free(sols);
            return NULL;
        }
        if (t == 1) {
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)x), s }, 2);
            sols = realloc(sols, (size_t)(ns + 1) * sizeof(Expr*));
            sols[ns++] = eq;
        } else {
            expr_free(s);
        }
    }

    Expr* out;
    if (ns == 0) { free(sols); out = expr_new_symbol(SYM_False); }
    else if (ns == 1) { out = sols[0]; free(sols); }
    else { out = expr_new_function(expr_new_symbol(SYM_Or), sols, (size_t)ns); free(sols); }
    return eval_and_free(out);
}
