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
    int s;
    if (a->denom) {
        /* Rational atom p/q REL 0.  A pole (q == 0) makes p/q undefined, so the
         * point is excluded from the solution set (false for every relation).
         * Off a pole, sign(p/q) = sign(p) * sign(q). */
        int sq = rru_poly_sign_at(a->denom, x, sample);
        if (sq == -2) return -1;
        if (sq == 0)  return 0;
        int sp = rru_poly_sign_at(a->poly, x, sample);
        if (sp == -2) return -1;
        s = sp * sq;
    } else {
        s = rru_poly_sign_at(a->poly, x, sample);
        if (s == -2) return -1;
    }
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
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

/* Gather the distinct real breakpoints of F in ascending order, deduped: the
 * roots of every numerator AND, for a rational atom p/q, the roots of q (the
 * poles), so a rational relation is decided on a diagram whose cells are the
 * intervals on which p/q is continuous and single-signed.  Returns an owned
 * array of length *m_out (caller frees the entries and the array); sets
 * *ok=false to decline (a non-polynomial numerator or denominator, an ELEM, or
 * an undecidable ordering). */
static Expr** collect_breakpoints(const RForm* F, const Expr* x, int* m_out, bool* ok) {
    *ok = true; *m_out = 0;
    Expr** polys = NULL; int np = 0, pcap = 0;
    for (int i = 0; *ok && i < F->n; i++) {
        RConj* c = F->c[i];
        for (int k = 0; *ok && k < c->n; k++) {
            RAtom* a = &c->a[k];
            if (a->rel == R_ELEM) { *ok = false; break; }
            /* The numerator, and the denominator when the atom is rational,
             * both contribute breakpoints (roots and poles). */
            Expr* cand[2]; int ncand = 0;
            cand[ncand++] = a->poly;
            if (a->denom) cand[ncand++] = a->denom;
            for (int ci = 0; *ok && ci < ncand; ci++) {
                Expr* p = cand[ci];
                bool dup = false;
                for (int j = 0; j < np; j++) if (expr_eq(polys[j], p)) { dup = true; break; }
                if (dup) continue;
                if (!rru_is_polynomial(p, x)) { *ok = false; break; }
                if (np == pcap) { pcap = pcap ? pcap * 2 : 8; polys = realloc(polys, (size_t)pcap * sizeof(Expr*)); }
                polys[np++] = p;            /* borrowed */
            }
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

    /* 6. Emit the union of satisfying cells (shared with reduce_realdiag). */
    Expr* result = rru_emit_sign_diagram(bp, m, truth, x);

    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp); free(truth);
    return result;
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

    /* Sample truth at every integer in a window [A, B] that strictly brackets
     * all breakpoints, so T is constant on each tail: T[A] is the truth of the
     * whole ray x <= A, T[B] of x >= B.  The set is a union of a possible left
     * ray, isolated middle points, and a possible right ray. */
    bool ap, bp_ok;
    double lo = rru_approx_double(bp[0], &ap);
    double hi = rru_approx_double(bp[m - 1], &bp_ok);
    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp);
    if (!ap || !bp_ok || !isfinite(lo) || !isfinite(hi)) return NULL;

    long A = (long)floor(lo) - 1;
    long B = (long)ceil(hi) + 1;
    if (B < A || B - A > 200000) return NULL;   /* absurd window -> decline */

    long win = B - A + 1;
    int* T = malloc((size_t)win * sizeof(int));
    for (long k = A; k <= B; k++) {
        Expr* s = expr_new_integer(k);
        int t = form_truth_at(F, x, s);
        expr_free(s);
        if (t == -1) { free(T); return NULL; }
        T[k - A] = t;
    }

    bool has_left  = (T[0] == 1);
    bool has_right = (T[win - 1] == 1);

    long lend = A - 1;                 /* last integer of the left ray */
    if (has_left)  { long k = A; while (k <= B && T[k - A] == 1) k++; lend = k - 1; }
    long rstart = B + 1;               /* first integer of the right ray */
    if (has_right) { long k = B; while (k >= A && T[k - A] == 1) k--; rstart = k + 1; }

    /* Both tails satisfied and their runs meet: the whole line qualifies. */
    if (has_left && has_right && lend >= rstart - 1) { free(T); return expr_new_symbol(SYM_True); }

    Expr** parts = NULL; int np = 0;
    if (has_left) {
        parts = realloc(parts, (size_t)(np + 1) * sizeof(Expr*));
        parts[np++] = expr_new_function(expr_new_symbol(SYM_LessEqual),
            (Expr*[]){ expr_copy((Expr*)x), expr_new_integer(lend) }, 2);
    }
    long mid_lo = has_left  ? lend + 1  : A;
    long mid_hi = has_right ? rstart - 1 : B;
    for (long k = mid_lo; k <= mid_hi; k++) {
        if (T[k - A] == 1) {
            parts = realloc(parts, (size_t)(np + 1) * sizeof(Expr*));
            parts[np++] = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)x), expr_new_integer(k) }, 2);
        }
    }
    if (has_right) {
        parts = realloc(parts, (size_t)(np + 1) * sizeof(Expr*));
        parts[np++] = expr_new_function(expr_new_symbol(SYM_GreaterEqual),
            (Expr*[]){ expr_copy((Expr*)x), expr_new_integer(rstart) }, 2);
    }
    free(T);

    Expr* out;
    if (np == 0)      { free(parts); out = expr_new_symbol(SYM_False); }
    else if (np == 1) { out = parts[0]; free(parts); }
    else              { out = expr_new_function(expr_new_symbol(SYM_Or), parts, (size_t)np); free(parts); }
    return eval_and_free(out);
}
