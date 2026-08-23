/*
 * reduce_univar.c
 *
 * Univariate real sign-diagram engine for `Reduce` (REDUCE_PLAN.md, Phase 2).
 * See reduce_univar.h for the algorithm.  Everything is done by reusing existing
 * builtins (Solve over Reals for root isolation, PolynomialQ / ReplaceAll / N)
 * plus flint_qqbar_compare as the exact real-algebraic sign oracle; rational
 * samples are signed exactly without FLINT.
 */
#include "reduce_univar.h"

#include "eval.h"
#include "sym_names.h"
#include "flint_qqbar.h"

#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ *
 *  Exact real sign oracle                                             *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Sign of a numeric literal; *is_num=false if `e` is not a plain number. */
static int number_sign(const Expr* e, bool* is_num) {
    *is_num = true;
    switch (e->type) {
        case EXPR_INTEGER: return (e->data.integer > 0) - (e->data.integer < 0);
        case EXPR_REAL:    return (e->data.real > 0.0) - (e->data.real < 0.0);
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint);
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_sgn(e->data.mpfr);
#endif
        case EXPR_FUNCTION:
            if (is_head(e, SYM_Rational) && e->data.function.arg_count == 2) {
                const Expr* n = e->data.function.args[0];
                if (n->type == EXPR_INTEGER) return (n->data.integer > 0) - (n->data.integer < 0);
                if (n->type == EXPR_BIGINT)  return mpz_sgn(n->data.bigint);
            }
            *is_num = false; return 0;
        default: *is_num = false; return 0;
    }
}

/* Sign of a constant real value: -1, 0, 1, or -2 (undecidable / FLINT off).
 * Plain numbers are signed exactly; anything else goes to the qqbar oracle. */
static int sign_of(const Expr* e) {
    bool isn; int s = number_sign(e, &isn);
    if (isn) return s;
    Expr* zero = expr_new_integer(0);
    int r = flint_qqbar_compare(e, zero);   /* sign(e - 0) */
    expr_free(zero);
    return r;
}

/* Sign of (a - b): -1, 0, 1, or -2. */
static int sign_compare(const Expr* a, const Expr* b) {
    Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Subtract),
        (Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
    int s = sign_of(diff);
    expr_free(diff);
    return s;
}

/* ------------------------------------------------------------------ *
 *  Sample-point selection                                             *
 * ------------------------------------------------------------------ */

static double approx_double(const Expr* e, bool* ok) {
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy((Expr*)e) }, 1));
    double d = 0.0; *ok = true;
    if (v->type == EXPR_REAL)         d = v->data.real;
    else if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
    else if (v->type == EXPR_BIGINT)  d = mpz_get_d(v->data.bigint);
#ifdef USE_MPFR
    else if (v->type == EXPR_MPFR)    d = mpfr_get_d(v->data.mpfr, MPFR_RNDN);
#endif
    else *ok = false;
    expr_free(v);
    return d;
}

/* num/den as a canonical Integer or Rational. */
static Expr* make_rational(long num, long den) {
    if (den == 1) return expr_new_integer(num);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(num),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_integer(den), expr_new_integer(-1) }, 2) }, 2));
}

/* A rational strictly between `lo` and `hi` (either NULL = unbounded), certified
 * exactly with sign_compare.  Returns owned Expr, or NULL if none is found or a
 * comparison is undecidable. */
static Expr* rational_between(const Expr* lo, const Expr* hi) {
    bool ok;
    if (!lo && !hi) return expr_new_integer(0);

    if (!lo) {                              /* (-inf, hi) */
        double b = approx_double(hi, &ok);
        if (!ok) return NULL;
        long s = (long)floor(b) - 1;
        for (int i = 0; i < 64; i++, s--) {
            Expr* c = expr_new_integer(s);
            int cmp = sign_compare(c, hi);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp < 0) return c;
            expr_free(c);
        }
        return NULL;
    }
    if (!hi) {                              /* (lo, +inf) */
        double a = approx_double(lo, &ok);
        if (!ok) return NULL;
        long s = (long)ceil(a) + 1;
        for (int i = 0; i < 64; i++, s++) {
            Expr* c = expr_new_integer(s);
            int cmp = sign_compare(c, lo);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp > 0) return c;
            expr_free(c);
        }
        return NULL;
    }

    double a = approx_double(lo, &ok); if (!ok) return NULL;
    double b = approx_double(hi, &ok); if (!ok) return NULL;
    double mid = 0.5 * (a + b);
    for (long D = 1; D <= (1L << 30); D *= 2) {
        long num = (long)llround(mid * (double)D);
        Expr* c = make_rational(num, D);
        int c1 = sign_compare(c, lo);       /* want > 0 */
        int c2 = sign_compare(c, hi);       /* want < 0 */
        if (c1 == -2 || c2 == -2) { expr_free(c); return NULL; }
        if (c1 > 0 && c2 < 0) return c;
        expr_free(c);
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Truth of the formula at a sample point                             *
 * ------------------------------------------------------------------ */

/* Sign of poly(sample): -1, 0, 1, or -2. */
static int poly_sign_at(const Expr* poly, const Expr* x, const Expr* sample) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)sample) }, 2);
    Expr* val = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rule }, 2));
    int s = sign_of(val);
    expr_free(val);
    return s;
}

/* Truth of one atom at the sample: 1 true, 0 false, -1 undecided. */
static int atom_truth_at(const RAtom* a, const Expr* x, const Expr* sample) {
    if (a->rel == R_ELEM) return -1;
    int s = poly_sign_at(a->poly, x, sample);
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
 *  Root collection                                                    *
 * ------------------------------------------------------------------ */

static bool poly_is_polynomial(const Expr* poly, const Expr* x) {
    Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_PolynomialQ),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)x) }, 2));
    bool r = (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True);
    expr_free(q);
    return r;
}

/* Append the real roots of `poly` (via Solve over Reals) to *arr.  Returns
 * false to bail on an unexpected Solve shape. */
static bool collect_roots(const Expr* poly, const Expr* x,
                          Expr*** arr, int* n, int* cap) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Solve),
        (Expr*[]){ eqn, expr_copy((Expr*)x), expr_new_symbol(SYM_Reals) }, 3);
    Expr* sols = eval_and_free(call);

    bool ok = is_head(sols, SYM_List);
    for (size_t i = 0; ok && i < sols->data.function.arg_count; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) { ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) { ok = false; break; }
        Expr* val = rule->data.function.args[1];
        if (is_head(val, SYM_ConditionalExpression)) { ok = false; break; }
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*)); }
        (*arr)[(*n)++] = expr_copy(val);
    }
    expr_free(sols);
    return ok;
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
            if (!poly_is_polynomial(a->poly, x)) { *ok = false; break; }
            if (np == pcap) { pcap = pcap ? pcap * 2 : 8; polys = realloc(polys, (size_t)pcap * sizeof(Expr*)); }
            polys[np++] = a->poly;          /* borrowed */
        }
    }
    Expr** roots = NULL; int nr = 0, rcap = 0;
    for (int j = 0; *ok && j < np; j++)
        if (!collect_roots(polys[j], x, &roots, &nr, &rcap)) *ok = false;
    free(polys);
    if (!*ok) { for (int i = 0; i < nr; i++) expr_free(roots[i]); free(roots); return NULL; }
    for (int i = 0; *ok && i < nr; i++) {
        int mn = i;
        for (int j = i + 1; j < nr; j++) {
            int c = sign_compare(roots[j], roots[mn]);
            if (c == -2) { *ok = false; break; }
            if (c < 0) mn = j;
        }
        if (*ok && mn != i) { Expr* t = roots[i]; roots[i] = roots[mn]; roots[mn] = t; }
    }
    int m = 0;
    Expr** bp = malloc((size_t)(nr + 1) * sizeof(Expr*));
    for (int i = 0; *ok && i < nr; i++) {
        if (m == 0) { bp[m++] = roots[i]; continue; }
        int c = sign_compare(bp[m - 1], roots[i]);
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
            sample = rational_between(lo, hi);
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
    double lo = approx_double(bp[0], &ap);
    double hi = approx_double(bp[m - 1], &bp_ok);
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
