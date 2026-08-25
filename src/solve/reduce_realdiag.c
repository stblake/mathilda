/*
 * reduce_realdiag.c
 *
 * General univariate real sign diagram for `Reduce` (Phase 9).  See
 * reduce_realdiag.h for the method.  This is the fallback the reals dispatch in
 * reduce.c runs after the exact polynomial reduce_univar() declines on a
 * radical / rational-pole / bounded-domain-transcendental atom.
 *
 * The 1-D cell decomposition and its emission are shared with reduce_univar via
 * reduce_real_util.c; what is new here is (a) breakpoints from Solve roots +
 * domain boundaries rather than PolynomialQ roots, (b) a per-sample DOMAIN GATE
 * driven by reduce_realfn's head table, and (c) a numeric-sign fallback so a
 * transcendental breakpoint (a multiple of Pi) does not force a decline.
 */
#include "reduce_realdiag.h"
#include "reduce_realfn.h"
#include "reduce_real_util.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <math.h>

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}
static bool is_true_sym(const Expr* e)  { return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True; }
static bool is_false_sym(const Expr* e) { return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False; }

static bool rd_contains_x(const Expr* e, const Expr* x) {
    if (!e || !x || x->type != EXPR_SYMBOL) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == x->data.symbol.name;
    if (e->type == EXPR_FUNCTION) {
        if (rd_contains_x(e->data.function.head, x)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (rd_contains_x(e->data.function.args[i], x)) return true;
    }
    return false;
}

/* Heads that make an atom PIECEWISE in x: Floor/Ceiling/Round/IntegerPart/Mod/
 * Max/Min/Abs.  This engine collects breakpoints from Solve roots, poles and
 * head-table domain boundaries only -- it has no way to locate the transition
 * points of these heads, so an atom that still carries one (with x under it)
 * after reduce_realfn preprocessing declined to rewrite it cannot be sampled
 * soundly.  It must DECLINE (return NULL -> unevaluated) rather than emit a
 * single-cell verdict from an arbitrary sample.  Radicals / Log / inverse-trig
 * are deliberately NOT here: the domain table and Solve handle those. */
static bool expr_has_opaque_piecewise(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_Floor || h == SYM_Ceiling || h == SYM_Round
            || h == SYM_IntegerPart || h == SYM_Mod
            || h == SYM_Max || h == SYM_Min || h == SYM_Abs) {
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (rd_contains_x(e->data.function.args[i], x)) return true;
        }
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_has_opaque_piecewise(e->data.function.args[i], x)) return true;
    return false;
}

/* poly with the single rule x -> sample applied, evaluated (adopts nothing). */
static Expr* subst_eval(const Expr* poly, const Expr* x, const Expr* sample) {
    Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)sample) }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
        (Expr*[]){ expr_copy((Expr*)poly), rule }, 2));
}

/* Current-call options, set by reduce_univar_general (the single entry).  The
 * static helpers below run only inside that call and never re-enter, so file
 * scope is a safe home for the per-call config (the numeric-fallback tolerance
 * from WorkingPrecision, and the ReduceOpts forwarded onto soft root
 * isolation).  Default: exact-first tolerance 1e-9, Solve defaults. */
static const ReduceOpts* g_rd_opts = NULL;
static double            g_rd_tol  = 1e-9;

/* Sign of poly at the sample: -1/0/1, or -2 undecidable.  Exact (qqbar) first,
 * then a numeric N fallback for transcendental constants (multiples of Pi),
 * at the WorkingPrecision tolerance g_rd_tol. */
static int gen_sign_at(const Expr* poly, const Expr* x, const Expr* sample) {
    Expr* val = subst_eval(poly, x, sample);
    int s = rru_sign_of(val);
    if (s == -2) {
        bool ok; double d = rru_approx_double(val, &ok);
        if (ok && isfinite(d)) {
            double tol = g_rd_tol * (1.0 + fabs(d));
            s = (d > tol) ? 1 : (d < -tol ? -1 : -2);
        }
    }
    expr_free(val);
    return s;
}

/* Compare two constant breakpoints: -1/0/1, or -2 undecidable.  Exact then
 * numeric (near-coincident distinct values return -2 -> the engine declines). */
static int cmp_bp(const Expr* a, const Expr* b) {
    int s = rru_sign_compare(a, b);
    if (s != -2) return s;
    bool oka, okb;
    double da = rru_approx_double(a, &oka);
    double db = rru_approx_double(b, &okb);
    if (!oka || !okb || !isfinite(da) || !isfinite(db)) return -2;
    double diff = da - db;
    double tol = g_rd_tol * (1.0 + fabs(da) + fabs(db));
    if (diff > tol) return 1;
    if (diff < -tol) return -1;
    return -2;
}

/* ------------------------------------------------------------------ *
 *  Free-parameter guard                                               *
 * ------------------------------------------------------------------ */

static void collect_value_syms(const Expr* e, const char* xname,
                               const char*** arr, int* n, int* cap) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* nm = e->data.symbol.name;
        if (nm == xname) return;
        for (int i = 0; i < *n; i++) if ((*arr)[i] == nm) return;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *arr = realloc(*arr, (size_t)*cap * sizeof(char*)); }
        (*arr)[(*n)++] = nm;
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        if (e->data.function.head->type == EXPR_FUNCTION)
            collect_value_syms(e->data.function.head, xname, arr, n, cap);
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_value_syms(e->data.function.args[i], xname, arr, n, cap);
    }
}

/* True iff any leaf symbol other than x is not a numeric constant (a real free
 * parameter, e.g. `a`).  Pi/E/EulerGamma/... are NumericQ and allowed. */
static bool has_free_parameter(const RForm* F, const Expr* x) {
    const char* xname = (x->type == EXPR_SYMBOL) ? x->data.symbol.name : NULL;
    const char** syms = NULL; int ns = 0, cap = 0;
    for (int i = 0; i < F->n; i++)
        for (int k = 0; k < F->c[i]->n; k++) {
            collect_value_syms(F->c[i]->a[k].poly, xname, &syms, &ns, &cap);
            if (F->c[i]->a[k].denom) collect_value_syms(F->c[i]->a[k].denom, xname, &syms, &ns, &cap);
        }
    bool freeparam = false;
    for (int i = 0; i < ns && !freeparam; i++) {
        Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_NumericQ),
            (Expr*[]){ expr_new_symbol(syms[i]) }, 1));
        if (!is_true_sym(q)) freeparam = true;
        expr_free(q);
    }
    free(syms);
    return freeparam;
}

/* ------------------------------------------------------------------ *
 *  Breakpoint collection (soft Solve)                                 *
 * ------------------------------------------------------------------ */

/* Append the real roots of `poly` to (*arr,*n,*cap) IFF Solve returns a clean
 * finite list of numeric-real constants; otherwise append nothing (an identity
 * region or a parametric/conditional answer is tolerated, not fatal). */
static void soft_roots(const Expr* poly, const Expr* x, Expr*** arr, int* n, int* cap) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[3] = { eqn, expr_copy((Expr*)x), expr_new_symbol(SYM_Reals) };
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, 3, g_rd_opts));
    if (!is_head(sols, SYM_List)) { expr_free(sols); return; }

    /* Keep each row that pins x to a numeric-real constant; silently skip the
     * rest (an empty identity row {{}}, a ConditionalExpression, a parametric or
     * complex value).  Collecting the clean rows individually recovers isolated
     * roots (e.g. x == 1/2 from a squared factor) even when the same Solve also
     * returns an interval it cannot express as a rule. */
    size_t nr = sols->data.function.arg_count;
    for (size_t i = 0; i < nr; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) continue;
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
        Expr* val = rule->data.function.args[1];
        if (is_head(val, SYM_ConditionalExpression)) continue;
        bool ok; double d = rru_approx_double(val, &ok);
        if (!ok || !isfinite(d)) continue;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*)); }
        (*arr)[(*n)++] = expr_copy(val);
    }
    expr_free(sols);
}

/* Collect equation-root breakpoints of an atom expression, descending its
 * polynomial factor structure: a Times roots each factor, a Power[base,int]
 * roots its base, everything else is rooted whole (soft).  This recovers an
 * isolated root of a polynomial factor (e.g. x==1/2 from (2x-1)^2) even when
 * Solve on the whole product returns the identity-region {{}} it cannot express
 * as rules. */
static void collect_factor_roots(const Expr* poly, const Expr* x,
                                 Expr*** arr, int* n, int* cap) {
    if (is_head(poly, SYM_Times)) {
        for (size_t i = 0; i < poly->data.function.arg_count; i++)
            collect_factor_roots(poly->data.function.args[i], x, arr, n, cap);
        return;
    }
    if (is_head(poly, SYM_Power) && poly->data.function.arg_count == 2
        && poly->data.function.args[1]->type == EXPR_INTEGER
        && poly->data.function.args[1]->data.integer > 0) {
        collect_factor_roots(poly->data.function.args[0], x, arr, n, cap);
        return;
    }
    soft_roots(poly, x, arr, n, cap);
}

/* Insertion-sort the breakpoints with cmp_bp, dropping structural duplicates.
 * Returns false (declining) if any comparison is undecidable. */
static bool sort_dedup_bp(Expr** bp, int* m) {
    for (int i = 1; i < *m; i++) {
        Expr* key = bp[i];
        int j = i - 1;
        bool dup = false;
        while (j >= 0) {
            int c = cmp_bp(bp[j], key);
            if (c == -2) return false;
            if (c == 0) { dup = true; break; }
            if (c < 0) break;                 /* bp[j] < key: position found */
            bp[j + 1] = bp[j];
            j--;
        }
        if (dup) {                            /* drop key */
            expr_free(key);
            for (int t = j + 1; t < *m - 1; t++) bp[t] = bp[t + 1];
            (*m)--;
            i--;
        } else {
            bp[j + 1] = key;
        }
    }
    return true;
}

/* A rational strictly between lo and hi (NULL == unbounded), certified with
 * cmp_bp (numeric-tolerant, so transcendental bounds work).  Owned, or NULL. */
static Expr* mk_ratio(long num, long den) {
    if (den == 1) return expr_new_integer(num);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(num),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_integer(den), expr_new_integer(-1) }, 2) }, 2));
}
static Expr* gen_sample_between(const Expr* lo, const Expr* hi) {
    bool ok;
    if (!lo && !hi) return expr_new_integer(0);
    if (!lo) {
        double b = rru_approx_double(hi, &ok); if (!ok) return NULL;
        long s = (long)floor(b) - 1;
        for (int i = 0; i < 64; i++, s--) {
            Expr* c = expr_new_integer(s);
            int cmp = cmp_bp(c, hi);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp < 0) return c;
            expr_free(c);
        }
        return NULL;
    }
    if (!hi) {
        double a = rru_approx_double(lo, &ok); if (!ok) return NULL;
        long s = (long)ceil(a) + 1;
        for (int i = 0; i < 64; i++, s++) {
            Expr* c = expr_new_integer(s);
            int cmp = cmp_bp(c, lo);
            if (cmp == -2) { expr_free(c); return NULL; }
            if (cmp > 0) return c;
            expr_free(c);
        }
        return NULL;
    }
    double a = rru_approx_double(lo, &ok); if (!ok) return NULL;
    double b = rru_approx_double(hi, &ok); if (!ok) return NULL;
    double mid = 0.5 * (a + b);
    for (long D = 1; D <= (1L << 30); D *= 2) {
        long num = (long)llround(mid * (double)D);
        Expr* c = mk_ratio(num, D);
        int c1 = cmp_bp(c, lo), c2 = cmp_bp(c, hi);
        if (c1 == -2 || c2 == -2) { expr_free(c); return NULL; }
        if (c1 > 0 && c2 < 0) return c;
        expr_free(c);
    }
    return NULL;
}

/* ------------------------------------------------------------------ *
 *  Per-sample truth (domain gate + pole gate + decide)                *
 * ------------------------------------------------------------------ */

/* Truth of one atom at the sample: 1 / 0 / -1 (undecided).  Assumes the domain
 * gate already passed. */
static int atom_truth_general(const RAtom* a, const Expr* x, const Expr* sample) {
    if (a->rel == R_ELEM) return -1;

    int sq = 1;
    if (a->denom) {
        sq = gen_sign_at(a->denom, x, sample);
        if (sq == -2) return -1;
        if (sq == 0) return 0;                 /* pole -> atom false */
    }

    if (a->rel == R_EQ || a->rel == R_NE) {
        const char* head = (a->rel == R_EQ) ? SYM_Equal : SYM_Unequal;
        Expr* rel = expr_new_function(expr_new_symbol(head),
            (Expr*[]){ expr_copy(a->poly), expr_new_integer(0) }, 2);
        Expr* rule = expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_copy((Expr*)x), expr_copy((Expr*)sample) }, 2);
        Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_ReplaceAll),
            (Expr*[]){ rel, rule }, 2));
        int r = is_true_sym(v) ? 1 : (is_false_sym(v) ? 0 : -1);
        expr_free(v);
        return r;
    }
    /* R_LT / R_LE : sign of poly * sign of denom. */
    int sp = gen_sign_at(a->poly, x, sample);
    if (sp == -2) return -1;
    int s = sp * sq;
    if (a->rel == R_LT) return s < 0;
    return s <= 0;   /* R_LE */
}

/* Truth of the whole DNF at the sample: 1 / 0 / -1.  Each conjunct's real-domain
 * constraints (`ccons[i]`, from that conjunct's radical/Log/inverse-trig atoms) gate
 * ONLY that conjunct -- a domain condition born in one Abs sign-branch must never
 * exclude the mutually-exclusive other branch (the bug that made `Sqrt[Abs[x]]<1`
 * collapse `x>=0 && x<=0` to `x==0`).  A conjunct is true only where its own domain
 * holds AND its atoms hold; the DNF is true iff some conjunct is. */
static int form_truth_general(const RForm* F, const Expr* x, const Expr* sample,
                              RDomCon* const* ccons, const int* cncon) {
    int any_undec = 0;
    for (int i = 0; i < F->n; i++) {
        RConj* c = F->c[i];
        int conj = 1;                          /* 1 true / 0 false / -1 undetermined */
        /* Domain gate for THIS conjunct.  A decidably-failed constraint excludes
         * the point soundly, so keep scanning past an undecidable one (a later
         * constraint may still fail decidably -- e.g. the inner x-1>=0 of a nested
         * radical fails where the outer radical's sign is non-real/undecidable). */
        bool dom_undec = false;
        for (int t = 0; t < cncon[i]; t++) {
            int s = gen_sign_at(ccons[i][t].poly, x, sample);
            if (s == -2) { dom_undec = true; continue; }
            bool okdom = ccons[i][t].strict ? (s > 0) : (s >= 0);
            if (!okdom) { conj = 0; break; }   /* out of this conjunct's domain */
        }
        if (conj == 1 && dom_undec) conj = -1; /* in-domain where decidable, else undetermined */
        if (conj == 1)
            for (int k = 0; k < c->n; k++) {
                int t = atom_truth_general(&c->a[k], x, sample);
                if (t == -1) { conj = -1; break; }
                if (t == 0)  { conj = 0;  break; }
            }
        if (conj == 1)  return 1;              /* a fully-true conjunct decides the DNF */
        if (conj == -1) any_undec = 1;
    }
    return any_undec ? -1 : 0;                 /* else: all false, or decline if any undetermined */
}

/* Free the per-conjunct domain-constraint arrays (polys + the row arrays). */
static void free_per_conj_cons(RDomCon** ccons, const int* cncon, int n) {
    if (!ccons) return;
    for (int i = 0; i < n; i++) {
        for (int t = 0; t < cncon[i]; t++) expr_free(ccons[i][t].poly);
        free(ccons[i]);
    }
    free(ccons);
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* reduce_univar_general(const RForm* F, const Expr* x, Expr** vars, int nv,
                            const ReduceOpts* opts) {
    (void)vars; (void)nv;
    /* Publish the per-call config for the static helpers below: WorkingPrecision
     * -> numeric-fallback tolerance (Infinity keeps the exact-first 1e-9), and
     * the ReduceOpts forwarded onto soft root isolation (Cubics / Quartics). */
    g_rd_opts = opts;
    double wp_digits;
    g_rd_tol = reduce_opts_wp_digits(opts, &wp_digits) ? pow(10.0, -wp_digits) : 1e-9;

    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);
    if (has_free_parameter(F, x)) return NULL;

    /* Soundness gate: an atom still carrying a piecewise head (Floor/Max/Abs/...)
     * in x -- one that preprocessing could not rewrite away -- has transition
     * points this sampler cannot see, so decline rather than risk a single-cell
     * verdict. */
    for (int i = 0; i < F->n; i++)
        for (int k = 0; k < F->c[i]->n; k++) {
            const RAtom* a = &F->c[i]->a[k];
            if (expr_has_opaque_piecewise(a->poly, x)) return NULL;
            if (a->denom && expr_has_opaque_piecewise(a->denom, x)) return NULL;
        }

    /* Domain constraints, collected PER CONJUNCT (ccons[i] holds the real-domain
     * constraints of conjunct i's atoms), so the gate in form_truth_general scopes
     * each to its own conjunct instead of ANDing them all globally. */
    size_t nc = (F->n > 0) ? (size_t)F->n : 1;   /* bounded count keeps -Walloc-size quiet */
    RDomCon** ccons = calloc(nc, sizeof(RDomCon*));
    int* cncon = calloc(nc, sizeof(int));
    int* cccap = calloc(nc, sizeof(int));
    for (int i = 0; i < F->n; i++)
        for (int k = 0; k < F->c[i]->n; k++) {
            const RAtom* a = &F->c[i]->a[k];
            if (a->rel == R_ELEM) {
                free_per_conj_cons(ccons, cncon, F->n); free(cncon); free(cccap);
                return NULL;
            }
            reduce_real_domain_collect(a->poly, x, &ccons[i], &cncon[i], &cccap[i]);
            if (a->denom) reduce_real_domain_collect(a->denom, x, &ccons[i], &cncon[i], &cccap[i]);
        }

    /* Breakpoints: equation roots, poles, and domain boundaries.  The domain
     * boundaries stay a UNION across all conjuncts -- extra breakpoints only refine
     * the cell decomposition and are always sound. */
    Expr** bp = NULL; int m = 0, cap = 0;
    for (int i = 0; i < F->n; i++)
        for (int k = 0; k < F->c[i]->n; k++) {
            const RAtom* a = &F->c[i]->a[k];
            collect_factor_roots(a->poly, x, &bp, &m, &cap);
            if (a->denom) soft_roots(a->denom, x, &bp, &m, &cap);
        }
    for (int i = 0; i < F->n; i++)
        for (int t = 0; t < cncon[i]; t++) soft_roots(ccons[i][t].poly, x, &bp, &m, &cap);

    bool ok = sort_dedup_bp(bp, &m);
    if (!ok) goto decline;

    {
        int ncells = 2 * m + 1;
        int* truth = malloc((size_t)ncells * sizeof(int));
        for (int idx = 0; idx < ncells && ok; idx++) {
            Expr* sample;
            if (idx % 2 == 1) {                       /* point cell */
                sample = expr_copy(bp[(idx + 1) / 2 - 1]);
            } else {                                  /* interval cell */
                int j = idx / 2;
                const Expr* lo = (j == 0) ? NULL : bp[j - 1];
                const Expr* hi = (j == m) ? NULL : bp[j];
                sample = gen_sample_between(lo, hi);
                if (!sample) { ok = false; break; }
            }
            int t = form_truth_general(F, x, sample, ccons, cncon);
            expr_free(sample);
            if (t == -1) { ok = false; break; }
            truth[idx] = t;
        }
        if (!ok) { free(truth); goto decline; }

        Expr* result = rru_emit_sign_diagram(bp, m, truth, x);
        free(truth);
        for (int i = 0; i < m; i++) expr_free(bp[i]);
        free(bp);
        free_per_conj_cons(ccons, cncon, F->n); free(cncon); free(cccap);
        return result;
    }

decline:
    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp);
    free_per_conj_cons(ccons, cncon, F->n); free(cncon); free(cccap);
    return NULL;
}
