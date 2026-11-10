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
#include "reduce_eq.h"
#include "zero_test.h"
#include "root_trans.h"

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

/* Residuals whose real roots Solve could not find (borrowed from the current
 * call's RForm; per-call scratch, reset by reduce_univar_general). */
static const Expr** g_rd_unrooted    = NULL;
static int          g_rd_nunrooted   = 0;
static int          g_rd_capunrooted = 0;

/* Certified roots of those residuals (isolate_unrooted): the breakpoint value
 * (owned copy) and the residual it is a root of (borrowed).  At such a sample
 * any polynomial carrying the residual as a factor is EXACTLY zero -- the one
 * fact a numeric sign evaluation at a Root[{f, x0}] can never establish. */
static Expr**       g_rd_regval = NULL;
static const Expr** g_rd_regres = NULL;
static int          g_rd_nreg = 0, g_rd_capreg = 0;

/* Kleene mask over g_rd_unrooted for the current sample: a residual whose
 * roots were NOT isolated in the sample's region is "unknown" there, and so
 * is every sign / atom that contains it.  NULL == all known. */
static bool*        g_rd_unknown = NULL;

/* r is poly itself or a Times / positive-Power factor of it (pointer). */
static bool rd_has_factor(const Expr* poly, const Expr* r) {
    if (poly == r) return true;
    if (is_head(poly, SYM_Times)) {
        for (size_t i = 0; i < poly->data.function.arg_count; i++)
            if (rd_has_factor(poly->data.function.args[i], r)) return true;
        return false;
    }
    if (is_head(poly, SYM_Power) && poly->data.function.arg_count == 2
        && poly->data.function.args[1]->type == EXPR_INTEGER
        && poly->data.function.args[1]->data.integer > 0)
        return rd_has_factor(poly->data.function.args[0], r);
    return false;
}

/* r occurs (by pointer) anywhere in e. */
static bool rd_contains_ptr(const Expr* e, const Expr* r) {
    if (e == r) return true;
    if (!e || e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rd_contains_ptr(e->data.function.args[i], r)) return true;
    return false;
}

/* poly carries a residual that is unknown at the current sample. */
static bool rd_poly_unknown(const Expr* poly) {
    if (!g_rd_unknown) return false;
    for (int u = 0; u < g_rd_nunrooted; u++)
        if (g_rd_unknown[u] && rd_contains_ptr(poly, g_rd_unrooted[u])) return true;
    return false;
}

/* The sample is a certified root of a residual factor of poly. */
static bool rd_registered_zero(const Expr* poly, const Expr* sample) {
    for (int i = 0; i < g_rd_nreg; i++)
        if (rd_has_factor(poly, g_rd_regres[i]) && expr_eq(g_rd_regval[i], sample))
            return true;
    return false;
}

/* Sign of poly at the sample: -1/0/1, or -2 undecidable.  Exact (qqbar) first,
 * then a numeric N fallback for transcendental constants (multiples of Pi),
 * at the WorkingPrecision tolerance g_rd_tol. */
static int gen_sign_at(const Expr* poly, const Expr* x, const Expr* sample) {
    if (rd_poly_unknown(poly)) return -2;
    if (rd_registered_zero(poly, sample)) return 0;
    Expr* val = subst_eval(poly, x, sample);
    int s = rru_sign_of(val);
    if (s == -2) {
        bool ok; double d = rru_approx_double(val, &ok);
        if (ok && isfinite(d)) {
            double tol = g_rd_tol * (1.0 + fabs(d));
            s = (d > tol) ? 1 : (d < -tol ? -1 : -2);
        }
    }
    /* Within tolerance of 0: the sample is typically an exact root of poly
     * (a breakpoint Solve returned, -ProductLog[-1/3] for E^x - 3 x); the
     * exact zero test confirms it. */
    if (s == -2 && zero_test_decide(val) == ZERO_TEST_TRUE) s = 0;
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

/* Class of N[poly /. x -> s]: 0 numerically zero, 1 non-real, 2 real nonzero,
 * -1 not a finite number. */
static int unrooted_class_at(const Expr* poly, const Expr* x, double s) {
    Expr* v = subst_eval(poly, x, expr_new_real(s));
    Expr* n = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ v }, 1));
    int cls = -1;
    double re = 0.0, im = 0.0;
    bool ok = true;
    if (n->type == EXPR_REAL) re = n->data.real;
    else if (n->type == EXPR_INTEGER) re = (double)n->data.integer;
    else if (is_head(n, SYM_Complex) && n->data.function.arg_count == 2) {
        bool o1, o2;
        re = rru_approx_double(n->data.function.args[0], &o1);
        im = rru_approx_double(n->data.function.args[1], &o2);
        ok = o1 && o2;
    } else ok = false;
    expr_free(n);
    if (ok && isfinite(re) && isfinite(im)) {
        double tol = 1e-9 * (1.0 + fabs(s));
        if (fabs(im) > tol)      cls = 1;
        else if (fabs(re) > tol) cls = 2;
        else                     cls = 0;
    }
    return cls;
}

/* Unrooted residual u must, on each open cell between consecutive
 * breakpoints, be numerically zero at all probes or non-real at all probes --
 * a cell where it is real and nonzero could hide a root (Sin[x] - x/2), and
 * then its roots must be isolated (isolate_unrooted) or the diagram declines. */
/* Probe class of unrooted residual u on interval cell j (0..m) of the
 * breakpoints: 0 numerically zero at all five probes, 1 non-real at all five,
 * -1 otherwise (real and nonzero somewhere, mixed, or not evaluable). */
static int unrooted_cell_class(const Expr* x, Expr** bp, int m, int j, int u) {
    bool ok_lo = true, ok_hi = true;
    double lo = (j == 0) ? 0.0 : rru_approx_double(bp[j - 1], &ok_lo);
    double hi = (j == m) ? 0.0 : rru_approx_double(bp[j], &ok_hi);
    if (!ok_lo || !ok_hi) return -1;
    double probes[5];
    if (j == 0 && j == m) {
        double p[5] = { -10.5, -1.25, 0.375, 1.75, 11.25 };
        for (int t = 0; t < 5; t++) probes[t] = p[t];
    } else if (j == 0) {
        double p[5] = { 0.125, 0.5, 1.5, 7.0, 40.0 };
        for (int t = 0; t < 5; t++) probes[t] = hi - p[t];
    } else if (j == m) {
        double p[5] = { 0.125, 0.5, 1.5, 7.0, 40.0 };
        for (int t = 0; t < 5; t++) probes[t] = lo + p[t];
    } else {
        for (int t = 0; t < 5; t++) probes[t] = lo + (hi - lo) * (t + 1) / 6.0;
    }
    int cls0 = -2;
    for (int t = 0; t < 5; t++) {
        int c = unrooted_class_at(g_rd_unrooted[u], x, probes[t]);
        if (c < 0 || c == 2) return -1;
        if (cls0 == -2) cls0 = c;
        else if (c != cls0) return -1;
    }
    return cls0;
}

/* Unrooted residual u must, on each open cell between consecutive
 * breakpoints, be numerically zero at all probes or non-real at all probes --
 * a cell where it is real and nonzero could hide a root (Sin[x] - x/2), and
 * then its roots must be isolated (isolate_unrooted) or the diagram declines. */
static bool unrooted_cells_ok(const Expr* x, Expr** bp, int m, int u) {
    for (int j = 0; j <= m; j++)
        if (unrooted_cell_class(x, bp, m, j, u) < 0) return false;
    return true;
}

/* Record that residual r vanishes at `sample` (owned copy taken). */
static void rd_register_zero(const Expr* sample, const Expr* r) {
    if (g_rd_nreg == g_rd_capreg) {
        g_rd_capreg = g_rd_capreg ? 2 * g_rd_capreg : 8;
        g_rd_regval = realloc(g_rd_regval, (size_t)g_rd_capreg * sizeof(Expr*));
        g_rd_regres = (const Expr**)realloc((void*)g_rd_regres,
                            (size_t)g_rd_capreg * sizeof(Expr*));
    }
    g_rd_regval[g_rd_nreg] = expr_copy((Expr*)sample);
    g_rd_regres[g_rd_nreg] = r;
    g_rd_nreg++;
}

/* True iff `e` carries an Element[...] (the C[k] in Integers of a family). */
static bool contains_element(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_head(e, SYM_Element)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_element(e->data.function.args[i])) return true;
    return false;
}

/* Append the real roots of `poly` to (*arr,*n,*cap).  Returns false -- the
 * caller must DECLINE -- when the real roots are not all known: Solve left the
 * equation unevaluated (Sin[x] - x/2 == 0), or a periodic family has infinitely
 * many real members (Sin[x] == 0) that no finite breakpoint list can carry.
 * Silently skipping either was a wrong answer, not a decline: with no
 * breakpoints the sign diagram is ONE cell decided at the single sample x = 0,
 * which made Reduce[Sin[x] > x/2, x, Reals], Reduce[x + Sin[x] > 0, x, Reals]
 * and Reduce[Sin[x] > 0, x, Reals] all return False.
 *
 * The rows are first collapsed to their real members exactly as the equation
 * path does (reduce_collapse_reals_families), so an exponential's imaginary-
 * period family contributes its one real member (Exp[x] == 2 -> Log[2]).  An
 * empty identity row {{}} is tolerated (it adds no breakpoint), and a
 * non-family ConditionalExpression contributes its value: an extra breakpoint
 * only refines the cells, so it is always sound. */
static bool soft_roots(const Expr* poly, const Expr* x, Expr*** arr, int* n, int* cap) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[3] = { eqn, expr_copy((Expr*)x), expr_new_symbol(SYM_Reals) };
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, 3, g_rd_opts));
    if (!is_head(sols, SYM_List)) {
        expr_free(sols);
        /* Solve cannot express an identity (ArcSin[x] + ArcCos[x] - Pi/2): a
         * residual that is identically zero has no isolated roots to miss. */
        if (zero_test_decide(poly) == ZERO_TEST_TRUE) return true;
        /* Roots unknown.  Record the residual: the driver accepts it only if,
         * on every cell of the known breakpoints, it is identically zero or
         * nowhere real (Log[x^2] - 2 Log[-x]: zero for x < 0, non-real for
         * x > 0), where it cannot change sign; else it declines. */
        if (g_rd_nunrooted == g_rd_capunrooted) {
            g_rd_capunrooted = g_rd_capunrooted ? 2 * g_rd_capunrooted : 4;
            g_rd_unrooted = (const Expr**)realloc((void*)g_rd_unrooted,
                                    (size_t)g_rd_capunrooted * sizeof(Expr*));
        }
        g_rd_unrooted[g_rd_nunrooted++] = poly;           /* borrowed */
        return true;
    }
    sols = reduce_collapse_reals_families(sols);
    if (!sols) return false;

    /* Keep each row that pins x to a numeric-real constant.  Collecting the
     * clean rows individually recovers isolated roots (e.g. x == 1/2 from a
     * squared factor) even when the same Solve also returns an interval it
     * cannot express as a rule. */
    bool ok_all = true;
    size_t nr = sols->data.function.arg_count;
    for (size_t i = 0; i < nr && ok_all; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) continue;
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) continue;
        Expr* val = rule->data.function.args[1];
        if (is_head(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            const Expr* cond = val->data.function.args[1];
            if (contains_element(cond)) {
                ok_all = false;              /* a real family: infinitely many roots */
                break;
            }
            val = val->data.function.args[0];
        }
        bool ok; double d = rru_approx_double(val, &ok);
        if (!ok || !isfinite(d)) continue;
        if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*)); }
        (*arr)[(*n)++] = expr_copy(val);
    }
    expr_free(sols);
    return ok_all;
}

/* Collect equation-root breakpoints of an atom expression, descending its
 * polynomial factor structure: a Times roots each factor, a Power[base,int]
 * roots its base, everything else is rooted whole (soft).  This recovers an
 * isolated root of a polynomial factor (e.g. x==1/2 from (2x-1)^2) even when
 * Solve on the whole product returns the identity-region {{}} it cannot express
 * as rules. */
static bool collect_factor_roots(const Expr* poly, const Expr* x,
                                 Expr*** arr, int* n, int* cap) {
    if (!rd_contains_x(poly, x)) return true;          /* constant factor: no roots */
    if (is_head(poly, SYM_Times)) {
        for (size_t i = 0; i < poly->data.function.arg_count; i++)
            if (!collect_factor_roots(poly->data.function.args[i], x, arr, n, cap))
                return false;
        return true;
    }
    if (is_head(poly, SYM_Power) && poly->data.function.arg_count == 2
        && poly->data.function.args[1]->type == EXPR_INTEGER
        && poly->data.function.args[1]->data.integer > 0) {
        return collect_factor_roots(poly->data.function.args[0], x, arr, n, cap);
    }
    return soft_roots(poly, x, arr, n, cap);
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
 *  Certified isolation of unrooted residuals                          *
 * ------------------------------------------------------------------ */

/* An exact rational isolation endpoint for breakpoint e: e itself when it is
 * already an Integer / Rational, else a dyadic rational just outside it
 * (dir -1: below, +1: above).  Widening only enlarges the searched span, and
 * every root found in the margin is still a genuine root -- an extra
 * breakpoint, which only refines the cells. */
static Expr* rd_outer_rational(const Expr* e, int dir) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || is_head(e, SYM_Rational))
        return expr_copy((Expr*)e);
    bool ok; double d = rru_approx_double(e, &ok);
    if (!ok || !isfinite(d) || fabs(d) > 1e12) return NULL;
    double scale = 1048576.0;                         /* 2^20 */
    double t = d * scale;
    double k = (dir < 0) ? floor(t) - 1.0 : ceil(t) + 1.0;
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer((int64_t)k),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_new_integer(1048576), expr_new_integer(-1) }, 2) }, 2));
}

/* Isolate, with certified arb arithmetic (root_trans_isolate), the real roots
 * of every unrooted residual that unrooted_cells_ok rejected, and append them
 * to the breakpoints.  The searched span is
 *   [-R, R]          when root_trans_dominance_bound proves f != 0 for |x| >= R
 *                    (Sin[x] - x/2: R = 3; x^2 - Cos[x]: R = 2);
 *   [bp_0, bp_{m-1}] otherwise -- the bounded cells between the known
 *                    breakpoints (x Sin[x] < 1 && -5 < x < 5), and the residual
 *                    is marked *unbounded[u]: on the two unbounded end cells its
 *                    roots stay unknown and the cell must be decided without it
 *                    (Kleene evaluation, see g_rd_unknown) or the call declines
 *                    (Sin[x] > 0 has infinitely many roots).
 * Each root is registered as a certified zero of its residual.  Returns false
 * (decline) when some residual's roots cannot be certified. */
static bool isolate_unrooted(const Expr* x, Expr*** bp, int* m, int* cap,
                             const bool* need, bool* unbounded) {
    int m0 = *m;
    Expr** nb = NULL; int nn = 0, ncap = 0;
    bool ok = true;
    for (int u = 0; ok && u < g_rd_nunrooted; u++) {
        unbounded[u] = false;
        if (!need[u]) continue;
        const Expr* f = g_rd_unrooted[u];
        Expr* R = root_trans_dominance_bound(f, x);
        Expr *lo = NULL, *hi = NULL;
        if (R) {
            if (!(R->type == EXPR_INTEGER && R->data.integer == 0)) {
                hi = R;
                lo = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(R) }, 2));
            } else expr_free(R);
        } else {
            unbounded[u] = true;
            if (m0 > 0) {
                lo = rd_outer_rational((*bp)[0], -1);
                hi = rd_outer_rational((*bp)[m0 - 1], +1);
                if (!lo || !hi) ok = false;
            }
        }
        if (ok && lo && hi) {
            int n0 = nn;
            ok = root_trans_isolate(f, x, lo, hi, &nb, &nn, &ncap);
            for (int i = n0; ok && i < nn; i++) {
                if (g_rd_nreg == g_rd_capreg) {
                    g_rd_capreg = g_rd_capreg ? 2 * g_rd_capreg : 8;
                    g_rd_regval = realloc(g_rd_regval, (size_t)g_rd_capreg * sizeof(Expr*));
                    g_rd_regres = (const Expr**)realloc((void*)g_rd_regres,
                                        (size_t)g_rd_capreg * sizeof(Expr*));
                }
                g_rd_regval[g_rd_nreg] = expr_copy(nb[i]);
                g_rd_regres[g_rd_nreg] = f;
                g_rd_nreg++;
            }
        }
        if (lo) expr_free(lo);
        if (hi) expr_free(hi);
    }
    for (int i = 0; i < nn; i++) {
        if (ok) {
            if (*m == *cap) { *cap = *cap ? 2 * *cap : 8; *bp = realloc(*bp, (size_t)*cap * sizeof(Expr*)); }
            (*bp)[(*m)++] = nb[i];
        } else expr_free(nb[i]);
    }
    free(nb);
    return ok;
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
        /* At a transcendental Root breakpoint (or with an unknown residual)
         * the truth is the certified sign, never a numeric Equal. */
        if (rd_poly_unknown(a->poly)) return -1;
        if (rd_registered_zero(a->poly, sample) || root_trans_is_object(sample)) {
            int s = gen_sign_at(a->poly, x, sample);
            if (s == -2) return -1;
            return (a->rel == R_EQ) ? (s == 0) : (s != 0);
        }
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
        /* Kleene AND: an undetermined atom does not stop the scan -- a later
         * decidably-false atom (x < 5 outside the region of x Sin[x] < 1)
         * still makes the conjunct false. */
        if (conj != 0)
            for (int k = 0; k < c->n; k++) {
                int t = atom_truth_general(&c->a[k], x, sample);
                if (t == -1) { conj = -1; continue; }
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
    bool ok = true;
    /* Per-call scratch (see soft_roots / isolate_unrooted).  Solve may re-enter
     * Reduce and so this function: save the caller's state, restore on exit. */
    const Expr** saved_unrooted = g_rd_unrooted;
    int saved_nunrooted = g_rd_nunrooted, saved_capunrooted = g_rd_capunrooted;
    Expr** saved_regval = g_rd_regval; const Expr** saved_regres = g_rd_regres;
    int saved_nreg = g_rd_nreg, saved_capreg = g_rd_capreg;
    bool* saved_unknown = g_rd_unknown;
    g_rd_unrooted = NULL; g_rd_nunrooted = 0; g_rd_capunrooted = 0;
    g_rd_regval = NULL; g_rd_regres = NULL; g_rd_nreg = 0; g_rd_capreg = 0;
    g_rd_unknown = NULL;
    bool* unbounded = NULL;               /* per residual: roots unknown off-span */
    bool* fallback = NULL;                /* per residual: decided cell by cell   */
    Expr *span_lo = NULL, *span_hi = NULL;/* known-breakpoint span (copies)      */
    Expr* result = NULL;
    for (int i = 0; ok && i < F->n; i++)
        for (int k = 0; ok && k < F->c[i]->n; k++) {
            const RAtom* a = &F->c[i]->a[k];
            ok = collect_factor_roots(a->poly, x, &bp, &m, &cap);
            if (ok && a->denom) ok = collect_factor_roots(a->denom, x, &bp, &m, &cap);
        }
    for (int i = 0; ok && i < F->n; i++)
        for (int t = 0; ok && t < cncon[i]; t++)
            ok = collect_factor_roots(ccons[i][t].poly, x, &bp, &m, &cap);
    if (!ok) goto out;

    ok = sort_dedup_bp(bp, &m);
    if (!ok) goto out;
    if (g_rd_nunrooted > 0) {
        /* A residual whose roots Solve could not find either cannot change
         * sign on any cell (unrooted_cells_ok), or its roots are isolated with
         * certified interval arithmetic (Sin[x] > x/2). */
        bool* need = calloc((size_t)g_rd_nunrooted, sizeof(bool));
        bool any = false;
        for (int u = 0; u < g_rd_nunrooted; u++)
            if (!unrooted_cells_ok(x, bp, m, u)) { need[u] = true; any = true; }
        if (any) {
            /* Keep the known breakpoints: isolation may decline part-way. */
            Expr** bp0 = malloc((size_t)(m > 0 ? m : 1) * sizeof(Expr*));
            int m0 = m;
            for (int i = 0; i < m; i++) bp0[i] = expr_copy(bp[i]);
            unbounded = calloc((size_t)g_rd_nunrooted, sizeof(bool));
            if (m > 0) { span_lo = expr_copy(bp[0]); span_hi = expr_copy(bp[m - 1]); }
            ok = isolate_unrooted(x, &bp, &m, &cap, need, unbounded);
            if (ok) ok = sort_dedup_bp(bp, &m);
            if (!ok) {
                /* Isolation declined -- typically a residual that VANISHES on
                 * a whole interval (ArcCos[1 - 2 x^2] - 2 ArcSin[x] on [0, 1]),
                 * where no root can be isolated.  Fall back to deciding each
                 * cell separately: where the residual is zero at every probe
                 * it is taken as zero; elsewhere it is unknown, and the cell
                 * must be decided without it (Kleene) or the call declines. */
                for (int i = 0; i < m; i++) expr_free(bp[i]);
                m = 0;
                for (int i = 0; i < m0; i++) {
                    if (m == cap) { cap = cap ? 2 * cap : 8; bp = realloc(bp, (size_t)cap * sizeof(Expr*)); }
                    bp[m++] = bp0[i];
                }
                m0 = 0;
                if (span_lo) { expr_free(span_lo); span_lo = NULL; }
                if (span_hi) { expr_free(span_hi); span_hi = NULL; }
                free(unbounded); unbounded = NULL;
                fallback = need; need = NULL;
                ok = true;
            }
            for (int i = 0; i < m0; i++) expr_free(bp0[i]);
            free(bp0);
        }
        free(need);
        if (!ok) goto out;
    }

    {
        int ncells = 2 * m + 1;
        int* truth = malloc((size_t)ncells * sizeof(int));
        bool* unknown = (unbounded || fallback)
                      ? calloc((size_t)g_rd_nunrooted, sizeof(bool)) : NULL;
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
            /* Residuals whose roots are unknown at this sample: those with no
             * dominance bound, outside the isolated span [span_lo, span_hi].
             * The cell's truth is constant there only if it is decided
             * WITHOUT them (Kleene); else form_truth_general returns -1. */
            g_rd_unknown = NULL;
            if (unknown && unbounded) {
                bool outside = true;
                if (span_lo) {
                    int c1 = cmp_bp(sample, span_lo), c2 = cmp_bp(sample, span_hi);
                    outside = (c1 == -2 || c2 == -2 || c1 < 0 || c2 > 0);
                }
                for (int u = 0; u < g_rd_nunrooted; u++)
                    unknown[u] = unbounded[u] && outside;
                g_rd_unknown = unknown;
            }
            int t;
            if (fallback) {
                if (idx % 2 == 0) {
                    /* Interval cell: zero where every probe is zero, else unknown. */
                    for (int u = 0; u < g_rd_nunrooted; u++) {
                        if (!fallback[u]) continue;
                        if (unrooted_cell_class(x, bp, m, idx / 2, u) == 0) {
                            rd_register_zero(sample, g_rd_unrooted[u]);
                            unknown[u] = false;
                        } else {
                            unknown[u] = true;
                        }
                    }
                    g_rd_unknown = unknown;
                    t = form_truth_general(F, x, sample, ccons, cncon);
                } else {
                    /* Breakpoint: decided without the residual if possible,
                     * else by its exact value there. */
                    for (int u = 0; u < g_rd_nunrooted; u++) unknown[u] = fallback[u];
                    g_rd_unknown = unknown;
                    t = form_truth_general(F, x, sample, ccons, cncon);
                    if (t == -1) {
                        g_rd_unknown = NULL;
                        t = form_truth_general(F, x, sample, ccons, cncon);
                    }
                }
            } else {
                t = form_truth_general(F, x, sample, ccons, cncon);
            }
            g_rd_unknown = NULL;
            expr_free(sample);
            if (t == -1) { ok = false; break; }
            truth[idx] = t;
        }
        free(unknown);
        if (ok) result = rru_emit_sign_diagram(bp, m, truth, x);
        free(truth);
    }

out:
    for (int i = 0; i < m; i++) expr_free(bp[i]);
    free(bp);
    free_per_conj_cons(ccons, cncon, F->n); free(cncon); free(cccap);
    free(unbounded);
    free(fallback);
    if (span_lo) expr_free(span_lo);
    if (span_hi) expr_free(span_hi);
    for (int i = 0; i < g_rd_nreg; i++) expr_free(g_rd_regval[i]);
    free(g_rd_regval); free((void*)g_rd_regres);
    free((void*)g_rd_unrooted);
    g_rd_unrooted = saved_unrooted;
    g_rd_nunrooted = saved_nunrooted; g_rd_capunrooted = saved_capunrooted;
    g_rd_regval = saved_regval; g_rd_regres = saved_regres;
    g_rd_nreg = saved_nreg; g_rd_capreg = saved_capreg;
    g_rd_unknown = saved_unknown;
    return ok ? result : NULL;
}
