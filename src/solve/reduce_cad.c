/*
 * reduce_cad.c
 *
 * Cylindrical Algebraic Decomposition over the reals for `Reduce`
 * (REDUCE_PLAN.md, Phase 6).  See reduce_cad.h for the shape of the method.
 *
 * This pass implements the 2-variable case (after pruning variables that appear
 * in no atom).  Variables are (vx, vy) in the given order; vy is eliminated
 * first (projected out) so the base decomposition is in vx and the output reads
 * in the given variable order.
 *
 * Everything numeric is delegated to the evaluator (Discriminant / Resultant /
 * Coefficient / Exponent / FactorList / Solve[..,Reals] / ReplaceAll /
 * PolynomialQ) and to the exact real-algebraic primitives in reduce_real_util.h
 * (the qqbar sign oracle, rational-sample selection, Solve-based real-root
 * isolation).  The single structural rule that keeps the method sound at nv==2:
 * factor the atom polynomials into a distinct-irreducible squarefree basis and
 * run projection + fibre isolation on that basis, so a projection resultant /
 * discriminant is never identically zero and a fibre polynomial can nullify only
 * on a 0-dimensional section (never inside an open base interval).
 *
 * Hard invariant: an undecidable sign/ordering (qqbar returning -2), a fibre we
 * cannot cleanly isolate/order/emit, an irrational section sample (deferred), or
 * an interval nullification all make the engine return NULL (Reduce stays
 * unevaluated).  Soundness over completeness.
 */
#include "reduce_cad.h"
#include "reduce_real_util.h"
#include "reduce_univar.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Small node builders (each CONSUMES its Expr* arguments)            *
 * ------------------------------------------------------------------ */

static Expr* mkfun1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mkfun2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mkfun3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == name;
}

/* A plain (possibly non-rational) number, i.e. not a symbolic expression. */
static bool is_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_FUNCTION && is_head(e, SYM_Rational)) return true;
    return false;
}

/* An exact rational number (Integer / BigInt / Rational). */
static bool is_rational_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_FUNCTION && is_head(e, SYM_Rational)
        && e->data.function.arg_count == 2) return true;
    return false;
}

static bool is_zero(const Expr* e) {
    return e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Does the interned symbol `name` occur anywhere in `e`? */
static bool contains_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (contains_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Evaluator round-trips                                              *
 * ------------------------------------------------------------------ */

/* ReplaceAll[p, v -> s]. */
static Expr* subst1(const Expr* p, const Expr* v, const Expr* s) {
    Expr* rule = mkfun2(SYM_Rule, expr_copy((Expr*)v), expr_copy((Expr*)s));
    return eval_and_free(mkfun2(SYM_ReplaceAll, expr_copy((Expr*)p), rule));
}

/* ReplaceAll[p, {vx -> sx, vy -> sy}]. */
static Expr* subst2(const Expr* p, const Expr* vx, const Expr* sx,
                    const Expr* vy, const Expr* sy) {
    Expr* r1 = mkfun2(SYM_Rule, expr_copy((Expr*)vx), expr_copy((Expr*)sx));
    Expr* r2 = mkfun2(SYM_Rule, expr_copy((Expr*)vy), expr_copy((Expr*)sy));
    Expr* rules = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ r1, r2 }, 2);
    return eval_and_free(mkfun2(SYM_ReplaceAll, expr_copy((Expr*)p), rules));
}

/* Exponent[p, v] as an int; -1 when the degree is not a plain integer. */
static int degree_in(const Expr* p, const Expr* v) {
    Expr* e = eval_and_free(mkfun2(SYM_Exponent, expr_copy((Expr*)p), expr_copy((Expr*)v)));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

/* PolynomialQ[p, {vx, vy}]. */
static bool is_poly2(const Expr* p, const Expr* vx, const Expr* vy) {
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_copy((Expr*)vx), expr_copy((Expr*)vy) }, 2);
    Expr* q = eval_and_free(mkfun2(SYM_PolynomialQ, expr_copy((Expr*)p), vlist));
    bool r = is_sym(q, SYM_True);
    expr_free(q);
    return r;
}

/* ------------------------------------------------------------------ *
 *  Owned Expr* dynamic array helpers                                  *
 * ------------------------------------------------------------------ */

static void arr_push(Expr*** a, int* n, int* cap, Expr* v) {
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *a = realloc(*a, (size_t)*cap * sizeof(Expr*)); }
    (*a)[(*n)++] = v;
}

/* Append the distinct non-constant irreducible factors of `poly` (via
 * FactorList) to the basis (*a, *n, *cap), deduping by structural equality.
 * Returns false to bail on an unexpected FactorList shape. */
static bool add_factors(const Expr* poly, Expr*** a, int* n, int* cap) {
    Expr* fl = eval_and_free(mkfun1(SYM_FactorList, expr_copy((Expr*)poly)));
    bool ok = is_head(fl, SYM_List);
    for (size_t i = 0; ok && i < fl->data.function.arg_count; i++) {
        Expr* pair = fl->data.function.args[i];
        if (!is_head(pair, SYM_List) || pair->data.function.arg_count < 1) { ok = false; break; }
        Expr* fac = pair->data.function.args[0];
        if (is_number(fac)) continue;                 /* numeric content */
        bool dup = false;
        for (int j = 0; j < *n; j++) if (expr_eq((*a)[j], fac)) { dup = true; break; }
        if (dup) continue;
        arr_push(a, n, cap, expr_copy(fac));
    }
    expr_free(fl);
    return ok;
}

/* Add a projection polynomial (factored into its non-constant factors) to P_x. */
static bool add_proj(const Expr* poly, Expr*** a, int* n, int* cap) {
    if (is_number(poly)) return true;                 /* constant: no breakpoints */
    return add_factors(poly, a, n, cap);
}

/* ------------------------------------------------------------------ *
 *  Order + dedup a root array (ascending, exact) with optional        *
 *  provenance.  roots[]/fac[] are parallel (fac may be NULL).         *
 * ------------------------------------------------------------------ */

/* On success: dropped duplicates are freed, survivors compacted to roots[0..*m),
 * *m set, returns true.  On an undecidable comparison: nothing is freed (every
 * roots[0..nr) pointer stays valid for the caller to free) and returns false. */
static bool order_dedup(Expr** roots, int* fac, int nr, int* m) {
    if (nr < 2) { *m = (nr > 0) ? nr : 0; return true; }
    for (int i = 0; i < nr; i++) {
        int mn = i;
        for (int j = i + 1; j < nr; j++) {
            int c = rru_sign_compare(roots[j], roots[mn]);
            if (c == -2) return false;
            if (c < 0) mn = j;
        }
        if (mn != i) {
            Expr* t = roots[i]; roots[i] = roots[mn]; roots[mn] = t;
            if (fac) { int tf = fac[i]; fac[i] = fac[mn]; fac[mn] = tf; }
        }
    }
    bool* drop = nr ? calloc((size_t)nr, sizeof(bool)) : NULL;
    for (int i = 1; i < nr; i++) {
        int c = rru_sign_compare(roots[i - 1], roots[i]);   /* adjacent (sorted) */
        if (c == -2) { free(drop); return false; }          /* nothing freed yet */
        if (c == 0) drop[i] = true;
    }
    int w = 0;
    for (int i = 0; i < nr; i++) {
        if (drop && drop[i]) { expr_free(roots[i]); continue; }
        roots[w] = roots[i]; if (fac) fac[w] = fac[i]; w++;
    }
    free(drop);
    *m = w;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Segment emission (interval / point / one-sided) -- symbolic-bound  *
 *  aware; mirrors reduce_univar's seg_to_expr.                        *
 * ------------------------------------------------------------------ */

static Expr* cad_seg(const Expr* lo, bool lo_open,
                     const Expr* hi, bool hi_open, const Expr* v) {
    bool has_lo = (lo != NULL), has_hi = (hi != NULL);
    if (!has_lo && !has_hi) return expr_new_symbol(SYM_True);
    if (!has_lo) {
        const char* op = hi_open ? SYM_Less : SYM_LessEqual;
        return mkfun2(op, expr_copy((Expr*)v), expr_copy((Expr*)hi));
    }
    if (!has_hi) {
        const char* op = lo_open ? SYM_Greater : SYM_GreaterEqual;
        return mkfun2(op, expr_copy((Expr*)v), expr_copy((Expr*)lo));
    }
    if (expr_eq(lo, hi)) {                             /* isolated point: v == b */
        return mkfun2(SYM_Equal, expr_copy((Expr*)v), expr_copy((Expr*)lo));
    }
    const char* op1 = lo_open ? SYM_Less : SYM_LessEqual;
    const char* op2 = hi_open ? SYM_Less : SYM_LessEqual;
    return expr_new_function(expr_new_symbol(SYM_Inequality),
        (Expr*[]){ expr_copy((Expr*)lo), expr_new_symbol(op1),
                   expr_copy((Expr*)v),  expr_new_symbol(op2),
                   expr_copy((Expr*)hi) }, 5);
}

/* ------------------------------------------------------------------ *
 *  Truth of the DNF formula at a full (vx->sx, vy->sy) sample point   *
 * ------------------------------------------------------------------ */

static int atom_truth(const RAtom* a, const Expr* vx, const Expr* sx,
                      const Expr* vy, const Expr* sy) {
    if (a->rel == R_ELEM) return -1;
    Expr* v = subst2(a->poly, vx, sx, vy, sy);
    int s = rru_sign_of(v);
    expr_free(v);
    if (s == -2) return -1;
    switch (a->rel) {
        case R_EQ: return s == 0;
        case R_NE: return s != 0;
        case R_LT: return s < 0;
        case R_LE: return s <= 0;
        default:   return -1;
    }
}

static int form_truth(const RForm* F, const Expr* vx, const Expr* sx,
                      const Expr* vy, const Expr* sy) {
    if (F->is_true) return 1;
    for (int i = 0; i < F->n; i++) {
        RConj* c = F->c[i];
        int conj = 1;
        for (int k = 0; k < c->n; k++) {
            int t = atom_truth(&c->a[k], vx, sx, vy, sy);
            if (t == -1) return -1;
            if (t == 0) { conj = 0; break; }
        }
        if (conj) return 1;
    }
    return 0;
}

/* Partial-CAD pruning: is every conjunction of F already falsified at vx=sx by
 * an atom that does not involve vy?  Best-effort (an undecidable pure-x atom is
 * treated as non-killing, so we never prune a cell that might be satisfiable). */
static bool xcell_dead(const RForm* F, const Expr* vx, const Expr* sx, const Expr* vy) {
    for (int c = 0; c < F->n; c++) {
        RConj* cj = F->c[c];
        bool killed = false;
        for (int k = 0; k < cj->n && !killed; k++) {
            RAtom* a = &cj->a[k];
            if (a->rel == R_ELEM) continue;
            if (contains_symbol(a->poly, vy->data.symbol.name)) continue;   /* not pure-x */
            int s = rru_poly_sign_at(a->poly, vx, sx);
            if (s == -2) continue;
            int t;
            switch (a->rel) {
                case R_EQ: t = (s == 0); break;
                case R_NE: t = (s != 0); break;
                case R_LT: t = (s < 0);  break;
                case R_LE: t = (s <= 0); break;
                default:   t = 1; break;
            }
            if (!t) killed = true;
        }
        if (!killed) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Symbolic fibre bounds                                              *
 * ------------------------------------------------------------------ */

/* The branch expression (a function of vx) of basis factor B[i] whose value at
 * vx=sx equals the numeric fibre root `rt`.  Solve[B[i]==0, vy] is cached per
 * factor.  Returns a BORROWED pointer into the cache, or NULL if the factor's
 * Solve is not a clean orderable list of branches, or no branch matches. */
static Expr* symbolic_branch(Expr** B, int i, Expr** cache, char* cache_bad,
                             const Expr* vx, const Expr* sx,
                             const Expr* vy, const Expr* rt) {
    if (cache_bad[i]) return NULL;
    if (!cache[i]) {
        Expr* eqn = mkfun2(SYM_Equal, expr_copy(B[i]), expr_new_integer(0));
        Expr* sol = eval_and_free(mkfun2(SYM_Solve, eqn, expr_copy((Expr*)vy)));
        if (!is_head(sol, SYM_List)) { expr_free(sol); cache_bad[i] = 1; return NULL; }
        cache[i] = sol;
    }
    Expr* sol = cache[i];
    for (size_t r = 0; r < sol->data.function.arg_count; r++) {
        Expr* row = sol->data.function.args[r];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) return NULL;
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) return NULL;
        Expr* branch = rule->data.function.args[1];
        if (is_head(branch, SYM_ConditionalExpression)) return NULL;
        Expr* bs = subst1(branch, vx, sx);
        int c = rru_sign_compare(bs, rt);
        expr_free(bs);
        if (c == 0) return branch;          /* matched */
        /* c == -2 (non-real / undecided) or +/-1: not this branch, keep looking */
    }
    return NULL;                            /* no clean match -> caller bails */
}

/* Emission bound for a fibre breakpoint: the numeric root (over a point base
 * cell, where vx is fixed) or the matched symbolic branch, Simplify'd for
 * readability (over an interval base cell).  NULL => caller bails. */
static Expr* bound_expr(bool point_xcell, Expr** B, int facj,
                        Expr** cache, char* cache_bad,
                        const Expr* vx, const Expr* sx,
                        const Expr* vy, const Expr* rt) {
    if (point_xcell) return expr_copy((Expr*)rt);
    Expr* br = symbolic_branch(B, facj, cache, cache_bad, vx, sx, vy, rt);
    if (!br) return NULL;
    return eval_and_free(mkfun1(SYM_Simplify, expr_copy(br)));
}

/* ------------------------------------------------------------------ *
 *  Structured y-region (a union of ordered segments)                  *
 * ------------------------------------------------------------------ *
 * A y-region is the satisfying set in vy at a fixed base sample.  Over an
 * interval base cell the segment bounds are SYMBOLIC (functions of vx, matched
 * branch expressions); over a point section they are numeric.  Keeping the
 * region structured (rather than emitting an Expr immediately) is what lets the
 * boundary-merge pass compare an interval's template, evaluated at a breakpoint,
 * against that breakpoint's own fibre region. */

typedef struct { Expr* lo; bool lo_open; Expr* hi; bool hi_open; } YSeg;
typedef struct { YSeg* seg; int n, cap; bool all_true; } YRegion;

static void yregion_init(YRegion* r) { r->seg = NULL; r->n = 0; r->cap = 0; r->all_true = false; }

static void yregion_free(YRegion* r) {
    for (int i = 0; i < r->n; i++) { if (r->seg[i].lo) expr_free(r->seg[i].lo); if (r->seg[i].hi) expr_free(r->seg[i].hi); }
    free(r->seg); yregion_init(r);
}

/* Push a segment, TAKING OWNERSHIP of lo/hi (NULL => unbounded on that side). */
static void yregion_push(YRegion* r, Expr* lo, bool lo_open, Expr* hi, bool hi_open) {
    if (r->n == r->cap) { r->cap = r->cap ? r->cap * 2 : 4; r->seg = realloc(r->seg, (size_t)r->cap * sizeof(YSeg)); }
    r->seg[r->n].lo = lo; r->seg[r->n].lo_open = lo_open;
    r->seg[r->n].hi = hi; r->seg[r->n].hi_open = hi_open; r->n++;
}

static bool yregion_empty(const YRegion* r) { return !r->all_true && r->n == 0; }

/* Emit a region as a vy-formula (True / a segment / an Or of segments), or NULL
 * for the empty region. */
static Expr* yregion_to_expr(const YRegion* r, const Expr* vy) {
    if (r->all_true) return expr_new_symbol(SYM_True);
    if (r->n == 0) return NULL;
    if (r->n == 1) return cad_seg(r->seg[0].lo, r->seg[0].lo_open, r->seg[0].hi, r->seg[0].hi_open, vy);
    Expr** parts = malloc((size_t)r->n * sizeof(Expr*));
    for (int i = 0; i < r->n; i++)
        parts[i] = cad_seg(r->seg[i].lo, r->seg[i].lo_open, r->seg[i].hi, r->seg[i].hi_open, vy);
    Expr* o = expr_new_function(expr_new_symbol(SYM_Or), parts, (size_t)r->n);
    free(parts);
    return o;
}

/* Substitute vx=pt into every symbolic bound of T, giving a numeric region. */
static void yregion_at_point(const YRegion* T, const Expr* vx, const Expr* pt, YRegion* out) {
    yregion_init(out);
    out->all_true = T->all_true;
    for (int i = 0; i < T->n; i++) {
        Expr* lo = T->seg[i].lo ? subst1(T->seg[i].lo, vx, pt) : NULL;
        Expr* hi = T->seg[i].hi ? subst1(T->seg[i].hi, vx, pt) : NULL;
        yregion_push(out, lo, T->seg[i].lo_open, hi, T->seg[i].hi_open);
    }
}

/* Membership of a numeric sample: 1 in, 0 out, -1 undecidable. */
static int region_member(const YRegion* r, const Expr* sy) {
    if (r->all_true) return 1;
    for (int i = 0; i < r->n; i++) {
        const YSeg* s = &r->seg[i];
        bool in = true;
        if (s->lo) { int c = rru_sign_compare(sy, s->lo); if (c == -2) return -1; if (s->lo_open ? (c <= 0) : (c < 0)) in = false; }
        if (in && s->hi) { int c = rru_sign_compare(sy, s->hi); if (c == -2) return -1; if (s->hi_open ? (c >= 0) : (c > 0)) in = false; }
        if (in) return 1;
    }
    return 0;
}

/* Do two numeric regions describe the same subset of the vy line?  Decided by
 * sampling every cell of their combined breakpoint decomposition -- sound: any
 * undecidable comparison returns false (so the merge conservatively declines). */
static bool regions_equal(const YRegion* A, const YRegion* B) {
    if (A->all_true || B->all_true) return A->all_true && B->all_true;
    Expr** pts = NULL; int np = 0, pcap = 0;
    for (int i = 0; i < A->n; i++) { if (A->seg[i].lo) arr_push(&pts, &np, &pcap, expr_copy(A->seg[i].lo)); if (A->seg[i].hi) arr_push(&pts, &np, &pcap, expr_copy(A->seg[i].hi)); }
    for (int i = 0; i < B->n; i++) { if (B->seg[i].lo) arr_push(&pts, &np, &pcap, expr_copy(B->seg[i].lo)); if (B->seg[i].hi) arr_push(&pts, &np, &pcap, expr_copy(B->seg[i].hi)); }
    int m = 0;
    if (np > 0 && !order_dedup(pts, NULL, np, &m)) { for (int i = 0; i < np; i++) { expr_free(pts[i]); } free(pts); return false; }
    bool eq = true;
    for (int idx = 0; eq && idx < 2 * m + 1; idx++) {
        Expr* sy;
        if (idx % 2 == 1) { sy = expr_copy(pts[(idx + 1) / 2 - 1]); }
        else { int j = idx / 2; const Expr* lo = (j == 0) ? NULL : pts[j - 1]; const Expr* hi = (j == m) ? NULL : pts[j]; sy = rru_rational_between(lo, hi); if (!sy) { eq = false; break; } }
        int a = region_member(A, sy), b = region_member(B, sy);
        expr_free(sy);
        if (a == -1 || b == -1 || a != b) eq = false;
    }
    for (int i = 0; i < m; i++) { expr_free(pts[i]); }
    free(pts);
    return eq;
}

static bool exprs_equal_or_null(const Expr* a, const Expr* b) {
    if (a == NULL && b == NULL) return true;
    if (a == NULL || b == NULL) return false;
    return expr_eq(a, b);
}

/* Two interval templates are the SAME symbolic description (so the breakpoint
 * between them is a candidate for removal): identical segment count, bounds
 * (structurally) and openness. */
static bool templates_equal(const YRegion* A, const YRegion* B) {
    if (A->all_true != B->all_true) return false;
    if (A->all_true) return true;
    if (A->n != B->n) return false;
    for (int i = 0; i < A->n; i++) {
        if (A->seg[i].lo_open != B->seg[i].lo_open || A->seg[i].hi_open != B->seg[i].hi_open) return false;
        if (!exprs_equal_or_null(A->seg[i].lo, B->seg[i].lo)) return false;
        if (!exprs_equal_or_null(A->seg[i].hi, B->seg[i].hi)) return false;
    }
    return true;
}

/* Is the breakpoint x=pt (whose own fibre region is `sec`) exactly the limit of
 * interval template `T`?  If so the interval's closed x-range may absorb it. */
static bool breakpoint_absorbable(const YRegion* T, const Expr* vx, const Expr* pt, const YRegion* sec) {
    if (yregion_empty(sec)) return false;
    YRegion at; yregion_at_point(T, vx, pt, &at);
    bool ok = regions_equal(&at, sec);
    yregion_free(&at);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Lift one base cell to its structured satisfying y-region           *
 * ------------------------------------------------------------------ */

/* Fills *out with the satisfying set in vy at base sample vx=sx (over an
 * interval cell the bounds are symbolic functions of vx; over a point section
 * they are numeric).  Returns false to bail (undecidable / un-emittable /
 * interval nullification), in which case *out is left empty.  An empty *out
 * (yregion_empty) means the fibre is entirely false. */
static bool lift_fiber(const RForm* F, const Expr* vx, const Expr* sx, bool is_point,
                       const Expr* vy, Expr** B, int nb,
                       Expr** cache, char* cache_bad, YRegion* out) {
    yregion_init(out);
    bool fail = false;

    /* Isolate the fibre's real roots (with provenance) from the basis. */
    Expr** roots = NULL; int* fac = NULL; int nr = 0, cap = 0;
    for (int i = 0; i < nb; i++) {
        Expr* fs = subst1(B[i], vx, sx);
        if (is_zero(fs)) {                   /* nullification */
            expr_free(fs);
            if (is_point) continue;          /* 0-dim section: sound to skip */
            fail = true; break;              /* positive-dim interval: bail */
        }
        int dy = degree_in(fs, vy);
        if (dy <= 0) { expr_free(fs); continue; }      /* constant in vy */
        if (!rru_collect_roots(fs, vy, &roots, &nr, &cap, &fac, i, NULL)) { expr_free(fs); fail = true; break; }
        expr_free(fs);
    }
    if (fail) { for (int i = 0; i < nr; i++) { expr_free(roots[i]); } free(roots); free(fac); return false; }

    int my = 0;
    if (nr > 0 && !order_dedup(roots, fac, nr, &my)) {
        for (int i = 0; i < nr; i++) { expr_free(roots[i]); } free(roots); free(fac);
        return false;
    }

    /* No fibre breakpoints: the whole vy line is one cell. */
    if (my == 0) {
        free(roots); free(fac);
        Expr* z = expr_new_integer(0);
        int t = form_truth(F, vx, sx, vy, z);
        expr_free(z);
        if (t == -1) return false;
        out->all_true = (t == 1);        /* t==0 => empty region */
        return true;
    }

    /* Truth at one sample per fibre cell (2my+1, alternating interval/point). */
    int ny = 2 * my + 1;
    int* yt = malloc((size_t)ny * sizeof(int));
    for (int idx = 0; idx < ny && !fail; idx++) {
        Expr* sy;
        if (idx % 2 == 1) {
            sy = expr_copy(roots[(idx + 1) / 2 - 1]);
        } else {
            int j = idx / 2;
            const Expr* lo = (j == 0) ? NULL : roots[j - 1];
            const Expr* hi = (j == my) ? NULL : roots[j];
            sy = rru_rational_between(lo, hi);
            if (!sy) { fail = true; break; }
        }
        int t = form_truth(F, vx, sx, vy, sy);
        expr_free(sy);
        if (t == -1) { fail = true; break; }
        yt[idx] = t;
    }
    if (fail) {
        free(yt);
        for (int i = 0; i < my; i++) { expr_free(roots[i]); } free(roots); free(fac);
        return false;
    }

    bool all_true = true, all_false = true;
    for (int i = 0; i < ny; i++) { if (yt[i]) all_false = false; else all_true = false; }

    if (all_false) {
        /* out stays empty */
    } else if (all_true) {
        out->all_true = true;
    } else {
        /* Merge maximal runs of true cells into segments with symbolic bounds. */
        bool active = false; Expr* lo = NULL; bool lo_open = false;
        for (int idx = 0; idx < ny && !fail; idx++) {
            bool is_int = (idx % 2 == 0);
            if (yt[idx] && !active) {
                active = true;
                if (is_int) {
                    int j = idx / 2;
                    if (j == 0) { lo = NULL; lo_open = false; }
                    else { lo = bound_expr(is_point, B, fac[j - 1], cache, cache_bad, vx, sx, vy, roots[j - 1]);
                           lo_open = true; if (!lo) fail = true; }
                } else {
                    int k = (idx + 1) / 2 - 1;
                    lo = bound_expr(is_point, B, fac[k], cache, cache_bad, vx, sx, vy, roots[k]);
                    lo_open = false; if (!lo) fail = true;
                }
            } else if (!yt[idx] && active) {
                active = false;
                Expr* hi; bool hi_open;
                if (is_int) { int j = idx / 2; hi = bound_expr(is_point, B, fac[j - 1], cache, cache_bad, vx, sx, vy, roots[j - 1]); hi_open = false; }
                else { int k = (idx + 1) / 2 - 1; hi = bound_expr(is_point, B, fac[k], cache, cache_bad, vx, sx, vy, roots[k]); hi_open = true; }
                if (!hi) { fail = true; if (lo) expr_free(lo); lo = NULL; break; }
                yregion_push(out, lo, lo_open, hi, hi_open);     /* takes ownership */
                lo = NULL;
            }
        }
        if (!fail && active) { yregion_push(out, lo, lo_open, NULL, false); lo = NULL; }
        if (fail && lo) expr_free(lo);
    }

    free(yt);
    for (int i = 0; i < my; i++) { expr_free(roots[i]); }
    free(roots); free(fac);
    if (fail) { yregion_free(out); return false; }
    return true;
}

/* ================================================================== *
 *  n-variable recursive CAD (Phase 6d)                                *
 * ------------------------------------------------------------------ *
 * A single recursive engine for nu>=3 effective variables.  It is the
 * level-indexed generalization of the 2-variable driver above: the base
 * decomposition + fibre lift become a recursive descent over a projection
 * stack, and the fibre lift (`lift_fiber`) becomes the parameter-generalized
 * leaf `cad_leaf`.  Every soundness decision reduces exactly to the 2-var
 * behaviour at nu==2 (the 2-var path is kept separate and untouched, so its
 * carefully-tuned boundary-merged output never changes).
 *
 * v1 scope (rational-fibre regime): breakpoints at every NON-innermost level
 * must be rational (given the rational assignment above), else the engine
 * declines -- an irrational section value would create an algebraic-coefficient
 * fibre (real-algebraic-coefficient isolation is Phase 6b, deferred).  Emission
 * is a correct DNF of cells (the innermost dimension merged into a YRegion, the
 * whole flattened by evaluate); the n-D boundary merge that closes an outer
 * range is layered on later (Stage B).  Interval nullification bails (6e).
 * ================================================================== */

/* A set of owned, distinct-irreducible polynomial factors -- one projection
 * level of the McCallum stack. */
typedef struct { Expr** p; int n, cap; } PolySet;

static void polyset_free(PolySet* s) {
    for (int i = 0; i < s->n; i++) expr_free(s->p[i]);
    free(s->p); s->p = NULL; s->n = s->cap = 0;
}

/* ReplaceAll[p, {vv[0]->asg[0], ..., vv[k-1]->asg[k-1]}] (k rules; k==0 copies). */
static Expr* subst_n(const Expr* p, Expr** vv, Expr** asg, int k) {
    if (k <= 0) return expr_copy((Expr*)p);
    Expr** rules = malloc((size_t)k * sizeof(Expr*));
    for (int i = 0; i < k; i++)
        rules[i] = mkfun2(SYM_Rule, expr_copy(vv[i]), expr_copy(asg[i]));
    Expr* rlist = expr_new_function(expr_new_symbol(SYM_List), rules, (size_t)k);
    free(rules);
    return eval_and_free(mkfun2(SYM_ReplaceAll, expr_copy((Expr*)p), rlist));
}

/* PolynomialQ[p, {vv[0..nu-1]}]. */
static bool is_poly_n(const Expr* p, Expr** vv, int nu) {
    Expr** vs = malloc((size_t)nu * sizeof(Expr*));
    for (int i = 0; i < nu; i++) vs[i] = expr_copy(vv[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)nu);
    free(vs);
    Expr* q = eval_and_free(mkfun2(SYM_PolynomialQ, expr_copy((Expr*)p), vlist));
    bool r = is_sym(q, SYM_True);
    expr_free(q);
    return r;
}

/* Truth of one atom at a full assignment vv[0..nu-1]->asg[0..nu-1]. */
static int atom_truth_n(const RAtom* a, Expr** vv, Expr** asg, int nu) {
    if (a->rel == R_ELEM) return -1;
    Expr* v = subst_n(a->poly, vv, asg, nu);
    int s = rru_sign_of(v);
    expr_free(v);
    if (s == -2) return -1;
    switch (a->rel) {
        case R_EQ: return s == 0;
        case R_NE: return s != 0;
        case R_LT: return s < 0;
        case R_LE: return s <= 0;
        default:   return -1;
    }
}

static int form_truth_n(const RForm* F, Expr** vv, Expr** asg, int nu) {
    if (F->is_true) return 1;
    for (int i = 0; i < F->n; i++) {
        RConj* c = F->c[i];
        int conj = 1;
        for (int k = 0; k < c->n; k++) {
            int t = atom_truth_n(&c->a[k], vv, asg, nu);
            if (t == -1) return -1;
            if (t == 0) { conj = 0; break; }
        }
        if (conj) return 1;
    }
    return 0;
}

/* Partial-CAD pruning at level i: is every conjunction already falsified by an
 * atom that involves none of the still-free variables vv[i+1..nu-1]?  Such an
 * atom is testable at asg[0..i] (i+1 fixed values).  Undecidable atoms are
 * treated as non-killing (never over-prune).  Generalizes xcell_dead. */
static bool cell_dead_n(const RForm* F, Expr** vv, Expr** asg, int i, int nu) {
    for (int c = 0; c < F->n; c++) {
        RConj* cj = F->c[c];
        bool killed = false;
        for (int k = 0; k < cj->n && !killed; k++) {
            RAtom* a = &cj->a[k];
            if (a->rel == R_ELEM) continue;
            bool testable = true;
            for (int j = i + 1; j < nu; j++)
                if (contains_symbol(a->poly, vv[j]->data.symbol.name)) { testable = false; break; }
            if (!testable) continue;
            Expr* v = subst_n(a->poly, vv, asg, i + 1);
            int s = rru_sign_of(v);
            expr_free(v);
            if (s == -2) continue;
            int t;
            switch (a->rel) {
                case R_EQ: t = (s == 0); break;
                case R_NE: t = (s != 0); break;
                case R_LT: t = (s < 0);  break;
                case R_LE: t = (s <= 0); break;
                default:   t = 1; break;
            }
            if (!t) killed = true;
        }
        if (!killed) return false;
    }
    return true;
}

/* One McCallum projection step: eliminate `ve` from every poly of `in`
 * (leading coefficient + discriminant per poly, pairwise resultant, all
 * factored), appending distinct x-factors into `out`.  Mirrors the 2-var
 * projection body.  Returns false to bail on an unexpected round-trip shape. */
static bool cad_project_out(const PolySet* in, const Expr* ve, PolySet* out) {
    for (int i = 0; i < in->n; i++) {
        int dv = degree_in(in->p[i], ve);
        if (dv <= 0) { if (!add_proj(in->p[i], &out->p, &out->n, &out->cap)) return false; continue; }
        Expr* disc = eval_and_free(mkfun2("Discriminant", expr_copy(in->p[i]), expr_copy((Expr*)ve)));
        bool okd = add_proj(disc, &out->p, &out->n, &out->cap); expr_free(disc); if (!okd) return false;
        Expr* lc = eval_and_free(mkfun3(SYM_Coefficient, expr_copy(in->p[i]), expr_copy((Expr*)ve), expr_new_integer(dv)));
        bool okl = add_proj(lc, &out->p, &out->n, &out->cap); expr_free(lc); if (!okl) return false;
    }
    for (int i = 0; i < in->n; i++) {
        if (degree_in(in->p[i], ve) <= 0) continue;
        for (int j = i + 1; j < in->n; j++) {
            if (degree_in(in->p[j], ve) <= 0) continue;
            Expr* r = eval_and_free(mkfun3(SYM_Resultant, expr_copy(in->p[i]), expr_copy(in->p[j]), expr_copy((Expr*)ve)));
            bool okr = add_proj(r, &out->p, &out->n, &out->cap); expr_free(r); if (!okr) return false;
        }
    }
    return true;
}

/* The branch expression of `polys->p[fac]` (a function of vv[0..i-1]) whose
 * value at the outer assignment asg[0..i-1] equals the numeric root `rt`, from
 * Solve[polys->p[fac]==0, vv[i]] (cached per factor).  BORROWED pointer, or NULL
 * if Solve is not a clean orderable branch list or no branch matches.
 * Generalizes symbolic_branch: solve var is vv[i], outer subst is i values. */
static Expr* symbolic_branch_lvl(PolySet* polys, int fac, Expr** cache, char* cache_bad,
                                 Expr** vv, Expr** asg, int i, const Expr* rt) {
    if (cache_bad[fac]) return NULL;
    if (!cache[fac]) {
        Expr* eqn = mkfun2(SYM_Equal, expr_copy(polys->p[fac]), expr_new_integer(0));
        Expr* sol = eval_and_free(mkfun2(SYM_Solve, eqn, expr_copy(vv[i])));
        if (!is_head(sol, SYM_List)) { expr_free(sol); cache_bad[fac] = 1; return NULL; }
        cache[fac] = sol;
    }
    Expr* sol = cache[fac];
    for (size_t r = 0; r < sol->data.function.arg_count; r++) {
        Expr* row = sol->data.function.args[r];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) return NULL;
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) return NULL;
        Expr* branch = rule->data.function.args[1];
        if (is_head(branch, SYM_ConditionalExpression)) return NULL;
        Expr* bs = subst_n(branch, vv, asg, i);
        int c = rru_sign_compare(bs, rt);
        expr_free(bs);
        if (c == 0) return branch;
    }
    return NULL;
}

/* Emission bound for a breakpoint at level i: the numeric root (when all outer
 * variables are fixed to sections, so the bound is a constant) or the matched
 * symbolic branch, Simplify'd.  NULL => caller bails.  Generalizes bound_expr. */
static Expr* bound_expr_lvl(bool no_outer_sector, PolySet* polys, int fac,
                            Expr** cache, char* cache_bad,
                            Expr** vv, Expr** asg, int i, const Expr* rt) {
    if (no_outer_sector) return expr_copy((Expr*)rt);
    Expr* br = symbolic_branch_lvl(polys, fac, cache, cache_bad, vv, asg, i, rt);
    if (!br) return NULL;
    return eval_and_free(mkfun1(SYM_Simplify, expr_copy(br)));
}

/* Lift the innermost fibre (level d-1) at outer assignment asg[0..d-2] to its
 * satisfying y-region in vv[d-1].  This is lift_fiber generalized to an outer
 * assignment vector; `no_outer_sector` is the generalized `is_point` (true iff
 * every outer variable is a section, so fibre bounds are numeric).  Returns
 * false to bail; *out is left empty then.  asg[d-1] is used as scratch. */
static bool cad_leaf(const RForm* F, Expr** vv, int d, Expr** asg, bool no_outer_sector,
                     PolySet* basis, Expr** cache, char* cache_bad, YRegion* out) {
    yregion_init(out);
    bool fail = false;
    const Expr* vy = vv[d - 1];

    Expr** roots = NULL; int* fac = NULL; int nr = 0, cap = 0;
    for (int i = 0; i < basis->n; i++) {
        /* A factor free of the fibre variable contributes no fibre root; skip it
         * (it may vanish at an outer section without that being a McCallum
         * nullification -- only a fibre factor collapsing is). */
        if (!contains_symbol(basis->p[i], vy->data.symbol.name)) continue;
        Expr* fs = subst_n(basis->p[i], vv, asg, d - 1);
        if (is_zero(fs)) {                   /* fibre factor nullified */
            expr_free(fs);
            if (no_outer_sector) continue;   /* 0-dim outer cell: sound to skip */
            fail = true; break;              /* positive-dim: bail (6e) */
        }
        int dy = degree_in(fs, vy);
        if (dy <= 0) { expr_free(fs); continue; }
        if (!rru_collect_roots(fs, vy, &roots, &nr, &cap, &fac, i, NULL)) { expr_free(fs); fail = true; break; }
        expr_free(fs);
    }
    if (fail) { for (int i = 0; i < nr; i++) { expr_free(roots[i]); } free(roots); free(fac); return false; }

    int my = 0;
    if (nr > 0 && !order_dedup(roots, fac, nr, &my)) {
        for (int i = 0; i < nr; i++) { expr_free(roots[i]); } free(roots); free(fac);
        return false;
    }

    if (my == 0) {
        free(roots); free(fac);
        asg[d - 1] = expr_new_integer(0);
        int t = form_truth_n(F, vv, asg, d);
        expr_free(asg[d - 1]); asg[d - 1] = NULL;
        if (t == -1) return false;
        out->all_true = (t == 1);
        return true;
    }

    int ny = 2 * my + 1;
    int* yt = malloc((size_t)ny * sizeof(int));
    for (int idx = 0; idx < ny && !fail; idx++) {
        Expr* sy;
        if (idx % 2 == 1) sy = expr_copy(roots[(idx + 1) / 2 - 1]);
        else {
            int j = idx / 2;
            const Expr* lo = (j == 0) ? NULL : roots[j - 1];
            const Expr* hi = (j == my) ? NULL : roots[j];
            sy = rru_rational_between(lo, hi);
            if (!sy) { fail = true; break; }
        }
        asg[d - 1] = sy;                     /* borrow for the truth evaluation */
        int t = form_truth_n(F, vv, asg, d);
        asg[d - 1] = NULL;
        expr_free(sy);
        if (t == -1) { fail = true; break; }
        yt[idx] = t;
    }
    if (fail) {
        free(yt);
        for (int i = 0; i < my; i++) { expr_free(roots[i]); } free(roots); free(fac);
        return false;
    }

    bool all_true = true, all_false = true;
    for (int i = 0; i < ny; i++) { if (yt[i]) all_false = false; else all_true = false; }

    if (all_false) {
        /* out stays empty */
    } else if (all_true) {
        out->all_true = true;
    } else {
        bool active = false; Expr* lo = NULL; bool lo_open = false;
        for (int idx = 0; idx < ny && !fail; idx++) {
            bool is_int = (idx % 2 == 0);
            if (yt[idx] && !active) {
                active = true;
                if (is_int) {
                    int j = idx / 2;
                    if (j == 0) { lo = NULL; lo_open = false; }
                    else { lo = bound_expr_lvl(no_outer_sector, basis, fac[j - 1], cache, cache_bad, vv, asg, d - 1, roots[j - 1]);
                           lo_open = true; if (!lo) fail = true; }
                } else {
                    int k = (idx + 1) / 2 - 1;
                    lo = bound_expr_lvl(no_outer_sector, basis, fac[k], cache, cache_bad, vv, asg, d - 1, roots[k]);
                    lo_open = false; if (!lo) fail = true;
                }
            } else if (!yt[idx] && active) {
                active = false;
                Expr* hi; bool hi_open;
                if (is_int) { int j = idx / 2; hi = bound_expr_lvl(no_outer_sector, basis, fac[j - 1], cache, cache_bad, vv, asg, d - 1, roots[j - 1]); hi_open = false; }
                else { int k = (idx + 1) / 2 - 1; hi = bound_expr_lvl(no_outer_sector, basis, fac[k], cache, cache_bad, vv, asg, d - 1, roots[k]); hi_open = true; }
                if (!hi) { fail = true; if (lo) expr_free(lo); lo = NULL; break; }
                yregion_push(out, lo, lo_open, hi, hi_open);
                lo = NULL;
            }
        }
        if (!fail && active) { yregion_push(out, lo, lo_open, NULL, false); lo = NULL; }
        if (fail && lo) expr_free(lo);
    }

    free(yt);
    for (int i = 0; i < my; i++) { expr_free(roots[i]); }
    free(roots); free(fac);
    if (fail) { yregion_free(out); return false; }
    return true;
}

/* ================================================================== *
 *  Stage B: recursive CAD region tree + n-D boundary merge            *
 * ------------------------------------------------------------------ *
 * The recursion builds a cell TREE rather than emitting directly, so the
 * 2-var boundary merge (which closes a non-strict outer range by absorbing a
 * section into the adjacent interval) generalizes to n variables.  Merge
 * decisions are decided by SAMPLING: two cells are "equal" iff their emitted
 * sub-formulas agree at a grid of sample points drawn from the cell structure;
 * any undecidable comparison leaves the (already-correct) unmerged form.
 * ================================================================== */

typedef struct CADRegion CADRegion;

/* One vv[level]-cell's satisfying set below it, in vv[level+1..d-1]. */
typedef struct {
    bool       empty;   /* nothing satisfies below this cell (dead / false)   */
    bool       leaf;    /* true: yr (innermost vv[d-1]); false: sub (nested)  */
    YRegion    yr;      /* valid when leaf                                    */
    CADRegion* sub;     /* valid when !leaf; owned                            */
    Expr*      sample;  /* numeric vv[level] value used to build this cell    */
} CADCell;

struct CADRegion {
    int      level, m;
    Expr**   bpsym;     /* [m] symbolic vv[level] breakpoint bounds (fn of vv[0..level-1]) */
    Expr**   bpnum;     /* [m] numeric vv[level] breakpoint values at the build context   */
    CADCell* iv;        /* [m+1] interval cells */
    CADCell* sec;       /* [m]   section cells  */
    bool     all_true;  /* every cell non-empty and fully satisfied */
    bool     all_false; /* no cell has any satisfying point         */
};

static void cad_region_free(CADRegion* R);

static void cadcell_free(CADCell* c) {
    if (c->sample) expr_free(c->sample);
    if (c->leaf) yregion_free(&c->yr);
    else if (c->sub) cad_region_free(c->sub);
}

static void cad_region_free(CADRegion* R) {
    if (!R) return;
    for (int k = 0; k < R->m; k++) {
        if (R->bpsym && R->bpsym[k]) expr_free(R->bpsym[k]);
        if (R->bpnum && R->bpnum[k]) expr_free(R->bpnum[k]);
    }
    free(R->bpsym); free(R->bpnum);
    if (R->iv) for (int j = 0; j <= R->m; j++) cadcell_free(&R->iv[j]);
    if (R->sec) for (int k = 0; k < R->m; k++) cadcell_free(&R->sec[k]);
    free(R->iv); free(R->sec);
    free(R);
}

/* Build the cell tree for vv[level..d-1] at the (rational) outer assignment
 * asg[0..level-1].  Returns NULL to bail (any of the Stage-A soundness declines).
 * asg[level] is used as scratch by the recursion. */
static CADRegion* cad_build(const RForm* F, Expr** vv, int d, int level, Expr** asg,
                            bool no_outer_sector, PolySet* pstack,
                            Expr** caches[], char* caches_bad[]) {
    PolySet* polys = &pstack[level];
    Expr** roots = NULL; int* fac = NULL; int nr = 0, cap = 0; bool fail = false;
    for (int j = 0; j < polys->n && !fail; j++) {
        if (!contains_symbol(polys->p[j], vv[level]->data.symbol.name)) continue;
        Expr* fs = subst_n(polys->p[j], vv, asg, level);
        if (is_zero(fs)) { expr_free(fs); if (no_outer_sector) continue; fail = true; break; }
        int dv = degree_in(fs, vv[level]);
        if (dv <= 0) { expr_free(fs); continue; }
        if (!rru_collect_roots(fs, vv[level], &roots, &nr, &cap, &fac, j, NULL)) { expr_free(fs); fail = true; break; }
        expr_free(fs);
    }
    if (fail) { for (int q = 0; q < nr; q++) { expr_free(roots[q]); } free(roots); free(fac); return NULL; }

    int m = 0;
    if (nr > 0 && !order_dedup(roots, fac, nr, &m)) {
        for (int q = 0; q < nr; q++) { expr_free(roots[q]); } free(roots); free(fac); return NULL;
    }
    if (level < d - 1)
        for (int q = 0; q < m; q++)
            if (!is_rational_number(roots[q])) { for (int r = 0; r < m; r++) { expr_free(roots[r]); } free(roots); free(fac); return NULL; }

    size_t mz = (m > 0) ? (size_t)m : 0;
    CADRegion* R = calloc(1, sizeof(CADRegion));
    R->level = level; R->m = m;
    R->bpnum = mz ? calloc(mz, sizeof(Expr*)) : NULL;
    R->bpsym = mz ? calloc(mz, sizeof(Expr*)) : NULL;
    R->iv = calloc(mz + 1, sizeof(CADCell));
    R->sec = mz ? calloc(mz, sizeof(CADCell)) : NULL;

    for (int k = 0; k < m && !fail; k++) {
        R->bpnum[k] = expr_copy(roots[k]);
        R->bpsym[k] = bound_expr_lvl(no_outer_sector, polys, fac[k], caches[level], caches_bad[level], vv, asg, level, roots[k]);
        if (!R->bpsym[k]) fail = true;
    }

    bool all_true = true, any = false;
    for (int idx = 0; idx <= 2 * m && !fail; idx++) {
        bool is_section = (idx % 2 == 1);
        CADCell* cell; Expr* s;
        if (is_section) { int k = (idx - 1) / 2; cell = &R->sec[k]; s = expr_copy(roots[k]); }
        else {
            int j = idx / 2;
            const Expr* lo = (j == 0) ? NULL : roots[j - 1];
            const Expr* hi = (j == m) ? NULL : roots[j];
            cell = &R->iv[j]; s = rru_rational_between(lo, hi);
            if (!s) { fail = true; break; }
        }
        cell->sample = s;
        asg[level] = s;                       /* borrow (cell owns s) */
        if (cell_dead_n(F, vv, asg, level, d)) { cell->empty = true; all_true = false; asg[level] = NULL; continue; }
        bool child_nos = no_outer_sector && is_section;
        if (level == d - 2) {
            cell->leaf = true;
            if (!cad_leaf(F, vv, d, asg, child_nos, &pstack[d - 1], caches[d - 1], caches_bad[d - 1], &cell->yr)) { fail = true; asg[level] = NULL; break; }
            if (yregion_empty(&cell->yr)) { cell->empty = true; all_true = false; }
            else { any = true; if (!cell->yr.all_true) all_true = false; }
        } else {
            cell->sub = cad_build(F, vv, d, level + 1, asg, child_nos, pstack, caches, caches_bad);
            if (!cell->sub) { fail = true; asg[level] = NULL; break; }
            if (cell->sub->all_false) { cell->empty = true; all_true = false; }
            else { any = true; if (!cell->sub->all_true) all_true = false; }
        }
        asg[level] = NULL;
    }

    for (int q = 0; q < m; q++) { expr_free(roots[q]); }
    free(roots); free(fac);
    if (fail) { cad_region_free(R); return NULL; }
    R->all_true = all_true; R->all_false = !any;
    return R;
}

/* ---- emission (raw + merged) ------------------------------------- */

static Expr* cad_region_expr(CADRegion* R, Expr** vv, int d, bool merge,
                             Expr** ctxvars, Expr** ctxvals, int nctx);

/* Emit one non-empty cell's subtree.  For a nested cell the outer context ctx
 * (used only by the merge) is extended by this cell's (vv[level], sample). */
static Expr* cad_cell_expr(CADCell* c, Expr** vv, int d, bool merge,
                           Expr** ctxvars, Expr** ctxvals, int nctx) {
    if (c->leaf) return yregion_to_expr(&c->yr, vv[d - 1]);
    return cad_region_expr(c->sub, vv, d, merge, ctxvars, ctxvals, nctx);
}

/* Truth of an emitted sub-formula E at the full numeric assignment asg[0..d-1]:
 * 1 (true), 0 (false), or -1 (did not decide). */
static int formula_truth_at(const Expr* E, Expr** vv, Expr** asg, int d) {
    Expr* s = subst_n(E, vv, asg, d);
    int r = is_sym(s, SYM_True) ? 1 : (is_sym(s, SYM_False) ? 0 : -1);
    expr_free(s);
    return r;
}

/* Sample the cell `c` (a subtree in vv[pos..d-1]) at the fixed outer assignment
 * asg[0..pos-1], and at each full sample point require the two emitted formulas
 * Ea, Eb to agree.  Sets *ok=false on any disagreement/undecidable; counts
 * decided points in *judged. */
static void cad_sample_cell(CADCell* c, int pos, Expr** vv, int d, Expr** asg,
                            const Expr* Ea, const Expr* Eb, int* judged, bool* ok) {
    if (!*ok) return;
    if (c->empty) {                            /* no substructure: one 0-filled probe */
        for (int t = pos; t < d; t++) asg[t] = expr_new_integer(0);
        int a = formula_truth_at(Ea, vv, asg, d), b = formula_truth_at(Eb, vv, asg, d);
        if (a == -1 || b == -1 || a != b) *ok = false; else (*judged)++;
        for (int t = pos; t < d; t++) { expr_free(asg[t]); asg[t] = NULL; }
        return;
    }
    if (c->leaf) {                             /* pos == d-1: sample vv[d-1] over the yr bounds */
        Expr** bv = NULL; int nb = 0, cb = 0;
        for (int i = 0; i < c->yr.n; i++) {
            if (c->yr.seg[i].lo) arr_push(&bv, &nb, &cb, subst_n(c->yr.seg[i].lo, vv, asg, d - 1));
            if (c->yr.seg[i].hi) arr_push(&bv, &nb, &cb, subst_n(c->yr.seg[i].hi, vv, asg, d - 1));
        }
        int mm = 0;
        if (nb > 0 && !order_dedup(bv, NULL, nb, &mm)) { for (int i = 0; i < nb; i++) expr_free(bv[i]); free(bv); *ok = false; return; }
        for (int idx = 0; *ok && idx <= 2 * mm; idx++) {
            Expr* v;
            if (mm == 0) v = expr_new_integer(0);
            else if (idx % 2 == 1) v = expr_copy(bv[(idx - 1) / 2]);
            else { int j = idx / 2; const Expr* lo = (j == 0) ? NULL : bv[j - 1]; const Expr* hi = (j == mm) ? NULL : bv[j]; v = rru_rational_between(lo, hi); if (!v) continue; }
            asg[d - 1] = v;
            int a = formula_truth_at(Ea, vv, asg, d), b = formula_truth_at(Eb, vv, asg, d);
            if (a == -1 || b == -1 || a != b) *ok = false; else (*judged)++;
            expr_free(v); asg[d - 1] = NULL;
            if (mm == 0) break;
        }
        for (int i = 0; i < mm; i++) expr_free(bv[i]);
        free(bv);
        return;
    }
    CADRegion* R = c->sub; int L = R->level;   /* == pos */
    for (int j = 0; *ok && j <= R->m; j++) {
        Expr* lo = (j == 0) ? NULL : subst_n(R->bpsym[j - 1], vv, asg, L);
        Expr* hi = (j == R->m) ? NULL : subst_n(R->bpsym[j], vv, asg, L);
        Expr* s = rru_rational_between(lo, hi);
        if (lo) expr_free(lo);
        if (hi) expr_free(hi);
        if (!s) continue;                      /* degenerate interval: no interior probe */
        asg[L] = s;
        cad_sample_cell(&R->iv[j], L + 1, vv, d, asg, Ea, Eb, judged, ok);
        expr_free(s); asg[L] = NULL;
    }
    for (int k = 0; *ok && k < R->m; k++) {
        Expr* s = subst_n(R->bpsym[k], vv, asg, L);
        asg[L] = s;
        cad_sample_cell(&R->sec[k], L + 1, vv, d, asg, Ea, Eb, judged, ok);
        expr_free(s); asg[L] = NULL;
    }
}

/* Do cells A and B describe the same set in vv[nctx..d-1], at the outer context
 * ctxvars/ctxvals[0..nctx-1]?  Decided by sampling A's and B's structure and
 * requiring their emitted formulas to agree.  Sound: any doubt returns false. */
static bool cad_cell_equiv(CADCell* A, CADCell* B, Expr** vv, int d,
                           Expr** ctxvars, Expr** ctxvals, int nctx) {
    Expr* Ea = A->empty ? expr_new_symbol(SYM_False) : cad_cell_expr(A, vv, d, false, NULL, NULL, 0);
    Expr* Eb = B->empty ? expr_new_symbol(SYM_False) : cad_cell_expr(B, vv, d, false, NULL, NULL, 0);
    Expr** asg = calloc((size_t)d, sizeof(Expr*));
    for (int i = 0; i < nctx; i++) asg[i] = ctxvals[i];   /* borrowed */
    (void)ctxvars;
    int judged = 0; bool ok = true;
    cad_sample_cell(B, nctx, vv, d, asg, Ea, Eb, &judged, &ok);
    if (ok) cad_sample_cell(A, nctx, vv, d, asg, Ea, Eb, &judged, &ok);
    for (int i = 0; i < nctx; i++) asg[i] = NULL;          /* release borrows */
    free(asg);
    expr_free(Ea); expr_free(Eb);
    return ok && judged > 0;
}

/* Structural sameness of two interval templates (so a run of intervals sharing
 * one description can be grouped across an absorbed breakpoint). */
static bool cad_templates_equal(CADCell* A, CADCell* B) {
    if (A->empty || B->empty) return A->empty && B->empty;
    if (A->leaf != B->leaf) return false;
    if (A->leaf) return templates_equal(&A->yr, &B->yr);
    CADRegion* RA = A->sub; CADRegion* RB = B->sub;
    if (RA->m != RB->m) return false;
    for (int k = 0; k < RA->m; k++) if (!exprs_equal_or_null(RA->bpsym[k], RB->bpsym[k])) return false;
    for (int j = 0; j <= RA->m; j++) if (!cad_templates_equal(&RA->iv[j], &RB->iv[j])) return false;
    for (int k = 0; k < RA->m; k++) if (!cad_templates_equal(&RA->sec[k], &RB->sec[k])) return false;
    return true;
}

/* Is section sec[k] the limit of interval iv[j] at breakpoint bpnum[k]?  I.e.
 * does iv[j]'s template, with vv[level] fixed to the breakpoint, equal sec[k]? */
static bool cad_absorbable(CADRegion* R, CADCell* ivcell, CADCell* seccell, int k,
                           Expr** vv, int d, Expr** ctxvars, Expr** ctxvals, int nctx) {
    if (seccell->empty) return false;
    int L = R->level;
    Expr** cv = malloc((size_t)(nctx + 1) * sizeof(Expr*));
    Expr** cx = malloc((size_t)(nctx + 1) * sizeof(Expr*));
    for (int i = 0; i < nctx; i++) { cv[i] = ctxvars[i]; cx[i] = ctxvals[i]; }
    cv[nctx] = vv[L]; cx[nctx] = R->bpnum[k];
    bool r = cad_cell_equiv(ivcell, seccell, vv, d, cv, cx, nctx + 1);
    free(cv); free(cx);
    return r;
}

/* Emit region R (a decomposition of vv[level]) as an Expr.  With merge, runs of
 * consecutive interval cells sharing one template absorb the sections between
 * (and at the ends) whose fibre equals the template's limit there, closing the
 * vv[level] range; without merge, every non-empty cell is emitted open. */
static Expr* cad_region_expr(CADRegion* R, Expr** vv, int d, bool merge,
                             Expr** ctxvars, Expr** ctxvals, int nctx) {
    if (R->all_true) return expr_new_symbol(SYM_True);
    if (R->all_false) return expr_new_symbol(SYM_False);
    int L = R->level, m = R->m;

    Expr** parts = NULL; int np = 0, cap = 0;
    bool* absorbed = (merge && m) ? calloc((size_t)m, sizeof(bool)) : NULL;

    for (int j = 0; j <= m; ) {
        if (R->iv[j].empty) { j++; continue; }
        int hi = j;
        if (merge) {
            while (hi < m && !R->iv[hi + 1].empty
                   && cad_templates_equal(&R->iv[j], &R->iv[hi + 1])
                   && cad_absorbable(R, &R->iv[j], &R->sec[hi], hi, vv, d, ctxvars, ctxvals, nctx)) {
                absorbed[hi] = true; hi++;
            }
        }
        bool left_closed = false, right_closed = false;
        if (merge) {
            left_closed  = (j > 0)   && !absorbed[j - 1] && cad_absorbable(R, &R->iv[j], &R->sec[j - 1], j - 1, vv, d, ctxvars, ctxvals, nctx);
            right_closed = (hi < m)   && cad_absorbable(R, &R->iv[j], &R->sec[hi], hi, vv, d, ctxvars, ctxvals, nctx);
            if (left_closed)  absorbed[j - 1] = true;
            if (right_closed) absorbed[hi]    = true;
        }
        const Expr* lo = (j == 0)  ? NULL : R->bpsym[j - 1];
        const Expr* hb = (hi == m) ? NULL : R->bpsym[hi];
        Expr* xcond = cad_seg(lo, !left_closed, hb, !right_closed, vv[L]);

        Expr** ecv = malloc((size_t)(nctx + 1) * sizeof(Expr*));
        Expr** ecx = malloc((size_t)(nctx + 1) * sizeof(Expr*));
        for (int i = 0; i < nctx; i++) { ecv[i] = ctxvars[i]; ecx[i] = ctxvals[i]; }
        ecv[nctx] = vv[L]; ecx[nctx] = R->iv[j].sample;
        Expr* child = cad_cell_expr(&R->iv[j], vv, d, merge, ecv, ecx, nctx + 1);
        free(ecv); free(ecx);

        if (is_sym(child, SYM_True)) { expr_free(child); arr_push(&parts, &np, &cap, xcond); }
        else arr_push(&parts, &np, &cap, expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ xcond, child }, 2));
        j = hi + 1;
    }

    for (int k = 0; k < m; k++) {
        if ((absorbed && absorbed[k]) || R->sec[k].empty) continue;
        Expr* xeq = mkfun2(SYM_Equal, expr_copy(vv[L]), expr_copy(R->bpsym[k]));
        Expr** ecv = malloc((size_t)(nctx + 1) * sizeof(Expr*));
        Expr** ecx = malloc((size_t)(nctx + 1) * sizeof(Expr*));
        for (int i = 0; i < nctx; i++) { ecv[i] = ctxvars[i]; ecx[i] = ctxvals[i]; }
        ecv[nctx] = vv[L]; ecx[nctx] = R->bpnum[k];
        Expr* child = cad_cell_expr(&R->sec[k], vv, d, merge, ecv, ecx, nctx + 1);
        free(ecv); free(ecx);
        if (is_sym(child, SYM_True)) { expr_free(child); arr_push(&parts, &np, &cap, xeq); }
        else arr_push(&parts, &np, &cap, expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ xeq, child }, 2));
    }

    free(absorbed);
    Expr* out;
    if (np == 0) out = expr_new_symbol(SYM_False);
    else if (np == 1) out = parts[0];
    else out = expr_new_function(expr_new_symbol(SYM_Or), parts, (size_t)np);
    free(parts);
    return out;
}

/* Driver for nu>=3 effective variables (vv[0..d-1] in the given order). */
static Expr* reduce_cad_nvar(const RForm* F, Expr** vv, int d) {
    PolySet* pstack = calloc((size_t)d, sizeof(PolySet));
    Expr** asg = calloc((size_t)d, sizeof(Expr*));
    Expr*** caches = calloc((size_t)d, sizeof(Expr**));
    char** caches_bad = calloc((size_t)d, sizeof(char*));
    CADRegion* root = NULL;
    Expr* out = NULL;

    /* Gate + basis (pstack[d-1]): every atom a d-variate polynomial relation,
     * factored into the distinct-irreducible squarefree basis. */
    for (int c = 0; c < F->n; c++) {
        RConj* cj = F->c[c];
        for (int k = 0; k < cj->n; k++) {
            RAtom* a = &cj->a[k];
            if (a->rel == R_ELEM || a->nonconst_denom) goto done;
            if (!is_poly_n(a->poly, vv, d)) goto done;
            if (!add_factors(a->poly, &pstack[d - 1].p, &pstack[d - 1].n, &pstack[d - 1].cap)) goto done;
        }
    }
    if (pstack[d - 1].n == 0) goto done;

    /* Iterated McCallum projection: pstack[k-1] eliminates vv[k] from pstack[k]. */
    for (int k = d - 1; k >= 1; k--)
        if (!cad_project_out(&pstack[k], vv[k], &pstack[k - 1])) goto done;

    for (int i = 0; i < d; i++) {
        caches[i]     = pstack[i].n ? calloc((size_t)pstack[i].n, sizeof(Expr*)) : NULL;
        caches_bad[i] = pstack[i].n ? calloc((size_t)pstack[i].n, sizeof(char)) : NULL;
    }

    root = cad_build(F, vv, d, 0, asg, true, pstack, caches, caches_bad);
    if (!root) goto done;
    out = cad_region_expr(root, vv, d, true, NULL, NULL, 0);
    out = eval_and_free(out);

done:
    cad_region_free(root);
    for (int i = 0; i < d; i++) {
        if (caches[i]) { for (int j = 0; j < pstack[i].n; j++) if (caches[i][j]) expr_free(caches[i][j]); free(caches[i]); }
        free(caches_bad[i]);
        polyset_free(&pstack[i]);
    }
    free(caches); free(caches_bad);
    free(pstack); free(asg);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Quantifier-elimination seam (REDUCE_PLAN.md, Phase 7, Case B)       *
 * ------------------------------------------------------------------ */

/* The Exists/ForAll verdict of one CAD cell over its bound-variable subtree.
 *   Exists (quant==0): the fibre is non-empty -- cad_build already set
 *     cell->empty EXACTLY when no bound point satisfies (cell_dead_n, an empty
 *     leaf yr, or a nested sub->all_false), so !empty is the Exists fold.
 *   ForAll (quant==1): every bound point satisfies -- the all_true roll-up
 *     (a leaf's yr.all_true, or a nested sub's all_true); an empty fibre fails. */
static int qe_cell_verdict(const CADCell* c, int quant) {
    if (quant == 0) return c->empty ? 0 : 1;              /* Exists */
    if (c->empty) return 0;                               /* ForAll, empty fibre */
    return c->leaf ? (c->yr.all_true ? 1 : 0)
                   : (c->sub->all_true ? 1 : 0);
}

/* Single-free-variable quantifier elimination.  Eliminate boundvars[0..nbound-1]
 * from the DNF `F` (a statement in freevar and the bound vars) under `quant`
 * (0 = Exists, 1 = ForAll), returning the quantifier-free description in freevar,
 * or NULL to decline.  freevar is the OUTERMOST CAD level (level 0) and the bound
 * vars are inner, so McCallum projection eliminates the bound vars first; the
 * free variable's 2m+1 cells then carry the per-cell quantifier verdict, which
 * the shared 1-D sign diagram merges into the answer.  F and every Expr* are
 * BORROWED; the returned Expr is freshly owned.  Mirrors reduce_cad_nvar's setup,
 * projection and single-exit teardown for leak-freedom. */
Expr* reduce_cad_qe(const RForm* F, Expr* freevar,
                    Expr** boundvars, int nbound, int quant) {
    /* Step 0: a constant statement decides directly (a literal-True/False fibre
     * is the same answer under either quantifier). */
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);
    if (nbound < 1) return NULL;                          /* front-end must strip */

    int d = 1 + nbound;
    Expr** vv = malloc((size_t)d * sizeof(Expr*));
    vv[0] = freevar;
    for (int i = 0; i < nbound; i++) vv[1 + i] = boundvars[i];

    PolySet* pstack = calloc((size_t)d, sizeof(PolySet));
    Expr** asg = calloc((size_t)d, sizeof(Expr*));
    Expr*** caches = calloc((size_t)d, sizeof(Expr**));
    char** caches_bad = calloc((size_t)d, sizeof(char*));
    CADRegion* root = NULL;
    int* truth = NULL;
    Expr* out = NULL;

    /* Gate + basis (pstack[d-1]): every atom a d-variate polynomial relation,
     * factored into the distinct-irreducible squarefree basis. */
    for (int c = 0; c < F->n; c++) {
        RConj* cj = F->c[c];
        for (int k = 0; k < cj->n; k++) {
            RAtom* a = &cj->a[k];
            if (a->rel == R_ELEM || a->nonconst_denom) goto done;
            if (!is_poly_n(a->poly, vv, d)) goto done;
            if (!add_factors(a->poly, &pstack[d - 1].p, &pstack[d - 1].n, &pstack[d - 1].cap)) goto done;
        }
    }
    if (pstack[d - 1].n == 0) goto done;

    /* Iterated McCallum projection: pstack[k-1] eliminates vv[k] from pstack[k]. */
    for (int k = d - 1; k >= 1; k--)
        if (!cad_project_out(&pstack[k], vv[k], &pstack[k - 1])) goto done;

    for (int i = 0; i < d; i++) {
        caches[i]     = pstack[i].n ? calloc((size_t)pstack[i].n, sizeof(Expr*)) : NULL;
        caches_bad[i] = pstack[i].n ? calloc((size_t)pstack[i].n, sizeof(char)) : NULL;
    }

    root = cad_build(F, vv, d, 0, asg, true, pstack, caches, caches_bad);
    if (!root) goto done;

    /* The free variable is level 0: its interval/section cells (2m+1, alternating)
     * carry the quantifier verdict over the bound subtree.  rru_emit_sign_diagram
     * COPIES the breakpoints, so it is called while `root` is still alive. */
    {
        int m = root->m;
        truth = malloc((size_t)(2 * m + 1) * sizeof(int));
        for (int j = 0; j <= m; j++) truth[2 * j]     = qe_cell_verdict(&root->iv[j], quant);
        for (int k = 0; k < m; k++)  truth[2 * k + 1] = qe_cell_verdict(&root->sec[k], quant);
        out = rru_emit_sign_diagram(root->bpsym, m, truth, freevar);
    }

done:
    cad_region_free(root);
    free(truth);
    for (int i = 0; i < d; i++) {
        if (caches[i]) { for (int j = 0; j < pstack[i].n; j++) if (caches[i][j]) expr_free(caches[i][j]); free(caches[i]); }
        free(caches_bad[i]);
        polyset_free(&pstack[i]);
    }
    free(caches); free(caches_bad);
    free(pstack); free(asg); free(vv);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* reduce_cad(const RForm* F, Expr** vars, int nv) {
    if (F->is_true) return expr_new_symbol(SYM_True);
    if (F->n == 0)  return expr_new_symbol(SYM_False);

    /* Effective variables: those that appear in some atom. */
    int* used = malloc((size_t)nv * sizeof(int)); int nu = 0;
    for (int i = 0; i < nv; i++) {
        const char* nm = vars[i]->data.symbol.name;
        bool appears = false;
        for (int c = 0; c < F->n && !appears; c++)
            for (int k = 0; k < F->c[c]->n && !appears; k++)
                if (contains_symbol(F->c[c]->a[k].poly, nm)) appears = true;
        if (appears) used[nu++] = i;
    }
    if (nu == 1) {                            /* really univariate -> sign diagram */
        /* CAD does not thread Reduce options; the delegated 1-D sign diagram
         * keeps Solve's defaults (Cubics / Quartics forwarding on the
         * multivariate-real path is a scoped future refinement). */
        Expr* r = reduce_univar(F, vars[used[0]], vars, nv, NULL);
        free(used);
        return r;
    }
    if (nu >= 3) {                            /* Phase 6d: n-variable recursive CAD */
        Expr** vv = malloc((size_t)nu * sizeof(Expr*));
        for (int i = 0; i < nu; i++) vv[i] = vars[used[i]];
        free(used);
        Expr* r = reduce_cad_nvar(F, vv, nu);
        free(vv);
        return r;
    }

    Expr* vx = vars[used[0]];
    Expr* vy = vars[used[1]];
    free(used);

    Expr* out = NULL;
    Expr** B = NULL;   int nb = 0,  bcap = 0;
    Expr** px = NULL;  int npx = 0, pcap = 0;
    Expr** bxr = NULL; int nbx = 0, bxcap = 0; int mx = 0;
    Expr** cache = NULL; char* cache_bad = NULL;
    Expr** xparts = NULL; int nxp = 0, xcap = 0;
    YRegion* R = NULL;    /* interval-cell regions R[0..mx] */
    YRegion* S = NULL;    /* section-cell regions S[0..mx-1] */
    bool* absorbed = NULL;

    /* Gate + basis: every atom is a bivariate polynomial relation; factor into
     * the distinct-irreducible squarefree basis B. */
    for (int c = 0; c < F->n; c++) {
        RConj* cj = F->c[c];
        for (int k = 0; k < cj->n; k++) {
            RAtom* a = &cj->a[k];
            if (a->rel == R_ELEM || a->nonconst_denom) goto done;
            if (!is_poly2(a->poly, vx, vy)) goto done;
            if (!add_factors(a->poly, &B, &nb, &bcap)) goto done;
        }
    }
    if (nb == 0) goto done;

    /* Projection (eliminate vy): leading coefficient + discriminant per factor,
     * pairwise resultant, all factored; keep the distinct x-factors. */
    for (int i = 0; i < nb; i++) {
        int dy = degree_in(B[i], vy);
        if (dy <= 0) { if (!add_proj(B[i], &px, &npx, &pcap)) goto done; continue; }
        Expr* disc = eval_and_free(mkfun2("Discriminant", expr_copy(B[i]), expr_copy(vy)));
        bool okd = add_proj(disc, &px, &npx, &pcap); expr_free(disc); if (!okd) goto done;
        Expr* lc = eval_and_free(mkfun3(SYM_Coefficient, expr_copy(B[i]), expr_copy(vy), expr_new_integer(dy)));
        bool okl = add_proj(lc, &px, &npx, &pcap); expr_free(lc); if (!okl) goto done;
    }
    for (int i = 0; i < nb; i++) {
        if (degree_in(B[i], vy) <= 0) continue;
        for (int j = i + 1; j < nb; j++) {
            if (degree_in(B[j], vy) <= 0) continue;
            Expr* r = eval_and_free(mkfun3(SYM_Resultant, expr_copy(B[i]), expr_copy(B[j]), expr_copy(vy)));
            bool okr = add_proj(r, &px, &npx, &pcap); expr_free(r); if (!okr) goto done;
        }
    }

    /* Base decomposition: real breakpoints of the projection in vx. */
    for (int i = 0; i < npx; i++)
        if (!rru_collect_roots(px[i], vx, &bxr, &nbx, &bxcap, NULL, 0, NULL)) {
            for (int q = 0; q < nbx; q++) { expr_free(bxr[q]); } free(bxr); bxr = NULL; nbx = 0;
            goto done;
        }
    if (!order_dedup(bxr, NULL, nbx, &mx)) {
        for (int q = 0; q < nbx; q++) expr_free(bxr[q]);
        free(bxr); bxr = NULL;
        goto done;
    }
    /* v1: every section sample must be rational (irrational-coefficient fibre
     * isolation is deferred).  Interval samples are rational by construction. */
    for (int i = 0; i < mx; i++) if (!is_rational_number(bxr[i])) goto done;

    cache = nb ? calloc((size_t)nb, sizeof(Expr*)) : NULL;
    cache_bad = nb ? calloc((size_t)nb, sizeof(char)) : NULL;

    /* Lift the base decomposition: interval regions R[0..mx] (symbolic bounds)
     * and point-section regions S[0..mx-1] (numeric bounds).  A cell falsified
     * by a purely-outer atom (partial CAD) or with an empty fibre stays empty. */
    R = calloc((size_t)(mx + 1), sizeof(YRegion));
    S = mx ? calloc((size_t)mx, sizeof(YRegion)) : NULL;
    absorbed = mx ? calloc((size_t)mx, sizeof(bool)) : NULL;
    for (int j = 0; j <= mx; j++) {
        const Expr* xlo = (j == 0) ? NULL : bxr[j - 1];
        const Expr* xhi = (j == mx) ? NULL : bxr[j];
        Expr* sx = rru_rational_between(xlo, xhi);
        if (!sx) goto done;
        if (!xcell_dead(F, vx, sx, vy) && !lift_fiber(F, vx, sx, false, vy, B, nb, cache, cache_bad, &R[j])) { expr_free(sx); goto done; }
        expr_free(sx);
    }
    for (int k = 0; k < mx; k++)
        if (!xcell_dead(F, vx, bxr[k], vy) && !lift_fiber(F, vx, bxr[k], true, vy, B, nb, cache, cache_bad, &S[k])) goto done;

    /* Every interval AND section fibre fully satisfied => the whole plane. */
    {
        bool all_true = true;
        for (int j = 0; j <= mx && all_true; j++) if (!R[j].all_true) all_true = false;
        for (int k = 0; k < mx && all_true; k++) if (!S[k].all_true) all_true = false;
        if (all_true) { out = expr_new_symbol(SYM_True); goto done; }
    }

    /* Merge maximal runs of consecutive non-empty interval cells that share one
     * symbolic template, absorbing an interior/boundary breakpoint whenever the
     * template's limit there equals that section's own fibre (closed x-range). */
    for (int j = 0; j <= mx; ) {
        if (yregion_empty(&R[j])) { j++; continue; }
        int hi = j;
        while (hi < mx && !yregion_empty(&R[hi + 1])
               && templates_equal(&R[j], &R[hi + 1])
               && breakpoint_absorbable(&R[j], vx, bxr[hi], &S[hi])) {
            absorbed[hi] = true; hi++;
        }
        bool left_closed  = (j > 0)    && !absorbed[j - 1] && breakpoint_absorbable(&R[j], vx, bxr[j - 1], &S[j - 1]);
        bool right_closed = (hi < mx)  && breakpoint_absorbable(&R[j], vx, bxr[hi], &S[hi]);
        if (left_closed)  absorbed[j - 1] = true;
        if (right_closed) absorbed[hi]    = true;

        const Expr* xlo = (j == 0)  ? NULL : bxr[j - 1];
        const Expr* xhi = (hi == mx)? NULL : bxr[hi];
        Expr* xcond = cad_seg(xlo, !left_closed, xhi, !right_closed, vx);
        Expr* yexpr = yregion_to_expr(&R[j], vy);
        if (is_sym(yexpr, SYM_True)) { expr_free(yexpr); arr_push(&xparts, &nxp, &xcap, xcond); }
        else arr_push(&xparts, &nxp, &xcap, expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ xcond, yexpr }, 2));
        j = hi + 1;
    }
    /* Standalone sections: breakpoints no run absorbed. */
    for (int k = 0; k < mx; k++) {
        if (absorbed[k] || yregion_empty(&S[k])) continue;
        Expr* xeq = mkfun2(SYM_Equal, expr_copy(vx), expr_copy(bxr[k]));
        Expr* yexpr = yregion_to_expr(&S[k], vy);
        if (is_sym(yexpr, SYM_True)) { expr_free(yexpr); arr_push(&xparts, &nxp, &xcap, xeq); }
        else arr_push(&xparts, &nxp, &xcap, expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ xeq, yexpr }, 2));
    }

    if (nxp == 0) out = expr_new_symbol(SYM_False);
    else if (nxp == 1) { out = xparts[0]; nxp = 0; }
    else { out = expr_new_function(expr_new_symbol(SYM_Or), xparts, (size_t)nxp); nxp = 0; }
    out = eval_and_free(out);   /* flatten And/Or, fuse Inequality chains */

done:
    for (int i = 0; i < nxp; i++) { expr_free(xparts[i]); } free(xparts);
    if (R) { for (int i = 0; i <= mx; i++) { yregion_free(&R[i]); } free(R); }
    if (S) { for (int i = 0; i < mx; i++) { yregion_free(&S[i]); } free(S); }
    free(absorbed);
    if (cache) { for (int i = 0; i < nb; i++) { if (cache[i]) expr_free(cache[i]); } free(cache); }
    free(cache_bad);
    for (int i = 0; i < mx; i++) { expr_free(bxr[i]); } free(bxr);
    for (int i = 0; i < npx; i++) { expr_free(px[i]); } free(px);
    for (int i = 0; i < nb; i++) { expr_free(B[i]); } free(B);
    return out;
}
