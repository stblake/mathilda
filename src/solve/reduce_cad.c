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
        if (!rru_collect_roots(fs, vy, &roots, &nr, &cap, &fac, i)) { expr_free(fs); fail = true; break; }
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
        Expr* r = reduce_univar(F, vars[used[0]], vars, nv);
        free(used);
        return r;
    }
    if (nu != 2) { free(used); return NULL; } /* nv>=3 effective: deferred (Phase 6d) */

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
        if (!rru_collect_roots(px[i], vx, &bxr, &nbx, &bxcap, NULL, 0)) {
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
