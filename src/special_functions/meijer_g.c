/*
 * meijer_g.c -- MeijerG, the Meijer G-function, for Mathilda.
 *
 *   MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z]
 *   MeijerG[{{a1..an},{an+1..ap}}, {{b1..bm},{bm+1..bq}}, z, r]
 *
 * defined by the Mellin-Barnes integral
 *
 *            r    /  prod_{j<=m} G(b_j+rs) prod_{j<=n} G(1-a_j-rs)
 *   G  =  ------  |  ------------------------------------------------ z^-s ds
 *         2 pi i  /  prod_{j>n} G(a_j+rs) prod_{j>m} G(1-b_j-rs)
 *
 * (G = Gamma; r = 1 is the ordinary G-function). Substituting t = r s shows
 * MeijerG[a, b, z, r] = MeijerG[a, b, w] with w = z^(1/r) on the principal
 * branch, which is how every path below treats r.
 *
 * Evaluation order in builtin_meijer_g:
 *
 *   1. validate: 3 or 4 arguments (else hdiv + argt), each parameter argument a
 *      pair of Lists (else hdiv); r == 1 drops to the 3-argument form;
 *   2. buffer z (NDArray / packed List): the machine kernel over the buffer;
 *      a plain List z threads element-wise (MeijerG is not Listable: its first
 *      two arguments ARE lists);
 *   3. canonical order: each of the four parameter groups sorted (the function
 *      is symmetric within a group), as the reference system does;
 *   3b. Gamma-pair cancellation: equal parameters in opposite groups drop out
 *      (mg_cancel_pairs); hdiv when that leaves m = n = 0;
 *   4. consistency: m = n = 0, or some a_k - b_j (k <= n, j <= m) a positive
 *      integer (the contour cannot separate the two pole families): hdiv;
 *   5. machine-precision call: sf_mg_eval (sf_machine_meijerg.c), the helper
 *      shared with the buffer and Compile[] surfaces, before any exact rewrite;
 *   6. arbitrary precision (any MPFR leaf, all leaves numeric): the arb Slater
 *      engine (flint_num_meijerg), at the minimum input precision;
 *   7. z = 0: 0, a finite value, or ComplexInfinity where the expansion decides;
 *   8. automatic reduction to simpler functions (see mg_try_reduce);
 *   9. otherwise unevaluated.
 *
 * Memory contract: the builtin takes ownership of res but must not free it.
 */

#include "meijer_g.h"
#include "sf_machine_meijerg.h"
#include "sf_scalar.h"          /* sf_machine_leaf, sf_machine_complex_result */
#include "flint_num_bridge.h"
#include "ndarray.h"
#include "ndarray_internal.h"   /* nd_parallel_for */
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "sym_names.h"
#include "arithmetic.h"         /* is_complex */
#include "message.h"
#include "common.h"             /* builtin_arg_error */
#include "print.h"              /* expr_to_string */
#include "loadmodule.h"         /* mathilda_load_module */
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/*  Expression builders (every argument is consumed)                   */
/* ------------------------------------------------------------------ */

static Expr* mk(const char* h, size_t n, Expr** a) {
    return expr_new_function(expr_new_symbol(h), a, n);
}
static Expr* mk1(const char* h, Expr* x) { Expr* a[1] = { x }; return mk(h, 1, a); }
static Expr* mk2(const char* h, Expr* x, Expr* y) { Expr* a[2] = { x, y }; return mk(h, 2, a); }
static Expr* mk3(const char* h, Expr* x, Expr* y, Expr* z) {
    Expr* a[3] = { x, y, z }; return mk(h, 3, a);
}
static Expr* I_(long v) { return expr_new_integer(v); }
static Expr* half(void) { return mk2(SYM_Rational, I_(1), I_(2)); }
static Expr* add(Expr* a, Expr* b) { return mk2(SYM_Plus, a, b); }
static Expr* mul(Expr* a, Expr* b) { return mk2(SYM_Times, a, b); }
static Expr* neg(Expr* a) { return mul(I_(-1), a); }
static Expr* sub(Expr* a, Expr* b) { return add(a, neg(b)); }
static Expr* pw(Expr* a, Expr* b) { return mk2(SYM_Power, a, b); }
static Expr* dvd(Expr* a, Expr* b) { return mul(a, pw(b, I_(-1))); }
static Expr* gam(Expr* a) { return mk1(SYM_Gamma, a); }
static Expr* ev(Expr* e) { return eval_and_free(e); }
static Expr* cp(const Expr* e) { return expr_copy((Expr*)e); }

static Expr* mk_list(Expr** items, size_t n) {
    return expr_new_function(expr_new_symbol(SYM_List), items, n);
}

static bool is_head(const Expr* e, const char* h) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == h;
}
static bool is_list(const Expr* e) { return is_head(e, SYM_List); }

static bool is_int_value(const Expr* e, long v) {
    return e && e->type == EXPR_INTEGER && e->data.integer == v;
}

/* An exact numeric leaf: Integer, BigInt, Rational, or Complex of those. */
static bool is_exact_number(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (is_head(e, SYM_Rational)) return true;
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2)
        return is_exact_number(e->data.function.args[0]) &&
               is_exact_number(e->data.function.args[1]);
    return false;
}

static bool is_inexact_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* A numeric leaf of any kind (exact or inexact, real or Complex). */
static bool is_number(const Expr* e) {
    if (!e) return false;
    if (is_exact_number(e) || is_inexact_real(e)) return true;
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2)
        return is_number(e->data.function.args[0]) && is_number(e->data.function.args[1]);
    return false;
}

static bool is_exact_zero(const Expr* e) { return is_int_value(e, 0); }

static bool is_any_zero(const Expr* e) {
    if (is_exact_zero(e)) return true;
    if (e && e->type == EXPR_REAL) return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    return false;
}

/* The value of a REAL numeric leaf as a double (exact or inexact). */
static bool real_value(const Expr* e, double* out) {
    bool ix = false;
    if (sf_machine_leaf(e, out, &ix)) return true;
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    return false;
}

/* Integer-valued: an Integer / BigInt, or an inexact real with an integral
 * value. *v receives the value when it fits a long. */
static bool integer_valued(const Expr* e, long* v) {
    if (e->type == EXPR_INTEGER) { *v = (long)e->data.integer; return true; }
    if (e->type == EXPR_BIGINT) { *v = LONG_MAX; return true; }
    if (is_inexact_real(e)) {
        double d;
        if (!real_value(e, &d) || d != nearbyint(d)) return false;
        *v = (fabs(d) < 1e15) ? (long)d : LONG_MAX;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  The parsed call                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    size_t n, p, m, q;
    const Expr** a;          /* a_1 .. a_p (borrowed) */
    const Expr** b;          /* b_1 .. b_q (borrowed) */
    const Expr* z;
    const Expr* r;           /* NULL: the 3-argument form */
} MgView;

static void mg_view_free(MgView* g) { free(g->a); free(g->b); g->a = g->b = NULL; }

/* {{x..},{y..}} -> the two inner lists, or false. */
static bool pair_of_lists(const Expr* e, const Expr** l1, const Expr** l2) {
    if (!is_list(e) || e->data.function.arg_count != 2) return false;
    *l1 = e->data.function.args[0];
    *l2 = e->data.function.args[1];
    return is_list(*l1) && is_list(*l2);
}

static bool mg_view(const Expr* g, MgView* v) {
    memset(v, 0, sizeof *v);
    if (g->type != EXPR_FUNCTION) return false;
    size_t na = g->data.function.arg_count;
    if (na != 3 && na != 4) return false;
    const Expr *a1, *a2, *b1, *b2;
    if (!pair_of_lists(g->data.function.args[0], &a1, &a2)) return false;
    if (!pair_of_lists(g->data.function.args[1], &b1, &b2)) return false;
    v->n = a1->data.function.arg_count;
    v->p = v->n + a2->data.function.arg_count;
    v->m = b1->data.function.arg_count;
    v->q = v->m + b2->data.function.arg_count;
    v->a = malloc(sizeof(Expr*) * (v->p ? v->p : 1));
    v->b = malloc(sizeof(Expr*) * (v->q ? v->q : 1));
    if (!v->a || !v->b) { mg_view_free(v); return false; }
    for (size_t j = 0; j < v->n; j++) v->a[j] = a1->data.function.args[j];
    for (size_t j = v->n; j < v->p; j++) v->a[j] = a2->data.function.args[j - v->n];
    for (size_t j = 0; j < v->m; j++) v->b[j] = b1->data.function.args[j];
    for (size_t j = v->m; j < v->q; j++) v->b[j] = b2->data.function.args[j - v->m];
    v->z = g->data.function.args[2];
    v->r = (na == 4) ? g->data.function.args[3] : NULL;
    return true;
}

/* MeijerG[{{a..},{a..}},{{b..},{b..}}, z (, r)] from a view, with z replaced by
 * `z` (consumed). Parameters are copied. */
static Expr* mg_rebuild(const MgView* g, Expr* z) {
    Expr** a1 = malloc(sizeof(Expr*) * (g->n ? g->n : 1));
    Expr** a2 = malloc(sizeof(Expr*) * (g->p - g->n ? g->p - g->n : 1));
    Expr** b1 = malloc(sizeof(Expr*) * (g->m ? g->m : 1));
    Expr** b2 = malloc(sizeof(Expr*) * (g->q - g->m ? g->q - g->m : 1));
    for (size_t j = 0; j < g->n; j++) a1[j] = cp(g->a[j]);
    for (size_t j = g->n; j < g->p; j++) a2[j - g->n] = cp(g->a[j]);
    for (size_t j = 0; j < g->m; j++) b1[j] = cp(g->b[j]);
    for (size_t j = g->m; j < g->q; j++) b2[j - g->m] = cp(g->b[j]);
    Expr* al = mk2(SYM_List, mk_list(a1, g->n), mk_list(a2, g->p - g->n));
    Expr* bl = mk2(SYM_List, mk_list(b1, g->m), mk_list(b2, g->q - g->m));
    free(a1); free(a2); free(b1); free(b2);
    if (g->r) return mk(SYM_MeijerG, 4, (Expr*[]){ al, bl, z, cp(g->r) });
    return mk3(SYM_MeijerG, al, bl, z);
}

/* ------------------------------------------------------------------ */
/*  Messages                                                           */
/* ------------------------------------------------------------------ */

static void msg_hdiv(const Expr* call) {
    char* s = expr_to_string((Expr*)call);
    mth_message("MeijerG", "hdiv", "%s does not exist. Arguments are not consistent.",
                s ? s : "MeijerG[...]");
    free(s);
}

/* ------------------------------------------------------------------ */
/*  Canonical order and consistency                                    */
/* ------------------------------------------------------------------ */

static int cmp_expr(const void* x, const void* y) {
    return expr_compare(*(const Expr* const*)x, *(const Expr* const*)y);
}

static bool group_sorted(const Expr** v, size_t lo, size_t hi) {
    for (size_t j = lo + 1; j < hi; j++)
        if (expr_compare(v[j - 1], v[j]) > 0) return false;
    return true;
}

/* NULL when already canonical; otherwise the sorted call (owned). */
static Expr* mg_canonical(const MgView* g) {
    if (group_sorted(g->a, 0, g->n) && group_sorted(g->a, g->n, g->p) &&
        group_sorted(g->b, 0, g->m) && group_sorted(g->b, g->m, g->q)) return NULL;
    MgView s = *g;
    s.a = malloc(sizeof(Expr*) * (g->p ? g->p : 1));
    s.b = malloc(sizeof(Expr*) * (g->q ? g->q : 1));
    memcpy(s.a, g->a, sizeof(Expr*) * g->p);
    memcpy(s.b, g->b, sizeof(Expr*) * g->q);
    qsort(s.a, g->n, sizeof(Expr*), cmp_expr);
    qsort(s.a + g->n, g->p - g->n, sizeof(Expr*), cmp_expr);
    qsort(s.b, g->m, sizeof(Expr*), cmp_expr);
    qsort(s.b + g->m, g->q - g->m, sizeof(Expr*), cmp_expr);
    Expr* out = mg_rebuild(&s, cp(g->z));
    free(s.a); free(s.b);
    return out;
}

/* Evaluated x - y (copies). */
static Expr* diff_of(const Expr* x, const Expr* y) { return ev(sub(cp(x), cp(y))); }

/* x - y an integer >= lo (exact, or an integer-valued inexact difference). */
static bool diff_int_at_least(const Expr* x, const Expr* y, long lo, long* v) {
    Expr* d = diff_of(x, y);
    bool hit = integer_valued(d, v) && *v >= lo;
    expr_free(d);
    return hit;
}

/* Some a_k - b_j (k <= n, j <= m) a positive integer: a pole of Gamma(b_j+s)
 * coincides with one of Gamma(1-a_k-s), and no contour separates them. */
static size_t mg_collisions(const MgView* g) {
    long v;
    size_t c = 0;
    for (size_t k = 0; k < g->n; k++)
        for (size_t j = 0; j < g->m; j++)
            if (diff_int_at_least(g->a[k], g->b[j], 1, &v)) c++;
    return c;
}
static bool mg_collision(const MgView* g) { return mg_collisions(g) > 0; }

/* A colliding parameter set can still be consistent when a numerator Gamma is
 * cancelled by a denominator one: for b_j - a_k = N a non-negative integer,
 *
 *   Gamma(b_j+s) / Gamma(a_k+s) = (-1)^N Gamma(1-a_k-s) / Gamma(1-b_j-s)
 *
 * (both are the polynomial (a_k+s)_N), so with j <= m, k > n the integrand --
 * and G -- is unchanged by moving b_j to the second b group and a_k to the
 * first a group, times (-1)^N; symmetrically for a_k (k <= n), b_j (j > m)
 * with b_j - a_k = N >= 0. Mathematica applies exactly this move:
 * MeijerG[{{31,4},{2,3}},{{5,3/2},{1/10}},3] is
 * -MeijerG[{{2,4,31},{3}},{{3/2},{1/10,5}},3]. Returns the rewritten call
 * (owned) after one move that strictly reduces the number of colliding pairs
 * (so moves can never cycle), or NULL when none does. */
static Expr* mg_move_pair(const MgView* g) {
    long N;
    size_t c0 = mg_collisions(g);   /* accept only a move that strictly reduces it */
    for (size_t j = 0; j < g->m; j++)
        for (size_t k = g->n; k < g->p; k++) {
            if (!diff_int_at_least(g->b[j], g->a[k], 0, &N)) continue;
            MgView s = *g;
            s.a = malloc(sizeof(Expr*) * (g->p ? g->p : 1));
            s.b = malloc(sizeof(Expr*) * (g->q ? g->q : 1));
            /* a_k joins the first a group (index n), b_j the second b group. */
            size_t t = 0;
            for (size_t i = 0; i < g->n; i++) s.a[t++] = g->a[i];
            s.a[t++] = g->a[k];
            for (size_t i = g->n; i < g->p; i++) if (i != k) s.a[t++] = g->a[i];
            t = 0;
            for (size_t i = 0; i < g->m; i++) if (i != j) s.b[t++] = g->b[i];
            s.b[t++] = g->b[j];
            for (size_t i = g->m; i < g->q; i++) s.b[t++] = g->b[i];
            s.n = g->n + 1; s.m = g->m - 1;
            if (mg_collisions(&s) >= c0) { free(s.a); free(s.b); continue; }
            Expr* out = mg_rebuild(&s, cp(g->z));
            free(s.a); free(s.b);
            return (N % 2) ? neg(out) : out;
        }
    for (size_t k = 0; k < g->n; k++)
        for (size_t j = g->m; j < g->q; j++) {
            if (!diff_int_at_least(g->b[j], g->a[k], 0, &N)) continue;
            MgView s = *g;
            s.a = malloc(sizeof(Expr*) * (g->p ? g->p : 1));
            s.b = malloc(sizeof(Expr*) * (g->q ? g->q : 1));
            size_t t = 0;
            for (size_t i = 0; i < g->n; i++) if (i != k) s.a[t++] = g->a[i];
            s.a[t++] = g->a[k];
            for (size_t i = g->n; i < g->p; i++) s.a[t++] = g->a[i];
            t = 0;
            for (size_t i = 0; i < g->m; i++) s.b[t++] = g->b[i];
            s.b[t++] = g->b[j];
            for (size_t i = g->m; i < g->q; i++) if (i != j) s.b[t++] = g->b[i];
            s.n = g->n - 1; s.m = g->m + 1;
            if (mg_collisions(&s) >= c0) { free(s.a); free(s.b); continue; }
            Expr* out = mg_rebuild(&s, cp(g->z));
            free(s.a); free(s.b);
            return (N % 2) ? neg(out) : out;
        }
    return NULL;
}

/* Gamma-pair cancellation. An a_k of the second a group (k > n) equal to a
 * b_j of the first b group (j <= m) contributes Gamma(b_j+s)/Gamma(a_k+s) = 1
 * to the Mellin-Barnes integrand, and an a_k of the first group (k <= n)
 * equal to a b_j of the second (j > m) contributes
 * Gamma(1-a_k-s)/Gamma(1-b_j-s) = 1: both pairs drop out, lowering
 * (m, p) or (n, q) by one, with the contour unaffected (one pole family
 * fewer to separate). MeijerG[{{}, {0}}, {{0, 0}, {}}, z] is G^{1,0}_{0,1}
 * = E^-z this way. Exact (structural) equality only.
 *
 * Returns the rewritten call (owned, all equal pairs removed at once), or
 * NULL when no pair cancels. A cancellation that would leave m = n = 0 has
 * no contour at all -- the call does not exist -- and sets *hdiv instead. */
static Expr* mg_cancel_pairs(const MgView* g, bool* hdiv) {
    *hdiv = false;
    bool* adrop = calloc(g->p ? g->p : 1, sizeof(bool));
    bool* bdrop = calloc(g->q ? g->q : 1, sizeof(bool));
    size_t dm = 0, dp = 0, dn = 0, dq = 0;
    for (size_t k = g->n; k < g->p; k++)
        for (size_t j = 0; j < g->m; j++)
            if (!bdrop[j] && expr_eq((Expr*)g->a[k], (Expr*)g->b[j])) {
                adrop[k] = bdrop[j] = true; dp++; dm++; break;
            }
    for (size_t k = 0; k < g->n; k++)
        for (size_t j = g->m; j < g->q; j++)
            if (!bdrop[j] && expr_eq((Expr*)g->a[k], (Expr*)g->b[j])) {
                adrop[k] = bdrop[j] = true; dn++; dq++; break;
            }
    Expr* out = NULL;
    if (dm + dn > 0) {
        if (g->m - dm == 0 && g->n - dn == 0) {
            *hdiv = true;
        } else {
            MgView s = *g;
            s.a = malloc(sizeof(Expr*) * (g->p ? g->p : 1));
            s.b = malloc(sizeof(Expr*) * (g->q ? g->q : 1));
            size_t t = 0;
            for (size_t k = 0; k < g->p; k++) if (!adrop[k]) s.a[t++] = g->a[k];
            t = 0;
            for (size_t j = 0; j < g->q; j++) if (!bdrop[j]) s.b[t++] = g->b[j];
            s.n = g->n - dn; s.p = g->p - dn - dp;
            s.m = g->m - dm; s.q = g->q - dm - dq;
            out = mg_rebuild(&s, cp(g->z));
            free(s.a); free(s.b);
        }
    }
    free(adrop); free(bdrop);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Machine precision                                                  */
/* ------------------------------------------------------------------ */

/* One leaf as a machine (re, im) pair; *inexact OR-ed. */
static bool mg_leaf(const Expr* e, double* re, double* im, bool* inexact) {
    bool ix = false;
    if (sf_machine_leaf(e, re, &ix)) { *im = 0.0; *inexact |= ix; return true; }
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2) {
        bool ia = false, ib = false;
        if (!sf_machine_leaf(e->data.function.args[0], re, &ia) ||
            !sf_machine_leaf(e->data.function.args[1], im, &ib)) return false;
        *inexact |= (ia || ib);
        return true;
    }
    return false;
}

/* The parameters (and r) as doubles. *inexact reports a machine Real among
 * them. False when p or q is too large or a leaf is not a machine scalar. */
static bool mg_machine_params(const MgView* g, double* ar, double* ai, double* br,
                              double* bi, double* r, bool* inexact) {
    if (g->p > SF_MG_MAXP || g->q > SF_MG_MAXP) return false;
    for (size_t j = 0; j < g->p; j++) if (!mg_leaf(g->a[j], &ar[j], &ai[j], inexact)) return false;
    for (size_t j = 0; j < g->q; j++) if (!mg_leaf(g->b[j], &br[j], &bi[j], inexact)) return false;
    *r = 1.0;
    if (g->r) {
        bool ix = false;
        if (!sf_machine_leaf(g->r, r, &ix)) return false;
        *inexact |= ix;
    }
    return true;
}

static Expr* mg_machine_result(double rr, double ri, bool is_c) {
    return is_c ? sf_machine_complex_result(rr, ri) : expr_new_real(rr);
}

/* The scalar machine gate: every leaf a machine scalar, at least one inexact. */
static Expr* mg_machine(const MgView* g) {
    double ar[SF_MG_MAXP], ai[SF_MG_MAXP], br[SF_MG_MAXP], bi[SF_MG_MAXP], r;
    bool inexact = false;
    if (!mg_machine_params(g, ar, ai, br, bi, &r, &inexact)) return NULL;
    double zr, zi;
    if (!mg_leaf(g->z, &zr, &zi, &inexact)) return NULL;
    if (!inexact) return NULL;                  /* exact call stays exact */
    /* A log-case plan (residue table) is memoised, so repeated scalar calls
     * with the same parameters pay its build once. */
    SfMgPlan scratch;
    const SfMgPlan* P = sf_mg_plan_get(ar, ai, g->n, g->p, br, bi, g->m, g->q, r,
                                       false, &scratch);
    if (!P) return NULL;
    double rr, ri; bool is_c;
    bool ok = sf_mg_eval(P, zr, zi, &rr, &ri, &is_c);
    if (P == &scratch) sf_mg_plan_free(&scratch);
    if (!ok) return NULL;
    return mg_machine_result(rr, ri, is_c);
}

/* Buffer path: z an NDArray / packed List. */
typedef struct {
    const SfMgPlan* P; const void* in; NDType dt; double* out; int any_c;
} MgNdCtx;

static bool mg_nd_chunk(void* vctx, size_t lo, size_t hi) {
    MgNdCtx* c = (MgNdCtx*)vctx;
    for (size_t i = lo; i < hi; i++) {
        double xr, xi, rr, ri; bool is_c;
        ndt_get(c->in, i, c->dt, &xr, &xi);
        if (!sf_mg_eval(c->P, xr, xi, &rr, &ri, &is_c)) return false;
        c->out[2 * i] = rr; c->out[2 * i + 1] = is_c ? ri : 0.0;
        if (is_c && ri != 0.0) c->any_c = 1;   /* benign race: only ever set to 1 */
    }
    return true;
}

static Expr* mg_ndarray(const MgView* g) {
    const Expr* z = g->z;
    NDType zdt = z->data.ndarray.dtype;
    if (zdt == NDT_INT64 || zdt == NDT_BOOL) return NULL;   /* exact data: element-wise */
    double ar[SF_MG_MAXP], ai[SF_MG_MAXP], br[SF_MG_MAXP], bi[SF_MG_MAXP], r;
    bool inexact = true;     /* a machine buffer makes every element call inexact */
    if (!mg_machine_params(g, ar, ai, br, bi, &r, &inexact)) return NULL;
    int any_c = 0;
    SfMgPlan* scratch = malloc(sizeof *scratch);
    if (!scratch) return NULL;
    const SfMgPlan* P = sf_mg_plan_get(ar, ai, g->n, g->p, br, bi, g->m, g->q, r,
                                       false, scratch);
    if (!P) { free(scratch); return NULL; }
    size_t sz = ndarray_size(z);
    double* buf = malloc(sizeof(double) * 2 * (sz ? sz : 1));
    bool ok = false;
    if (buf) {
        MgNdCtx ctx = { P, z->data.ndarray.data, zdt, buf, ndt_is_complex(zdt) ? 1 : 0 };
        ok = nd_parallel_for(sz, mg_nd_chunk, &ctx);
        any_c = ctx.any_c;
    }
    if (P == scratch) sf_mg_plan_free(scratch);
    free(scratch);
    if (!buf) return NULL;
    if (!ok) { free(buf); return NULL; }
    NDType odt = any_c ? ((zdt == NDT_FLOAT32 || zdt == NDT_COMPLEX32) ? NDT_COMPLEX32
                                                                           : NDT_COMPLEX64)
                           : zdt;
    void* out = malloc(ndt_elem_size(odt) * (sz ? sz : 1));
    if (!out) { free(buf); return NULL; }
    for (size_t i = 0; i < sz; i++) ndt_set(out, i, odt, buf[2 * i], buf[2 * i + 1]);
    free(buf);
    return expr_new_ndarray_like(z, z->data.ndarray.rank, z->data.ndarray.dims, out, odt);
}

/* ------------------------------------------------------------------ */
/*  Arbitrary precision                                                */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
static Expr* demote_to_machine(Expr* e) {
    double d;
    if (e->type == EXPR_MPFR) { d = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return expr_new_real(d); }
    if (e->type == EXPR_REAL) return expr_copy(e);
    Expr *re = NULL, *im = NULL;
    if (is_complex(e, &re, &im)) {
        double a, b;
        if (real_value(re, &a) && real_value(im, &b)) return sf_machine_complex_result(a, b);
    }
    return NULL;
}
#endif

static bool leaf_has_machine_real(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2)
        return leaf_has_machine_real(e->data.function.args[0]) ||
               leaf_has_machine_real(e->data.function.args[1]);
    return false;
}

static bool leaf_inexact(const Expr* e) {
    if (is_inexact_real(e)) return true;
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2)
        return leaf_inexact(e->data.function.args[0]) || leaf_inexact(e->data.function.args[1]);
    return false;
}

static bool mg_all_numeric_inexact(const MgView* g) {
    bool inexact = false;
    for (size_t j = 0; j < g->p; j++) { if (!is_number(g->a[j])) return false; inexact |= leaf_inexact(g->a[j]); }
    for (size_t j = 0; j < g->q; j++) { if (!is_number(g->b[j])) return false; inexact |= leaf_inexact(g->b[j]); }
    if (!is_number(g->z) || (g->r && !is_number(g->r))) return false;
    inexact |= leaf_inexact(g->z);
    if (g->r) inexact |= leaf_inexact(g->r);
    return inexact;
}

/* Every leaf numeric and at least one inexact: the arb engine. A machine Real
 * anywhere makes the answer machine precision (precision contagion). */
static Expr* mg_arbitrary(const MgView* g) {
    bool inexact = false, machine = false;
    for (size_t j = 0; j < g->p; j++) {
        if (!is_number(g->a[j])) return NULL;
        inexact |= leaf_inexact(g->a[j]); machine |= leaf_has_machine_real(g->a[j]);
    }
    for (size_t j = 0; j < g->q; j++) {
        if (!is_number(g->b[j])) return NULL;
        inexact |= leaf_inexact(g->b[j]); machine |= leaf_has_machine_real(g->b[j]);
    }
    if (!is_number(g->z) || (g->r && !is_number(g->r))) return NULL;
    inexact |= leaf_inexact(g->z); machine |= leaf_has_machine_real(g->z);
    if (g->r) { inexact |= leaf_inexact(g->r); machine |= leaf_has_machine_real(g->r); }
    if (!inexact || is_any_zero(g->z)) return NULL;
    Expr* v = flint_num_meijerg(g->a, g->n, g->p, g->b, g->m, g->q, g->z, g->r);
#ifdef USE_MPFR
    if (v && machine) { Expr* d = demote_to_machine(v); expr_free(v); v = d; }
#else
    (void)machine;
#endif
    return v;
}

/* ------------------------------------------------------------------ */
/*  Slater expansions as expressions                                   */
/* ------------------------------------------------------------------ */

/* One expansion of G^{m,n}_{p,q}: in t = w (z-series) or t = 1/w (dual),
 * w = z^(1/r). A and B are owned, evaluated parameter copies (for the dual,
 * A = 1 - b, B = 1 - a, with n,m and p,q swapped). */
typedef struct {
    size_t n, p, m, q;
    Expr** A; Expr** B;
    bool dual;
    const Expr* z; const Expr* r;
} MgExp;

static void mg_exp_free(MgExp* E) {
    if (E->A) for (size_t j = 0; j < E->p; j++) expr_free(E->A[j]);
    if (E->B) for (size_t j = 0; j < E->q; j++) expr_free(E->B[j]);
    free(E->A); free(E->B);
    E->A = E->B = NULL;
}

static void mg_exp_init(const MgView* g, bool dual, MgExp* E) {
    memset(E, 0, sizeof *E);
    E->dual = dual; E->z = g->z; E->r = g->r;
    if (!dual) {
        E->n = g->n; E->p = g->p; E->m = g->m; E->q = g->q;
        E->A = malloc(sizeof(Expr*) * (E->p ? E->p : 1));
        E->B = malloc(sizeof(Expr*) * (E->q ? E->q : 1));
        for (size_t j = 0; j < E->p; j++) E->A[j] = cp(g->a[j]);
        for (size_t j = 0; j < E->q; j++) E->B[j] = cp(g->b[j]);
    } else {
        E->n = g->m; E->p = g->q; E->m = g->n; E->q = g->p;
        E->A = malloc(sizeof(Expr*) * (E->p ? E->p : 1));
        E->B = malloc(sizeof(Expr*) * (E->q ? E->q : 1));
        for (size_t j = 0; j < E->p; j++) E->A[j] = ev(sub(I_(1), cp(g->b[j])));
        for (size_t j = 0; j < E->q; j++) E->B[j] = ev(sub(I_(1), cp(g->a[j])));
    }
}

/* t^e as a power of z: Power[z, +-e/r]. Consumes e. */
static Expr* tpow(const MgExp* E, Expr* e) {
    Expr* ex = E->dual ? neg(e) : e;
    if (E->r) ex = dvd(ex, cp(E->r));
    return pw(cp(E->z), ex);
}

/* x = (-1)^{p-m-n} t. */
static int mg_sign(const MgExp* E) {
    return (((long)E->p - (long)E->m - (long)E->n) % 2 != 0) ? -1 : 1;
}
static Expr* mg_x(const MgExp* E) {
    Expr* t = tpow(E, I_(1));
    return mg_sign(E) < 0 ? neg(t) : t;
}

/* An exact non-positive integer: -N, N >= 0. */
static bool exact_nonpos_int(const Expr* e, long* N) {
    if (e->type == EXPR_INTEGER && e->data.integer <= 0) { *N = (long)(-e->data.integer); return true; }
    return false;
}

/* Two of B_0..B_{m-1} provably an integer apart (exact difference). */
static bool mg_log_case(const MgExp* E) {
    for (size_t h = 0; h < E->m; h++)
        for (size_t j = h + 1; j < E->m; j++) {
            Expr* d = diff_of(E->B[j], E->B[h]);
            long v;
            bool hit = integer_valued(d, &v);
            expr_free(d);
            if (hit) return true;
        }
    return false;
}

/* Term pieces: Cnum, the t-exponent, the pFq parameters, and which lower
 * parameters are regularized (those from j > m). Owned by the caller. */
typedef struct {
    Expr* C;            /* numeric/symbolic coefficient */
    Expr* expo;         /* exponent of t */
    Expr** u; size_t pu;
    Expr** l; bool* reg; size_t ql;
    int sgn;            /* sign of x = sgn * t */
} MgTerm;

static void mg_term_free(MgTerm* T) {
    expr_free(T->C); expr_free(T->expo);
    for (size_t k = 0; k < T->pu; k++) expr_free(T->u[k]);
    for (size_t k = 0; k < T->ql; k++) expr_free(T->l[k]);
    free(T->u); free(T->l); free(T->reg);
    memset(T, 0, sizeof *T);
}

/* Build term h. Regularized lower parameters that are exact non-positive
 * integers are shifted out: R(u; -N, l'; x) = prod (u)_{N+1} x^{N+1}
 * R(u+N+1; N+2, l'+N+1; x). */
static void mg_term(const MgExp* E, size_t h, MgTerm* T) {
    memset(T, 0, sizeof *T);
    size_t p = E->p, q = E->q, m = E->m, n = E->n;
    Expr* C = I_(1);
    for (size_t j = 0; j < m; j++)
        if (j != h) C = mul(C, gam(sub(cp(E->B[j]), cp(E->B[h]))));
    for (size_t j = 0; j < n; j++)
        C = mul(C, gam(sub(add(I_(1), cp(E->B[h])), cp(E->A[j]))));
    for (size_t j = n; j < p; j++)
        C = dvd(C, gam(sub(cp(E->A[j]), cp(E->B[h]))));
    T->pu = p; T->ql = q ? q - 1 : 0;
    T->u = malloc(sizeof(Expr*) * (p ? p : 1));
    T->l = malloc(sizeof(Expr*) * (q ? q : 1));
    T->reg = malloc(sizeof(bool) * (q ? q : 1));
    for (size_t k = 0; k < p; k++) T->u[k] = ev(sub(add(I_(1), cp(E->B[h])), cp(E->A[k])));
    size_t nl = 0;
    for (size_t j = 0; j < q; j++) {
        if (j == h) continue;
        T->l[nl] = ev(sub(add(I_(1), cp(E->B[h])), cp(E->B[j])));
        T->reg[nl] = (j >= m);
        nl++;
    }
    T->sgn = mg_sign(E);
    Expr* expo = cp(E->B[h]);
    for (int guard = 0; guard < 64; guard++) {
        size_t hit = (size_t)-1; long N = 0;
        for (size_t k = 0; k < nl; k++)
            if (T->reg[k] && exact_nonpos_int(T->l[k], &N)) { hit = k; break; }
        if (hit == (size_t)-1) break;
        for (size_t k = 0; k < p; k++) {
            Expr* pf = mk2(SYM_Pochhammer, cp(T->u[k]), I_(N + 1));
            C = mul(C, pf);
            T->u[k] = ev(add(T->u[k], I_(N + 1)));
        }
        for (size_t k = 0; k < nl; k++) {
            if (k == hit) { expr_free(T->l[k]); T->l[k] = I_(N + 2); }
            else T->l[k] = ev(add(T->l[k], I_(N + 1)));
        }
        if (T->sgn < 0 && ((N + 1) % 2 != 0)) C = neg(C);
        expo = add(expo, I_(N + 1));
    }
    T->C = ev(C);
    T->expo = ev(expo);
}

/* The regularized pFq of the term, as HypergeometricPFQ[u, l, x] divided by
 * Gamma of the regularized lower parameters. */
static Expr* mg_pfq_expr(const MgTerm* T, Expr* x) {
    Expr** u = malloc(sizeof(Expr*) * (T->pu ? T->pu : 1));
    Expr** l = malloc(sizeof(Expr*) * (T->ql ? T->ql : 1));
    for (size_t k = 0; k < T->pu; k++) u[k] = cp(T->u[k]);
    for (size_t k = 0; k < T->ql; k++) l[k] = cp(T->l[k]);
    Expr* f = mk3(SYM_HypergeometricPFQ, mk_list(u, T->pu), mk_list(l, T->ql), x);
    free(u); free(l);
    for (size_t k = 0; k < T->ql; k++)
        if (T->reg[k]) f = dvd(f, gam(cp(T->l[k])));
    return f;
}

static bool eq_rational(const Expr* e, long num, long den) {
    Expr* t = (den == 1) ? I_(num) : ev(mk2(SYM_Rational, I_(num), I_(den)));
    bool r = expr_eq(e, t);
    expr_free(t);
    return r;
}

/* The term as an expression, `named` rewriting the pFq to a classical
 * function where the shape has one. */
static Expr* mg_term_expr(const MgExp* E, const MgTerm* T, bool named) {
    Expr* extra = NULL;          /* extra exponent of t (Bessel normalisation) */
    Expr* F = NULL;
    /* Gamma(l) for each non-regularized lower: the named forms below are the
     * REGULARIZED functions. */
    Expr* nonreg = I_(1);
    for (size_t k = 0; k < T->ql; k++)
        if (!T->reg[k]) nonreg = mul(nonreg, gam(cp(T->l[k])));
    if (named) {
        if (T->pu == 0 && T->ql == 0) {
            F = mk1(SYM_Exp, mg_x(E));
        } else if (T->pu == 1 && T->ql == 0) {
            F = pw(sub(I_(1), mg_x(E)), neg(cp(T->u[0])));
        } else if (T->pu == 0 && T->ql == 1) {
            /* 0F1~(;c;-+t) = t^{(1-c)/2} J/I_{c-1}(2 t^{1/2}). */
            const char* bes = (T->sgn < 0) ? SYM_BesselJ : SYM_BesselI;
            extra = mul(half(), sub(I_(1), cp(T->l[0])));
            F = mk2(bes, sub(cp(T->l[0]), I_(1)), mul(I_(2), tpow(E, half())));
        } else if (T->pu == 1 && T->ql == 1) {
            Expr* d = ev(sub(cp(T->l[0]), mul(I_(2), cp(T->u[0]))));
            bool kummer = is_exact_zero(d);
            expr_free(d);
            if (kummer) {
                /* 1F1~(u;2u;x) = e^{x/2} G(u+1/2)/G(2u) (x/4)^{1/2-u} I_{u-1/2}(x/2). */
                Expr* x = mg_x(E);
                Expr* uh = add(cp(T->u[0]), half());
                F = mul(mk1(SYM_Exp, mul(half(), cp(x))),
                        dvd(gam(cp(uh)), gam(mul(I_(2), cp(T->u[0])))));
                F = mul(F, pw(mul(mk2(SYM_Rational, I_(1), I_(4)), cp(x)),
                              sub(half(), cp(T->u[0]))));
                F = mul(F, mk2(SYM_BesselI, sub(cp(T->u[0]), half()), mul(half(), x)));
                expr_free(uh);
            }
        } else if (T->pu == 2 && T->ql == 1 && is_int_value(T->l[0], 1)) {
            bool h1 = eq_rational(T->u[0], 1, 2), h2 = eq_rational(T->u[1], 1, 2);
            bool m1 = eq_rational(T->u[0], -1, 2), m2 = eq_rational(T->u[1], -1, 2);
            Expr* two_pi = dvd(I_(2), expr_new_symbol(SYM_Pi));
            if (h1 && h2)                    F = mul(two_pi, mk1(SYM_EllipticK, mg_x(E)));
            else if ((h1 && m2) || (m1 && h2)) F = mul(two_pi, mk1(SYM_EllipticE, mg_x(E)));
            else expr_free(two_pi);
        }
    }
    Expr* expo = cp(T->expo);
    if (F) {
        F = mul(nonreg, F);
        if (extra) expo = add(expo, extra);
    } else {
        expr_free(nonreg);
        F = mg_pfq_expr(T, mg_x(E));
    }
    return mul(mul(cp(T->C), tpow(E, expo)), F);
}

/* Sum of all terms of expansion E (no log case). */
static Expr* mg_sum_expr(const MgExp* E, bool named) {
    if (E->m == 0) return I_(0);
    Expr** t = malloc(sizeof(Expr*) * E->m);
    for (size_t h = 0; h < E->m; h++) {
        MgTerm T;
        mg_term(E, h, &T);
        t[h] = mg_term_expr(E, &T, named);
        mg_term_free(&T);
    }
    Expr* s = (E->m == 1) ? t[0] : mk(SYM_Plus, E->m, t);
    free(t);
    return s;
}

/* Which expansion represents G everywhere (as an analytic function of z off
 * the cut): 0 = in w, 1 = in 1/w, -1 = none (p == q with m + n == p, where the
 * two expansions are not continuations of each other). For p == q the one with
 * fewer terms (ties: in w). */
/* -1, 0, +1 as |z| is below, on, or above 1 for a numeric z; 2 when z is
 * not a number (or the comparison does not decide). */
static int mg_abs_cmp_one(const Expr* z) {
    if (!is_number(z)) return 2;
    Expr* az = ev(mk1(SYM_Abs, cp(z)));
    int c = 2;
    if (is_int_value(az, 1)) c = 0;
    else {
        Expr* lt = ev(mk2(SYM_Less, cp(az), I_(1)));
        Expr* gt = ev(mk2(SYM_Greater, cp(az), I_(1)));
        if (lt->type == EXPR_SYMBOL && lt->data.symbol.name == SYM_True) c = -1;
        else if (gt->type == EXPR_SYMBOL && gt->data.symbol.name == SYM_True) c = 1;
        expr_free(lt); expr_free(gt);
    }
    expr_free(az);
    return c;
}

static int mg_pick(const MgView* g) {
    if (g->p < g->q) return 0;
    if (g->p > g->q) return 1;
    /* A numeric z picks the expansion that CONVERGES there, |z| < 1 in z and
     * |z| > 1 in 1/z: the other one needs every term's pF(p-1) continued past
     * its circle (and, for real z beyond it, onto a branch cut that only the
     * sum of the terms does not have). For m + n <= p the two expansions are
     * not continuations of each other but each IS G in its own disk (the
     * Mellin-Barnes contour closes left for |z| < 1 and right for |z| > 1),
     * so the convergent one is still the value: MeijerG[{{1},{1}},{{1/2},{0}},4]
     * = Sqrt[Pi]. On |z| = 1 that case picks neither (see mg_unit_ok). */
    int c = mg_abs_cmp_one(g->z);
    if (c == -1) return 0;
    if (c == 1) return 1;
    if (g->m + g->n <= g->p) return -1;
    if (c == 0) return 1;       /* callers gate |z| = 1 themselves */
    return (g->m <= g->n) ? 0 : 1;
}

/* Re(sum a - sum b) for exact numeric parameters; false otherwise. */
static bool mg_param_excess(const MgView* g, double* out) {
    Expr* s = I_(0);
    for (size_t j = 0; j < g->p; j++) s = add(s, cp(g->a[j]));
    for (size_t j = 0; j < g->q; j++) s = sub(s, cp(g->b[j]));
    s = ev(mk1(SYM_Re, s));
    bool ok = is_exact_number(s) && real_value(s, out);
    expr_free(s);
    return ok;
}

/* p == q and z exactly 1: may a Slater expansion be summed AT z = 1? Each
 * term is a pF(p-1) at x = (-1)^{p-m-n} (the dual expansion has the same
 * parity), whose parameter excess sum(lower) - sum(upper) works out to
 * Re(sum a - sum b) - 1 for every term of either expansion. The series
 * converges at x = 1 for excess > 0 and at x = -1 for excess > -1; Abel's
 * theorem then makes the sum the radial limit, which is G(1): G is continuous
 * at z = 1 when m + n > p (arg 0 is inside its sector), and for m + n = p the
 * value on the circle is that inner limit. Only z == 1 is admitted: at other
 * points of the circle the term powers z^{b_h} sit on or near their cuts. */
static bool mg_unit_ok(const MgView* g) {
    if (g->p != g->q || !is_int_value(g->z, 1) || g->r) return false;
    double ex;
    if (!mg_param_excess(g, &ex)) return false;
    ex -= 1.0;
    long par = ((long)g->p - (long)g->m - (long)g->n) % 2;
    return par == 0 ? (ex > 1e-12) : (ex > -1.0 + 1e-12);
}

Expr* meijerg_slater_expand(const Expr* call, bool named) {
    MgView g;
    if (!mg_view(call, &g)) return NULL;
    Expr* out = NULL;
    int which = mg_pick(&g);
    if (which >= 0) {
        MgExp E;
        mg_exp_init(&g, which == 1, &E);
        if (!mg_log_case(&E)) out = ev(mg_sum_expr(&E, named));
        mg_exp_free(&E);
    }
    mg_view_free(&g);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Automatic reduction                                                */
/* ------------------------------------------------------------------ */

static Expr* mg_named_shape(const MgView* g, bool dual);

/* The Slater sum of expansion `dual` (named), or NULL in the log case. */
static Expr* mg_expansion_or_null(const MgView* g, bool dual) {
    MgExp E;
    mg_exp_init(g, dual, &E);
    Expr* out = mg_log_case(&E) ? NULL : ev(mg_sum_expr(&E, true));
    mg_exp_free(&E);
    return out;
}

/* G^{2,2}_{2,2}(w | a1, a2; b1, b2) in 1 - w (PBM vol. 3, 8.4.49.22 type):
 *   prod_{i,j} Gamma(1 - a_i + b_j) / Gamma(2 - a1 - a2 + b1 + b2)
 *     * w^b1 2F1(1 - a1 + b1, 1 - a2 + b1; 2 - a1 - a2 + b1 + b2; 1 - w).
 * Both sides have the cut w <= 0. The 2F1 is analytic at w = 1, so this
 * closes the logarithmic case b1 - b2 in Z (where the Slater sum does not
 * exist) and the unit argument (the two-G theorem's K_nu K_mu integrals
 * land there: G(1 | 1/2, 1/2; 0, 0)/4 = Pi^2/4). Checked numerically
 * against the Mellin-Barnes integral at w = 7/10, 13/10, 1/2, 4. Only for
 * FunctionExpand: the automatic MeijerG evaluation keeps the G. Declines
 * where a Gamma of the prefactor has a pole (the G itself is then not
 * defined by the contour). Evaluated and owned, or NULL. */
static Expr* mg_g2222(const MgView* g) {
    if (g->m != 2 || g->n != 2 || g->p != 2 || g->q != 2 || g->r) return NULL;
    const Expr *a1 = g->a[0], *a2 = g->a[1], *b1 = g->b[0], *b2 = g->b[1];
    Expr* g_args[4] = {
        ev(add(sub(I_(1), cp(a1)), cp(b1))), ev(add(sub(I_(1), cp(a1)), cp(b2))),
        ev(add(sub(I_(1), cp(a2)), cp(b1))), ev(add(sub(I_(1), cp(a2)), cp(b2))) };
    Expr* c = ev(add(sub(sub(I_(2), cp(a1)), cp(a2)), add(cp(b1), cp(b2))));
    long N;
    bool pole = exact_nonpos_int(c, &N);
    for (int i = 0; i < 4 && !pole; i++) pole = exact_nonpos_int(g_args[i], &N);
    if (pole) {
        for (int i = 0; i < 4; i++) expr_free(g_args[i]);
        expr_free(c);
        return NULL;
    }
    Expr* pref = dvd(mul(mul(gam(cp(g_args[0])), gam(g_args[1])),
                         mul(gam(cp(g_args[2])), gam(g_args[3]))),
                     gam(cp(c)));
    Expr* f = mk(SYM_Hypergeometric2F1, 4,
                 (Expr*[]){ g_args[0], g_args[2], c, sub(I_(1), cp(g->z)) });
    return ev(mul(pref, mul(pw(cp(g->z), cp(b1)), f)));
}

Expr* meijerg_function_expand(const Expr* call) {
    MgView g;
    if (!mg_view(call, &g)) return NULL;
    Expr* out = mg_named_shape(&g, false);
    if (!out && !g.r) out = mg_named_shape(&g, true);
    if (!out) {
        int which = mg_pick(&g);
        /* p == q at z = 1: the Slater terms are pF(p-1) at +-1 and may each
         * diverge where G is finite -- G^{2,2}_{2,2}(1 | 1/3, 2/3; 0, 1/3)
         * gave ComplexInfinity for 5.5009. Sum them only when every term
         * converges there (mg_unit_ok), as the automatic reduction does. */
        if (g.p == g.q && is_int_value(g.z, 1) && !mg_unit_ok(&g)) which = -2;
        if (which >= 0) {
            out = mg_expansion_or_null(&g, which == 1);
        } else if (!is_number(g.z)) {
            /* p == q, m + n <= p, symbolic z: G is one expansion inside the
             * unit circle and the other outside -- a Piecewise, the inner one
             * also owning z == 1 when its terms converge there. */
            Expr* in = mg_expansion_or_null(&g, false);
            Expr* outer = in ? mg_expansion_or_null(&g, true) : NULL;
            if (in && outer) {
                Expr* inner_cond = mk2(SYM_Less, mk1(SYM_Abs, cp(g.z)), I_(1));
                double ex;
                long par = ((long)g.p - (long)g.m - (long)g.n) % 2;
                if (!g.r && mg_param_excess(&g, &ex) &&
                    (par == 0 ? (ex - 1.0 > 1e-12) : (ex - 1.0 > -1.0 + 1e-12)))
                    inner_cond = mk2(SYM_Or, mk2(SYM_Equal, cp(g.z), I_(1)), inner_cond);
                Expr* outer_cond = mk2(SYM_Greater, mk1(SYM_Abs, cp(g.z)), I_(1));
                Expr* pairs = mk2(SYM_List, mk2(SYM_List, in, inner_cond),
                                  mk2(SYM_List, outer, outer_cond));
                out = ev(mk2(SYM_Piecewise, pairs, expr_new_symbol(SYM_Indeterminate)));
            } else {
                expr_free(in); expr_free(outer);
            }
        }
    }
    /* the Slater route declined (log case, or unit argument) */
    if (!out) out = mg_g2222(&g);
    mg_view_free(&g);
    return out;
}

/* The policy, calibrated against Mathematica 14:
 *   - G^{2,0}_{0,2}(;b,c)        -> 2 w^{(b+c)/2} K_{b-c}(2 sqrt w)    (always)
 *   - G^{2,1}_{1,2}(a;b,c)       -> G(1-a+b) G(1-a+c) w^b U(1-a+b, 1+b-c, w)
 *   - one Slater term (m = 1 in w, or n = 1 in 1/w): always, the pFq named
 *     where possible (power, Exp, Bessel, 1F1(u;2u), EllipticK/E);
 *   - two Slater terms whose pFq is a 0F1, 1F1 or 2F1, when the two exponent
 *     parameters are not an integer apart (no log case);
 *   - anything else stays MeijerG. */
/* The named shapes that are not Slater sums (their expansions are log cases
 * for integer-spaced parameters), read in t = w (dual false) or, through
 * G^{m,n}_{p,q}(w | a; b) = G^{n,m}_{q,p}(1/w | 1-b; 1-a), in t = 1/w (dual
 * true):
 *   G^{2,0}_{0,2}(t | ;b,c)  = 2 t^{(b+c)/2} K_{b-c}(2 sqrt t)
 *   G^{2,1}_{1,2}(t | a;b,c) = G(1-a+b) G(1-a+c) t^b U(1-a+b, 1+b-c, t)
 * Evaluated and owned, or NULL. */
static Expr* mg_named_shape(const MgView* g, bool dual) {
    MgExp E; mg_exp_init(g, dual, &E);
    Expr* r = NULL;
    if (E.p == 0 && E.q == 2 && E.m == 2) {
        Expr* s = add(cp(E.B[0]), cp(E.B[1]));
        r = mul(mul(I_(2), tpow(&E, mul(half(), s))),
                mk2(SYM_BesselK, sub(cp(E.B[0]), cp(E.B[1])),
                    mul(I_(2), tpow(&E, half()))));
    } else if (E.p == 1 && E.n == 1 && E.q == 2 && E.m == 2) {
        const Expr *a = E.A[0], *b = E.B[0], *c = E.B[1];
        r = mul(mul(gam(sub(add(I_(1), cp(b)), cp(a))), gam(sub(add(I_(1), cp(c)), cp(a)))),
                mul(tpow(&E, cp(b)),
                    mk3(SYM_HypergeometricU, sub(add(I_(1), cp(b)), cp(a)),
                        sub(add(I_(1), cp(b)), cp(c)), tpow(&E, I_(1)))));
    }
    mg_exp_free(&E);
    return r ? ev(r) : NULL;
}

static Expr* mg_try_reduce(const MgView* g) {
    /* Named special shapes, in w only: the dual spellings (in 1/w) are left
     * to FunctionExpand, as the reference system does --
     * MeijerG[{{0,0},{}},{{},{}},z] stays, FunctionExpand gives BesselK. */
    Expr* named = mg_named_shape(g, false);
    if (named) return named;
    int which = mg_pick(g);
    /* p == q on the unit circle: the Slater terms' pFq (a 2F1 for p = 2) sit on
     * their circle of convergence and may diverge individually where G does
     * not. Keep MeijerG, whose numerics handle it -- except at z = 1 when
     * every term converges there (mg_unit_ok): the sum is then G(1), from the
     * expansion in z when m + n <= p (the inner limit), else the usual one.
     * MeijerG[{{1/2},{1}},{{0,1/2},{}},1] = Pi/2. */
    if (g->p == g->q && is_exact_number(g->z) && mg_abs_cmp_one(g->z) == 0) {
        if (!mg_unit_ok(g)) return NULL;
        if (which < 0) which = 0;
    }
    if (which < 0) return NULL;
    MgExp E;
    mg_exp_init(g, which == 1, &E);
    Expr* out = NULL;
    size_t pu = E.p, ql = E.q ? E.q - 1 : 0;
    bool named_shape = (pu == 0 && ql == 1) || (pu == 1 && ql == 1) || (pu == 2 && ql == 1);
    if (E.m == 0) {
        out = I_(0);
    } else if (E.m == 1 || (E.m == 2 && named_shape && !mg_log_case(&E))) {
        out = ev(mg_sum_expr(&E, true));
    }
    mg_exp_free(&E);
    return out;
}

/* ------------------------------------------------------------------ */
/*  z = 0                                                              */
/* ------------------------------------------------------------------ */

/* Real part of an exact numeric leaf, as a sign: -1, 0, +1; false if not an
 * exact (or inexact) number. */
static bool re_sign(const Expr* e, int* s) {
    const Expr* re = e;
    if (is_head(e, SYM_Complex) && e->data.function.arg_count == 2) re = e->data.function.args[0];
    double d;
    if (!real_value(re, &d)) return false;
    *s = (d > 0) - (d < 0);
    return true;
}

static Expr* mg_at_zero(const MgView* g) {
    if (g->p > g->q) {
        if (g->n == 1) return expr_new_symbol(SYM_ComplexInfinity);
        return NULL;
    }
    bool inexact = !is_exact_zero(g->z);
    Expr* out = NULL;
    if (g->m == 0) {
        out = I_(0);
    } else {
        int minsign = 2;
        for (size_t h = 0; h < g->m; h++) {
            int s;
            if (!re_sign(g->b[h], &s)) return NULL;     /* symbolic: undecided */
            if (s < minsign) minsign = s;
        }
        MgExp E;
        mg_exp_init(g, false, &E);
        if (minsign > 0) {
            out = I_(0);
        } else if (mg_log_case(&E)) {
            out = NULL;                                  /* Log terms: leave it */
        } else {
            /* Terms with Re b_h < 0 blow up (unless their coefficient is zero);
             * terms with b_h = 0 give their coefficient; Re b_h > 0 give 0. */
            Expr* acc = I_(0);
            bool blowup = false, undecided = false;
            for (size_t h = 0; h < g->m && !undecided; h++) {
                int s = 0; (void)re_sign(g->b[h], &s);
                if (s > 0) continue;
                MgTerm T; mg_term(&E, h, &T);
                Expr* c = cp(T.C);
                for (size_t k = 0; k < T.ql; k++)
                    if (T.reg[k]) c = dvd(c, gam(cp(T.l[k])));
                c = ev(c);
                bool zero_c = is_exact_zero(c);
                if (!zero_c) {
                    if (s < 0) blowup = true;
                    else if (is_exact_zero(T.expo)) acc = add(acc, cp(c));
                    else undecided = true;   /* Re 0, Im != 0: oscillates */
                }
                expr_free(c);
                mg_term_free(&T);
            }
            if (undecided) { expr_free(acc); out = NULL; }
            else if (blowup) { expr_free(acc); out = expr_new_symbol(SYM_ComplexInfinity); }
            else out = ev(acc);
        }
        mg_exp_free(&E);
    }
    if (out && inexact && !is_head(out, SYM_DirectedInfinity) &&
        !(out->type == EXPR_SYMBOL && out->data.symbol.name == SYM_ComplexInfinity))
        out = ev(mk1(SYM_N, out));
    return out;
}

/* ------------------------------------------------------------------ */
/*  Series support                                                     */
/* ------------------------------------------------------------------ */

/* sum_h t^{e_h} sum_{k<=order} c_{h,k} x^k for expansion E, with x = sgn t,
 * coefficients from the regularized series  R = sum_k prod (u)_k / (k! prod
 * Gamma(l+k)) x^k  (finite at every lower-parameter pole). */
static Expr* mg_truncated(const MgExp* E, long order) {
    if (E->m == 0) return I_(0);
    size_t nt = E->m;
    Expr** terms = malloc(sizeof(Expr*) * nt);
    for (size_t h = 0; h < nt; h++) {
        MgTerm T;
        mg_term(E, h, &T);
        /* Undo the integer shift's bookkeeping: the truncated series is built
         * directly from the regularized coefficients, so the shifted term is
         * already in T (a leading t^{N+1} folded into expo). */
        Expr** ks = malloc(sizeof(Expr*) * (size_t)(order + 1));
        for (long k = 0; k <= order; k++) {
            Expr* c = cp(T.C);
            for (size_t i = 0; i < T.pu; i++) c = mul(c, mk2(SYM_Pochhammer, cp(T.u[i]), I_(k)));
            for (size_t i = 0; i < T.ql; i++) {
                if (T.reg[i]) c = dvd(c, gam(add(cp(T.l[i]), I_(k))));
                else          c = dvd(c, mk2(SYM_Pochhammer, cp(T.l[i]), I_(k)));
            }
            c = dvd(c, mk1(SYM_Factorial, I_(k)));
            if (T.sgn < 0 && (k % 2)) c = neg(c);
            ks[k] = mul(c, tpow(E, add(cp(T.expo), I_(k))));
        }
        terms[h] = mk(SYM_Plus, (size_t)(order + 1), ks);
        free(ks);
        mg_term_free(&T);
    }
    Expr* s = (nt == 1) ? terms[0] : mk(SYM_Plus, nt, terms);
    free(terms);
    return s;
}

/* ---- The log case: the eps -> 0 limit, coefficient by coefficient ----
 *
 * With B_j -> B_j + (j+1) eps (j < m) every Slater coefficient becomes a
 * Laurent series in eps whose poles cancel across the terms of one cluster of
 * integer-spaced B's; the limit is the sum of the eps^0 coefficients. Each
 * coefficient is a product of four kinds of factor, all expanded here exactly
 * (Series cannot expand Gamma about its poles):
 *   Gamma(al + be eps)   regular:  Gamma(al) exp(sum_j psi^(j-1)(al) (be eps)^j / j!)
 *                        al = -N:  Gamma(1 + d) / (d prod_{i=1..N} (d - i)), d = be eps
 *   linear (al + be eps) and its reciprocal (Pochhammer factors)
 *   z^((h+1) eps) = exp((h+1) eps Log[z]).
 * A truncated Laurent series sum_{i<len} c[i] eps^(v+i). */
typedef struct { int v; int len; Expr** c; } LSer;

static LSer ls_new(int v, int len) {
    LSer r; r.v = v; r.len = len;
    r.c = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
    for (int i = 0; i < len; i++) r.c[i] = I_(0);
    return r;
}
static void ls_free(LSer* a) {
    for (int i = 0; i < a->len; i++) expr_free(a->c[i]);
    free(a->c); a->c = NULL; a->len = 0;
}
static void ls_set(LSer* a, int i, Expr* e) {
    if (i < 0 || i >= a->len) { expr_free(e); return; }
    expr_free(a->c[i]); a->c[i] = e;
}
/* a * b, both consumed; result truncated to `len` terms. */
static LSer ls_mul(LSer a, LSer b, int len) {
    LSer r = ls_new(a.v + b.v, len);
    for (int i = 0; i < len; i++) {
        Expr* acc = I_(0);
        for (int j = 0; j <= i; j++)
            if (j < a.len && i - j < b.len)
                acc = add(acc, mul(cp(a.c[j]), cp(b.c[i - j])));
        ls_set(&r, i, ev(acc));
    }
    ls_free(&a); ls_free(&b);
    return r;
}
/* exp(sum_{j>=1} L[j] eps^j) for a valuation-0 series L with L[0] = 0. */
static LSer ls_exp0(const Expr* const* L, int len) {
    LSer r = ls_new(0, len);
    ls_set(&r, 0, I_(1));
    for (int nn = 1; nn < len; nn++) {
        Expr* acc = I_(0);
        for (int j = 1; j <= nn; j++)
            acc = add(acc, mul(mul(I_(j), cp(L[j])), cp(r.c[nn - j])));
        ls_set(&r, nn, ev(dvd(acc, I_(nn))));
    }
    return r;
}
/* (al + be eps)^(+1 or -1). al, be borrowed. */
static LSer ls_linear(const Expr* al, const Expr* be, bool recip, int len) {
    if (is_exact_zero(al)) {
        LSer r = ls_new(recip ? -1 : 1, len);
        ls_set(&r, 0, recip ? ev(pw(cp(be), I_(-1))) : cp(be));
        return r;
    }
    LSer r = ls_new(0, len);
    if (!recip) {
        ls_set(&r, 0, cp(al));
        if (len > 1) ls_set(&r, 1, cp(be));
        return r;
    }
    /* 1/(al (1 + (be/al) eps)) = sum (-be/al)^i / al */
    Expr* t = ev(neg(dvd(cp(be), cp(al))));
    Expr* inv = ev(pw(cp(al), I_(-1)));
    for (int i = 0; i < len; i++) ls_set(&r, i, ev(mul(cp(inv), pw(cp(t), I_(i)))));
    expr_free(t); expr_free(inv);
    return r;
}
/* Gamma(al + be eps)^(+1 or -1). */
static LSer ls_gamma(const Expr* al, const Expr* be, bool inv, int len) {
    long N;
    if (exact_nonpos_int(al, &N)) {
        /* Gamma(-N + d) = Gamma(1 + d) / (d prod_{i=1..N} (d - i)). */
        Expr* one = I_(1);
        LSer r = ls_gamma(one, be, inv, len);
        expr_free(one);
        Expr* zero = I_(0);
        r = ls_mul(r, ls_linear(zero, be, !inv, len), len);
        expr_free(zero);
        for (long i = 1; i <= N; i++) {
            Expr* mi = I_(-i);
            r = ls_mul(r, ls_linear(mi, be, !inv, len), len);
            expr_free(mi);
        }
        return r;
    }
    Expr** L = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
    L[0] = I_(0);
    for (int j = 1; j < len; j++) {
        Expr* t = mul(mk2(SYM_PolyGamma, I_(j - 1), cp(al)),
                      dvd(pw(cp(be), I_(j)), mk1(SYM_Factorial, I_(j))));
        L[j] = ev(inv ? neg(t) : t);
    }
    LSer e = ls_exp0((const Expr* const*)L, len);
    for (int j = 0; j < len; j++) expr_free(L[j]);
    free(L);
    LSer g = ls_new(0, 1);
    ls_set(&g, 0, ev(inv ? pw(gam(cp(al)), I_(-1)) : gam(cp(al))));
    return ls_mul(g, e, len);
}
/* exp(be eps L) with L = Log[z]. */
static LSer ls_zpow(const Expr* be, const Expr* logz, int len) {
    Expr** L = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
    L[0] = I_(0);
    for (int j = 1; j < len; j++) L[j] = j == 1 ? ev(mul(cp(be), cp(logz))) : I_(0);
    LSer e = ls_exp0((const Expr* const*)L, len);
    for (int j = 0; j < len; j++) expr_free(L[j]);
    free(L);
    return e;
}

/* Valuation of Gamma(al + .)^(+-1) and of a linear (al + .)^(+-1). */
static int val_gamma(const Expr* al, bool inv) { long N; return exact_nonpos_int(al, &N) ? (inv ? 1 : -1) : 0; }
static int val_linear(const Expr* al, bool recip) { return is_exact_zero(al) ? (recip ? -1 : 1) : 0; }

/* The coefficient of z^{B_h + k} in the log case (eps^0 part), or NULL when
 * it is zero. `logz` is Log[z] (z^(1/r) for the generalized form: logz is then
 * Log[z]/r). */
static Expr* mg_log_coef(const MgExp* E, size_t h, long k, const Expr* logz) {
    size_t p = E->p, q = E->q, m = E->m, n = E->n;
    /* Collect factors as (kind, al, be, flag): kind 0 = Gamma, 1 = linear. */
    typedef struct { int kind; Expr* al; Expr* be; bool inv; } Fac;
    size_t cap = 4 * (p + q + 2) + (size_t)(k + 1) * (p + q + 2);
    Fac* F = malloc(sizeof(Fac) * cap);
    size_t nf = 0;
    long eh = (long)h + 1;
    for (size_t j = 0; j < m; j++) if (j != h) {
        F[nf].kind = 0; F[nf].al = ev(sub(cp(E->B[j]), cp(E->B[h])));
        F[nf].be = I_((long)j + 1 - eh); F[nf].inv = false; nf++;
    }
    for (size_t j = 0; j < n; j++) {
        F[nf].kind = 0; F[nf].al = ev(sub(add(I_(1), cp(E->B[h])), cp(E->A[j])));
        F[nf].be = I_(eh); F[nf].inv = false; nf++;
    }
    for (size_t j = n; j < p; j++) {
        F[nf].kind = 0; F[nf].al = ev(sub(cp(E->A[j]), cp(E->B[h])));
        F[nf].be = I_(-eh); F[nf].inv = true; nf++;
    }
    for (size_t i = 0; i < p; i++)
        for (long t = 0; t < k; t++) {
            F[nf].kind = 1; F[nf].al = ev(add(sub(add(I_(1), cp(E->B[h])), cp(E->A[i])), I_(t)));
            F[nf].be = I_(eh); F[nf].inv = false; nf++;
        }
    for (size_t j = 0; j < q; j++) {
        if (j == h) continue;
        if (j < m) {
            for (long t = 0; t < k; t++) {
                F[nf].kind = 1;
                F[nf].al = ev(add(sub(add(I_(1), cp(E->B[h])), cp(E->B[j])), I_(t)));
                F[nf].be = I_(eh - ((long)j + 1)); F[nf].inv = true; nf++;
            }
        } else {
            F[nf].kind = 0;
            F[nf].al = ev(add(sub(add(I_(1), cp(E->B[h])), cp(E->B[j])), I_(k)));
            F[nf].be = I_(eh); F[nf].inv = true; nf++;
        }
    }
    int V = 0;
    for (size_t i = 0; i < nf; i++)
        V += F[i].kind == 0 ? val_gamma(F[i].al, F[i].inv) : val_linear(F[i].al, F[i].inv);
    Expr* out = NULL;
    if (V <= 0) {
        int len = 1 - V;
        LSer acc = ls_new(0, len);
        ls_set(&acc, 0, I_(1));
        for (size_t i = 0; i < nf; i++) {
            LSer f = F[i].kind == 0 ? ls_gamma(F[i].al, F[i].be, F[i].inv, len)
                                    : ls_linear(F[i].al, F[i].be, F[i].inv, len);
            /* Re-base f to its own valuation's len terms (already). */
            acc = ls_mul(acc, f, len);
        }
        Expr* beh = I_(eh);
        acc = ls_mul(acc, ls_zpow(beh, logz, len), len);
        expr_free(beh);
        /* eps^0 sits at index -acc.v (acc.v == V by construction). */
        int idx = -acc.v;
        Expr* c0 = (idx >= 0 && idx < acc.len) ? cp(acc.c[idx]) : I_(0);
        ls_free(&acc);
        c0 = dvd(c0, mk1(SYM_Factorial, I_(k)));
        if (mg_sign(E) < 0 && (k % 2)) c0 = neg(c0);
        out = ev(c0);
        if (is_exact_zero(out)) { expr_free(out); out = NULL; }
    }
    for (size_t i = 0; i < nf; i++) { expr_free(F[i].al); expr_free(F[i].be); }
    free(F);
    return out;
}

/* The log-case truncated expansion. For each term h, powers k <= order - s_h
 * where s_h = B_h - (smallest B of its integer-spaced cluster): only those
 * absolute powers are complete, i.e. receive every cluster member's pole part,
 * so only their eps^0 coefficients are the limit. */
static Expr* mg_truncated_log(const MgExp* E, long order) {
    if (E->m == 0) return I_(0);
    /* log w = Log[z] / r, and +-: the dual expansion is in 1/w. */
    Expr* logz = mk1(SYM_Log, cp(E->z));
    if (E->r) logz = dvd(logz, cp(E->r));
    if (E->dual) logz = neg(logz);
    logz = ev(logz);
    size_t cap = E->m * (size_t)(order + 1) + 1;
    Expr** terms = malloc(sizeof(Expr*) * cap);
    size_t nt = 0;
    for (size_t h = 0; h < E->m; h++) {
        long s_h = 0;
        for (size_t j = 0; j < E->m; j++) {
            Expr* d = diff_of(E->B[h], E->B[j]);
            long v;
            if (integer_valued(d, &v) && v > s_h) s_h = v;
            expr_free(d);
        }
        for (long k = 0; k + s_h <= order; k++) {
            Expr* c = mg_log_coef(E, h, k, logz);
            if (!c) continue;
            terms[nt++] = mul(c, tpow(E, add(cp(E->B[h]), I_(k))));
        }
    }
    expr_free(logz);
    Expr* s = nt == 0 ? I_(0) : (nt == 1 ? terms[0] : mk(SYM_Plus, nt, terms));
    free(terms);
    return s;
}

Expr* meijerg_series_terms(const Expr* call, bool at_infinity, long order) {
    MgView g;
    if (!mg_view(call, &g)) return NULL;
    if (order < 0) order = 0;
    if (order > 40) order = 40;
    Expr* out = NULL;
    /* At 0 the expansion in w (p <= q: for p == q it is the |z| < 1 branch,
     * the right local expansion even when m + n == p); at Infinity the dual. */
    bool ok = at_infinity ? (g.p >= g.q) : (g.p <= g.q);
    if (ok) {
        MgExp E;
        mg_exp_init(&g, at_infinity, &E);
        out = ev(mg_log_case(&E) ? mg_truncated_log(&E, order) : mg_truncated(&E, order));
        mg_exp_free(&E);
    }
    mg_view_free(&g);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Entry point                                                        */
/* ------------------------------------------------------------------ */

Expr* builtin_meijer_g(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t na = res->data.function.arg_count;
    if (na != 3 && na != 4) {
        msg_hdiv(res);
        return builtin_arg_error("MeijerG", na, 3, 4);
    }
    MgView g;
    if (!mg_view(res, &g)) { msg_hdiv(res); return NULL; }

    Expr* out = NULL;

    /* (1) r == 1: the ordinary G-function. */
    if (g.r && is_int_value(g.r, 1)) {
        MgView g3 = g; g3.r = NULL;
        out = mg_rebuild(&g3, cp(g.z));
        goto done;
    }

    /* (2) A machine buffer: the machine kernel over it. On decline the buffer
     * must not reach the scalar steps below, which would leave it unevaluated. */
    if (is_ndarray(g.z)) {
        out = mg_ndarray(&g);
        if (!out) out = ndarray_delist_and_reeval(res);
        goto done;
    }
    /* (2b) Thread over a List third argument. */
    if (is_list(g.z)) {
        size_t n = g.z->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t i = 0; i < n; i++)
            items[i] = ev(mg_rebuild(&g, cp(g.z->data.function.args[i])));
        out = mk_list(items, n);
        free(items);
        goto done;
    }

    /* (3) Canonical parameter order. */
    out = mg_canonical(&g);
    if (out) goto done;

    /* (3b) Gamma-pair cancellation (equal a_k, b_j in opposite groups). */
    {
        bool hdiv = false;
        out = mg_cancel_pairs(&g, &hdiv);
        if (hdiv) msg_hdiv(res);
        if (out || hdiv) goto done;
    }

    /* (4) Consistency: m = n = 0 has no contour at all; a pole collision is
     * resolved by cancelling a Gamma pair when one exists, else hdiv. */
    if (g.m == 0 && g.n == 0) { msg_hdiv(res); goto done; }
    if (mg_collision(&g)) {
        out = mg_move_pair(&g);
        if (!out) msg_hdiv(res);
        goto done;
    }

    /* (5) Machine precision, (6) arbitrary precision. */
    out = mg_machine(&g);
    if (out) goto done;
    out = mg_arbitrary(&g);
    if (out) goto done;

    /* (7) z = 0. */
    if (is_any_zero(g.z)) { out = mg_at_zero(&g); goto done; }

    /* (8) Reduction to simpler functions -- not for an all-numeric inexact
     * call whose numerics declined (a pole, or no expansion implemented): the
     * rewrite would only trade MeijerG for unevaluable HypergeometricPFQs. */
    if (mg_all_numeric_inexact(&g)) goto done;
    out = mg_try_reduce(&g);

done:
    mg_view_free(&g);
    return out;
}

/* ------------------------------------------------------------------------ */
/* MeijerGReduce[expr, x]                                                   */
/* ------------------------------------------------------------------------ */

/* The rewrite is Roach`MeijerGReduce in src/internal/RoachIntegrate.m: it
 * reuses the function -> G table (Roach`GForm1, PBM vol. 3 sec. 8.4) that the
 * Method -> "Roach" integrator decomposes integrands with, so the two cannot
 * drift apart. The module is loaded lazily, once per symbol-table generation,
 * exactly as try_roach (src/calculus/integrate.c) loads it.
 *
 * Returns the rewritten expression, or NULL (unevaluated) when no factor of
 * expr has a known G form or the second argument is not a symbol. */
static unsigned mgr_load_gen = 0;
static bool mgr_load_ok = false;

static Expr* builtin_meijer_g_reduce(Expr* res) {
    size_t na = res->data.function.arg_count;
    if (na != 2) return builtin_arg_error("MeijerGReduce", na, 2, 2);
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    if (mgr_load_gen != symtab_generation()) {
        mgr_load_gen = symtab_generation();
        mth_msg_suppress_push();
        mgr_load_ok = mathilda_load_module("RoachIntegrate.m") != 0;
        mth_msg_suppress_pop();
        if (!mgr_load_ok)
            mth_message("MeijerGReduce", "nofile",
                "cannot locate src/internal/RoachIntegrate.m on disk.");
    }
    if (!mgr_load_ok) return NULL;
    Expr* call = expr_new_function(expr_new_symbol("Roach`MeijerGReduce"),
        (Expr*[]){ expr_copy(res->data.function.args[0]), expr_copy(x) }, 2);
    /* the table's pattern conditions run Simplify/Refine probes whose
     * messages are internal; they must reach neither the user nor Check[] */
    unsigned long fired0 = mth_msg_fired_count();
    mth_msg_suppress_push();
    Expr* r = evaluate(call);
    mth_msg_suppress_pop();
    mth_msg_fired_restore(fired0);
    expr_free(call);
    if (!r) return NULL;
    if ((r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_DollarFailed) ||
        (r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL &&
         strcmp(r->data.function.head->data.symbol.name, "Roach`MeijerGReduce") == 0)) {
        expr_free(r);
        return NULL;
    }
    return r;
}

void meijer_g_init(void) {
    symtab_add_builtin("MeijerG", builtin_meijer_g);
    symtab_get_def("MeijerG")->attributes |= ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_add_builtin("MeijerGReduce", builtin_meijer_g_reduce);
    symtab_get_def("MeijerGReduce")->attributes |= ATTR_PROTECTED;
}
