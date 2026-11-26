/*
 * simp_sf_hyper.c -- special-function Simplify descriptors: the
 * hypergeometric functions 2F1, 1F1, 0F1 (all evaluated to
 * HypergeometricPFQ) and HypergeometricU.
 *
 * Derivations in docs/design/simp_special_functions.md (section 10); every
 * operator below was checked numerically at real and complex parameters.
 *
 * Lattice access.  The evaluator rewrites 2F1/1F1/0F1 to
 * HypergeometricPFQ[{a1..ap}, {b1..bq}, z], whose shift parameters sit inside
 * lists.  The canon hook maps a PFQ of shape (2,1), (1,1), (0,1) onto a
 * private head with the parameters as direct arguments --
 * SfEngine$F21[a, b, c, z], SfEngine$F11[a, b, z], SfEngine$F01[b, z] --
 * with the upper parameters ordered by their CLASS (symbolic part, then
 * rational part mod 1), never by their offset, so 2F1(a, b; c) and
 * 2F1(1 + b, a; c) land in one family.  Any other PFQ has both lists put in
 * that same canonical order.  The finalize hook maps the private heads back.
 *
 * Shift matrices from raising operators.  With theta = z d/dz, each
 * function F satisfies a second-order ODE  theta^2 F = alpha F + beta theta F
 * and each parameter direction d a raising relation
 * F(p + e_d) = (u_d + v_d theta) F(p).  In the engine's basis
 * (F, F_0) with F_0 = F(p + e_0) the up matrix for d follows mechanically:
 * G = (u_d + v_d theta) F, G_0 = (u_0' + v_0' theta) G (u_0', v_0' at p + e_d),
 * theta acting on coefficients by z d/dz and on theta F through the ODE; then
 * theta F = (F_0 - u_0 F)/v_0 converts to the basis.  The data (all checked):
 *
 *   2F1(a, b; c; z):  alpha = a b z/(1 - z),  beta = (z (a + b) - c + 1)/(1 - z)
 *     a: (1, 1/a)   b: (1, 1/b)   [(theta + a) F = a F(a+1), from the series]
 *     c: u = c (c - a - b)/((c - a)(c - b)),  v = c (1 - z)/((c - a)(c - b) z)
 *        [Gauss: c(c-1)(z-1) F(c-1) + c(c-1-(2c-a-b-1) z) F + (c-a)(c-b) z F(c+1)
 *         = 0 with (theta + c - 1) F = (c - 1) F(c-1)]
 *   1F1(a; b; z):  alpha = a z,  beta = z - b + 1
 *     a: (1, 1/a)   b: u = b/(b - a),  v = -b/((b - a) z)
 *        [b(b-1) M(b-1) + b(1-b-z) M + z(b-a) M(b+1) = 0, (theta + b - 1) M = (b-1) M(b-1)]
 *   0F1(; b; z):   alpha = z,  beta = 1 - b
 *     b: (0, b/z)   [d/dz 0F1(; b; z) = 0F1(; b + 1; z)/b]
 *   U(a, b, z), direction a only, three-term:
 *     U(a-1) + (b - 2a - z) U(a) + a (a - b + 1) U(a+1) = 0
 *
 * Transformations (canon, when the partner atom is present):
 *   Kummer   1F1(a; b; z) = E^z 1F1(b - a; b; -z)              (entire; all z)
 *            U(a, b, z) = z^(1-b) U(a - b + 1, 2 - b, z)        (all z checked)
 *   Euler    2F1(a, b; c; z) = (1 - z)^(c-a-b) 2F1(c - a, c - b; c; z)
 *   Pfaff    2F1(a, b; c; z) = (1 - z)^(-a) 2F1(a, c - b; c; z/(z - 1))
 *            -- both only for z < 1 proved: on the cut (1, Inf) the branch
 *            conventions could not be checked numerically.
 */

#include "hypergeopfq.h"   /* hyper_named_to_pfq */
#include <stdlib.h>
#include <string.h>
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "simp.h"
#include "simp_internal.h"
#include "simp_sf_engine.h"

static const char* SYM_SfF21 = NULL;
static const char* SYM_SfF11 = NULL;
static const char* SYM_SfF01 = NULL;

static void hyper_init(void) {
    Expr* a = sf_sym("SfEngine$F21"); SYM_SfF21 = a->data.symbol.name; expr_free(a);
    Expr* b = sf_sym("SfEngine$F11"); SYM_SfF11 = b->data.symbol.name; expr_free(b);
    Expr* c = sf_sym("SfEngine$F01"); SYM_SfF01 = c->data.symbol.name; expr_free(c);
}

/* ---- theta-basis operator algebra ---------------------------------------- */

/* u F + v theta F, owned coefficients */
typedef struct { Expr* u; Expr* v; } TV;

static void tv_free(TV* t) { expr_free(t->u); expr_free(t->v); }

/* z d/dz of a coefficient */
static Expr* zD(const Expr* c, const Expr* z) {
    return sf_times(sf_cp(z), eval_and_free(sf_call("D", (Expr*[]){ sf_cp(c), sf_cp(z) }, 2)));
}

/* theta (u F + v theta F) = (z u' + v alpha) F + (u + z v' + v beta) theta F */
static TV tv_theta(const TV* x, const Expr* alpha, const Expr* beta, const Expr* z) {
    TV r;
    r.u = sf_plus(zD(x->u, z), sf_times(sf_cp(x->v), sf_cp(alpha)));
    r.v = sf_plus(sf_plus(sf_cp(x->u), zD(x->v, z)), sf_times(sf_cp(x->v), sf_cp(beta)));
    return r;
}

/* (r0 + r1 theta) applied to x */
static TV tv_apply(const Expr* r0, const Expr* r1, const TV* x,
                   const Expr* alpha, const Expr* beta, const Expr* z) {
    TV th = tv_theta(x, alpha, beta, z);
    TV r;
    r.u = sf_plus(sf_times(sf_cp(r0), sf_cp(x->u)), sf_times(sf_cp(r1), sf_cp(th.u)));
    r.v = sf_plus(sf_times(sf_cp(r0), sf_cp(x->v)), sf_times(sf_cp(r1), sf_cp(th.v)));
    tv_free(&th);
    return r;
}

/* u F + v theta F in the basis (F, F_0), F_0 = (b0 + b1 theta) F:
 * theta F = (F_0 - b0 F)/b1  =>  (u - v b0/b1) F + (v/b1) F_0.
 * Writes the two coefficients into row[0], row[1]. */
static void tv_to_basis(const TV* x, const Expr* b0, const Expr* b1, Expr** row) {
    row[0] = sf_minus(sf_cp(x->u), sf_div(sf_times(sf_cp(x->v), sf_cp(b0)), sf_cp(b1)));
    row[1] = sf_div(sf_cp(x->v), sf_cp(b1));
}

/* Raising data of one family at a point: ODE and per-direction operators. */
typedef struct {
    Expr* alpha; Expr* beta;           /* theta^2 F = alpha F + beta theta F */
    Expr* r0[3]; Expr* r1[3];          /* direction d: F(p+e_d) = (r0 + r1 theta) F */
    int   nd;
} Raise;

static void raise_free(Raise* R) {
    expr_free(R->alpha); expr_free(R->beta);
    for (int d = 0; d < R->nd; d++) { expr_free(R->r0[d]); expr_free(R->r1[d]); }
}

typedef void (*RaiseFn)(Expr* const* a, Raise* R);

static Expr* subst1(const Expr* e, const Expr* from, const Expr* to) {
    if (expr_eq(e, from)) return sf_cp(to);
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = subst1(e->data.function.args[i], from, to);
    Expr* r = expr_new_function(sf_cp(e->data.function.head), a, n);
    free(a);
    return r;
}

/* The up matrix for direction d at point args a (`shifted` = args at
 * p + e_d), rows 0..1 into M.  theta = z d/dz differentiates the
 * coefficients in the function's own argument, which may be an expression
 * (2 z): the matrix is built with a private symbol in its place and the
 * argument substituted back afterwards.  `zi` is the argument's index. */
static bool up_from_raise_sym(RaiseFn fill, Expr* const* a0, Expr* const* shifted0, size_t n,
                              size_t zi, int d, Expr** M);

static bool up_from_raise(RaiseFn fill, Expr* const* a, Expr* const* shifted, int d,
                          const Expr* z, Expr** M) {
    Raise R, Rs;
    fill(a, &R);
    if (d >= R.nd) { raise_free(&R); return false; }
    fill(shifted, &Rs);
    TV one = { sf_int(1), sf_int(0) };
    /* G = F(p + e_d), G_0 = F(p + e_d + e_0) = (r0' + r1' theta) G */
    TV G = tv_apply(R.r0[d], R.r1[d], &one, R.alpha, R.beta, z);
    TV G0 = tv_apply(Rs.r0[0], Rs.r1[0], &G, R.alpha, R.beta, z);
    tv_to_basis(&G, R.r0[0], R.r1[0], &M[0]);
    tv_to_basis(&G0, R.r0[0], R.r1[0], &M[3]);
    M[2] = sf_int(0); M[5] = sf_int(0);
    tv_free(&one); tv_free(&G); tv_free(&G0);
    raise_free(&R); raise_free(&Rs);
    return true;
}

static bool up_from_raise_sym(RaiseFn fill, Expr* const* a0, Expr* const* shifted0, size_t n,
                              size_t zi, int d, Expr** M) {
    Expr* u = sf_sym("SfEngine$theta");
    Expr** a = malloc(sizeof(Expr*) * n);
    Expr** s = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        a[i] = (i == zi) ? sf_cp(u) : sf_cp(a0[i]);
        s[i] = (i == zi) ? sf_cp(u) : sf_cp(shifted0[i]);
    }
    bool ok = up_from_raise(fill, a, s, d, u, M);
    if (ok)
        for (int i = 0; i < 6; i++) {
            Expr* t = eval_and_free(subst1(M[i], u, a0[zi]));
            expr_free(M[i]);
            M[i] = t;
        }
    for (size_t i = 0; i < n; i++) { expr_free(a[i]); expr_free(s[i]); }
    free(a); free(s);
    expr_free(u);
    return ok;
}

/* args a with shift argument `idx` increased by one (owned array) */
static Expr** args_bumped(Expr* const* a, size_t n, size_t idx) {
    Expr** s = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++)
        s[i] = (i == idx) ? sf_plus(sf_cp(a[i]), sf_int(1)) : sf_cp(a[i]);
    return s;
}

static void args_free(Expr** a, size_t n) { for (size_t i = 0; i < n; i++) expr_free(a[i]); free(a); }

/* ---- 2F1 ------------------------------------------------------------------ */

static void raise_f21(Expr* const* p, Raise* R) {
    const Expr *a = p[0], *b = p[1], *c = p[2], *z = p[3];
    Expr* omz = sf_minus(sf_int(1), sf_cp(z));
    R->alpha = sf_div(sf_times(sf_times(sf_cp(a), sf_cp(b)), sf_cp(z)), sf_cp(omz));
    R->beta = sf_div(sf_plus(sf_minus(sf_times(sf_cp(z), sf_plus(sf_cp(a), sf_cp(b))), sf_cp(c)), sf_int(1)), sf_cp(omz));
    R->nd = 3;
    R->r0[0] = sf_int(1); R->r1[0] = sf_pow(sf_cp(a), sf_int(-1));
    R->r0[1] = sf_int(1); R->r1[1] = sf_pow(sf_cp(b), sf_int(-1));
    Expr* den = sf_times(sf_minus(sf_cp(c), sf_cp(a)), sf_minus(sf_cp(c), sf_cp(b)));
    R->r0[2] = sf_div(sf_times(sf_cp(c), sf_minus(sf_minus(sf_cp(c), sf_cp(a)), sf_cp(b))), sf_cp(den));
    R->r1[2] = sf_div(sf_times(sf_cp(c), omz), sf_times(den, sf_cp(z)));
}

static bool up_f21(Expr* const* a, int d, Expr** M) {
    if (d < 0 || d > 2) return false;
    Expr** s = args_bumped(a, 4, (size_t)d);
    bool ok = up_from_raise_sym(raise_f21, a, s, 4, 3, d, M);
    args_free(s, 4);
    return ok;
}

/* ---- 1F1 ------------------------------------------------------------------ */

static void raise_f11(Expr* const* p, Raise* R) {
    const Expr *a = p[0], *b = p[1], *z = p[2];
    R->alpha = sf_times(sf_cp(a), sf_cp(z));
    R->beta = sf_plus(sf_minus(sf_cp(z), sf_cp(b)), sf_int(1));
    R->nd = 2;
    R->r0[0] = sf_int(1); R->r1[0] = sf_pow(sf_cp(a), sf_int(-1));
    Expr* bma = sf_minus(sf_cp(b), sf_cp(a));
    R->r0[1] = sf_div(sf_cp(b), sf_cp(bma));
    R->r1[1] = sf_times(sf_int(-1), sf_div(sf_cp(b), sf_times(bma, sf_cp(z))));
}

static bool up_f11(Expr* const* a, int d, Expr** M) {
    if (d < 0 || d > 1) return false;
    Expr** s = args_bumped(a, 3, (size_t)d);
    bool ok = up_from_raise_sym(raise_f11, a, s, 3, 2, d, M);
    args_free(s, 3);
    return ok;
}

/* ---- 0F1 ------------------------------------------------------------------ */

static void raise_f01(Expr* const* p, Raise* R) {
    const Expr *b = p[0], *z = p[1];
    R->alpha = sf_cp(z);
    R->beta = sf_minus(sf_int(1), sf_cp(b));
    R->nd = 1;
    R->r0[0] = sf_int(0); R->r1[0] = sf_div(sf_cp(b), sf_cp(z));
}

static bool up_f01(Expr* const* a, int d, Expr** M) {
    if (d != 0) return false;
    Expr** s = args_bumped(a, 2, 0);
    bool ok = up_from_raise_sym(raise_f01, a, s, 2, 1, d, M);
    args_free(s, 2);
    return ok;
}

/* ---- U: three-term in a ----------------------------------------------------- */

static bool up_u(Expr* const* p, int d, Expr** M) {
    if (d != 0) return false;
    const Expr *a = p[0], *b = p[1], *z = p[2];
    /* U(a+2) = -[U(a) + (b - 2(a+1) - z) U(a+1)] / ((a+1)(a - b + 2)) */
    Expr* den = sf_times(sf_plus(sf_cp(a), sf_int(1)), sf_plus(sf_minus(sf_cp(a), sf_cp(b)), sf_int(2)));
    Expr* A = sf_times(sf_int(-1), sf_pow(sf_cp(den), sf_int(-1)));
    Expr* B = sf_times(sf_int(-1), sf_div(sf_minus(sf_minus(sf_cp(b), sf_times(sf_int(2), sf_plus(sf_cp(a), sf_int(1)))), sf_cp(z)), den));
    M[0] = sf_int(0); M[1] = sf_int(1); M[2] = sf_int(0);
    M[3] = A;         M[4] = B;         M[5] = sf_int(0);
    return true;
}

/* ---- canonical parameter order ----------------------------------------------- */

/* Class key of a parameter for ordering: symbolic part, then den, then
 * class; ties broken by the rational offset, then structurally. */
static int param_cmp(const Expr* x, const Expr* y) {
    Expr *sx = NULL, *sy = NULL; int64_t nx, dx, ny, dy;
    bool okx = sf_split(x, &sx, &nx, &dx), oky = sf_split(y, &sy, &ny, &dy);
    int r = 0;
    if (okx && oky) {
        if (!sx && sy) r = -1;
        else if (sx && !sy) r = 1;
        else if (sx && sy) r = expr_compare(sx, sy);
        if (r == 0 && dx != dy) r = dx < dy ? -1 : 1;
        if (r == 0) {
            int64_t cx = ((nx % dx) + dx) % dx, cy = ((ny % dy) + dy) % dy;
            if (cx != cy) r = cx < cy ? -1 : 1;
        }
        if (r == 0 && nx != ny) r = nx < ny ? -1 : 1;
    }
    if (r == 0) r = expr_compare(x, y);
    if (sx) expr_free(sx);
    if (sy) expr_free(sy);
    return r;
}

static void sort_params(Expr** v, size_t n) {
    for (size_t i = 1; i < n; i++)                         /* lists are short */
        for (size_t j = i; j > 0 && param_cmp(v[j - 1], v[j]) > 0; j--) {
            Expr* t = v[j]; v[j] = v[j - 1]; v[j - 1] = t;
        }
}

static bool is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

static Expr* canon_pfq(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)root; (void)ctx;
    const Expr* up = atom->data.function.args[0];
    const Expr* lo = atom->data.function.args[1];
    const Expr* z = atom->data.function.args[2];
    if (!is_list(up) || !is_list(lo)) return NULL;
    size_t p = up->data.function.arg_count, q = lo->data.function.arg_count;
    Expr** u = malloc(sizeof(Expr*) * (p ? p : 1));
    Expr** l = malloc(sizeof(Expr*) * (q ? q : 1));
    for (size_t i = 0; i < p; i++) u[i] = sf_cp(up->data.function.args[i]);
    for (size_t i = 0; i < q; i++) l[i] = sf_cp(lo->data.function.args[i]);
    sort_params(u, p);
    sort_params(l, q);
    Expr* out;
    if (p == 2 && q == 1)
        out = sf_call(SYM_SfF21, (Expr*[]){ u[0], u[1], l[0], sf_cp(z) }, 4);
    else if (p == 1 && q == 1)
        out = sf_call(SYM_SfF11, (Expr*[]){ u[0], l[0], sf_cp(z) }, 3);
    else if (p == 0 && q == 1)
        out = sf_call(SYM_SfF01, (Expr*[]){ l[0], sf_cp(z) }, 2);
    else {
        Expr* ul = sf_call(SYM_List, u, p);
        Expr* ll = sf_call(SYM_List, l, q);
        out = sf_call(SYM_HypergeometricPFQ, (Expr*[]){ ul, ll, sf_cp(z) }, 3);
        if (expr_eq(out, atom)) { expr_free(out); out = NULL; }
        free(u); free(l);
        return out;
    }
    free(u); free(l);
    return out;
}

/* Hypergeometric0F1 / 1F1 / 2F1: the PFQ call they stand for, canonicalised
 * as PFQ is (into the private lattice heads). */
static Expr* canon_named(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    Expr* pfq = hyper_named_to_pfq(atom);
    if (!pfq) return NULL;
    Expr* out = canon_pfq(pfq, root, ctx);
    if (!out) return pfq;
    expr_free(pfq);
    return out;
}

/* ---- transformations (partner present) --------------------------------------- */

static bool contains_expr(const Expr* e, const Expr* sub) {
    if (expr_eq(e, sub)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_expr(e->data.function.args[i], sub)) return true;
    return false;
}

/* Is the partner `cand` (an F21/F11/U atom, built unevaluated) in root,
 * after putting it in canonical parameter order? */
static Expr* canonical_atom(Expr* cand) {
    Expr* ev = eval_and_free(cand);
    return ev;
}

static Expr* canon_f11(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    const Expr *a = atom->data.function.args[0], *b = atom->data.function.args[1], *z = atom->data.function.args[2];
    /* Kummer: rewrite the side whose argument reads negative */
    if (!sf_reads_negative(z)) return NULL;
    Expr* bma = sf_mk1("Expand", sf_minus(sf_cp(b), sf_cp(a)));
    Expr* mz = sf_mk1("Expand", sf_times(sf_int(-1), sf_cp(z)));
    Expr* partner = canonical_atom(sf_call(SYM_SfF11, (Expr*[]){ sf_cp(bma), sf_cp(b), sf_cp(mz) }, 3));
    bool present = contains_expr(root, partner);
    expr_free(partner);
    if (!present) { expr_free(bma); expr_free(mz); return NULL; }
    /* 1F1(a; b; z) = E^z 1F1(b - a; b; -z) */
    return sf_times(sf_pow(sf_sym(SYM_E), sf_cp(z)), sf_call(SYM_SfF11, (Expr*[]){ bma, sf_cp(b), mz }, 3));
}

/* 2F1 Euler (same z, partner F21[c-a, c-b, c, z]) for z < 1 proved. */
static bool proved_lt1(const AssumeCtx* ctx, const Expr* z) {
    if (!ctx) return false;
    Expr* pred = sf_call(SYM_Less, (Expr*[]){ sf_cp(z), sf_int(1) }, 2);
    int r = assume_reduce_entails(ctx, pred);
    expr_free(pred);
    return r == 1;
}

static Expr* f21_atom(Expr* a, Expr* b, Expr* c, Expr* z) {
    if (param_cmp(a, b) > 0) { Expr* t = a; a = b; b = t; }
    return sf_call(SYM_SfF21, (Expr*[]){ a, b, c, z }, 4);
}

static bool rewrite_this_side(const Expr* mine, const Expr* partner) {
    size_t x = score_with_func(mine, NULL), y = score_with_func(partner, NULL);
    if (x != y) return x > y;
    return expr_compare(mine, partner) > 0;
}

static Expr* canon_f21(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    const Expr *a = atom->data.function.args[0], *b = atom->data.function.args[1];
    const Expr *c = atom->data.function.args[2], *z = atom->data.function.args[3];
    if (!proved_lt1(ctx, z)) return NULL;
    /* Euler: partner F21[c - a, c - b, c, z] */
    Expr* cma = sf_mk1("Expand", sf_minus(sf_cp(c), sf_cp(a)));
    Expr* cmb = sf_mk1("Expand", sf_minus(sf_cp(c), sf_cp(b)));
    Expr* partner = f21_atom(sf_cp(cma), sf_cp(cmb), sf_cp(c), sf_cp(z));
    bool present = contains_expr(root, partner) && rewrite_this_side(atom, partner);
    if (present) {
        Expr* pre = sf_pow(sf_minus(sf_int(1), sf_cp(z)),
                           sf_mk1("Expand", sf_minus(sf_minus(sf_cp(c), sf_cp(a)), sf_cp(b))));
        expr_free(cma); expr_free(cmb);
        return sf_times(pre, partner);
    }
    expr_free(partner);
    expr_free(cma); expr_free(cmb);
    /* Pfaff: partner F21[a, c - b, c, z/(z - 1)] (either upper parameter) */
    for (int side = 0; side < 2; side++) {
        const Expr* keep = side ? b : a;
        const Expr* other = side ? a : b;
        Expr* w = eval_and_free(sf_call("Together", (Expr*[]){ sf_div(sf_cp(z), sf_minus(sf_cp(z), sf_int(1))) }, 1));
        Expr* cmo = sf_mk1("Expand", sf_minus(sf_cp(c), sf_cp(other)));
        Expr* part = f21_atom(sf_cp(keep), cmo, sf_cp(c), w);
        if (contains_expr(root, part) && rewrite_this_side(atom, part)) {
            Expr* pre = sf_pow(sf_minus(sf_int(1), sf_cp(z)), sf_times(sf_int(-1), sf_cp(keep)));
            return sf_times(pre, part);
        }
        expr_free(part);
    }
    return NULL;
}

/* Kummer for U: U(a, b, z) = z^(1-b) U(a - b + 1, 2 - b, z), partner present */
static Expr* canon_u(const Expr* atom, const Expr* root, const AssumeCtx* ctx) {
    (void)ctx;
    const Expr *a = atom->data.function.args[0], *b = atom->data.function.args[1], *z = atom->data.function.args[2];
    Expr* a2 = sf_mk1("Expand", sf_plus(sf_minus(sf_cp(a), sf_cp(b)), sf_int(1)));
    Expr* b2 = sf_mk1("Expand", sf_minus(sf_int(2), sf_cp(b)));
    Expr* partner = sf_call(SYM_HypergeometricU, (Expr*[]){ a2, b2, sf_cp(z) }, 3);
    if (contains_expr(root, partner) && rewrite_this_side(atom, partner)) {
        Expr* pre = sf_pow(sf_cp(z), sf_minus(sf_int(1), sf_cp(b)));
        return sf_times(pre, partner);
    }
    expr_free(partner);
    return NULL;
}

/* ---- finalize: private heads back to the named heads ---------------------------
 * Hypergeometric2F1 / 1F1 / 0F1 evaluate through HypergeometricPFQ and keep their
 * head when it cannot reduce them -- the form the reference writes these in. */

static Expr* hyper_finalize_walk(const Expr* e, bool* changed) {
    if (e->type != EXPR_FUNCTION) return sf_cp(e);
    size_t n = e->data.function.arg_count;
    Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) a[i] = hyper_finalize_walk(e->data.function.args[i], changed);
    const Expr* h = e->data.function.head;
    const char* hn = h->type == EXPR_SYMBOL ? h->data.symbol.name : NULL;
    Expr* r;
    if (hn == SYM_SfF21 && n == 4) {
        r = sf_call(SYM_Hypergeometric2F1, (Expr*[]){ a[0], a[1], a[2], a[3] }, 4);
        *changed = true;
    } else if (hn == SYM_SfF11 && n == 3) {
        r = sf_call(SYM_Hypergeometric1F1, (Expr*[]){ a[0], a[1], a[2] }, 3);
        *changed = true;
    } else if (hn == SYM_SfF01 && n == 2) {
        r = sf_call(SYM_Hypergeometric0F1, (Expr*[]){ a[0], a[1] }, 2);
        *changed = true;
    } else {
        r = expr_new_function(sf_cp(h), a, n);
    }
    free(a);
    return r;
}

static Expr* hyper_finalize(Expr* e) {
    if (!sf_contains_head(e, SYM_SfF21) && !sf_contains_head(e, SYM_SfF11) && !sf_contains_head(e, SYM_SfF01))
        return e;
    bool changed = false;
    Expr* r = hyper_finalize_walk(e, &changed);
    expr_free(e);
    return eval_and_free(r);
}

/* ---- module table ---------------------------------------------------------- */

static const SfDesc HYPER_DESCS[] = {
    { &SYM_HypergeometricPFQ, 3, 0, 1, { 0 },       NULL,   NULL, canon_pfq, NULL, NULL },
    { &SYM_Hypergeometric2F1, 4, 0, 1, { 0 },       NULL,   NULL, canon_named, NULL, NULL },
    { &SYM_Hypergeometric1F1, 3, 0, 1, { 0 },       NULL,   NULL, canon_named, NULL, NULL },
    { &SYM_Hypergeometric0F1, 2, 0, 1, { 0 },       NULL,   NULL, canon_named, NULL, NULL },
    { &SYM_SfF21,             4, 3, 2, { 0, 1, 2 }, up_f21, NULL, canon_f21, NULL, NULL },
    { &SYM_SfF11,             3, 2, 2, { 0, 1 },    up_f11, NULL, canon_f11, NULL, NULL },
    { &SYM_SfF01,             2, 1, 2, { 0 },       up_f01, NULL, NULL,      NULL, NULL },
    { &SYM_HypergeometricU,   3, 1, 2, { 0 },       up_u,   NULL, canon_u,   NULL, NULL },
};

const SfModule SF_MODULE_HYPER = {
    HYPER_DESCS, sizeof HYPER_DESCS / sizeof HYPER_DESCS[0], NULL, 0, NULL, 0,
    hyper_init, hyper_finalize,
};
