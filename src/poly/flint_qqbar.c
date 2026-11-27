/* flint_qqbar.c — exact algebraic-number canonicalisation for RootReduce.
 *
 * See flint_qqbar.h for the contract. The engine is FLINT's `qqbar` (exact
 * real/complex algebraic numbers). A constant algebraic-number Expr is folded
 * bottom-up into a single qqbar, whose minimal polynomial, degree, rationality,
 * quadratic-radical form and root index then produce the WL-faithful canonical
 * Expr: a rational, a quadratic radical, or Root[Function[minpoly&], k].
 */

#include "flint_qqbar.h"
#include "linkedlist.h"

#include "sym_names.h"
#include "eval.h"

#include <stdlib.h>
#include <string.h>

#ifdef USE_FLINT

#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_vec.h>
#include <flint/qqbar.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>

/* The maximal-order (Round 2 / Pohst-Zassenhaus) engine and its integral-basis
 * accessors, used by flint_qqbar_integral_basis below. numberfield_internal.h
 * guards its own FLINT includes and exposes nf_ok_basis / nf_ok_denom. */
#include "numberfield.h"
#include "numberfield_internal.h"

/* Intermediate qqbar degree above which we give up and fall back (identity /
 * the parametric engine). WL's degree-21 examples stay comfortably under. */
#define QQBAR_DEGREE_CAP 120
/* Set by to_qqbar when it declines because an intermediate degree passed the
 * cap (as opposed to a non-algebraic input): such a decline is "undecided",
 * never a "not algebraic". Reset by the caller that reads it. */
static int g_qqbar_cap_hit = 0;
/* Max distinct algebraic generators collected for the NumberField method. */
#define QQBAR_MAX_ATOMS  12

/* ------------------------------------------------------------------ */
/*  Small Expr / FLINT helpers                                         */
/* ------------------------------------------------------------------ */

static int head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

static Expr* expr_from_fmpz(const fmpz_t z) {
    mpz_t m; mpz_init(m); fmpz_get_mpz(m, z);
    Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(m));
    mpz_clear(m);
    return r;
}

static Expr* expr_from_fmpq(const fmpq_t q) {
    if (fmpz_is_one(fmpq_denref(q))) return expr_from_fmpz(fmpq_numref(q));
    Expr* args[2] = { expr_from_fmpz(fmpq_numref(q)), expr_from_fmpz(fmpq_denref(q)) };
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

/* Fill an fmpz from an integer-like Expr; 1 on success. */
static int fmpz_from_intlike(const Expr* e, fmpz_t out) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, (slong)e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

/* Fill an fmpq from an integer-like or Rational[p,q] Expr; 1 on success. */
static int fmpq_from_expr(const Expr* e, fmpq_t out) {
    fmpz_t z; fmpz_init(z);
    if (fmpz_from_intlike(e, z)) {
        fmpz_set(fmpq_numref(out), z); fmpz_one(fmpq_denref(out));
        fmpz_clear(z); return 1;
    }
    fmpz_clear(z);
    if (head_is(e, "Rational") && e->data.function.arg_count == 2) {
        fmpz_t n, d; fmpz_init(n); fmpz_init(d);
        int ok = fmpz_from_intlike(e->data.function.args[0], n) &&
                 fmpz_from_intlike(e->data.function.args[1], d);
        if (ok) fmpq_set_fmpz_frac(out, n, d);
        fmpz_clear(n); fmpz_clear(d);
        return ok;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Univariate integer poly Expr (Root body) -> fmpz_poly              */
/* ------------------------------------------------------------------ */

/* Treat `Slot[1]` (slot form) or the symbol named `var` (2-arg Function form)
 * as the polynomial variable. Only integer coefficients are accepted; a
 * rational coefficient or a foreign symbol makes the build fail (return 0),
 * which just leaves the Root object unconverted. */
static int is_poly_var(const Expr* e, const char* var) {
    if (var && e->type == EXPR_SYMBOL) return strcmp(e->data.symbol.name, var) == 0;
    return head_is(e, "Slot") && e->data.function.arg_count == 1 &&
           e->data.function.args[0]->type == EXPR_INTEGER &&
           e->data.function.args[0]->data.integer == 1;
}

static int build_fmpz_poly(const Expr* e, const char* var, fmpz_poly_t out) {
    if (!e) return 0;
    fmpz_t z; fmpz_init(z);
    if (fmpz_from_intlike(e, z)) { fmpz_poly_set_fmpz(out, z); fmpz_clear(z); return 1; }
    fmpz_clear(z);
    if (is_poly_var(e, var)) {
        fmpz_poly_zero(out);                 /* out may be a reused temp */
        fmpz_poly_set_coeff_si(out, 1, 1);
        return 1;
    }
    if (e->type != EXPR_FUNCTION) return 0;

    size_t n = e->data.function.arg_count;
    if (head_is(e, "Plus")) {
        fmpz_poly_zero(out);
        fmpz_poly_t t; fmpz_poly_init(t);
        for (size_t i = 0; i < n; i++) {
            if (!build_fmpz_poly(e->data.function.args[i], var, t)) { fmpz_poly_clear(t); return 0; }
            fmpz_poly_add(out, out, t);
        }
        fmpz_poly_clear(t);
        return 1;
    }
    if (head_is(e, "Times")) {
        fmpz_poly_set_si(out, 1);
        fmpz_poly_t t; fmpz_poly_init(t);
        for (size_t i = 0; i < n; i++) {
            if (!build_fmpz_poly(e->data.function.args[i], var, t)) { fmpz_poly_clear(t); return 0; }
            fmpz_poly_mul(out, out, t);
        }
        fmpz_poly_clear(t);
        return 1;
    }
    if (head_is(e, "Power") && n == 2) {
        const Expr* base = e->data.function.args[0];
        const Expr* exp  = e->data.function.args[1];
        if (exp->type != EXPR_INTEGER || exp->data.integer < 0) return 0;
        fmpz_poly_t b; fmpz_poly_init(b);
        if (!build_fmpz_poly(base, var, b)) { fmpz_poly_clear(b); return 0; }
        fmpz_poly_pow(out, b, (ulong)exp->data.integer);
        fmpz_poly_clear(b);
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  WL-faithful root ordering                                          */
/* ------------------------------------------------------------------ */

/* Root ordering, matching Mathilda's own Root[] canonical order EXACTLY (see
 * root_canonical_cmp_mpfr in root_numeric.c and the Root docstring): real roots
 * first, ascending by real part; then non-real roots by real part ascending,
 * then |Im| ascending, then the more-negative-Im member of a conjugate pair
 * first.  This MUST agree with the numeric Root builtin: RootReduce emits
 * Root[minpoly, k] with the index k this order assigns to the algebraic value,
 * and `N[Root[minpoly, k]]` re-derives the k-th root with the mpfr order above.
 * The previous code ordered non-real roots by FLINT's qqbar_cmp_root_order,
 * which disagrees for complex roots — so RootReduce of a degree>=3 (or
 * quartic-via-quadratic) complex algebraic number emitted a Root pointing at the
 * WRONG conjugate/sibling root (e.g. RootReduce[Sqrt[-1+I]] came back as its own
 * negative). */
static int wl_cmp_qq(const qqbar_t a, const qqbar_t b) {
    int ra = qqbar_is_real(a), rb = qqbar_is_real(b);
    if (ra != rb) return ra ? -1 : 1;
    if (ra) return qqbar_cmp_re(a, b);
    /* both non-real.  A conjugate pair has equal Re and |Im| -- a cheap test
     * (same minimal polynomial, conjugate enclosure) that spares qqbar_cmp_re
     * an exact proof of the equal real parts: the more-negative Im first. */
    {
        qqbar_t cb; qqbar_init(cb);
        qqbar_conj(cb, b);
        int conj = qqbar_equal(a, cb);
        qqbar_clear(cb);
        if (conj) {
            int sa = qqbar_sgn_im(a), sb = qqbar_sgn_im(b);
            return sa < sb ? -1 : (sa > sb ? 1 : 0);
        }
    }
    /* Re ascending */
    int c = qqbar_cmp_re(a, b);
    if (c != 0) return c;
    /* |Im| ascending.  Decided first on rigorous enclosures at rising
     * precision: two roots with equal real parts have equal |Im| only as a
     * conjugate pair, and the exact qqbar_im (a composed-polynomial
     * factorisation) inside this O(d^2) sort made the Root index of 2^(1/32)
     * cost minutes (audit NumberFieldIntegralBasis#1).  The exact path below
     * remains for the pairs the enclosures cannot separate. */
    for (slong prec = 64; prec <= 1024; prec *= 2) {
        acb_t za, zb; acb_init(za); acb_init(zb);
        arb_t ma, mb; arb_init(ma); arb_init(mb);
        qqbar_get_acb(za, a, prec);
        qqbar_get_acb(zb, b, prec);
        arb_abs(ma, acb_imagref(za));
        arb_abs(mb, acb_imagref(zb));
        int r = arb_lt(ma, mb) ? -1 : (arb_gt(ma, mb) ? 1 : 0);
        arb_clear(ma); arb_clear(mb); acb_clear(za); acb_clear(zb);
        if (r != 0) return r;
    }
    qqbar_t ia, ib; qqbar_init(ia); qqbar_init(ib);
    qqbar_im(ia, a); qqbar_im(ib, b);
    qqbar_abs(ia, ia); qqbar_abs(ib, ib);       /* real qqbars: |Im| */
    c = qqbar_cmp_re(ia, ib);                    /* compares their (real) values */
    qqbar_clear(ia); qqbar_clear(ib);
    if (c != 0) return c;
    /* tie: more-negative Im first */
    int sa = qqbar_sgn_im(a), sb = qqbar_sgn_im(b);
    if (sa < sb) return -1;
    if (sa > sb) return  1;
    return 0;
}

/* Selection-sort an index permutation of `roots` into WL order (n <= cap). */
static void wl_sort_indices(qqbar_srcptr roots, slong n, slong* idx) {
    for (slong i = 0; i < n; i++) idx[i] = i;
    for (slong i = 0; i < n; i++) {
        slong best = i;
        for (slong j = i + 1; j < n; j++)
            if (wl_cmp_qq(roots + idx[j], roots + idx[best]) < 0) best = j;
        slong t = idx[i]; idx[i] = idx[best]; idx[best] = t;
    }
}

/* ---- wl_root_index: sorted-roots memo, keyed by minimal polynomial --------
 * wl_root_index isolates a minpoly's d roots and O(d^2)-sorts them (each
 * comparison an arb-ball qqbar_cmp_re + a qqbar_im that factors a polynomial)
 * on EVERY call -- and qqbar_to_expr calls it once per degree>=3 conversion,
 * i.e. after every AlgebraicNumber Plus/Times in a field ansatz.  The WL sort
 * order depends only on the minpoly, so cache the sorted roots per minpoly and
 * reduce a repeat call to a linear qqbar_equal scan.  This collapses the split-
 * specials ansatz's thousands of full isolations+sorts over ONE field to one.
 * The map is a pure function of the minpoly (session-independent, never stale);
 * bounded, reset when it approaches full, and leaked at exit like the qqbar
 * cache above. */
#define WLROOTS_CACHE_SIZE 128u          /* power of two; open addressing */
typedef struct { fmpz_poly_t poly; qqbar_ptr roots; slong d; uint64_t h; int used; } WLRootsSlot;
static WLRootsSlot g_wlroots[WLROOTS_CACHE_SIZE];
static size_t g_wlroots_count = 0;

static uint64_t fmpz_poly_hashval(const fmpz_poly_struct* p) {
    uint64_t h = 1469598103934665603ULL;          /* FNV-1a offset basis */
    slong n = fmpz_poly_length(p);
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i < n; i++) {
        fmpz_poly_get_coeff_fmpz(c, p, i);
        uint64_t v = (uint64_t)fmpz_fdiv_ui(c, 4294967291u);   /* coeff mod a prime */
        if (fmpz_sgn(c) < 0) v ^= 0x9e3779b97f4a7c15ULL;       /* separate sign */
        h = (h ^ v) * 1099511628211ULL;
    }
    fmpz_clear(c);
    return h;
}

static void wlroots_reset(void) {
    for (unsigned i = 0; i < WLROOTS_CACHE_SIZE; i++) {
        if (g_wlroots[i].used) {
            fmpz_poly_clear(g_wlroots[i].poly);
            _qqbar_vec_clear(g_wlroots[i].roots, g_wlroots[i].d);
            g_wlroots[i].used = 0;
        }
    }
    g_wlroots_count = 0;
}

/* Return the cached roots of minpoly `p` (degree d) in WL sort order, isolating
 * and sorting on a miss; NULL if the table could not cache it (fall back). */
static qqbar_srcptr wlroots_get(const fmpz_poly_struct* p, slong d, uint64_t h) {
    unsigned mask = WLROOTS_CACHE_SIZE - 1u;
    unsigned i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < WLROOTS_CACHE_SIZE; probe++) {
        WLRootsSlot* s = &g_wlroots[i];
        if (!s->used) break;                                   /* empty => miss */
        if (s->h == h && s->d == d && fmpz_poly_equal(s->poly, p)) return s->roots;
        i = (i + 1u) & mask;
    }
    if (g_wlroots_count * 4u >= WLROOTS_CACHE_SIZE * 3u) wlroots_reset();   /* >75% */
    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, p, QQBAR_ROOTS_IRREDUCIBLE);
    slong* idx = malloc(sizeof(slong) * (size_t)d);
    wl_sort_indices(roots, d, idx);
    qqbar_ptr sorted = _qqbar_vec_init(d);
    for (slong j = 0; j < d; j++) qqbar_set(sorted + j, roots + idx[j]);
    free(idx);
    _qqbar_vec_clear(roots, d);
    i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < WLROOTS_CACHE_SIZE; probe++) {
        WLRootsSlot* s = &g_wlroots[i];
        if (!s->used) {
            fmpz_poly_init(s->poly); fmpz_poly_set(s->poly, p);
            s->roots = sorted; s->d = d; s->h = h; s->used = 1;   /* publish last */
            g_wlroots_count++;
            return s->roots;
        }
        i = (i + 1u) & mask;
    }
    _qqbar_vec_clear(sorted, d);   /* table full even after reset: cannot cache */
    return NULL;
}

/* 1-based WL index of x among the roots of its minimal polynomial. */
static slong wl_root_index(const qqbar_t x) {
    slong d = qqbar_degree(x);
    if (d <= 0) return 1;
    const fmpz_poly_struct* P = QQBAR_POLY(x);
    qqbar_srcptr sorted = wlroots_get(P, d, fmpz_poly_hashval(P));
    if (sorted) {
        for (slong j = 0; j < d; j++)
            if (qqbar_equal(sorted + j, x)) return j + 1;
        return 1;                                  /* x must be one of them */
    }
    /* fallback (cache unavailable): the original isolate-and-sort per call */
    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, P, QQBAR_ROOTS_IRREDUCIBLE);
    slong* idx = malloc(sizeof(slong) * (size_t)d);
    wl_sort_indices(roots, d, idx);
    slong k = 1;
    for (slong j = 0; j < d; j++)
        if (qqbar_equal(roots + idx[j], x)) { k = j + 1; break; }
    free(idx);
    _qqbar_vec_clear(roots, d);
    return k;
}

/* ------------------------------------------------------------------ */
/*  Expr -> qqbar (recursive field arithmetic)                         */
/* ------------------------------------------------------------------ */

/* True iff e is a negative Rational[p, q] (the exponent of a fractional power). */
static int is_negative_rational_expr(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2) return 0;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Rational) return 0;
    const Expr* p = e->data.function.args[0];
    if (p->type == EXPR_INTEGER) return p->data.integer < 0;
    if (p->type == EXPR_BIGINT) return mpz_sgn(p->data.bigint) < 0;
    return 0;
}

static int to_qqbar(const Expr* e, qqbar_t out);

/* Root[Function[...], k] -> the k-th root (WL order) of the body polynomial. */
/* ------------------------------------------------------------------ */
/*  Memo cache: Root[...] Expr -> qqbar.                                */
/*                                                                      */
/*  A Root object's qqbar value is a pure, session-independent constant */
/*  (fixed by its minimal polynomial and root index), so this map never */
/*  goes stale and needs no invalidation.  Isolating a Root's roots     */
/*  (qqbar_roots_fmpz_poly) is the dominant cost of the algebraic-number */
/*  engine whenever the same handful of generators are re-converted many */
/*  times -- RootReduce zero-tests, AlgebraicNumber add/mul/pow (each    */
/*  converts its generator), NullSpace/RowReduce ZeroTests -- so caching */
/*  it collapses hundreds of isolations to one per distinct generator.   */
/* ------------------------------------------------------------------ */
#define QQBAR_CACHE_SIZE 4096u          /* power of two; open addressing */
typedef struct { Expr* key; uint64_t h; qqbar_t val; int used; } QQBarCacheSlot;
static QQBarCacheSlot g_qqbar_cache[QQBAR_CACHE_SIZE];
static size_t g_qqbar_cache_count = 0;

static void flint_qqbar_cache_reset(void) {
    for (unsigned i = 0; i < QQBAR_CACHE_SIZE; i++) {
        if (g_qqbar_cache[i].used) {
            expr_free(g_qqbar_cache[i].key);
            qqbar_clear(g_qqbar_cache[i].val);
            g_qqbar_cache[i].key = NULL;
            g_qqbar_cache[i].used = 0;
        }
    }
    g_qqbar_cache_count = 0;
}

/* On hit, copy the cached value into `out` and return 1; else 0. */
static int qqbar_cache_get(const Expr* e, uint64_t h, qqbar_t out) {
    unsigned mask = QQBAR_CACHE_SIZE - 1u;
    unsigned i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < QQBAR_CACHE_SIZE; probe++) {
        QQBarCacheSlot* s = &g_qqbar_cache[i];
        if (!s->used) return 0;                     /* empty slot => absent */
        if (s->h == h && expr_eq(s->key, e)) { qqbar_set(out, s->val); return 1; }
        i = (i + 1u) & mask;
    }
    return 0;
}

/* Insert (e -> val); bounded, resetting the table when it approaches full
 * (distinct generators per session are few).  These are raw FLINT/Expr ops
 * with no evaluation checkpoint, so no async TimeConstrained longjmp can
 * interrupt mid-insert; `used` is still published last, defensively. */
static void qqbar_cache_put(const Expr* e, uint64_t h, const qqbar_t val) {
    if (g_qqbar_cache_count * 4u >= QQBAR_CACHE_SIZE * 3u)   /* >75% full */
        flint_qqbar_cache_reset();
    unsigned mask = QQBAR_CACHE_SIZE - 1u;
    unsigned i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < QQBAR_CACHE_SIZE; probe++) {
        QQBarCacheSlot* s = &g_qqbar_cache[i];
        if (!s->used) {
            s->key = expr_copy((Expr*)e);           /* refcount bump keeps it alive */
            s->h = h;
            qqbar_init(s->val);
            qqbar_set(s->val, val);
            s->used = 1;                             /* publish last */
            g_qqbar_cache_count++;
            return;
        }
        if (s->h == h && expr_eq(s->key, e)) return;   /* already present */
        i = (i + 1u) & mask;
    }
}

static int root_object_to_qqbar(const Expr* e, qqbar_t out) {
    uint64_t h = expr_hash(e);
    if (qqbar_cache_get(e, h, out)) return 1;
    size_t n = e->data.function.arg_count;
    if (n < 2) return 0;
    const Expr* fn = e->data.function.args[0];
    const Expr* ke = e->data.function.args[1];
    if (!head_is(fn, "Function") || ke->type != EXPR_INTEGER) return 0;

    const char* var = NULL;
    const Expr* body;
    if (fn->data.function.arg_count == 2 &&
        fn->data.function.args[0]->type == EXPR_SYMBOL) {
        var  = fn->data.function.args[0]->data.symbol.name;   /* Function[t, body] */
        body = fn->data.function.args[1];
    } else if (fn->data.function.arg_count == 1) {
        body = fn->data.function.args[0];                /* Function[body(Slot[1])] */
    } else {
        return 0;
    }

    fmpz_poly_t P; fmpz_poly_init(P);
    if (!build_fmpz_poly(body, var, P)) { fmpz_poly_clear(P); return 0; }
    slong d = fmpz_poly_degree(P);
    slong k = (slong)ke->data.integer;
    if (d < 1 || k < 1 || k > d) { fmpz_poly_clear(P); return 0; }

    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, P, 0);
    slong* idx = malloc(sizeof(slong) * (size_t)d);
    wl_sort_indices(roots, d, idx);
    qqbar_set(out, roots + idx[k - 1]);
    free(idx);
    _qqbar_vec_clear(roots, d);
    fmpz_poly_clear(P);
    qqbar_cache_put(e, h, out);
    return 1;
}

/* If `exp` has the form (I Pi r) for a rational r — the canonical
 * Times[Complex[0, r], Pi] — set `out` to r and return 1; else return 0. Used to
 * recognise E^(I Pi r) = exp(i pi r) as a root of unity. */
static int exp_arg_i_pi_rational(const Expr* exp, fmpq_t out) {
    if (!head_is(exp, "Times") || exp->data.function.arg_count != 2) return 0;
    const Expr* a0 = exp->data.function.args[0];
    const Expr* a1 = exp->data.function.args[1];
    const Expr* co = NULL;
    if (a0->type == EXPR_SYMBOL && a0->data.symbol.name == SYM_Pi)      co = a1;
    else if (a1->type == EXPR_SYMBOL && a1->data.symbol.name == SYM_Pi) co = a0;
    if (!co || !head_is(co, "Complex") || co->data.function.arg_count != 2) return 0;
    const Expr* reP = co->data.function.args[0];
    if (!(reP->type == EXPR_INTEGER && reP->data.integer == 0)) return 0;   /* purely imaginary */
    return fmpq_from_expr(co->data.function.args[1], out);
}

/* If `arg` is a rational multiple r of Pi — Pi itself (r = 1) or
 * Times[rational..., Pi] — set *out = r and return 1; else 0.  Used to
 * evaluate Sin/Cos/Tan/Cot/Sec/Csc at a rational multiple of Pi through the
 * root of unity zeta = exp(i pi r), extending the exact closed-form coverage
 * past the src/trig.c tables (which stop at denominator 12) to every rational
 * angle: Sin[Pi/7], Cos[3 Pi/7], Tan[Pi/14], ... become algebraic. */
static int trig_arg_pi_rational(const Expr* arg, fmpq_t out) {
    if (!arg) return 0;
    if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_Pi) {
        fmpq_set_si(out, 1, 1);
        return 1;
    }
    if (!head_is(arg, "Times")) return 0;
    size_t n = arg->data.function.arg_count;
    int pi_seen = 0, ok = 1;
    fmpq_t coeff, c; fmpq_init(coeff); fmpq_init(c);
    fmpq_set_si(coeff, 1, 1);
    for (size_t i = 0; i < n; i++) {
        const Expr* f = arg->data.function.args[i];
        if (f->type == EXPR_SYMBOL && f->data.symbol.name == SYM_Pi) {
            if (pi_seen) { ok = 0; break; }   /* Pi^2 etc. is not a rational angle */
            pi_seen = 1;
            continue;
        }
        if (!fmpq_from_expr(f, c)) { ok = 0; break; }   /* non-rational factor */
        fmpq_mul(coeff, coeff, c);
    }
    if (ok && pi_seen) fmpq_set(out, coeff);
    else ok = 0;
    fmpq_clear(coeff); fmpq_clear(c);
    return ok;
}

/* Is `e` a circular trig head Sin/Cos/Tan/Cot/Sec/Csc of a single argument? */
static int is_circular_trig_head(const Expr* e, size_t n) {
    return n == 1 && (head_is(e, "Sin") || head_is(e, "Cos") ||
                      head_is(e, "Tan") || head_is(e, "Cot") ||
                      head_is(e, "Sec") || head_is(e, "Csc"));
}

static int to_qqbar(const Expr* e, qqbar_t out) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) { qqbar_set_si(out, (slong)e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT) {
        fmpz_t z; fmpz_init(z); fmpz_set_mpz(z, e->data.bigint);
        qqbar_set_fmpz(out, z); fmpz_clear(z); return 1;
    }
    if (e->type == EXPR_SYMBOL) {
        /* GoldenRatio = (1 + Sqrt[5])/2 is the sole algebraic named constant;
         * every other bare symbol (Pi, E, EulerGamma, ...) is non-algebraic. */
        if (e->data.symbol.name == SYM_GoldenRatio) {
            qqbar_sqrt_ui(out, 5);        /* Sqrt[5]        */
            qqbar_add_ui(out, out, 1);    /* 1 + Sqrt[5]    */
            qqbar_div_ui(out, out, 2);    /* (1 + Sqrt[5])/2 */
            return 1;
        }
        return 0;
    }
    if (e->type != EXPR_FUNCTION) return 0;   /* real / string */

    size_t n = e->data.function.arg_count;

    if (head_is(e, "Rational") && n == 2) {
        fmpq_t q; fmpq_init(q);
        int ok = fmpq_from_expr(e, q);
        if (ok) qqbar_set_fmpq(out, q);
        fmpq_clear(q);
        return ok;
    }
    if (head_is(e, "Complex") && n == 2) {
        qqbar_t re, im, ii; qqbar_init(re); qqbar_init(im); qqbar_init(ii);
        int ok = to_qqbar(e->data.function.args[0], re) &&
                 to_qqbar(e->data.function.args[1], im);
        if (ok) { qqbar_i(ii); qqbar_mul(im, im, ii); qqbar_add(out, re, im); }
        qqbar_clear(re); qqbar_clear(im); qqbar_clear(ii);
        return ok;
    }
    if (head_is(e, "Plus") || head_is(e, "Times")) {
        int is_plus = head_is(e, "Plus");
        qqbar_t acc, t; qqbar_init(acc); qqbar_init(t);
        if (is_plus) qqbar_set_si(acc, 0); else qqbar_set_si(acc, 1);
        int ok = 1;
        for (size_t i = 0; i < n && ok; i++) {
            if (!to_qqbar(e->data.function.args[i], t)) { ok = 0; break; }
            if (is_plus) qqbar_add(acc, acc, t); else qqbar_mul(acc, acc, t);
            if (qqbar_degree(acc) > QQBAR_DEGREE_CAP) { ok = 0; g_qqbar_cap_hit = 1; break; }
        }
        if (ok) qqbar_set(out, acc);
        qqbar_clear(acc); qqbar_clear(t);
        return ok;
    }
    if (head_is(e, "Power") && n == 2) {
        const Expr* be = e->data.function.args[0];
        const Expr* xe = e->data.function.args[1];
        /* E^(I Pi r), r rational, is the root of unity exp(i pi r) — algebraic. */
        if (be->type == EXPR_SYMBOL && be->data.symbol.name == SYM_E) {
            fmpq_t r; fmpq_init(r);
            int ok = 0;
            if (exp_arg_i_pi_rational(xe, r)) {
                const fmpz* pn = fmpq_numref(r);
                const fmpz* pd = fmpq_denref(r);
                if (fmpz_fits_si(pn) && fmpz_fits_si(pd)) {
                    qqbar_exp_pi_i(out, fmpz_get_si(pn), (ulong)fmpz_get_si(pd));
                    ok = 1;
                }
            }
            fmpq_clear(r);
            return ok;   /* E with any other exponent is not algebraic */
        }
        qqbar_t base; qqbar_init(base);
        int ok = to_qqbar(be, base);
        /* A base that is zero -- possibly only after simplification, like
         * Sqrt[2]^2 - 2 -- raised to a negative power is not an algebraic
         * number (ComplexInfinity): decline before FLINT throws
         * "qqbar_pow_fmpq: division by zero". */
        if (ok && qqbar_is_zero(base) &&
            ((xe->type == EXPR_INTEGER && xe->data.integer < 0) ||
             (xe->type != EXPR_INTEGER && is_negative_rational_expr(xe))))
            ok = 0;
        if (ok) {
            if (xe->type == EXPR_INTEGER) {
                qqbar_pow_si(out, base, (slong)xe->data.integer);
            } else {
                fmpq_t ef; fmpq_init(ef);
                if (fmpq_from_expr(xe, ef)) {
                    /* x^(p/q) as (PRINCIPAL q-th root of x)^p.  FLINT's
                     * qqbar_pow_fmpq does NOT use the principal analytic branch
                     * for a fractional power -- it can return a different root of
                     * the value's minimal polynomial, so RootReduce would emit an
                     * algebraic number NOT equal to its input (e.g.
                     * Sqrt[-1/2 - I Sqrt[3]/2] came back as its own negative,
                     * -0.5+0.866 I instead of 0.5-0.866 I).  qqbar_root_ui IS the
                     * principal n-th root (argument in (-pi/n, pi/n]), matching
                     * Mathematica's Power convention and Mathilda's numeric N[].
                     * Raising by the integer numerator p afterwards is branch-
                     * unambiguous, so (x^(1/q))^p = the principal x^(p/q). */
                    const fmpz* nump = fmpq_numref(ef);
                    const fmpz* denp = fmpq_denref(ef);
                    if (fmpz_fits_si(nump) && fmpz_fits_si(denp)) {
                        ulong qden = (ulong)fmpz_get_si(denp);   /* denom > 0 */
                        slong pnum = fmpz_get_si(nump);
                        qqbar_root_ui(out, base, qden);          /* principal root */
                        qqbar_pow_si(out, out, pnum);
                    } else ok = 0;
                } else ok = 0;
                fmpq_clear(ef);
            }
        }
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) { ok = 0; g_qqbar_cap_hit = 1; }
        qqbar_clear(base);
        return ok;
    }
    if (head_is(e, "Sqrt") && n == 1) {
        qqbar_t base; qqbar_init(base);
        int ok = to_qqbar(e->data.function.args[0], base);
        if (ok) qqbar_sqrt(out, base);
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) { ok = 0; g_qqbar_cap_hit = 1; }
        qqbar_clear(base);
        return ok;
    }
    if ((head_is(e, "Re") || head_is(e, "Im") || head_is(e, "Abs") ||
         head_is(e, "Conjugate")) && n == 1) {
        /* Re/Im/Abs/Conjugate of a constant algebraic number is algebraic. */
        qqbar_t z; qqbar_init(z);
        int ok = to_qqbar(e->data.function.args[0], z);
        if (ok) {
            if      (head_is(e, "Re"))  qqbar_re(out, z);
            else if (head_is(e, "Im"))  qqbar_im(out, z);
            else if (head_is(e, "Abs")) qqbar_abs(out, z);
            else                        qqbar_conj(out, z);
        }
        qqbar_clear(z);
        return ok;
    }
    if (head_is(e, "AlgebraicNumber") && n == 2) {
        /* AlgebraicNumber[gen, {c0..cm}] = sum ci gen^i, via Horner. */
        const Expr* gen = e->data.function.args[0];
        const Expr* cl  = e->data.function.args[1];
        if (!head_is(cl, "List")) return 0;
        qqbar_t g; qqbar_init(g);
        if (!to_qqbar(gen, g)) { qqbar_clear(g); return 0; }
        qqbar_t acc; qqbar_init(acc); qqbar_set_si(acc, 0);
        fmpq_t c; fmpq_init(c);
        int ok = 1;
        size_t m = cl->data.function.arg_count;
        for (size_t i = m; i-- > 0; ) {          /* highest coeff first */
            if (!fmpq_from_expr(cl->data.function.args[i], c)) { ok = 0; break; }
            qqbar_mul(acc, acc, g);
            qqbar_add_fmpq(acc, acc, c);
            if (qqbar_degree(acc) > QQBAR_DEGREE_CAP) { ok = 0; g_qqbar_cap_hit = 1; break; }
        }
        fmpq_clear(c);
        if (ok) qqbar_set(out, acc);
        qqbar_clear(acc); qqbar_clear(g);
        return ok;
    }
    if (is_circular_trig_head(e, n)) {
        /* Sin/Cos/... of a rational multiple r of Pi, via zeta = exp(i pi r):
         *   Cos[r Pi] = (zeta + zeta^-1)/2,  Sin[r Pi] = (zeta - zeta^-1)/(2 i),
         * and Tan/Cot/Sec/Csc as the corresponding ratios (declining at a pole,
         * e.g. Tan[Pi/2], where the needed Cos/Sin vanishes). */
        fmpq_t r; fmpq_init(r);
        int ok = 0;
        if (trig_arg_pi_rational(e->data.function.args[0], r)) {
            const fmpz* pn = fmpq_numref(r);
            const fmpz* pd = fmpq_denref(r);
            if (fmpz_fits_si(pn) && fmpz_fits_si(pd)) {
                slong p = fmpz_get_si(pn);
                ulong q = (ulong)fmpz_get_si(pd);
                qqbar_t zeta, zinv, s, c, two_i;
                qqbar_init(zeta); qqbar_init(zinv);
                qqbar_init(s); qqbar_init(c); qqbar_init(two_i);
                qqbar_exp_pi_i(zeta, p, q);               /* exp( i pi r) */
                qqbar_exp_pi_i(zinv, -p, q);              /* exp(-i pi r) */
                qqbar_add(c, zeta, zinv); qqbar_div_ui(c, c, 2);   /* Cos */
                qqbar_sub(s, zeta, zinv);
                qqbar_i(two_i); qqbar_mul_si(two_i, two_i, 2);     /* 2 i   */
                qqbar_div(s, s, two_i);                            /* Sin   */
                ok = 1;
                if      (head_is(e, "Cos")) qqbar_set(out, c);
                else if (head_is(e, "Sin")) qqbar_set(out, s);
                else if (head_is(e, "Tan")) { if (qqbar_is_zero(c)) ok = 0; else qqbar_div(out, s, c); }
                else if (head_is(e, "Cot")) { if (qqbar_is_zero(s)) ok = 0; else qqbar_div(out, c, s); }
                else if (head_is(e, "Sec")) { if (qqbar_is_zero(c)) ok = 0; else qqbar_inv(out, c); }
                else                        { if (qqbar_is_zero(s)) ok = 0; else qqbar_inv(out, s); } /* Csc */
                if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) { ok = 0; g_qqbar_cap_hit = 1; }
                qqbar_clear(zeta); qqbar_clear(zinv);
                qqbar_clear(s); qqbar_clear(c); qqbar_clear(two_i);
            }
        }
        fmpq_clear(r);
        return ok;
    }
    if (head_is(e, "Root")) return root_object_to_qqbar(e, out);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  qqbar -> canonical Expr (rational / quadratic radical / Root)      */
/* ------------------------------------------------------------------ */

/* Build the minimal polynomial of x as an Expr in Slot[1]. */
static Expr* minpoly_slot_expr(const qqbar_t x) {
    const fmpz_poly_struct* P = QQBAR_POLY(x);
    slong len = fmpz_poly_length(P);
    Expr** terms = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
    size_t nt = 0;
    fmpz_t c; fmpz_init(c);
    for (slong i = 0; i < len; i++) {
        fmpz_poly_get_coeff_fmpz(c, P, i);
        if (fmpz_is_zero(c)) continue;
        Expr* mono;
        if (i == 0) {
            mono = NULL;
        } else {
            Expr* slot = expr_new_function(expr_new_symbol(SYM_Slot),
                            (Expr*[]){ expr_new_integer(1) }, 1);
            mono = (i == 1) ? slot
                 : expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ slot, expr_new_integer((int64_t)i) }, 2);
        }
        Expr* ce = expr_from_fmpz(c);
        Expr* term;
        if (!mono)                          term = ce;
        else if (fmpz_is_one(c))            { term = mono; expr_free(ce); }
        else                                term = expr_new_function(expr_new_symbol(SYM_Times),
                                                        (Expr*[]){ ce, mono }, 2);
        terms[nt++] = term;
    }
    fmpz_clear(c);
    Expr* body = (nt == 0) ? expr_new_integer(0)
               : (nt == 1) ? terms[0]
               : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    free(terms);
    return expr_new_function(expr_new_symbol(SYM_Function), (Expr*[]){ body }, 1);
}

/* ---- qqbar_to_expr Root-object memo, keyed by (minimal polynomial, index) ---- */
#define Q2E_CACHE_SIZE 256u          /* power of two; open addressing */
typedef struct { fmpz_poly_t poly; slong k; uint64_t h; Expr* val; int used; } Q2ESlot;
static Q2ESlot g_q2e_cache[Q2E_CACHE_SIZE];
static size_t g_q2e_count = 0;

static void q2e_cache_reset(void) {
    for (unsigned i = 0; i < Q2E_CACHE_SIZE; i++) {
        if (g_q2e_cache[i].used) {
            fmpz_poly_clear(g_q2e_cache[i].poly);
            expr_free(g_q2e_cache[i].val);
            g_q2e_cache[i].used = 0;
        }
    }
    g_q2e_count = 0;
}

/* On hit, return a fresh copy of the cached Root Expr; else NULL. */
static Expr* q2e_cache_get(const fmpz_poly_struct* p, slong k, uint64_t h) {
    unsigned mask = Q2E_CACHE_SIZE - 1u, i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < Q2E_CACHE_SIZE; probe++) {
        Q2ESlot* s = &g_q2e_cache[i];
        if (!s->used) return NULL;
        if (s->h == h && s->k == k && fmpz_poly_equal(s->poly, p)) return expr_copy(s->val);
        i = (i + 1u) & mask;
    }
    return NULL;
}

/* Insert (minpoly, k) -> a shared copy of `root`; bounded, reset when near full. */
static void q2e_cache_put(const fmpz_poly_struct* p, slong k, uint64_t h, const Expr* root) {
    if (g_q2e_count * 4u >= Q2E_CACHE_SIZE * 3u) q2e_cache_reset();
    unsigned mask = Q2E_CACHE_SIZE - 1u, i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < Q2E_CACHE_SIZE; probe++) {
        Q2ESlot* s = &g_q2e_cache[i];
        if (!s->used) {
            fmpz_poly_init(s->poly); fmpz_poly_set(s->poly, p);
            s->k = k; s->h = h; s->val = expr_copy((Expr*)root); s->used = 1;
            g_q2e_count++;
            return;
        }
        if (s->h == h && s->k == k && fmpz_poly_equal(s->poly, p)) return;
        i = (i + 1u) & mask;
    }
}

static Expr* qqbar_to_expr(const qqbar_t x) {
    if (qqbar_is_rational(x)) {
        fmpq_t q; fmpq_init(q); qqbar_get_fmpq(q, x);
        Expr* r = expr_from_fmpq(q); fmpq_clear(q);
        return r;
    }
    slong d = qqbar_degree(x);
    if (d == 2) {
        fmpz_t a, b, c, q; fmpz_init(a); fmpz_init(b); fmpz_init(c); fmpz_init(q);
        qqbar_get_quadratic(a, b, c, q, x, 0);        /* x = (a + b sqrt(c)) / q */
        Expr* half[2] = { expr_from_fmpz(c),
                          expr_new_function(expr_new_symbol(SYM_Rational),
                              (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) };
        Expr* sqrtc = expr_new_function(expr_new_symbol(SYM_Power), half, 2);
        Expr* bterm = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ expr_from_fmpz(b), sqrtc }, 2);
        Expr* num = expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ expr_from_fmpz(a), bterm }, 2);
        Expr* invq = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_from_fmpz(q), expr_new_integer(-1) }, 2);
        Expr* frac = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ num, invq }, 2);
        fmpz_clear(a); fmpz_clear(b); fmpz_clear(c); fmpz_clear(q);
        return eval_and_free(frac);
    }
    /* degree >= 3: Root[Function[minpoly&], k] (held).  Memoised by (minpoly, k):
     * poly_to_algnum rebuilds this for the FIXED field generator on every
     * AlgebraicNumber arithmetic result (`minpoly_slot_expr` re-allocates the whole
     * minimal-polynomial expression each time), so caching collapses thousands of
     * identical rebuilds -- the #1 cost of the Q(theta) ansatz assembly -- to one
     * per field.  A Root object is a pure, session-independent constant, so the
     * cache never goes stale; it shares the checked-in wl-roots cache reset. */
    slong k = wl_root_index(x);
    {
        uint64_t h = fmpz_poly_hashval(QQBAR_POLY(x)) ^ (uint64_t)(k * 0x9e3779b97f4a7c15ULL);
        Expr* cached = q2e_cache_get(QQBAR_POLY(x), k, h);
        if (cached) return cached;
        Expr* fn = minpoly_slot_expr(x);
        Expr* root = expr_new_function(expr_new_symbol(SYM_Root),
                       (Expr*[]){ fn, expr_new_integer((int64_t)k) }, 2);
        q2e_cache_put(QQBAR_POLY(x), k, h, root);
        return expr_copy(root);
    }
}

/* ------------------------------------------------------------------ */
/*  Method -> "NumberField": common-number-field re-expression         */
/* ------------------------------------------------------------------ */

/* Express x in Q(alpha) at escalating working precision.  64 bits resolves
 * membership for the low-degree fields the tower and dsolve build in the common
 * case, so a single 64-bit attempt is tried first (unchanged from before).  For
 * a HIGH-degree conjugate-root compositum (P8's quartic splitting field, degree
 * 8+), 64 bits is too low and the membership test spuriously fails, declining a
 * valid field; only there do we escalate 64 -> 256 -> 1024 -> 4096.  The
 * `qqbar_degree(alpha) > 6` gate keeps the low-degree path byte-for-byte
 * unchanged — it never escalates — which is what a prior ungated version got
 * wrong (it regressed dsolve). */
static int qqbar_express_in_field_esc(fmpq_poly_t res, const qqbar_t alpha,
                                      const qqbar_t x) {
    if (qqbar_express_in_field(res, alpha, x, 100000, 0, 64)) return 1;
    if (qqbar_degree(alpha) <= 6) return 0;   /* common path: no escalation */
    static const slong precs[3] = { 256, 1024, 4096 };
    for (int i = 0; i < 3; i++)
        if (qqbar_express_in_field(res, alpha, x, 100000, 0, precs[i])) return 1;
    return 0;
}

/* As qqbar_express_in_field_esc but escalates the working precision for EVERY
 * degree, not only degree > 6.  Used by the ToNumberField builtin, where the
 * caller has explicitly named the generator and a valid membership must not be
 * declined merely because the relation's coefficients need > 64 bits to resolve
 * (e.g. expressing Root[x^5-x-1] in Q(Root[2869 x^5 + ...]), MATHILDA_DIVERGENCES
 * A14).  The dsolve-facing `in_field` path keeps the degree<=6 gate, which the
 * ungated escalation regressed. */
static int qqbar_express_in_field_esc_all(fmpq_poly_t res, const qqbar_t alpha,
                                          const qqbar_t x) {
    static const slong precs[4] = { 64, 256, 1024, 4096 };
    for (int i = 0; i < 4; i++)
        if (qqbar_express_in_field(res, alpha, x, 100000, 0, precs[i])) return 1;
    return 0;
}

/* True if x is expressible in Q(alpha). */
static int in_field(const qqbar_t x, const qqbar_t alpha) {
    fmpq_poly_t f; fmpq_poly_init(f);
    int ok = qqbar_express_in_field_esc(f, alpha, x);
    fmpq_poly_clear(f);
    return ok;
}

/* Collect the distinct algebraic generators (radicals, roots of unity, the
 * imaginary unit, Root objects) of `e` into `atoms`. Returns 0 on overflow /
 * conversion failure. Integers, rationals, Plus, Times and integer powers are
 * descended, not treated as atoms. */
static int collect_atoms(const Expr* e, qqbar_ptr atoms, int* na) {
    if (!e) return 1;
    if (e->type != EXPR_FUNCTION) return 1;   /* int/bigint/symbol: no atom */
    size_t n = e->data.function.arg_count;

    int is_generator = 0;
    if (head_is(e, "Power") && n == 2) {
        const Expr* xe = e->data.function.args[1];
        if (xe->type == EXPR_INTEGER)             /* integer power: descend base */
            return collect_atoms(e->data.function.args[0], atoms, na);
        is_generator = 1;                          /* fractional power: a radical */
    } else if (head_is(e, "Sqrt") || head_is(e, "Root")) {
        is_generator = 1;
    } else if (head_is(e, "Complex")) {
        is_generator = 1;                          /* the imaginary unit */
    } else if (head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!collect_atoms(e->data.function.args[i], atoms, na)) return 0;
        return 1;
    } else if (head_is(e, "AlgebraicNumber") && n == 2) {
        /* An AlgebraicNumber lies in Q(generator): its atoms are the generator's
         * (the coefficients are rational, contributing none). */
        return collect_atoms(e->data.function.args[0], atoms, na);
    } else if (head_is(e, "Rational")) {
        return 1;
    } else {
        return 0;
    }

    if (is_generator) {
        if (*na >= QQBAR_MAX_ATOMS) return 0;
        qqbar_t v; qqbar_init(v);
        if (!to_qqbar(e, v)) { qqbar_clear(v); return 0; }
        if (qqbar_is_rational(v)) { qqbar_clear(v); return 1; }  /* not a real generator */
        for (int i = 0; i < *na; i++)
            if (qqbar_equal(atoms + i, v)) { qqbar_clear(v); return 1; }  /* dup */
        qqbar_set(atoms + *na, v);
        (*na)++;
        qqbar_clear(v);
    }
    return 1;
}

/* Evaluate the rational polynomial f at alpha via qqbar Horner. */
static void eval_poly_at(qqbar_t out, const fmpq_poly_t f, const qqbar_t alpha) {
    slong d = fmpq_poly_degree(f);
    qqbar_t acc; qqbar_init(acc); qqbar_set_si(acc, 0);
    fmpq_t c; fmpq_init(c);
    for (slong i = d; i >= 0; i--) {
        qqbar_mul(acc, acc, alpha);
        fmpq_poly_get_coeff_fmpq(c, f, i);
        qqbar_add_fmpq(acc, acc, c);
    }
    fmpq_clear(c);
    qqbar_set(out, acc);
    qqbar_clear(acc);
}

/* Re-derive `direct` through a single primitive element alpha of the field
 * generated by the atoms of `e` (the "AlgebraicNumber objects in a common
 * number field" route). Writes the reconstructed value to `out` and returns 1
 * only if it exactly matches `direct`; otherwise returns 0 and the caller uses
 * the direct value. Distinct computation, identical (verified) result. */
static int number_field_value(const Expr* e, const qqbar_t direct, qqbar_t out) {
    qqbar_ptr atoms = _qqbar_vec_init(QQBAR_MAX_ATOMS);
    int na = 0;
    int ok = collect_atoms(e, atoms, &na);
    if (ok && na == 0) { qqbar_set(out, direct); _qqbar_vec_clear(atoms, QQBAR_MAX_ATOMS); return 1; }

    qqbar_t alpha, cand, ca; qqbar_init(alpha); qqbar_init(cand); qqbar_init(ca);
    if (ok) qqbar_set(alpha, atoms + 0);
    for (int i = 1; i < na && ok; i++) {
        int found = 0;
        for (slong c = 1; c <= 8 && !found; c++) {
            qqbar_mul_si(ca, atoms + i, c);
            qqbar_add(cand, alpha, ca);
            if (in_field(atoms + i, cand) && in_field(alpha, cand)) {
                qqbar_set(alpha, cand); found = 1;
            }
        }
        if (!found) ok = 0;
    }

    if (ok) {
        fmpq_poly_t f; fmpq_poly_init(f);
        if (qqbar_express_in_field(f, alpha, direct, 100000, 0, 64)) {
            eval_poly_at(out, f, alpha);
            ok = qqbar_equal(out, direct);
        } else ok = 0;
        fmpq_poly_clear(f);
    }

    qqbar_clear(alpha); qqbar_clear(cand); qqbar_clear(ca);
    _qqbar_vec_clear(atoms, QQBAR_MAX_ATOMS);
    return ok;
}

/* ------------------------------------------------------------------ */
/*  AlgebraicNumber / ToNumberField internals                          */
/* ------------------------------------------------------------------ */

/* Algebraic-integer generator of Q(alpha): phi = lc * alpha, where lc is the
 * (positive) leading coefficient of alpha's primitive integer minimal
 * polynomial. `lc_out` (caller-initialised) receives lc. */
static void algint_generator(const qqbar_t alpha, qqbar_t phi_out, fmpz_t lc_out) {
    slong d = qqbar_degree(alpha);
    fmpz_poly_get_coeff_fmpz(lc_out, QQBAR_POLY(alpha), d);   /* lc > 0 */
    qqbar_t lcq; qqbar_init(lcq); qqbar_set_fmpz(lcq, lc_out);
    qqbar_mul(phi_out, alpha, lcq);
    qqbar_clear(lcq);
}

/* Load a List[c0,..] of integer/rational Exprs into an fmpq_poly. 1 on success. */
static int coeffs_to_fmpq_poly(const Expr* list, fmpq_poly_t out) {
    if (!head_is(list, "List")) return 0;
    fmpq_t c; fmpq_init(c);
    int ok = 1;
    size_t m = list->data.function.arg_count;
    for (size_t i = 0; i < m; i++) {
        if (!fmpq_from_expr(list->data.function.args[i], c)) { ok = 0; break; }
        fmpq_poly_set_coeff_fmpq(out, (slong)i, c);
    }
    fmpq_clear(c);
    return ok;
}

/* Build the canonical AlgebraicNumber value from a power-basis polynomial `p`
 * (already reduced mod the minimal polynomial of `phi`, degree < n) over the
 * algebraic-integer generator `phi` of degree n. Returns a plain Integer/Rational
 * when the value is rational (n == 1, or every higher coefficient is zero),
 * otherwise AlgebraicNumber[qqbar_to_expr(phi), {d0..d_{n-1}}] with the list
 * padded to n. Consumes nothing; returns a fresh owned Expr. */
/* ---- field-generator cache: theta expr -> (minimal polynomial, degree, canonical?)
 * algnum_binop/algnum_rational reduce a product mod theta's minimal polynomial;
 * recomputing it via to_qqbar(theta) on every field add/mul/pow was a large slice of the
 * Q(theta) ansatz cost (theta is fixed per field).  Keyed by the generator expr --
 * the minimal polynomial is a pure function of it, so the cache never goes stale.
 * `canon` records whether the generator IS the canonical (integer) generator qqbar
 * would emit, so the AlgebraicNumber builtin can skip re-canonicalising an
 * already-canonical value. */
#define GENFD_CACHE_SIZE 128u
typedef struct { Expr* key; uint64_t h; fmpz_poly_t minpoly; slong n; int canon; int used; } GenFDSlot;
static GenFDSlot g_genfd_cache[GENFD_CACHE_SIZE];
static size_t g_genfd_count = 0;

static void genfd_cache_reset(void) {
    for (unsigned i = 0; i < GENFD_CACHE_SIZE; i++) {
        if (g_genfd_cache[i].used) {
            expr_free(g_genfd_cache[i].key);
            fmpz_poly_clear(g_genfd_cache[i].minpoly);
            g_genfd_cache[i].used = 0;
        }
    }
    g_genfd_count = 0;
}

/* Release every qqbar-side memo (see flint_qqbar.h).  Called once at exit by
 * flint_bridge's exit handler, BEFORE flint_cleanup_master tears down FLINT's
 * thread-local pools: the cached qqbar/fmpz_poly values hold fmpz bigints, and
 * `leaks`/valgrind cannot follow FLINT's tagged fmpz pointers, so a reachable
 * cache entry was reported as an orphaned mpz page plus its limbs (155 "leaks"
 * after one DSolve whose Simplify calls filled the cache). */
void flint_qqbar_release_caches(void) {
    flint_qqbar_cache_reset();
    q2e_cache_reset();
    wlroots_reset();
    genfd_cache_reset();
}

/* On success set *Mout (a reference into the cache, valid until the next reset),
 * *nout (degree) and *canon (1 iff gen is the canonical integer generator); return
 * 1.  Return 0 iff gen is not a constant algebraic number. */
static int genfd_get(const Expr* gen, const fmpz_poly_struct** Mout, slong* nout, int* canon) {
    uint64_t h = expr_hash(gen);
    unsigned mask = GENFD_CACHE_SIZE - 1u, i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < GENFD_CACHE_SIZE; probe++) {
        GenFDSlot* s = &g_genfd_cache[i];
        if (!s->used) break;
        if (s->h == h && expr_eq(s->key, gen)) {
            *Mout = s->minpoly; *nout = s->n; if (canon) *canon = s->canon; return 1;
        }
        i = (i + 1u) & mask;
    }
    qqbar_t phi; qqbar_init(phi);
    if (!to_qqbar(gen, phi)) { qqbar_clear(phi); return 0; }
    slong n = qqbar_degree(phi);
    /* Canonical iff gen is EXACTLY the generator the AlgebraicNumber builtin would
     * emit for this value: an algebraic integer (minpoly monic, so algint_generator
     * does not rescale) whose printed form equals gen.  Only then can the builtin's
     * fast path return AlgebraicNumber[gen, coeffs] unchanged (same generator AND
     * same power basis for the coeffs). */
    fmpz_t lead; fmpz_init(lead);
    fmpz_poly_get_coeff_fmpz(lead, QQBAR_POLY(phi), n);
    int is_alg_int = fmpz_is_one(lead);
    fmpz_clear(lead);
    Expr* cg = is_alg_int ? qqbar_to_expr(phi) : NULL;
    int canonical = is_alg_int && cg && expr_eq(cg, gen);
    if (cg) expr_free(cg);
    if (g_genfd_count * 4u >= GENFD_CACHE_SIZE * 3u) genfd_cache_reset();
    i = (unsigned)h & mask;
    for (unsigned probe = 0; probe < GENFD_CACHE_SIZE; probe++) {
        GenFDSlot* s = &g_genfd_cache[i];
        if (!s->used) {
            s->key = expr_copy((Expr*)gen); s->h = h;
            fmpz_poly_init(s->minpoly); fmpz_poly_set(s->minpoly, QQBAR_POLY(phi));
            s->n = n; s->canon = canonical; s->used = 1; g_genfd_count++;
            *Mout = s->minpoly; *nout = n; if (canon) *canon = canonical;
            qqbar_clear(phi);
            return 1;
        }
        i = (i + 1u) & mask;
    }
    qqbar_clear(phi);
    return 0;   /* table full even after reset: shouldn't happen */
}

/* poly_to_algnum reusing a KNOWN generator expr `gen` (canonical, fixed per field)
 * instead of rebuilding it via qqbar_to_expr(phi) -- the hot path for field
 * arithmetic, where theta never changes.  A rational result (all higher coords 0)
 * still collapses to the plain number. */
static Expr* poly_to_algnum_gen(const Expr* gen, const fmpq_poly_t p, slong n) {
    Expr** dl = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
    int all_higher_zero = 1;
    fmpq_t di; fmpq_init(di);
    for (slong i = 0; i < n; i++) {
        fmpq_poly_get_coeff_fmpq(di, p, i);
        dl[i] = expr_from_fmpq(di);
        if (i >= 1 && !fmpq_is_zero(di)) all_higher_zero = 0;
    }
    fmpq_clear(di);
    Expr* result;
    if (n <= 1 || all_higher_zero) {
        for (slong i = 1; i < n; i++) expr_free(dl[i]);
        result = dl[0]; free(dl);
    } else {
        Expr* clist = expr_new_function(expr_new_symbol(SYM_List), dl, (size_t)n);
        free(dl);
        result = expr_new_function(expr_new_symbol(SYM_AlgebraicNumber),
                     (Expr*[]){ expr_copy((Expr*)gen), clist }, 2);
    }
    return result;
}

/* ------------------------------------------------------------------ */
/*  Native number-field matrix RREF (the ansatz-solve fast path)        */
/* ------------------------------------------------------------------ */

/* Convert one matrix entry to an nf_elem over the field generated by `gen`.
 * Accepts an integer/rational, or AlgebraicNumber[gen, {coords}] whose
 * generator is expr_eq to `gen`.  Returns 1 on success, 0 (leaving the matrix
 * to the classical path) for anything else. */
static int nf_entry_from_expr(const Expr* e, const Expr* gen, nf_elem_t out, const nf_t nf) {
    fmpq_t q; fmpq_init(q);
    if (fmpq_from_expr(e, q)) { nf_elem_set_fmpq(out, q, nf); fmpq_clear(q); return 1; }
    fmpq_clear(q);
    if (head_is(e, "AlgebraicNumber") && e->type == EXPR_FUNCTION
        && e->data.function.arg_count == 2) {
        if (!expr_eq(e->data.function.args[0], gen)) return 0;
        fmpq_poly_t p; fmpq_poly_init(p);
        int ok = coeffs_to_fmpq_poly(e->data.function.args[1], p);
        if (ok) nf_elem_set_fmpq_poly(out, p, nf);
        fmpq_poly_clear(p);
        return ok;
    }
    return 0;
}

/* Reduced row echelon form of an r×c matrix whose entries are rationals and
 * AlgebraicNumber[theta, ...] over ONE common theta, computed exactly with
 * FLINT antic nf_elem arithmetic (Gauss-Jordan over the field).  The RREF is
 * canonical, so the result is identical to the classical AlgebraicNumber
 * RowReduce (any Method) -- but the field arithmetic is native nf_elem, not
 * per-element evaluator dispatch, which is the ansatz-solve hot path.
 *
 * Returns a List-of-Lists Expr (caller owns), or NULL to fall through to the
 * classical path: no AlgebraicNumber entry (pure-rational -> the fmpq_mat
 * path), a second distinct generator, a non-field entry, or no FLINT.  Does not
 * take ownership of `flat`. */
Expr* flint_qqbar_nf_mat_rref(Expr** flat, int r, int c) {
    if (!flat || r <= 0 || c <= 0) return NULL;

    /* The field generator: the first AlgebraicNumber entry. */
    const Expr* gen = NULL;
    for (long i = 0; i < (long)r * c; i++) {
        const Expr* e = flat[i];
        if (e && e->type == EXPR_FUNCTION && head_is(e, "AlgebraicNumber")
            && e->data.function.arg_count == 2) { gen = e->data.function.args[0]; break; }
    }
    if (!gen) return NULL;   /* pure-rational (or empty): leave to the fmpq path */

    const fmpz_poly_struct* Mz; slong n;
    if (!genfd_get(gen, &Mz, &n, NULL) || n < 1) return NULL;

    fmpq_poly_t Mq; fmpq_poly_init(Mq); fmpq_poly_set_fmpz_poly(Mq, Mz);
    nf_t nf; nf_init(nf, Mq); fmpq_poly_clear(Mq);

    long total = (long)r * c;
    nf_elem_t* A = malloc(sizeof(nf_elem_t) * (size_t)total);
    for (long i = 0; i < total; i++) nf_elem_init(A[i], nf);

    int ok = 1;
    for (long i = 0; i < total && ok; i++)
        ok = nf_entry_from_expr(flat[i], gen, A[i], nf);

    Expr* out = NULL;
    if (ok) {
        nf_elem_t piv, factor, t;
        nf_elem_init(piv, nf); nf_elem_init(factor, nf); nf_elem_init(t, nf);
        int prow_count = 0;
        for (int col = 0; col < c && prow_count < r; col++) {
            int prow = -1;
            for (int i = prow_count; i < r; i++)
                if (!nf_elem_is_zero(A[(long)i * c + col], nf)) { prow = i; break; }
            if (prow < 0) continue;
            if (prow != prow_count)
                for (int j = 0; j < c; j++)
                    nf_elem_swap(A[(long)prow_count * c + j], A[(long)prow * c + j], nf);
            /* scale the pivot row so the pivot becomes 1 */
            nf_elem_inv(piv, A[(long)prow_count * c + col], nf);
            for (int j = 0; j < c; j++)
                nf_elem_mul(A[(long)prow_count * c + j], A[(long)prow_count * c + j], piv, nf);
            /* eliminate the pivot column from every other row */
            for (int i = 0; i < r; i++) {
                if (i == prow_count) continue;
                if (nf_elem_is_zero(A[(long)i * c + col], nf)) continue;
                nf_elem_set(factor, A[(long)i * c + col], nf);
                for (int j = 0; j < c; j++) {
                    nf_elem_mul(t, factor, A[(long)prow_count * c + j], nf);
                    nf_elem_sub(A[(long)i * c + j], A[(long)i * c + j], t, nf);
                }
            }
            prow_count++;
        }
        nf_elem_clear(piv, nf); nf_elem_clear(factor, nf); nf_elem_clear(t, nf);

        /* read the RREF back as canonical AlgebraicNumber / rational Exprs */
        Expr** rows = malloc(sizeof(Expr*) * (size_t)r);
        fmpq_poly_t pr2; fmpq_poly_init(pr2);
        for (int i = 0; i < r; i++) {
            Expr** el = malloc(sizeof(Expr*) * (size_t)c);
            for (int j = 0; j < c; j++) {
                nf_elem_get_fmpq_poly(pr2, A[(long)i * c + j], nf);
                el[j] = poly_to_algnum_gen(gen, pr2, n);
            }
            rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
            free(el);
        }
        fmpq_poly_clear(pr2);
        out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)r);
        free(rows);
    }

    for (long i = 0; i < total; i++) nf_elem_clear(A[i], nf);
    free(A);
    nf_clear(nf);
    return out;
}

static Expr* poly_to_algnum(const qqbar_t phi, const fmpq_poly_t p, slong n) {
    Expr** dl = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
    int all_higher_zero = 1;
    fmpq_t di; fmpq_init(di);
    for (slong i = 0; i < n; i++) {
        fmpq_poly_get_coeff_fmpq(di, p, i);
        dl[i] = expr_from_fmpq(di);
        if (i >= 1 && !fmpq_is_zero(di)) all_higher_zero = 0;
    }
    fmpq_clear(di);

    Expr* result;
    if (n <= 1 || all_higher_zero) {
        for (slong i = 1; i < n; i++) expr_free(dl[i]);
        result = dl[0];                 /* the rational value d0 */
        free(dl);
    } else {
        Expr* g = qqbar_to_expr(phi);
        Expr* clist = expr_new_function(expr_new_symbol(SYM_List), dl, (size_t)n);
        free(dl);
        result = expr_new_function(expr_new_symbol(SYM_AlgebraicNumber),
                     (Expr*[]){ g, clist }, 2);
    }
    return result;
}

/* poly_to_algnum for a run of values over one generator phi: the generator's
 * expression is built at most once, and only when a value needs it.  Building
 * it orders phi's conjugates, which is the expensive part -- for 2^(1/32)
 * (audit NumberFieldIntegralBasis#1) it cost most of 20 s when redone per row,
 * and a run whose values are all rational must not pay for it at all. */
static Expr* poly_to_algnum_lazy(const qqbar_t phi, Expr** gen,
                                 const fmpq_poly_t p, slong n) {
    if (n <= 1 || fmpq_poly_degree(p) <= 0)
        return poly_to_algnum(phi, p, n);       /* rational: no generator */
    if (!*gen) *gen = qqbar_to_expr(phi);
    return poly_to_algnum_gen(*gen, p, n);
}

/* Same-field binary op on two AlgebraicNumbers with expr_eq generators.
 * op 0 = add, 1 = mul. Returns a fresh combined value, or NULL to decline. */
static Expr* algnum_binop(const Expr* a, const Expr* b, int op) {
    if (!head_is(a, "AlgebraicNumber") || !head_is(b, "AlgebraicNumber")) return NULL;
    if (a->data.function.arg_count != 2 || b->data.function.arg_count != 2) return NULL;
    const Expr* g = a->data.function.args[0];
    if (!expr_eq(g, b->data.function.args[0])) return NULL;

    /* theta's minimal polynomial + degree from the generator cache -- no to_qqbar
     * per operation.  The result reuses the SAME generator expr g (theta is fixed
     * per field), so poly_to_algnum_gen never rebuilds it via qqbar_to_expr. */
    const fmpz_poly_struct* Mz; slong n;
    if (!genfd_get(g, &Mz, &n, NULL)) return NULL;

    fmpq_poly_t M, pa, pb, r;
    fmpq_poly_init(M); fmpq_poly_init(pa); fmpq_poly_init(pb); fmpq_poly_init(r);
    fmpq_poly_set_fmpz_poly(M, Mz);
    Expr* result = NULL;
    if (coeffs_to_fmpq_poly(a->data.function.args[1], pa) &&
        coeffs_to_fmpq_poly(b->data.function.args[1], pb)) {
        if (op == 0) { fmpq_poly_add(r, pa, pb); }
        else         { fmpq_poly_mul(r, pa, pb); fmpq_poly_rem(r, r, M); }
        result = poly_to_algnum_gen(g, r, n);
    }
    fmpq_poly_clear(M); fmpq_poly_clear(pa); fmpq_poly_clear(pb); fmpq_poly_clear(r);
    return result;
}

/* ------------------------------------------------------------------ */
/*  Public API                                                         */
/* ------------------------------------------------------------------ */

int flint_qqbar_is_constant_algebraic(const Expr* e) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_BIGINT:  return 1;
        case EXPR_SYMBOL:  return e->data.symbol.name == SYM_GoldenRatio; /* sole algebraic constant */
        case EXPR_FUNCTION: break;
        default:           return 0;   /* real / string / other free symbol */
    }
    size_t n = e->data.function.arg_count;
    if (head_is(e, "Rational")) return 1;
    if (head_is(e, "Root"))     return 1;   /* opaque; converter validates */
    if (head_is(e, "Complex") || head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!flint_qqbar_is_constant_algebraic(e->data.function.args[i])) return 0;
        return 1;
    }
    if (head_is(e, "Sqrt") && n == 1)
        return flint_qqbar_is_constant_algebraic(e->data.function.args[0]);
    if (head_is(e, "Power") && n == 2) {
        const Expr* be = e->data.function.args[0];
        const Expr* xe = e->data.function.args[1];
        if (be->type == EXPR_SYMBOL && be->data.symbol.name == SYM_E) {
            fmpq_t r; fmpq_init(r);            /* E^(I Pi rational) is a root of unity */
            int ok = exp_arg_i_pi_rational(xe, r);
            fmpq_clear(r);
            return ok;
        }
        int exp_ok = xe->type == EXPR_INTEGER || head_is(xe, "Rational") ||
                     xe->type == EXPR_BIGINT;
        return exp_ok && flint_qqbar_is_constant_algebraic(be);
    }
    if ((head_is(e, "Re") || head_is(e, "Im") || head_is(e, "Abs") ||
         head_is(e, "Conjugate")) && n == 1)
        return flint_qqbar_is_constant_algebraic(e->data.function.args[0]);
    if (is_circular_trig_head(e, n)) {
        /* Sin/Cos/Tan/Cot/Sec/Csc of a rational multiple of Pi is algebraic
         * (a root of unity combination).  Optimistic: a pole (Tan[Pi/2]) is
         * caught by to_qqbar declining, so flint_qqbar_canonical returns NULL. */
        fmpq_t r; fmpq_init(r);
        int ok = trig_arg_pi_rational(e->data.function.args[0], r);
        fmpq_clear(r);
        return ok;
    }
    if (head_is(e, "AlgebraicNumber") && n == 2) {
        /* AlgebraicNumber[gen, {rationals}]: gen constant-algebraic, coeffs rational. */
        if (!flint_qqbar_is_constant_algebraic(e->data.function.args[0])) return 0;
        const Expr* cl = e->data.function.args[1];
        if (!head_is(cl, "List")) return 0;
        for (size_t i = 0; i < cl->data.function.arg_count; i++) {
            const Expr* c = cl->data.function.args[i];
            if (!(c->type == EXPR_INTEGER || c->type == EXPR_BIGINT ||
                  head_is(c, "Rational")))
                return 0;
        }
        return 1;
    }
    return 0;
}

Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod method) {
    if (!flint_qqbar_is_constant_algebraic(e)) return NULL;
    qqbar_t val; qqbar_init(val);
    if (!to_qqbar(e, val)) { qqbar_clear(val); return NULL; }

    Expr* out;
    if (method == QQBAR_METHOD_NUMBERFIELD) {
        qqbar_t v2; qqbar_init(v2);
        out = qqbar_to_expr(number_field_value(e, val, v2) ? v2 : val);
        qqbar_clear(v2);
    } else {
        out = qqbar_to_expr(val);
    }
    qqbar_clear(val);
    return out;
}

Expr* flint_qqbar_reduce_coeffs(const Expr* e, QQBarMethod method) {
    if (!e) return NULL;
    /* Atomic number leaf: already canonical, avoid rebuilding it via qqbar. A
     * bare free symbol is not constant-algebraic, so flint_qqbar_canonical below
     * declines it and we fall through to the identity copy. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL)
        return expr_copy((Expr*)e);
    /* Maximal constant-algebraic subexpression: fold to one canonical number. */
    Expr* q = flint_qqbar_canonical(e, method);
    if (q) return q;
    /* Otherwise recurse into the (free-variable-bearing) structure. Non-function
     * atoms (a symbol, a string) are returned as-is. */
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = flint_qqbar_reduce_coeffs(e->data.function.args[i], method);
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return out;
}

/* The exact conversions of the operands of the comparisons, memoised.  A sign
 * diagram or a decomposition compares the same few algebraic numbers (its
 * breakpoints, -Sqrt[Root[...]] and the like) against each other over and
 * over, and each conversion of a compound one is a root extraction or a
 * resultant plus a factorisation.  Direct-mapped by the structural hash; a
 * conversion is a pure function of the expression, so an entry never goes
 * stale.  Bare numbers are not worth an entry. */
#define QQ_MEMO_SIZE 256
typedef struct { uint64_t h; Expr* e; qqbar_struct q; int ok; } QqMemo;
static QqMemo g_qq_memo[QQ_MEMO_SIZE];

static int to_qqbar_memo(const Expr* e, qqbar_t out) {
    if (e->type != EXPR_FUNCTION || head_is(e, "Rational")) return to_qqbar(e, out);
    uint64_t h = expr_hash(e);
    QqMemo* m = &g_qq_memo[h % QQ_MEMO_SIZE];
    if (m->e && m->h == h && expr_eq(m->e, e)) {
        if (m->ok) qqbar_set(out, &m->q);
        return m->ok;
    }
    int ok = to_qqbar(e, out);
    if (m->e) { expr_free(m->e); qqbar_clear(&m->q); }
    m->h = h; m->e = expr_copy((Expr*)e); m->ok = ok;
    qqbar_init(&m->q);
    if (ok) qqbar_set(&m->q, out);
    return ok;
}

int flint_qqbar_equal(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -1;
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -1;
    if (to_qqbar_memo(a, qa) && to_qqbar_memo(b, qb)) r = qqbar_equal(qa, qb) ? 1 : 0;
    qqbar_clear(qa); qqbar_clear(qb);
    return r;
}

/* A certified real enclosure of a constant algebraic expression, computed in
 * ball arithmetic over the tree: 1 when e is PROVABLY real (rational leaves,
 * sums and products of real parts, integer powers, roots and square roots of a
 * strictly positive ball, a real qqbar for a Root[...] / AlgebraicNumber[...]
 * leaf, whose exact form is cheap), 0 when realness cannot be
 * certified this way (the caller then takes the exact path).  No minimal
 * polynomial of a sum or product is ever formed. */
static void qq_arb_pow_si(arb_t r, const arb_t b, slong k, slong prec) {
    arb_pow_ui(r, b, (ulong)(k < 0 ? -k : k), prec);
    if (k < 0) arb_inv(r, r, prec);
}

static int qq_arb_real(const Expr* e, arb_t out, slong prec) {
    if (e->type == EXPR_INTEGER) { arb_set_si(out, (slong)e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT) {
        fmpz_t z; fmpz_init(z); fmpz_set_mpz(z, e->data.bigint);
        arb_set_fmpz(out, z); fmpz_clear(z); return 1;
    }
    if (e->type != EXPR_FUNCTION && e->type != EXPR_SYMBOL) return 0;
    size_t n = e->type == EXPR_FUNCTION ? e->data.function.arg_count : 0;
    if (head_is(e, "Rational") && n == 2) {
        fmpq_t q; fmpq_init(q);
        int ok = fmpq_from_expr(e, q);
        if (ok) arb_set_fmpq(out, q, prec);
        fmpq_clear(q);
        return ok;
    }
    if ((head_is(e, "Plus") || head_is(e, "Times")) && n > 0) {
        int plus = head_is(e, "Plus");
        arb_t t; arb_init(t);
        int ok = qq_arb_real(e->data.function.args[0], out, prec);
        for (size_t i = 1; i < n && ok; i++) {
            ok = qq_arb_real(e->data.function.args[i], t, prec);
            if (ok) { if (plus) arb_add(out, out, t, prec); else arb_mul(out, out, t, prec); }
        }
        arb_clear(t);
        return ok;
    }
    if ((head_is(e, "Power") && n == 2) || (head_is(e, "Sqrt") && n == 1)) {
        const Expr* be = e->data.function.args[0];
        if (be->type == EXPR_SYMBOL) return 0;                 /* E^..., not here */
        arb_t b; arb_init(b);
        int ok = qq_arb_real(be, b, prec);
        if (ok && head_is(e, "Sqrt")) {
            ok = arb_is_positive(b);
            if (ok) arb_sqrt(out, b, prec);
        } else if (ok) {
            const Expr* xe = e->data.function.args[1];
            fmpq_t ef; fmpq_init(ef);
            if (xe->type == EXPR_INTEGER) {
                ok = xe->data.integer >= 0 || !arb_contains_zero(b);
                if (ok) qq_arb_pow_si(out, b, (slong)xe->data.integer, prec);
            } else if (fmpq_from_expr(xe, ef) && fmpz_fits_si(fmpq_numref(ef)) &&
                       fmpz_fits_si(fmpq_denref(ef))) {
                /* the principal root of a positive real is the real root */
                ok = arb_is_positive(b);
                if (ok) {
                    arb_root_ui(out, b, (ulong)fmpz_get_si(fmpq_denref(ef)), prec);
                    qq_arb_pow_si(out, out, fmpz_get_si(fmpq_numref(ef)), prec);
                }
            } else ok = 0;
            fmpq_clear(ef);
        }
        arb_clear(b);
        return ok;
    }
    /* a Root[...] / AlgebraicNumber[...] leaf (or GoldenRatio): its exact form,
     * which is cheap (Root conversions are cached), when it is real.  Anything
     * else -- Re, Im, Abs, Complex -- is left to the exact path rather than
     * converted twice. */
    if (!(head_is(e, "Root") || head_is(e, "AlgebraicNumber") ||
          (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_GoldenRatio)))
        return 0;
    qqbar_t q; qqbar_init(q);
    int ok = to_qqbar(e, q) && qqbar_is_real(q);
    if (ok) qqbar_get_arb(out, q, prec);
    qqbar_clear(q);
    return ok;
}

int flint_qqbar_compare(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -2;
    /* Ball arithmetic first: two certified real enclosures that are disjoint
     * decide the order without the exact minimal polynomials of a and b (the
     * cost of every CAD sign test on nested radicals and Root objects).  Ties
     * and unprovable realness take the exact path below. */
    {
        arb_t xa, xb; arb_init(xa); arb_init(xb);
        int r = -2;
        if (qq_arb_real(a, xa, 128) && qq_arb_real(b, xb, 128)) {
            if (arb_lt(xa, xb)) r = -1;
            else if (arb_gt(xa, xb)) r = 1;
            else if (arb_is_exact(xa) && arb_is_exact(xb) && arb_equal(xa, xb)) r = 0;
        }
        arb_clear(xa); arb_clear(xb);
        if (r != -2) return r;
    }
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -2;
    if (to_qqbar_memo(a, qa) && to_qqbar_memo(b, qb) && qqbar_is_real(qa) && qqbar_is_real(qb))
        r = qqbar_cmp_re(qa, qb);
    qqbar_clear(qa); qqbar_clear(qb);
    return r;
}

int flint_qqbar_is_real(const Expr* e) {
    if (!flint_qqbar_is_constant_algebraic(e)) return -1;
    qqbar_t q; qqbar_init(q);
    int r = -1;
    if (to_qqbar(e, q)) r = qqbar_is_real(q) ? 1 : 0;
    qqbar_clear(q);
    return r;
}

/* See flint_qqbar.h.  qqbar_sgn_re / qqbar_sgn_im decide the sign exactly:
 * a zero part is certified by the minimal polynomial, never by an enclosure
 * that merely contains 0.  Memoised (to_qqbar_memo): N asks for the same
 * radicand once per precision of its adaptive climb. */
int flint_qqbar_part_sign(const Expr* e, int imag) {
    if (!flint_qqbar_is_constant_algebraic(e)) return -2;
    qqbar_t q; qqbar_init(q);
    int r = -2;
    if (to_qqbar_memo(e, q)) r = imag ? qqbar_sgn_im(q) : qqbar_sgn_re(q);
    qqbar_clear(q);
    return r;
}

Expr* flint_qqbar_algebraic_number(const Expr* gen, const Expr* coeffs) {
    if (!gen || !coeffs || !head_is(coeffs, "List")) return NULL;

    /* Canonical fast path.  The evaluator re-canonicalises every AlgebraicNumber it
     * sees; during the Q(theta) ansatz assembly that is the single largest cost
     * (to_qqbar + reduce + poly_to_algnum + qqbar_to_expr, per re-evaluation).  When
     * the generator is already the canonical integer generator (`canon` from the
     * cache) and the coeffs are already reduced -- rational and length <= deg minpoly
     * -- the slow path below returns exactly AlgebraicNumber[gen, coeffs] (padded to
     * n), or the plain rational when only the constant term is nonzero.  Return that
     * directly.  Anything else (non-canonical generator, an unreduced or non-rational
     * coeff) falls through to the exact slow path. */
    {
        const fmpz_poly_struct* Mz; slong n; int canon;
        if (genfd_get(gen, &Mz, &n, &canon) && canon) {
            size_t m = coeffs->data.function.arg_count;
            if ((slong)m <= n) {
                int all_rat = 1, higher_nonzero = 0;
                for (size_t i = 0; i < m; i++) {
                    const Expr* c = coeffs->data.function.args[i];
                    if (!(c->type == EXPR_INTEGER || c->type == EXPR_BIGINT || head_is(c, "Rational"))) { all_rat = 0; break; }
                    if (i >= 1 && !(c->type == EXPR_INTEGER && c->data.integer == 0)) higher_nonzero = 1;
                }
                if (all_rat) {
                    if (!higher_nonzero)
                        return m >= 1 ? expr_copy((Expr*)coeffs->data.function.args[0]) : expr_new_integer(0);
                    if ((slong)m == n)
                        return expr_new_function(expr_new_symbol(SYM_AlgebraicNumber),
                                   (Expr*[]){ expr_copy((Expr*)gen), expr_copy((Expr*)coeffs) }, 2);
                    /* m < n: pad the coeff list to length n (as poly_to_algnum does). */
                    Expr** dl = malloc(sizeof(Expr*) * (size_t)n);
                    for (slong i = 0; i < n; i++)
                        dl[i] = (size_t)i < m ? expr_copy((Expr*)coeffs->data.function.args[i]) : expr_new_integer(0);
                    Expr* clist = expr_new_function(expr_new_symbol(SYM_List), dl, (size_t)n);
                    free(dl);
                    return expr_new_function(expr_new_symbol(SYM_AlgebraicNumber),
                               (Expr*[]){ expr_copy((Expr*)gen), clist }, 2);
                }
            }
        }
    }

    qqbar_t alpha; qqbar_init(alpha);
    if (!to_qqbar(gen, alpha)) { qqbar_clear(alpha); return NULL; }

    fmpz_t lc; fmpz_init(lc);
    qqbar_t phi; qqbar_init(phi);
    algint_generator(alpha, phi, lc);        /* phi = lc*alpha, algebraic integer */
    slong n = qqbar_degree(phi);

    /* p(x) = sum_i (c_i / lc^i) x^i. */
    fmpq_poly_t p, M; fmpq_poly_init(p); fmpq_poly_init(M);
    fmpq_poly_set_fmpz_poly(M, QQBAR_POLY(phi));   /* monic minimal polynomial */
    fmpq_t c, lcpow, coef; fmpq_init(c); fmpq_init(lcpow); fmpq_init(coef);
    fmpq_one(lcpow);
    int ok = 1;
    size_t m = coeffs->data.function.arg_count;
    for (size_t i = 0; i < m; i++) {
        if (!fmpq_from_expr(coeffs->data.function.args[i], c)) { ok = 0; break; }
        fmpq_div(coef, c, lcpow);
        fmpq_poly_set_coeff_fmpq(p, (slong)i, coef);
        fmpq_mul_fmpz(lcpow, lcpow, lc);
    }
    fmpq_clear(c); fmpq_clear(lcpow); fmpq_clear(coef);

    Expr* result = NULL;
    if (ok) {
        fmpq_poly_rem(p, p, M);              /* reduce to power basis of phi */
        result = poly_to_algnum(phi, p, n);
    }
    fmpq_poly_clear(p); fmpq_poly_clear(M);
    fmpz_clear(lc); qqbar_clear(phi); qqbar_clear(alpha);
    return result;
}

/* Express x = res(alpha) exactly, escalating the working precision for a
 * high-degree compositum (see qqbar_express_in_field_esc). */
/* Membership used by flint_qqbar_to_number_field_common AFTER it has chosen a
 * primitive element by degree.  At that point every atom provably lies in
 * Q(alpha), so a failure can only be insufficient working precision and
 * escalating for every degree is right -- the same reasoning that justifies
 * _esc_all for an explicitly named generator (A14).  The shared `in_field` keeps
 * its degree<=6 fast path, since the dsolve-facing callers ask about membership
 * they have no reason to expect and the ungated escalation regressed them. */
static int qqbar_express_in_field_retry(fmpq_poly_t res, const qqbar_t alpha,
                                        const qqbar_t x) {
    return qqbar_express_in_field_esc_all(res, alpha, x);
}

Expr* flint_qqbar_to_number_field(const Expr* a, const Expr* theta) {
    if (!a || !theta) return NULL;
    qqbar_t av, tv; qqbar_init(av); qqbar_init(tv);
    if (!to_qqbar(a, av) || !to_qqbar(theta, tv)) {
        qqbar_clear(av); qqbar_clear(tv); return NULL;
    }
    fmpz_t lc; fmpz_init(lc);
    qqbar_t phi; qqbar_init(phi);
    algint_generator(tv, phi, lc);
    fmpz_clear(lc);
    slong n = qqbar_degree(phi);

    fmpq_poly_t f; fmpq_poly_init(f);
    Expr* result = NULL;
    if (qqbar_express_in_field_esc_all(f, phi, av))   /* escalate for named generator (A14) */
        result = poly_to_algnum(phi, f, n);   /* NULL stays NULL when a not in Q(theta) */
    fmpq_poly_clear(f);
    qqbar_clear(phi); qqbar_clear(av); qqbar_clear(tv);
    return result;
}

Expr* flint_qqbar_to_number_field_self(const Expr* x) {
    return flint_qqbar_to_number_field(x, x);
}

/* ------------------------------------------------------------------ */
/*  Fast common field of radicals of rationals and roots of unity      */
/* ------------------------------------------------------------------ */
/*
 * The common field of a set of atoms that are each a root of unity
 * ((-1)^(p/q), I, Complex[0, +-1]) or a radical r^(p/q) of a rational r is
 * a quotient of the explicit product algebra
 *
 *     A = Q[z]/(Phi_N(z)) (x) Q[g_1]/(g_1^n_1 - l_1) (x) ... ,
 *
 * z = exp(2 Pi I/N), g_i = l_i^(1/n_i) > 0 for the distinct PRIMES l_i of the
 * radicands (an atom r^(p/q) with r < 0 is |r|^(p/q) (-1)^(p/q), and |r|^(p/q)
 * is a product of prime radicals).  Every atom is a monomial in the generators
 * of A, so this is the whole computation, done by exact linear algebra in the
 * product basis of A -- no algebraic-number arithmetic:
 *
 *   1. theta = sum c_i gen_i.  Its powers 1, theta, .., theta^(D-1) are a basis
 *      of A (D = dim A) for all but finitely many c; then A = Q[theta] and the
 *      minimal polynomial m_A of theta over A is read off theta^D.
 *   2. A need not be a field: the atoms may be dependent (Sqrt[3] lies in
 *      Q(I, (-1)^(1/3)), (-1)^(2/3) in Q((-1)^(1/3)), Sqrt[2] in Q((-1)^(1/4))).
 *      A is a product of fields, m_A the product of their minimal polynomials,
 *      and the field the atoms actually generate is the image of A under the
 *      evaluation at the principal values -- the factor f of m_A that vanishes
 *      at the numerical value of theta.  That factor is found with ball
 *      arithmetic, precision raised until exactly one factor qualifies.
 *   3. Each atom's coordinates in the theta basis of A, reduced modulo f, are
 *      its image in K = Q[x]/(f) = Q(atoms).
 *
 * This replaces the qqbar primitive-element search for the common case of the
 * algebraic extension engines (fg_field_images, algfield): for
 * Q(I, Sqrt[3], 2^(1/3)) the qqbar search composed and factored sixteen
 * resultants per atom and then expressed every atom by LLL, ~0.25-0.7 s per
 * field build; this is a 12 x 12 rational matrix inverse.  Returns the same
 * List[AlgebraicNumber[phi, {..}] | rational, ..] as the qqbar search (phi
 * = theta, an algebraic integer), or NULL to decline: an atom of another
 * shape, a radicand that does not factor cheaply, or an algebra above
 * KF_MAX_DIM.
 */
#define KF_MAX_DIM    64
#define KF_MAX_PRIMES 8

typedef struct {
    fmpq_t u;                      /* root-of-unity part: (-1)^u            */
    fmpq_t c;                      /* rational factor                       */
    slong  np;
    fmpz   ell[KF_MAX_PRIMES];     /* primes of |radicand| ...               */
    fmpq   a[KF_MAX_PRIMES];       /* ... and their rational exponents       */
} KfAtom;

static void kf_vec_zero(fmpq* v, slong n) { for (slong i = 0; i < n; i++) fmpq_zero(v + i); }
static void kf_vec_set(fmpq* v, const fmpq* w, slong n) { for (slong i = 0; i < n; i++) fmpq_set(v + i, w + i); }

static ulong kf_lcm(ulong a, ulong b) {
    ulong x = a, y = b;
    while (y) { ulong r = x % y; x = y; y = r; }
    return x ? a / x * b : 0;
}

static void kf_atom_init(KfAtom* t) {
    fmpq_init(t->u); fmpq_init(t->c); fmpq_one(t->c); t->np = 0;
    for (slong i = 0; i < KF_MAX_PRIMES; i++) { fmpz_init(t->ell + i); fmpq_init(t->a + i); }
}
static void kf_atom_clear(KfAtom* t) {
    fmpq_clear(t->u); fmpq_clear(t->c);
    for (slong i = 0; i < KF_MAX_PRIMES; i++) { fmpz_clear(t->ell + i); fmpq_clear(t->a + i); }
}

/* Add ell^(e * x) (x rational) to the atom's prime list.  0 on overflow. */
static int kf_atom_add_prime(KfAtom* t, const fmpz_t ell, slong e, const fmpq_t x) {
    fmpq_t ex; fmpq_init(ex);
    fmpq_mul_si(ex, x, e);
    slong i = 0;
    for (; i < t->np; i++) if (fmpz_equal(t->ell + i, ell)) break;
    int ok = 1;
    if (i == t->np) {
        if (t->np == KF_MAX_PRIMES) ok = 0;
        else { fmpz_set(t->ell + i, ell); fmpq_zero(t->a + i); t->np++; }
    }
    if (ok) fmpq_add(t->a + i, t->a + i, ex);
    fmpq_clear(ex);
    return ok;
}

/* Factor |z| (z an integer, |z| < 2^64) into t's prime list with exponent
 * sign * x.  0 when z is too large to factor cheaply. */
static int kf_atom_add_int(KfAtom* t, const fmpz_t z, int sign, const fmpq_t x) {
    if (fmpz_bits(z) > 64) return 0;
    if (fmpz_is_pm1(z)) return 1;
    fmpz_factor_t f; fmpz_factor_init(f);
    fmpz_factor(f, z);
    int ok = 1;
    for (slong i = 0; i < f->num && ok; i++)
        ok = kf_atom_add_prime(t, f->p + i, sign * (slong)f->exp[i], x);
    fmpz_factor_clear(f);
    return ok;
}

/* Read one atom.  0 declines. */
static int kf_atom_parse(const Expr* e, KfAtom* t) {
    /* the imaginary unit (Complex[0, +-1]) */
    if (head_is(e, "Complex") && e->data.function.arg_count == 2) {
        const Expr* re = e->data.function.args[0];
        const Expr* im = e->data.function.args[1];
        if (re->type != EXPR_INTEGER || re->data.integer != 0 || im->type != EXPR_INTEGER ||
            (im->data.integer != 1 && im->data.integer != -1)) return 0;
        fmpq_set_si(t->u, im->data.integer == 1 ? 1 : 3, 2);
        return 1;
    }
    const Expr* base = NULL; fmpq_t x; fmpq_init(x);
    int ok = 0;
    if (head_is(e, "Sqrt") && e->data.function.arg_count == 1) {
        base = e->data.function.args[0]; fmpq_set_si(x, 1, 2); ok = 1;
    } else if (head_is(e, "Power") && e->data.function.arg_count == 2 &&
               head_is(e->data.function.args[1], "Rational")) {
        base = e->data.function.args[0];
        ok = fmpq_from_expr(e->data.function.args[1], x);
    }
    fmpq_t r; fmpq_init(r);
    if (ok) ok = fmpq_from_expr(base, r) && !fmpq_is_zero(r);
    if (ok) {
        if (fmpq_sgn(r) < 0) {                   /* principal branch: |r|^x (-1)^x */
            fmpq_neg(r, r);
            fmpq_set(t->u, x);
        }
        ok = kf_atom_add_int(t, fmpq_numref(r), 1, x) &&
             kf_atom_add_int(t, fmpq_denref(r), -1, x);
    }
    fmpq_clear(r); fmpq_clear(x);
    return ok;
}

/* The product algebra A: generator 0 is z (when phi > 1), then the primes. */
typedef struct {
    ulong  N; slong phi;           /* z = exp(2 Pi I/N), phi = deg Phi_N      */
    fmpz_poly_t Phi;
    slong  np;
    fmpz   ell[KF_MAX_PRIMES];
    slong  n[KF_MAX_PRIMES], stride[KF_MAX_PRIMES];
    slong  D;
} KfAlg;

/* out += c * gen * in  (gen -1 = z, else a prime index).  `in` has length D. */
static void kf_mul_gen(fmpq* out, const fmpq* in, const KfAlg* A, slong gen, const fmpq_t c) {
    fmpq_t s; fmpq_init(s);
    for (slong idx = 0; idx < A->D; idx++) {
        if (fmpq_is_zero(in + idx)) continue;
        fmpq_mul(s, in + idx, c);
        if (gen < 0) {
            slong j = idx % A->phi;
            if (j + 1 < A->phi) fmpq_add(out + idx + 1, out + idx + 1, s);
            else {                               /* z^phi = -sum_k Phi_k z^k */
                fmpq_t t; fmpq_init(t);
                for (slong k = 0; k < A->phi; k++) {
                    fmpq_mul_fmpz(t, s, A->Phi->coeffs + k);
                    fmpq_sub(out + idx - j + k, out + idx - j + k, t);
                }
                fmpq_clear(t);
            }
        } else {
            slong st = A->stride[gen], n = A->n[gen];
            slong tt = (idx / st) % n;
            if (tt + 1 < n) fmpq_add(out + idx + st, out + idx + st, s);
            else {                               /* g^n = ell */
                fmpq_mul_fmpz(s, s, A->ell + gen);
                fmpq_add(out + idx - (n - 1) * st, out + idx - (n - 1) * st, s);
            }
        }
    }
    fmpq_clear(s);
}

/* v <- gen^k v (k >= 0). */
static void kf_pow_gen(fmpq* v, const KfAlg* A, slong gen, ulong k) {
    fmpq* w = _fmpq_vec_init(A->D);
    fmpq_t one; fmpq_init(one); fmpq_one(one);
    for (ulong i = 0; i < k; i++) {
        kf_vec_zero(w, A->D);
        kf_mul_gen(w, v, A, gen, one);
        kf_vec_set(v, w, A->D);
    }
    fmpq_clear(one);
    _fmpq_vec_clear(w, A->D);
}

/* The atom as a vector in the product basis of A. */
static void kf_atom_vec(fmpq* v, const KfAtom* t, const KfAlg* A) {
    kf_vec_zero(v, A->D);
    fmpq_t c; fmpq_init(c); fmpq_set(c, t->c);
    /* (-1)^u = z^(u N / 2) */
    fmpq_t k; fmpq_init(k);
    fmpq_mul_ui(k, t->u, A->N);
    fmpq_div_2exp(k, k, 1);
    ulong zk = fmpz_fdiv_ui(fmpq_numref(k), A->N);   /* k is an integer here */
    /* each prime: ell^(a) = ell^s g^tt with a n = s n + tt */
    slong tts[KF_MAX_PRIMES];
    for (slong i = 0; i < A->np; i++) tts[i] = 0;
    for (slong j = 0; j < t->np; j++) {
        if (fmpz_is_zero(t->ell + j)) continue;  /* a rational power, folded into c */
        slong gi = 0;
        while (gi < A->np && !fmpz_equal(A->ell + gi, t->ell + j)) gi++;
        if (gi == A->np) continue;
        fmpq_t m; fmpq_init(m);
        fmpq_mul_si(m, t->a + j, A->n[gi]);     /* an integer by construction */
        fmpz_t s, r, nz; fmpz_init(s); fmpz_init(r); fmpz_init(nz);
        fmpz_set_si(nz, A->n[gi]);
        fmpz_fdiv_qr(s, r, fmpq_numref(m), nz);
        fmpz_clear(nz);
        tts[gi] = (slong)fmpz_get_si(r);
        fmpz_t pw; fmpz_init(pw);
        slong sv = fmpz_get_si(s);
        fmpz_pow_ui(pw, A->ell + gi, (ulong)(sv < 0 ? -sv : sv));
        if (sv >= 0) fmpq_mul_fmpz(c, c, pw); else fmpq_div_fmpz(c, c, pw);
        fmpz_clear(pw); fmpz_clear(s); fmpz_clear(r); fmpq_clear(m);
    }
    if (A->phi > 1) {
        fmpq_set(v + 0, c);
        kf_pow_gen(v, A, -1, zk);
    } else {
        /* N <= 2: z = +-1 is rational */
        if (A->N == 2 && (zk & 1)) fmpq_neg(c, c);
        fmpq_set(v + 0, c);
    }
    for (slong i = 0; i < A->np; i++) kf_pow_gen(v, A, i, (ulong)tts[i]);
    fmpq_clear(k); fmpq_clear(c);
}

/* (sum_k p_k x^k)(z) for a rational polynomial p and a ball z. */
static void kf_acb_eval_fmpq_poly(acb_t res, const fmpq_poly_t p, const acb_t z, slong prec) {
    acb_zero(res);
    fmpq_t c; fmpq_init(c);
    arb_t cr; arb_init(cr);
    for (slong k = fmpq_poly_degree(p); k >= 0; k--) {
        acb_mul(res, res, z, prec);
        fmpq_poly_get_coeff_fmpq(c, p, k);
        arb_set_fmpq(cr, c, prec);
        arb_add(acb_realref(res), acb_realref(res), cr, prec);
    }
    arb_clear(cr); fmpq_clear(c);
}

/* The subfield K = Q(atoms) of L = Q[x]/(fL) = Q(phiL), and the atoms' images
 * in it.  `aimg[i]` is atom i as an element of L.  The generators of the
 * product algebra may generate more than the atoms do ((-2)^(1/3) alone lives
 * in Q(2^(1/3), (-1)^(1/3)) but generates a cubic field), and the field the
 * caller asks for is Q(atoms) itself -- its consumers write answers back in the
 * product basis of the atoms, which spans only that.  A primitive element of K
 * is theta = sum c_i atom_i for all but finitely many c: it is one exactly
 * when every atom lies in Q(theta), i.e. is a Q-combination of 1, theta, ..,
 * theta^(m-1), m the degree of theta's minimal polynomial.  K's generator is
 * then the algebraic integer phi = Dn theta (Dn clearing the denominators of
 * that minimal polynomial), the form the qqbar path emits.  Returns the List
 * of images, or NULL. */
static Expr* kf_subfield_images(const fmpq_poly_struct* aimg, size_t na,
                                const fmpq_poly_t fL, const qqbar_t phiL) {
    const slong dL = fmpq_poly_degree(fL);
    Expr* result = NULL;
    fmpq_poly_t th, pw, mK, col;
    fmpq_poly_init(th); fmpq_poly_init(pw); fmpq_poly_init(mK); fmpq_poly_init(col);
    fmpq_mat_t M, R, B, Y;
    fmpq_mat_init(M, dL, dL); fmpq_mat_init(R, dL, (slong)na);
    fmpq_mat_init(B, 1, 1); fmpq_mat_init(Y, 1, 1);
    for (size_t i = 0; i < na; i++)
        for (slong r = 0; r < dL; r++)
            fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(R, r, (slong)i), aimg + i, r);
    fmpq_t c; fmpq_init(c);
    slong m = 0;
    int okK = 0;
    for (int trial = 0; trial < 8 && !okK; trial++) {
        fmpq_poly_zero(th);
        slong idx = 0;
        for (size_t i = 0; i < na; i++) {
            if (fmpq_poly_degree(aimg + i) < 1) continue;     /* a rational atom */
            fmpq_set_si(c, 1 + idx * (trial + 1), 1);
            fmpq_poly_scalar_mul_fmpq(pw, aimg + i, c);
            fmpq_poly_add(th, th, pw);
            idx++;
        }
        if (idx == 0) { m = 1; okK = 1; break; }             /* K = Q */
        /* multiplication by theta on L, and theta's minimal polynomial */
        fmpq_poly_one(pw);
        for (slong j = 0; j < dL; j++) {
            fmpq_poly_mul(col, th, pw);
            fmpq_poly_rem(col, col, fL);
            for (slong r = 0; r < dL; r++)
                fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(M, r, j), col, r);
            fmpq_poly_shift_left(pw, pw, 1);
        }
        fmpq_mat_minpoly(mK, M);
        m = fmpq_poly_degree(mK);
        if (m < 1) continue;
        /* B = [1, theta, .., theta^(m-1)] in L's power basis */
        fmpq_mat_clear(B); fmpq_mat_init(B, dL, m);
        fmpq_mat_clear(Y); fmpq_mat_init(Y, m, (slong)na);
        fmpq_poly_one(pw);
        for (slong j = 0; j < m; j++) {
            for (slong r = 0; r < dL; r++)
                fmpq_poly_get_coeff_fmpq(fmpq_mat_entry(B, r, j), pw, r);
            fmpq_poly_mul(pw, pw, th);
            fmpq_poly_rem(pw, pw, fL);
        }
        okK = fmpq_mat_can_solve(Y, B, R);
    }
    if (okK && m == 1) {
        /* every atom is rational */
        Expr** items = calloc(na, sizeof(Expr*));
        for (size_t i = 0; i < na; i++) {
            fmpq_poly_get_coeff_fmpq(c, aimg + i, 0);
            items[i] = expr_from_fmpq(c);
        }
        result = expr_new_function(expr_new_symbol(SYM_List), items, na);
        free(items);
    } else if (okK && m == dL) {
        /* the atoms generate all of L: keep its generator phiL, whose images
         * are already in hand */
        Expr* gen = qqbar_to_expr(phiL);
        if (gen && head_is(gen, "Root")) qqbar_cache_put(gen, expr_hash(gen), phiL);
        Expr** items = calloc(na, sizeof(Expr*));
        for (size_t i = 0; i < na; i++) items[i] = poly_to_algnum_lazy(phiL, &gen, aimg + i, dL);
        if (gen) expr_free(gen);
        result = expr_new_function(expr_new_symbol(SYM_List), items, na);
        free(items);
    } else if (okK) {
        /* phi = Dn theta:  x^m + sum a_k x^k  ->  x^m + sum a_k Dn^(m-k) x^k */
        fmpz_t Dn, pk; fmpz_init(Dn); fmpz_init(pk);
        fmpz_one(Dn);
        for (slong k = 0; k < m; k++) {
            fmpq_poly_get_coeff_fmpq(c, mK, k);
            fmpz_lcm(Dn, Dn, fmpq_denref(c));
        }
        fmpz_poly_t fphi; fmpz_poly_init(fphi);
        fmpz_poly_set_coeff_ui(fphi, m, 1);
        for (slong k = 0; k < m; k++) {
            fmpq_poly_get_coeff_fmpq(c, mK, k);
            fmpz_pow_ui(pk, Dn, (ulong)(m - k));
            fmpq_mul_fmpz(c, c, pk);                      /* an integer */
            fmpz_poly_set_coeff_fmpz(fphi, k, fmpq_numref(c));
        }
        qqbar_t phi; qqbar_init(phi);
        int have = 0;
        acb_t z, zt; acb_init(z); acb_init(zt);
        for (slong prec = 64; prec <= 4096 && !have; prec *= 2) {
            qqbar_get_acb(z, phiL, prec);
            kf_acb_eval_fmpq_poly(zt, th, z, prec);
            acb_mul_fmpz(zt, zt, Dn, prec);
            acb_t enc; acb_init(enc);
            if (_qqbar_validate_uniqueness(enc, fphi, zt, 2 * prec)) {
                fmpz_poly_set(QQBAR_POLY(phi), fphi);
                acb_set(QQBAR_ENCLOSURE(phi), enc);
                have = 1;
            }
            acb_clear(enc);
        }
        acb_clear(z); acb_clear(zt);
        if (have) {
            Expr* gen = NULL;
            if (m >= 2) {
                /* phi's surface form, and the qqbar it denotes recorded against
                 * it: every consumer reads the field back through to_qqbar(gen),
                 * which would otherwise isolate and re-sort phi's conjugates a
                 * second time */
                gen = qqbar_to_expr(phi);
                if (gen && head_is(gen, "Root")) qqbar_cache_put(gen, expr_hash(gen), phi);
            }
            Expr** items = calloc(na, sizeof(Expr*));
            fmpq_poly_t img; fmpq_poly_init(img);
            for (size_t i = 0; i < na; i++) {
                /* atom = sum_j y_j theta^j = sum_j (y_j / Dn^j) phi^j */
                fmpq_poly_zero(img);
                fmpz_one(pk);
                for (slong j = 0; j < m; j++) {
                    fmpq_div_fmpz(c, fmpq_mat_entry(Y, j, (slong)i), pk);
                    fmpq_poly_set_coeff_fmpq(img, j, c);
                    fmpz_mul(pk, pk, Dn);
                }
                items[i] = poly_to_algnum_lazy(phi, &gen, img, m);
            }
            fmpq_poly_clear(img);
            if (gen) expr_free(gen);
            result = expr_new_function(expr_new_symbol(SYM_List), items, na);
            free(items);
        }
        qqbar_clear(phi);
        fmpz_poly_clear(fphi);
        fmpz_clear(Dn); fmpz_clear(pk);
    }
    fmpq_clear(c);
    fmpq_mat_clear(M); fmpq_mat_clear(R); fmpq_mat_clear(B); fmpq_mat_clear(Y);
    fmpq_poly_clear(th); fmpq_poly_clear(pw); fmpq_poly_clear(mK); fmpq_poly_clear(col);
    return result;
}

static Expr* kf_number_field_common(const Expr* const* as, size_t na) {
    if (na == 0 || na > 64) return NULL;
    KfAtom* at = malloc(sizeof(KfAtom) * na);
    for (size_t i = 0; i < na; i++) kf_atom_init(at + i);
    int ok = 1;
    for (size_t i = 0; i < na && ok; i++) ok = kf_atom_parse(as[i], at + i);

    KfAlg A; memset(&A, 0, sizeof A);
    fmpz_poly_init(A.Phi);
    for (slong i = 0; i < KF_MAX_PRIMES; i++) fmpz_init(A.ell + i);
    A.N = 1;
    /* N = lcm of the denominators of u/2; n_ell = lcm of the exponent denominators */
    for (size_t i = 0; i < na && ok; i++) {
        fmpq_t h; fmpq_init(h);
        fmpq_div_2exp(h, at[i].u, 1);
        if (!fmpz_fits_si(fmpq_denref(h)) || fmpz_get_si(fmpq_denref(h)) > 1024) ok = 0;
        else A.N = (ulong)kf_lcm(A.N, fmpz_get_ui(fmpq_denref(h)));
        fmpq_clear(h);
        if (A.N > 1024) ok = 0;
        for (slong j = 0; j < at[i].np && ok; j++) {
            if (fmpz_is_one(fmpq_denref(at[i].a + j))) continue;  /* rational power */
            if (!fmpz_fits_si(fmpq_denref(at[i].a + j))) { ok = 0; break; }
            slong d = fmpz_get_si(fmpq_denref(at[i].a + j));
            slong gi = 0;
            while (gi < A.np && !fmpz_equal(A.ell + gi, at[i].ell + j)) gi++;
            if (gi == A.np) {
                if (A.np == KF_MAX_PRIMES) { ok = 0; break; }
                fmpz_set(A.ell + gi, at[i].ell + j); A.n[gi] = 1; A.np++;
            }
            A.n[gi] = (slong)kf_lcm((ulong)A.n[gi], (ulong)d);
            if (A.n[gi] > KF_MAX_DIM) ok = 0;
        }
    }
    /* a prime whose exponents are all integers is rational: fold it into c */
    for (size_t i = 0; i < na && ok; i++) {
        for (slong j = 0; j < at[i].np; j++) {
            slong gi = 0;
            while (gi < A.np && !fmpz_equal(A.ell + gi, at[i].ell + j)) gi++;
            if (gi < A.np) continue;
            /* every exponent of this prime is an integer */
            fmpz_t pw; fmpz_init(pw);
            slong e = fmpz_get_si(fmpq_numref(at[i].a + j));
            fmpz_pow_ui(pw, at[i].ell + j, (ulong)(e < 0 ? -e : e));
            if (e >= 0) fmpq_mul_fmpz(at[i].c, at[i].c, pw); else fmpq_div_fmpz(at[i].c, at[i].c, pw);
            fmpz_clear(pw);
            fmpz_zero(at[i].ell + j);            /* dropped: never matches a prime */
        }
    }
    Expr* result = NULL;
    if (ok) {
        A.phi = 1;
        if (A.N > 2) {
            fmpz_poly_cyclotomic(A.Phi, A.N);
            A.phi = fmpz_poly_degree(A.Phi);
        }
        A.D = A.phi;
        for (slong i = 0; i < A.np && ok; i++) {
            A.stride[i] = A.D;
            A.D *= A.n[i];
            if (A.D > KF_MAX_DIM) ok = 0;
        }
        if (A.D < 2) ok = 0;                     /* nothing algebraic: let qqbar answer */
    }
    if (ok) {
        const slong D = A.D;
        slong ngen = A.np + (A.phi > 1 ? 1 : 0);
        fmpq_mat_t P, Pinv; fmpq_mat_init(P, D, D); fmpq_mat_init(Pinv, D, D);
        fmpq* v = _fmpq_vec_init(D);
        fmpq* w = _fmpq_vec_init(D);
        fmpq* u = _fmpq_vec_init(D);
        slong* cs = malloc(sizeof(slong) * (size_t)ngen);
#define KF_GEN(g) ((A.phi > 1) ? (g) - 1 : (g))      /* z is generator 0 */
        int found = 0;
        for (int trial = 0; trial < 8 && !found; trial++) {
            for (slong g = 0; g < ngen; g++) cs[g] = 1 + g * (trial + 1);
            kf_vec_zero(v, D); fmpq_one(v + 0);
            for (slong j = 0; j < D; j++) {
                for (slong r = 0; r < D; r++) fmpq_set(fmpq_mat_entry(P, r, j), v + r);
                kf_vec_zero(w, D);
                fmpq_t cg; fmpq_init(cg);
                for (slong g = 0; g < ngen; g++) {
                    fmpq_set_si(cg, cs[g], 1);
                    kf_mul_gen(w, v, &A, KF_GEN(g), cg);
                }
                /* + the pairwise products gen_g gen_h: a theta linear in the
                 * generators has conjugates with equal real parts beyond the
                 * complex-conjugate pairs (Re(2^(1/3) w) = Re(2^(1/3) w^2)),
                 * and each such tie costs the root ordering behind theta's
                 * Root index an exact qqbar comparison */
                fmpq_one(cg);
                for (slong g = 0; g < ngen; g++)
                    for (slong h = g + 1; h < ngen; h++) {
                        kf_vec_zero(u, D);
                        kf_mul_gen(u, v, &A, KF_GEN(h), cg);
                        kf_mul_gen(w, u, &A, KF_GEN(g), cg);
                    }
                fmpq_clear(cg);
                kf_vec_set(v, w, D);
            }
            found = fmpq_mat_inv(Pinv, P);
        }
        /* v = theta^D; m_A(x) = x^D - sum y_j x^j, y = Pinv v */
        fmpz_poly_t mA; fmpz_poly_init(mA);
        acb_t th; acb_init(th);
        fmpq_poly_t yq; fmpq_poly_init(yq);
        if (found) {
            fmpq_poly_zero(yq);
            fmpq_t s, t; fmpq_init(s); fmpq_init(t);
            for (slong r = 0; r < D; r++) {
                fmpq_zero(s);
                for (slong j = 0; j < D; j++) {
                    fmpq_mul(t, fmpq_mat_entry(Pinv, r, j), v + j);
                    fmpq_add(s, s, t);
                }
                fmpq_neg(s, s);
                fmpq_poly_set_coeff_fmpq(yq, r, s);
            }
            fmpq_one(s);
            fmpq_poly_set_coeff_fmpq(yq, D, s);
            fmpq_clear(s); fmpq_clear(t);
            /* theta is an algebraic integer: m_A is monic over Z */
            if (!fmpz_is_one(fmpq_poly_denref(yq))) found = 0;
            else fmpq_poly_get_numerator(mA, yq);
        }
        /* pick the factor of m_A that vanishes at the principal value */
        qqbar_t phi; qqbar_init(phi);
        fmpz_poly_t f; fmpz_poly_init(f);
        if (found) {
            fmpz_poly_factor_t fac; fmpz_poly_factor_init(fac);
            fmpz_poly_factor(fac, mA);
            int chosen = 0;
            for (slong prec = 64; prec <= 4096 && !chosen; prec *= 2) {
                acb_zero(th);
                acb_ptr gv = _acb_vec_init(ngen);
                acb_t g; acb_init(g);
                for (slong gidx = 0; gidx < ngen; gidx++) {
                    slong gen = KF_GEN(gidx);
                    if (gen < 0) {
                        acb_set_si(gv + gidx, 2); acb_div_ui(gv + gidx, gv + gidx, A.N, prec);
                        acb_exp_pi_i(gv + gidx, gv + gidx, prec);
                    } else {
                        arb_t rr; arb_init(rr);
                        arb_set_fmpz(rr, A.ell + gen);
                        arb_root_ui(rr, rr, (ulong)A.n[gen], prec);
                        acb_set_arb(gv + gidx, rr);
                        arb_clear(rr);
                    }
                    acb_mul_si(g, gv + gidx, cs[gidx], prec);
                    acb_add(th, th, g, prec);
                }
                for (slong gg = 0; gg < ngen; gg++)
                    for (slong h = gg + 1; h < ngen; h++) {
                        acb_mul(g, gv + gg, gv + h, prec);
                        acb_add(th, th, g, prec);
                    }
                acb_clear(g);
                _acb_vec_clear(gv, ngen);
                /* candidates: factors whose value at theta contains 0 */
                slong cand = -1, ncand = 0;
                acb_t val; acb_init(val);
                for (slong k = 0; k < fac->num; k++) {
                    arb_fmpz_poly_evaluate_acb(val, fac->p + k, th, prec);
                    if (acb_contains_zero(val)) { cand = k; ncand++; }
                }
                acb_clear(val);
                if (ncand != 1) continue;
                /* phi = the root of that factor that theta approximates: certify
                 * a ball around theta isolating exactly one root (no ordering of
                 * the conjugates, which is what qqbar_roots_fmpz_poly costs) */
                acb_t enc; acb_init(enc);
                if (_qqbar_validate_uniqueness(enc, fac->p + cand, th, 2 * prec)) {
                    fmpz_poly_set(QQBAR_POLY(phi), fac->p + cand);
                    acb_set(QQBAR_ENCLOSURE(phi), enc);
                    fmpz_poly_set(f, fac->p + cand);
                    chosen = 1;
                }
                acb_clear(enc);
            }
            fmpz_poly_factor_clear(fac);
            found = chosen;
        }
        if (found) {
            /* the atoms' images in L = Q[x]/(f) = Q(generators of A) */
            fmpq_poly_t fq; fmpq_poly_init(fq);
            fmpq_poly_set_fmpz_poly(fq, f);
            fmpq_poly_struct* aimg = flint_malloc(sizeof(fmpq_poly_struct) * na);
            fmpq_t s, t; fmpq_init(s); fmpq_init(t);
            for (size_t i = 0; i < na; i++) {
                fmpq_poly_init(aimg + i);
                kf_atom_vec(w, at + i, &A);
                for (slong r = 0; r < D; r++) {
                    fmpq_zero(s);
                    for (slong j = 0; j < D; j++) {
                        if (fmpq_is_zero(w + j)) continue;
                        fmpq_mul(t, fmpq_mat_entry(Pinv, r, j), w + j);
                        fmpq_add(s, s, t);
                    }
                    fmpq_poly_set_coeff_fmpq(aimg + i, r, s);
                }
                fmpq_poly_rem(aimg + i, aimg + i, fq);
            }
            result = kf_subfield_images(aimg, na, fq, phi);
            for (size_t i = 0; i < na; i++) fmpq_poly_clear(aimg + i);
            flint_free(aimg);
            fmpq_clear(s); fmpq_clear(t);
            fmpq_poly_clear(fq);
        }
        fmpz_poly_clear(f);
        qqbar_clear(phi);
        fmpq_poly_clear(yq);
        acb_clear(th);
        fmpz_poly_clear(mA);
        free(cs);
        _fmpq_vec_clear(v, D); _fmpq_vec_clear(w, D); _fmpq_vec_clear(u, D);
#undef KF_GEN
        fmpq_mat_clear(P); fmpq_mat_clear(Pinv);
    }
    fmpz_poly_clear(A.Phi);
    for (slong i = 0; i < KF_MAX_PRIMES; i++) fmpz_clear(A.ell + i);
    for (size_t i = 0; i < na; i++) kf_atom_clear(at + i);
    free(at);
    return result;
}

Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n,
                                         int smallest) {
    /* radicals of rationals and roots of unity: the exact product-algebra
     * construction above, no qqbar search */
    Expr* fast = kf_number_field_common(as, n);
    if (fast) return fast;
    return flint_qqbar_to_number_field_common_fn(NULL, as, n, smallest, false);
}

Expr* flint_qqbar_to_number_field_common_fn(Expr* (*get)(const void*, size_t, bool*),
                                            const void* ctx, size_t n, int smallest, bool out_ll) {
    (void)smallest;                            /* Automatic and All coincide here */
    if (!ctx || n == 0) return NULL;

    qqbar_ptr vals = _qqbar_vec_init((slong)n);
    int ok = 1;
    for (size_t i = 0; i < n && ok; i++) {
        /* each number read on demand: from the array, or through `get` */
        bool own = false;
        const Expr* a = get ? get(ctx, i, &own) : ((const Expr* const*)ctx)[i];
        if (!to_qqbar(a, vals + i)) ok = 0;
        if (own) expr_free((Expr*)a);
    }

    /* Primitive element of Q(alpha, b), chosen by DEGREE rather than by trial
     * membership.
     *
     * alpha + c*b lies in Q(alpha, b) for every c, so Q(alpha + c*b) is always a
     * subfield of the compositum and equals it exactly when the degrees agree.
     * The largest deg(alpha + c*b) over a range of c is therefore the compositum
     * degree (all but finitely many c give a primitive element), so the best c
     * follows from degrees alone -- cheaply, and with no membership test that
     * could fail for a reason other than non-membership.
     *
     * That last point was a real defect, not a tidy-up.  The previous form took
     * the first c for which in_field(b, cand) && in_field(alpha, cand) held; those
     * run at 64 bits of working precision, and qqbar_express_in_field_esc
     * deliberately does not escalate when the generator has degree <= 6 (its
     * "common path" gate).  A relation needing more than 64 bits was therefore
     * read as non-membership, every c was rejected, and the whole compositum
     * failed -- for degree 6 SPECIFICALLY, since degrees 2-5 resolve inside 64
     * bits and degrees 7+ are allowed to escalate.  ToNumberField[{2^(1/6),
     * 2^(1/3)}] came back unevaluated and with it every multivariate gcd over a
     * degree-6 radical field, which silently answered 1. */
    qqbar_t alpha, cand, ca; qqbar_init(alpha); qqbar_init(cand); qqbar_init(ca);
    if (ok) qqbar_set(alpha, vals + 0);
    for (size_t i = 1; i < n && ok; i++) {
        /* [Q(alpha, b) : Q] <= deg(alpha) deg(b), so a candidate of exactly that
         * degree is already the compositum and the scan can stop: the remaining
         * candidates cannot beat it.  Without this every atom paid all sixteen
         * qqbar additions -- each a composed resultant of degree deg(alpha)deg(b)
         * that is then factored over Z -- although c = 1 is primitive almost
         * always; for Q(I, Sqrt[3], 2^(1/3)) that was ~0.25 s per field build, and
         * a Together over that field builds several. */
        slong bestdeg = 0, bestc = 0;
        const slong dmax = qqbar_degree(alpha) * qqbar_degree(vals + i);
        for (slong c = 1; c <= 16; c++) {
            qqbar_mul_si(ca, vals + i, c);
            qqbar_add(cand, alpha, ca);
            slong d = qqbar_degree(cand);
            if (d > bestdeg && d <= QQBAR_DEGREE_CAP) { bestdeg = d; bestc = c; }
            if (bestdeg == dmax) break;
        }
        if (bestc == 0) { ok = 0; break; }
        qqbar_mul_si(ca, vals + i, bestc);
        qqbar_add(cand, alpha, ca);
        qqbar_set(alpha, cand);
    }
    qqbar_clear(cand); qqbar_clear(ca);

    Expr* result = NULL;
    if (ok) {
        fmpz_t lc; fmpz_init(lc);
        qqbar_t phi; qqbar_init(phi);
        algint_generator(alpha, phi, lc);
        fmpz_clear(lc);
        slong nn = qqbar_degree(phi);

        /* calloc, not malloc: the loop below stops early on the first
         * express-in-field miss (`&& good`), so items[j+1 .. n-1] are never
         * assigned.  The error-path cleanup frees every non-NULL entry, so an
         * uninitialised tail would free wild pointers -- a heap corruption that
         * surfaced as a SEGV under ASan and as non-deterministic wrong numerics
         * (a later Root/AlgebraicNumber reading corrupted memory) in batch runs. */
        Expr** items = calloc(n, sizeof(Expr*));
        int good = 1;
        Expr* gen = NULL;                        /* built once, on demand */
        for (size_t i = 0; i < n && good; i++) {
            fmpq_poly_t f; fmpq_poly_init(f);
            if (qqbar_express_in_field_retry(f, phi, vals + i))
                items[i] = poly_to_algnum_lazy(phi, &gen, f, nn);
            else { items[i] = NULL; good = 0; }
            fmpq_poly_clear(f);
        }
        if (gen) expr_free(gen);
        if (good) {
            if (out_ll) {                      /* straight into a LinkedList */
                LinkedList* ll = ll_new();
                for (size_t t = 0; t < n; t++) ll_push_back(ll, items[t]);
                result = expr_new_linkedlist(ll);
            } else
            result = expr_new_function(expr_new_symbol(SYM_List), items, n);
        } else {
            for (size_t i = 0; i < n; i++) if (items[i]) expr_free(items[i]);
        }
        free(items);
        qqbar_clear(phi);
    }
    qqbar_clear(alpha);
    _qqbar_vec_clear(vals, (slong)n);
    return result;
}

/* NumberFieldIntegralBasis[a]: a Z-module basis of the ring of integers O_K of
 * K = Q(a).  Builds the algebraic-integer generator phi of Q(a) (monic minimal
 * polynomial f), hands f to the number-field layer -- nf_field_create runs
 * Dedekind at every ramified prime and, where Z[phi] is not maximal, enlarges it
 * to O_K by Round 2 (Pohst-Zassenhaus) -- then renders each basis row
 * omega_i = (1/D) * sum_j W[i][j] phi^j back as a canonical AlgebraicNumber (or a
 * plain Integer/Rational when the row is rational).  Returns a fresh List, or
 * NULL to decline: `a` is not a constant algebraic number, its degree exceeds the
 * cap, or the field layer cannot certify O_K (disc will not factor, or a prime is
 * too large for a single-word modulus). */
Expr* flint_qqbar_integral_basis(const Expr* a) {
    if (!a) return NULL;

    qqbar_t alpha; qqbar_init(alpha);
    if (!to_qqbar(a, alpha)) { qqbar_clear(alpha); return NULL; }

    fmpz_t lc; fmpz_init(lc);
    qqbar_t phi; qqbar_init(phi);
    algint_generator(alpha, phi, lc);        /* phi = lc*alpha, an algebraic integer */
    fmpz_clear(lc);
    slong n = qqbar_degree(phi);

    /* Q(a) = Q: the ring of integers is Z, with integral basis {1}. */
    if (n <= 1) {
        qqbar_clear(phi); qqbar_clear(alpha);
        Expr* one = expr_new_integer(1);
        return expr_new_function(expr_new_symbol(SYM_List), &one, 1);
    }
    if (n > QQBAR_DEGREE_CAP) { qqbar_clear(phi); qqbar_clear(alpha); return NULL; }

    /* Monic integer defining polynomial f[0..n] of phi (coeffs[n] == 1). */
    mpz_t* coeffs = malloc(sizeof(mpz_t) * (size_t)(n + 1));
    for (slong i = 0; i <= n; i++) {
        mpz_init(coeffs[i]);
        fmpz_t ci; fmpz_init(ci);
        fmpz_poly_get_coeff_fmpz(ci, QQBAR_POLY(phi), i);
        fmpz_get_mpz(coeffs[i], ci);
        fmpz_clear(ci);
    }

    NumberField* K = nf_field_create((const mpz_t*)coeffs, (int)n);
    for (slong i = 0; i <= n; i++) mpz_clear(coeffs[i]);
    free(coeffs);
    if (!K) { qqbar_clear(phi); qqbar_clear(alpha); return NULL; }

    /* The Round-2 basis is a correct Z-basis of the numerator lattice L (so
     * (1/D)*L = O_K), but its rows are arbitrarily ordered and unreduced.  Put
     * it in the standard integral-basis presentation: the lower-triangular
     * Hermite normal form in the theta-power basis, so omega_0 = 1 and omega_k
     * has degree exactly k.  FLINT's fmpz_mat_hnf gives the UPPER-triangular
     * HNF; the lower-triangular one is its 180-degree rotation applied to the
     * column-reversed lattice, i.e. numerator[i][j] = H[n-1-i][n-1-j] where
     * H = HNF(W with columns reversed).  (Monogenic W = I stays {1, theta, ..}.) */
    const mpz_t* W = nf_ok_basis(K);
    mpz_t D; mpz_init(D); nf_ok_denom(K, D);
    fmpz_t Dz; fmpz_init(Dz); fmpz_set_mpz(Dz, D);

    fmpz_mat_t A, H;
    fmpz_mat_init(A, (slong)n, (slong)n);
    fmpz_mat_init(H, (slong)n, (slong)n);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++)
            fmpz_set_mpz(fmpz_mat_entry(A, i, j), W[i * n + (n - 1 - j)]);
    fmpz_mat_hnf(H, A);

    Expr** items = malloc(sizeof(Expr*) * (size_t)n);
    fmpq_t coef; fmpq_init(coef);
    Expr* gen = NULL;                 /* see poly_to_algnum_lazy */
    for (slong i = 0; i < n; i++) {
        fmpq_poly_t row; fmpq_poly_init(row);
        for (slong j = 0; j < n; j++) {
            fmpq_set_fmpz_frac(coef, fmpz_mat_entry(H, n - 1 - i, n - 1 - j), Dz);
            fmpq_poly_set_coeff_fmpq(row, j, coef);
        }
        items[i] = poly_to_algnum_lazy(phi, &gen, row, n); /* rational row -> Integer/Rational */
        fmpq_poly_clear(row);
    }
    if (gen) expr_free(gen);
    fmpq_clear(coef);
    fmpz_mat_clear(A); fmpz_mat_clear(H);
    fmpz_clear(Dz); mpz_clear(D);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);

    nf_field_free(K);
    qqbar_clear(phi); qqbar_clear(alpha);
    return result;
}

/* AlgebraicIntegerQ[x]: 1 if x is an algebraic integer, 0 if it is not, -1 when
 * undecided (FLINT compiled out).  x is an algebraic integer iff its primitive
 * integer minimal polynomial is monic -- i.e. the (positive) leading coefficient
 * of QQBAR_POLY(x) is 1.  A rational p/q has minimal polynomial q*x - p, so only
 * integers (q == 1) qualify; anything that is not a constant algebraic number (a
 * free symbol, Pi, ...) fails the to_qqbar conversion and is not an algebraic
 * integer -> 0. */
/* Is `e` an algebraic integer by ring closure alone, with no minimal
 * polynomial: integers, Gaussian integers, GoldenRatio (x^2 - x - 1), and the
 * sums, products and positive-integer powers of such -- the algebraic integers
 * are a ring -- plus b^(p/q) with p > 0 for such b (a root of the monic
 * y^q - b^p over the algebraic integers, hence integral). A sufficient test:
 * 0 means "not shown", never "not an algebraic integer". It decides the
 * degree-200 2^(1/200) and the degree-128 sum of seven prime square roots,
 * which pass to_qqbar's degree cap. */
static int alg_int_by_closure(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return 1;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_GoldenRatio;
    if (e->type != EXPR_FUNCTION) return 0;
    size_t n = e->data.function.arg_count;
    if (head_is(e, "Complex") && n == 2)
        return alg_int_by_closure(e->data.function.args[0]) &&
               alg_int_by_closure(e->data.function.args[1]);
    if (head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!alg_int_by_closure(e->data.function.args[i])) return 0;
        return 1;
    }
    if (head_is(e, "Sqrt") && n == 1) return alg_int_by_closure(e->data.function.args[0]);
    if (head_is(e, "Power") && n == 2) {
        const Expr* xe = e->data.function.args[1];
        int positive = 0;
        if (xe->type == EXPR_INTEGER) positive = xe->data.integer > 0;
        else if (head_is(xe, "Rational") && xe->data.function.arg_count == 2) {
            const Expr* p = xe->data.function.args[0];
            positive = (p->type == EXPR_INTEGER && p->data.integer > 0) ||
                       (p->type == EXPR_BIGINT && mpz_sgn(p->data.bigint) > 0);
        }
        return positive && alg_int_by_closure(e->data.function.args[0]);
    }
    return 0;
}

int flint_qqbar_algebraic_integer_q(const Expr* x) {
    if (!x) return 0;
    if (alg_int_by_closure(x)) return 1;
    qqbar_t v; qqbar_init(v);
    g_qqbar_cap_hit = 0;
    if (!to_qqbar(x, v)) {
        /* Not algebraic -> False; past the degree cap -> undecided. */
        int capped = g_qqbar_cap_hit;
        g_qqbar_cap_hit = 0;
        qqbar_clear(v);
        return capped ? -1 : 0;
    }
    slong d = qqbar_degree(v);
    fmpz_t lead; fmpz_init(lead);
    fmpz_poly_get_coeff_fmpz(lead, QQBAR_POLY(v), d);
    int r = fmpz_is_one(lead) ? 1 : 0;
    fmpz_clear(lead);
    qqbar_clear(v);
    return r;
}

/* AlgebraicNumberDenominator[x]: the smallest positive integer d such that d*x is
 * an algebraic integer.  Returns 1 and sets *out to that d (Integer / BigInt) on
 * success; 0 if x is not a constant algebraic number; -1 if FLINT is compiled out.
 *
 * This is NOT qqbar_denominator (the leading coefficient a_n of the primitive
 * integer minimal polynomial p(x) = sum_i a_i x^i): a_n is only an upper bound.
 * Writing the monic minimal polynomial of d*x, its x^i coefficient is
 * (a_i/a_n) d^{n-i}, so d*x is an algebraic integer iff for every i < n the
 * denominator q_i of a_i/a_n divides d^{n-i}.  Since q_i | a_n, the minimal d
 * divides a_n, and per prime P | a_n:
 *     v_P(d) = max_{i<n} ceil( v_P(q_i) / (n-i) ),  v_P(q_i)=max(0, v_P(a_n)-v_P(a_i))
 * (a zero coefficient contributes no constraint).  One factorisation of a_n --
 * small in practice -- then a valuation scan gives the exact answer.  Example:
 * 1/5 + Sqrt[2] has p = 25 x^2 - 10 x - 49, a_n = 25, yet d = 5. */
int flint_qqbar_algebraic_number_denominator(const Expr* x, Expr** out) {
    if (!x || !out) return 0;
    qqbar_t v; qqbar_init(v);
    if (!to_qqbar(x, v)) { qqbar_clear(v); return 0; }
    slong n = qqbar_degree(v);
    fmpz_t an; fmpz_init(an);
    fmpz_poly_get_coeff_fmpz(an, QQBAR_POLY(v), n);   /* positive leading coeff */

    fmpz_t d; fmpz_init_set_ui(d, 1);
    if (!fmpz_is_one(an)) {
        fmpz_factor_t fac; fmpz_factor_init(fac);
        fmpz_factor(fac, an);                          /* an > 0: primes are positive */
        fmpz_t ai, tmp, pk;
        fmpz_init(ai); fmpz_init(tmp); fmpz_init(pk);
        for (slong j = 0; j < fac->num; j++) {
            const fmpz* P = fac->p + j;
            slong vP_an = (slong) fac->exp[j];
            slong best = 0;
            for (slong i = 0; i < n; i++) {            /* i = 0 .. n-1 (i = n gives 1) */
                fmpz_poly_get_coeff_fmpz(ai, QQBAR_POLY(v), i);
                if (fmpz_is_zero(ai)) continue;        /* a_i/a_n = 0: no constraint */
                slong vP_ai = fmpz_remove(tmp, ai, P); /* valuation of a_i at P */
                slong vP_qi = vP_an - vP_ai;
                if (vP_qi <= 0) continue;              /* a_i/a_n already P-integral */
                slong ei   = n - i;
                slong need = (vP_qi + ei - 1) / ei;    /* ceil(vP_qi / ei) */
                if (need > best) best = need;
            }
            if (best > 0) {
                fmpz_pow_ui(pk, P, (ulong) best);
                fmpz_mul(d, d, pk);
            }
        }
        fmpz_clear(ai); fmpz_clear(tmp); fmpz_clear(pk);
        fmpz_factor_clear(fac);
    }

    *out = expr_from_fmpz(d);
    fmpz_clear(d); fmpz_clear(an); qqbar_clear(v);
    return 1;
}

/* Absolute field norm N_{Q(a)/Q}(a) = product of the roots of a's minimal
 * polynomial.  From the primitive integer minimal polynomial
 * P(x) = c_n x^n + ... + c_0 (content 1, c_n > 0), the monic-over-Q constant
 * term is c_0/c_n and the product of the n roots is (-1)^n c_0/c_n. */
static void qqbar_abs_norm(const qqbar_t v, fmpq_t out) {
    slong n = qqbar_degree(v);
    fmpz_t c0, cn; fmpz_init(c0); fmpz_init(cn);
    fmpz_poly_get_coeff_fmpz(c0, QQBAR_POLY(v), 0);
    fmpz_poly_get_coeff_fmpz(cn, QQBAR_POLY(v), n);   /* positive leading coeff */
    fmpq_set_fmpz_frac(out, c0, cn);                  /* c_0 / c_n, reduced */
    if (n & 1) fmpq_neg(out, out);                    /* (-1)^n */
    fmpz_clear(c0); fmpz_clear(cn);
}

/* AlgebraicNumberNorm[a] (theta == NULL): the absolute norm N_{Q(a)/Q}(a).
 * AlgebraicNumberNorm[a, Extension -> theta] (theta != NULL): the relative norm
 * N_{Q(theta)/Q}(a), defined when a lies in K = Q(theta).  By transitivity of
 * the norm in the tower Q <= Q(a) <= K,
 *     N_{K/Q}(a) = N_{Q(a)/Q}(a)^{[K:Q(a)]} = (absolute norm)^{n/d},
 * with n = [K:Q] = deg minpoly(theta), d = [Q(a):Q] = deg minpoly(a); a in K
 * forces d | n by the tower law.  Membership is decided by qqbar_express_in_field.
 * Returns: 1 with *out set (a fresh owned Integer/Rational); 0 when a or theta is
 * not a constant algebraic number; 2 when theta is given but a is not in Q(theta);
 * -1 when FLINT is compiled out (the #else stub). */
int flint_qqbar_algebraic_number_norm(const Expr* a, const Expr* theta, Expr** out) {
    if (!a || !out) return 0;
    qqbar_t av; qqbar_init(av);
    if (!to_qqbar(a, av)) { qqbar_clear(av); return 0; }

    fmpq_t norm; fmpq_init(norm);
    qqbar_abs_norm(av, norm);
    int rc = 1;

    if (theta) {                              /* relative norm over Q(theta) */
        qqbar_t tv; qqbar_init(tv);
        if (!to_qqbar(theta, tv)) {
            rc = 0;                           /* theta not a constant algebraic number */
        } else {
            slong d = qqbar_degree(av);
            fmpz_t lc; fmpz_init(lc);
            qqbar_t phi; qqbar_init(phi);
            algint_generator(tv, phi, lc);    /* Q(phi) = Q(theta) */
            fmpz_clear(lc);
            slong n = qqbar_degree(phi);
            fmpq_poly_t f; fmpq_poly_init(f);
            if (qqbar_express_in_field_retry(f, phi, av)) {
                fmpq_pow_si(norm, norm, n / d);  /* (absolute norm)^{[K:Q(a)]} */
            } else {
                rc = 2;                          /* a is not an element of Q(theta) */
            }
            fmpq_poly_clear(f);
            qqbar_clear(phi);
        }
        qqbar_clear(tv);
    }

    if (rc == 1) *out = expr_from_fmpq(norm);
    fmpq_clear(norm);
    qqbar_clear(av);
    return rc;
}

/* Absolute field trace Tr_{Q(a)/Q}(a) = sum of the roots of a's minimal
 * polynomial.  From the primitive integer minimal polynomial
 * P(x) = c_n x^n + ... + c_0 (content 1, c_n > 0), the monic-over-Q coefficient
 * of x^{n-1} is c_{n-1}/c_n and the sum of the n roots is -c_{n-1}/c_n.  For a
 * degree-1 minimal polynomial (an integer or rational) this reads c_0/c_1 and
 * gives a itself. */
static void qqbar_abs_trace(const qqbar_t v, fmpq_t out) {
    slong n = qqbar_degree(v);
    fmpz_t c1, cn; fmpz_init(c1); fmpz_init(cn);
    fmpz_poly_get_coeff_fmpz(c1, QQBAR_POLY(v), n - 1); /* coeff of x^{n-1} */
    fmpz_poly_get_coeff_fmpz(cn, QQBAR_POLY(v), n);     /* positive leading coeff */
    fmpq_set_fmpz_frac(out, c1, cn);                    /* c_{n-1} / c_n, reduced */
    fmpq_neg(out, out);                                 /* sum of roots = -c_{n-1}/c_n */
    fmpz_clear(c1); fmpz_clear(cn);
}

/* AlgebraicNumberTrace[a] (theta == NULL): the absolute trace Tr_{Q(a)/Q}(a).
 * AlgebraicNumberTrace[a, Extension -> theta] (theta != NULL): the relative trace
 * Tr_{Q(theta)/Q}(a), defined when a lies in K = Q(theta).  By transitivity of
 * the trace in the tower Q <= Q(a) <= K,
 *     Tr_{K/Q}(a) = [K:Q(a)] * Tr_{Q(a)/Q}(a) = (n/d) * (absolute trace),
 * with n = [K:Q] = deg minpoly(theta), d = [Q(a):Q] = deg minpoly(a); a in K
 * forces d | n by the tower law.  Membership is decided by qqbar_express_in_field.
 * (Contrast the norm, which is multiplicative and raises to the power n/d; the
 * trace is additive and scales by n/d.)  Returns: 1 with *out set (a fresh owned
 * Integer/Rational); 0 when a or theta is not a constant algebraic number; 2 when
 * theta is given but a is not in Q(theta); -1 when FLINT is compiled out. */
int flint_qqbar_algebraic_number_trace(const Expr* a, const Expr* theta, Expr** out) {
    if (!a || !out) return 0;
    qqbar_t av; qqbar_init(av);
    if (!to_qqbar(a, av)) { qqbar_clear(av); return 0; }

    fmpq_t trace; fmpq_init(trace);
    qqbar_abs_trace(av, trace);
    int rc = 1;

    if (theta) {                              /* relative trace over Q(theta) */
        qqbar_t tv; qqbar_init(tv);
        if (!to_qqbar(theta, tv)) {
            rc = 0;                           /* theta not a constant algebraic number */
        } else {
            slong d = qqbar_degree(av);
            fmpz_t lc; fmpz_init(lc);
            qqbar_t phi; qqbar_init(phi);
            algint_generator(tv, phi, lc);    /* Q(phi) = Q(theta) */
            fmpz_clear(lc);
            slong n = qqbar_degree(phi);
            fmpq_poly_t f; fmpq_poly_init(f);
            if (qqbar_express_in_field_retry(f, phi, av)) {
                fmpq_mul_si(trace, trace, n / d);  /* (absolute trace) * [K:Q(a)] */
            } else {
                rc = 2;                            /* a is not an element of Q(theta) */
            }
            fmpq_poly_clear(f);
            qqbar_clear(phi);
        }
        qqbar_clear(tv);
    }

    if (rc == 1) *out = expr_from_fmpq(trace);
    fmpq_clear(trace);
    qqbar_clear(av);
    return rc;
}

Expr* flint_qqbar_algnum_add(const Expr* a, const Expr* b) {
    return algnum_binop(a, b, 0);
}

Expr* flint_qqbar_algnum_mul(const Expr* a, const Expr* b) {
    return algnum_binop(a, b, 1);
}

Expr* flint_qqbar_algnum_pow(const Expr* a, long p) {
    if (!head_is(a, "AlgebraicNumber") || a->data.function.arg_count != 2) return NULL;
    const Expr* g = a->data.function.args[0];
    const fmpz_poly_struct* Mz; slong n;             /* minpoly + degree from the cache */
    if (!genfd_get(g, &Mz, &n, NULL)) return NULL;

    fmpq_poly_t M, base, acc, tmp;
    fmpq_poly_init(M); fmpq_poly_init(base); fmpq_poly_init(acc); fmpq_poly_init(tmp);
    fmpq_poly_set_fmpz_poly(M, Mz);

    Expr* result = NULL;
    if (coeffs_to_fmpq_poly(a->data.function.args[1], base)) {
        int ok = 1;
        unsigned long e = (p < 0) ? (unsigned long)(-(long long)p) : (unsigned long)p;
        if (p < 0) {                        /* invert base modulo M first */
            fmpq_poly_t d, s, t; fmpq_poly_init(d); fmpq_poly_init(s); fmpq_poly_init(t);
            fmpq_poly_xgcd(d, s, t, base, M);   /* d = s*base + t*M */
            if (fmpq_poly_degree(d) == 0 && !fmpq_poly_is_zero(d)) {
                fmpq_t inv; fmpq_init(inv);
                fmpq_poly_get_coeff_fmpq(inv, d, 0);   /* d is a nonzero constant */
                fmpq_inv(inv, inv);
                fmpq_poly_scalar_mul_fmpq(base, s, inv);   /* base = base^{-1} mod M */
                fmpq_clear(inv);
            } else ok = 0;                    /* non-invertible (should not happen) */
            fmpq_poly_clear(d); fmpq_poly_clear(s); fmpq_poly_clear(t);
        }
        if (ok) {
            fmpq_poly_set_ui(acc, 1);         /* acc = 1 */
            while (e > 0) {                   /* square-and-multiply mod M */
                if (e & 1UL) { fmpq_poly_mul(tmp, acc, base); fmpq_poly_rem(acc, tmp, M); }
                e >>= 1;
                if (e > 0) { fmpq_poly_mul(tmp, base, base); fmpq_poly_rem(base, tmp, M); }
            }
            result = poly_to_algnum_gen(g, acc, n);
        }
    }
    fmpq_poly_clear(M); fmpq_poly_clear(base); fmpq_poly_clear(acc); fmpq_poly_clear(tmp);
    return result;
}

/* a (AlgebraicNumber) combined with a rational scalar r: op 0 = a + r, 1 = a*r. */
static Expr* algnum_rational(const Expr* a, const Expr* r, int op) {
    if (!head_is(a, "AlgebraicNumber") || a->data.function.arg_count != 2) return NULL;
    fmpq_t rq; fmpq_init(rq);
    if (!fmpq_from_expr(r, rq)) { fmpq_clear(rq); return NULL; }
    const Expr* g = a->data.function.args[0];
    const fmpz_poly_struct* Mz; slong n;              /* degree from the cache; reuse g */
    if (!genfd_get(g, &Mz, &n, NULL)) { fmpq_clear(rq); return NULL; }
    fmpq_poly_t p; fmpq_poly_init(p);
    Expr* result = NULL;
    if (coeffs_to_fmpq_poly(a->data.function.args[1], p)) {
        if (op == 0) {                        /* add r to the constant coefficient */
            fmpq_t c0; fmpq_init(c0);
            fmpq_poly_get_coeff_fmpq(c0, p, 0);
            fmpq_add(c0, c0, rq);
            fmpq_poly_set_coeff_fmpq(p, 0, c0);
            fmpq_clear(c0);
        } else {                              /* scale every coefficient by r */
            fmpq_poly_scalar_mul_fmpq(p, p, rq);
        }
        result = poly_to_algnum_gen(g, p, n);
    }
    fmpq_poly_clear(p); fmpq_clear(rq);
    return result;
}

Expr* flint_qqbar_algnum_add_rational(const Expr* a, const Expr* r) {
    return algnum_rational(a, r, 0);
}

Expr* flint_qqbar_algnum_scale_rational(const Expr* a, const Expr* r) {
    return algnum_rational(a, r, 1);
}

/* Minimal polynomial of the algebraic generator `gen`, returned as a List of
 * its integer coefficients low-to-high, {c0, c1, ..., cn} (the primitive
 * integer minimal polynomial served from the generator cache), or NULL if gen
 * is not a recognised algebraic generator.  The field-coefficient Expand in
 * flint_bridge.c uses this to build the reduction modulus M(tau) without
 * duplicating the qqbar minpoly machinery; keeping it here reuses the same
 * genfd cache the AlgebraicNumber arithmetic already populates. */
Expr* flint_qqbar_gen_minpoly_coeffs(const Expr* gen) {
    if (!gen) return NULL;
    const fmpz_poly_struct* Mz; slong n;
    if (!genfd_get(gen, &Mz, &n, NULL)) return NULL;
    if (n < 1) return NULL;
    Expr** cs = malloc(sizeof(Expr*) * (size_t)(n + 1));
    if (!cs) return NULL;
    fmpz_t z; fmpz_init(z);
    for (slong i = 0; i <= n; i++) {
        fmpz_poly_get_coeff_fmpz(z, Mz, i);
        cs[i] = expr_from_fmpz(z);
    }
    fmpz_clear(z);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), cs, (size_t)(n + 1));
    free(cs);
    return list;
}

/* ================================================================== */
/*  Algebraic-entry matrices: exact linear algebra over the number    */
/*  field generated by the entries' atoms                             */
/* ================================================================== */
/*
 * The classical exact paths (the fraction-free elimination in linsolve.c /
 * inv.c, the cofactor expansion in det.c) treat an algebraic constant such
 * as Sqrt[3], I or Root[-1 - 4 #1 + #1^3 &, 1] as an opaque indeterminate.
 * Nothing reduces it by its minimal polynomial, so every elimination step
 * grows the expressions, and each exact division, GCD and zero test goes
 * back through Together / Cancel / PolynomialGCD over a multivariate ring
 * in those indeterminates: a 12x12 LinearSolve over Q(I Sqrt[3]) did not
 * finish in 40 s.
 *
 * Here the input is mapped ONCE into a number field K = Q(phi):
 *   1. collect the distinct algebraic atoms of every entry (radicals
 *      b^(1/q) of a rational b, I, Root objects, any other
 *      constant-algebraic leaf);
 *   2. take a primitive element theta = a_1 + c_2 a_2 + ... chosen by
 *      degree (the flint_qqbar_to_number_field_common construction), phi
 *      its algebraic-integer multiple, and express every atom in Q(phi)
 *      exactly (precision-escalated);
 *   3. evaluate each entry as an nf_elem (Plus, Times, integer Power,
 *      radical powers, inverses), so all arithmetic is native antic nf_elem
 *      with every value reduced mod the minimal polynomial, and a zero test
 *      is exact and O(1);
 *   4. render results in a monomial basis of the ATOMS themselves
 *      (a_1^e_1 ... a_k^e_k, e_i < deg a_i, chosen greedily by ascending
 *      total degree until they span K), so an answer over Q(I, Sqrt[3])
 *      reads a + b I + c Sqrt[3] + d I Sqrt[3] -- never a primitive element
 *      the user did not write.
 * Every entry point declines (NULL / 0) outside that domain: a free symbol,
 * an inexact number, no algebraic atom at all (the pure-rational fmpq paths
 * own that case), too many atoms, a field above ALGF_MAX_DEGREE, or a
 * division by zero while mapping an entry (the classical path then reports
 * it in its own terms).
 */

#define ALGF_MAX_ATOMS      8
#define ALGF_MAX_DEGREE     64
#define ALGF_MAX_CANDIDATES 4096
#define ALGF_MAX_POW        100000

typedef struct {
    Expr* atoms[ALGF_MAX_ATOMS];   /* normalised atom expressions (owned) */
    int na;
    nf_t nf;
    int nf_ready;
    slong D;                       /* [K : Q] */
    nf_elem_struct* aval;          /* the na atom values in K */
    Expr** basis;                  /* D monomial expressions in the atoms */
    fmpq_mat_t binv;               /* power-basis coords -> monomial coords */
    int binv_ready;
} AlgF;

static void algf_clear(AlgF* F) {
    for (int i = 0; i < F->na; i++) expr_free(F->atoms[i]);
    if (F->nf_ready) {
        if (F->aval) {
            for (int i = 0; i < F->na; i++) nf_elem_clear(F->aval + i, F->nf);
            free(F->aval);
        }
        nf_clear(F->nf);
    }
    if (F->basis) {
        for (slong i = 0; i < F->D; i++) if (F->basis[i]) expr_free(F->basis[i]);
        free(F->basis);
    }
    if (F->binv_ready) fmpq_mat_clear(F->binv);
    memset(F, 0, sizeof(*F));
}

static int algf_is_rat(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || head_is(e, "Rational"));
}

static Expr* algf_unit_i(void) {
    return expr_new_function(expr_new_symbol(SYM_Complex),
                             (Expr*[]){ expr_new_integer(0), expr_new_integer(1) }, 2);
}

/* b^(1/q) for a rational base b and an integer q > 1. */
static Expr* algf_radical(const Expr* b, const fmpz_t q) {
    fmpq_t r; fmpq_init(r);
    fmpz_one(fmpq_numref(r)); fmpz_set(fmpq_denref(r), q);
    Expr* ex = expr_from_fmpq(r);
    fmpq_clear(r);
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy((Expr*)b), ex }, 2);
}

/* Index of atom `a` (borrowed), or -1. */
static int algf_find(const AlgF* F, const Expr* a) {
    for (int i = 0; i < F->na; i++) if (expr_eq(F->atoms[i], a)) return i;
    return -1;
}

/* Register an atom (consumes `a`).  0 when the atom table is full. */
static int algf_add(AlgF* F, Expr* a) {
    if (algf_find(F, a) >= 0) { expr_free(a); return 1; }
    if (F->na >= ALGF_MAX_ATOMS) { expr_free(a); return 0; }
    F->atoms[F->na++] = a;
    return 1;
}

/* Pass 1: validate `e` and collect its atoms. */
static int algf_collect(AlgF* F, const Expr* e) {
    if (!e) return 0;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return 1;
        case EXPR_SYMBOL:
            /* Every symbol declines -- GoldenRatio included: the classical
             * paths keep the reference's forms for it
             * (RowReduce[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}] is
             * {{1, 1/GoldenRatio}, {0, 0}}), where a basis rendering would
             * read -1 + GoldenRatio. */
            return 0;
        case EXPR_FUNCTION: break;
        default: return 0;      /* Real, MPFR, String, ... */
    }
    size_t n = e->data.function.arg_count;
    if (head_is(e, "Rational")) return n == 2;
    if (head_is(e, "Complex") && n == 2) {
        if (!algf_is_rat(e->data.function.args[0]) || !algf_is_rat(e->data.function.args[1]))
            return 0;
        return algf_add(F, algf_unit_i());
    }
    if (head_is(e, "Plus") || head_is(e, "Times")) {
        for (size_t i = 0; i < n; i++)
            if (!algf_collect(F, e->data.function.args[i])) return 0;
        return 1;
    }
    if (head_is(e, "Power") && n == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* xe = e->data.function.args[1];
        if (xe->type == EXPR_INTEGER) return algf_collect(F, b);
        if (algf_is_rat(b) && head_is(xe, "Rational")) {
            fmpq_t r; fmpq_init(r);
            int ok = fmpq_from_expr(xe, r);
            if (ok) ok = algf_add(F, algf_radical(b, fmpq_denref(r)));
            fmpq_clear(r);
            return ok;
        }
    }
    if (!flint_qqbar_is_constant_algebraic(e)) return 0;
    return algf_add(F, expr_copy((Expr*)e));
}

/* Pass 2: build K, the atom values and the rendering basis. */
static int algf_build(AlgF* F) {
    if (F->na == 0) return 0;
    /* MATHILDA_NO_ALG_FIELD=1 keeps every caller on its classical path (for
     * differential testing, as MATHILDA_NO_NF_RREF does for the
     * AlgebraicNumber RREF). */
    if (getenv("MATHILDA_NO_ALG_FIELD")) return 0;
    int na = F->na;
    qqbar_ptr av = _qqbar_vec_init(na);
    int ok = 1;
    for (int i = 0; i < na && ok; i++) ok = to_qqbar(F->atoms[i], av + i);

    /* primitive element by degree (see flint_qqbar_to_number_field_common) */
    qqbar_t theta, cand, ca; qqbar_init(theta); qqbar_init(cand); qqbar_init(ca);
    if (ok) qqbar_set(theta, av + 0);
    for (int i = 1; i < na && ok; i++) {
        if (qqbar_is_rational(av + i)) continue;
        slong bestdeg = 0, bestc = 0;
        for (slong c = 1; c <= 16; c++) {
            qqbar_mul_si(ca, av + i, c);
            qqbar_add(cand, theta, ca);
            slong d = qqbar_degree(cand);
            if (d > bestdeg && d <= ALGF_MAX_DEGREE) { bestdeg = d; bestc = c; }
        }
        if (bestc == 0) { ok = 0; break; }
        qqbar_mul_si(ca, av + i, bestc);
        qqbar_add(theta, theta, ca);
    }
    qqbar_clear(cand); qqbar_clear(ca);

    qqbar_t phi; qqbar_init(phi);
    if (ok) {
        fmpz_t lc; fmpz_init(lc);
        algint_generator(theta, phi, lc);
        fmpz_clear(lc);
        F->D = qqbar_degree(phi);
        if (F->D < 2 || F->D > ALGF_MAX_DEGREE) ok = 0;
    }
    if (ok) {
        fmpq_poly_t M; fmpq_poly_init(M);
        fmpq_poly_set_fmpz_poly(M, QQBAR_POLY(phi));
        nf_init(F->nf, M);
        fmpq_poly_clear(M);
        F->nf_ready = 1;
        F->aval = malloc(sizeof(nf_elem_struct) * (size_t)na);
        for (int i = 0; i < na; i++) nf_elem_init(F->aval + i, F->nf);
        fmpq_poly_t f; fmpq_poly_init(f);
        for (int i = 0; i < na && ok; i++) {
            if (qqbar_express_in_field_esc_all(f, phi, av + i))
                nf_elem_set_fmpq_poly(F->aval + i, f, F->nf);
            else ok = 0;
        }
        fmpq_poly_clear(f);
    }

    /* Rendering basis: monomials in the atoms by ascending total degree, each
     * kept when independent of those already chosen (incremental echelon). */
    if (ok) {
        slong D = F->D;
        slong degs[ALGF_MAX_ATOMS];
        slong ncand = 1, maxtot = 0;
        for (int i = 0; i < na && ok; i++) {
            degs[i] = qqbar_degree(av + i);
            maxtot += degs[i] - 1;
            ncand *= degs[i];
            if (ncand > ALGF_MAX_CANDIDATES) ok = 0;
        }
        F->basis = calloc((size_t)D, sizeof(Expr*));
        fmpq_mat_t B; fmpq_mat_init(B, D, D);   /* chosen monomials, power-basis coords */
        fmpq_mat_t E; fmpq_mat_init(E, D, D);   /* their echelon reduction */
        slong* pc = malloc(sizeof(slong) * (size_t)D);
        fmpq* v = _fmpq_vec_init(D);
        fmpq_t t, u; fmpq_init(t); fmpq_init(u);
        fmpq_poly_t pp; fmpq_poly_init(pp);
        nf_elem_t mv, pw; nf_elem_init(mv, F->nf); nf_elem_init(pw, F->nf);
        slong got = 0;
        int ev[ALGF_MAX_ATOMS];
        for (slong tot = 0; ok && tot <= maxtot && got < D; tot++) {
            for (slong idx = 0; idx < ncand && got < D; idx++) {
                slong rem = idx, s = 0;
                for (int i = na - 1; i >= 0; i--) {
                    ev[i] = (int)(rem % degs[i]); rem /= degs[i]; s += ev[i];
                }
                if (s != tot) continue;
                nf_elem_one(mv, F->nf);
                for (int i = 0; i < na; i++) if (ev[i]) {
                    nf_elem_pow(pw, F->aval + i, (ulong)ev[i], F->nf);
                    nf_elem_mul(mv, mv, pw, F->nf);
                }
                nf_elem_get_fmpq_poly(pp, mv, F->nf);
                for (slong j = 0; j < D; j++) {
                    fmpq_poly_get_coeff_fmpq(v + j, pp, j);
                    fmpq_set(fmpq_mat_entry(B, got, j), v + j);
                }
                for (slong r = 0; r < got; r++) {
                    if (fmpq_is_zero(v + pc[r])) continue;
                    fmpq_set(t, v + pc[r]);
                    for (slong j = 0; j < D; j++) {
                        fmpq_mul(u, t, fmpq_mat_entry(E, r, j));
                        fmpq_sub(v + j, v + j, u);
                    }
                }
                slong p = -1;
                for (slong j = 0; j < D; j++) if (!fmpq_is_zero(v + j)) { p = j; break; }
                if (p < 0) continue;            /* dependent: try the next monomial */
                fmpq_inv(t, v + p);
                for (slong j = 0; j < D; j++) fmpq_mul(fmpq_mat_entry(E, got, j), v + j, t);
                pc[got] = p;
                Expr* fac[ALGF_MAX_ATOMS]; size_t nfac = 0;
                for (int i = 0; i < na; i++) if (ev[i]) {
                    fac[nfac++] = ev[i] == 1 ? expr_copy(F->atoms[i])
                        : expr_new_function(expr_new_symbol(SYM_Power),
                              (Expr*[]){ expr_copy(F->atoms[i]), expr_new_integer(ev[i]) }, 2);
                }
                Expr* mono = nfac == 0 ? expr_new_integer(1)
                           : nfac == 1 ? fac[0]
                           : expr_new_function(expr_new_symbol(SYM_Times), fac, nfac);
                F->basis[got] = eval_and_free(mono);
                got++;
            }
        }
        if (ok && got == D) {
            fmpq_mat_init(F->binv, D, D);
            F->binv_ready = 1;
            if (!fmpq_mat_inv(F->binv, B)) ok = 0;
        } else ok = 0;
        nf_elem_clear(mv, F->nf); nf_elem_clear(pw, F->nf);
        fmpq_poly_clear(pp); fmpq_clear(t); fmpq_clear(u);
        _fmpq_vec_clear(v, D);
        free(pc);
        fmpq_mat_clear(B); fmpq_mat_clear(E);
    }
    qqbar_clear(phi); qqbar_clear(theta);
    _qqbar_vec_clear(av, na);
    return ok;
}

/* Pass 3: the value in K of a constant expression whose atoms were
 * collected.  0 = outside the domain, or a division by zero. */
static int algf_value(AlgF* F, const Expr* e, nf_elem_t out) {
    fmpq_t q; fmpq_init(q);
    int isq = fmpq_from_expr(e, q);
    if (isq) nf_elem_set_fmpq(out, q, F->nf);
    fmpq_clear(q);
    if (isq) return 1;
    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        if (head_is(e, "Complex") && n == 2) {
            Expr* iu = algf_unit_i();
            int k = algf_find(F, iu);
            expr_free(iu);
            if (k < 0) return 0;
            fmpq_t re, im; fmpq_init(re); fmpq_init(im);
            int ok = fmpq_from_expr(e->data.function.args[0], re) &&
                     fmpq_from_expr(e->data.function.args[1], im);
            if (ok) {
                nf_elem_scalar_mul_fmpq(out, F->aval + k, im, F->nf);
                nf_elem_add_fmpq(out, out, re, F->nf);
            }
            fmpq_clear(re); fmpq_clear(im);
            return ok;
        }
        if ((head_is(e, "Plus") || head_is(e, "Times")) && n >= 1) {
            int plus = head_is(e, "Plus");
            nf_elem_t t; nf_elem_init(t, F->nf);
            int ok = algf_value(F, e->data.function.args[0], out);
            for (size_t i = 1; i < n && ok; i++) {
                ok = algf_value(F, e->data.function.args[i], t);
                if (ok) {
                    if (plus) nf_elem_add(out, out, t, F->nf);
                    else nf_elem_mul(out, out, t, F->nf);
                }
            }
            nf_elem_clear(t, F->nf);
            return ok;
        }
        if (head_is(e, "Power") && n == 2) {
            const Expr* b = e->data.function.args[0];
            const Expr* xe = e->data.function.args[1];
            slong p = 0;
            int ok = 0, generic = 0;
            nf_elem_t base, inv; nf_elem_init(base, F->nf); nf_elem_init(inv, F->nf);
            if (xe->type == EXPR_INTEGER) {
                p = (slong)xe->data.integer;
                ok = algf_value(F, b, base);
            } else if (algf_is_rat(b) && head_is(xe, "Rational")) {
                /* b^(p/q) = (b^(1/q))^p, the principal branch (as in to_qqbar) */
                fmpq_t r; fmpq_init(r);
                if (fmpq_from_expr(xe, r) && fmpz_fits_si(fmpq_numref(r))) {
                    Expr* a = algf_radical(b, fmpq_denref(r));
                    int k = algf_find(F, a);
                    expr_free(a);
                    if (k >= 0) {
                        nf_elem_set(base, F->aval + k, F->nf);
                        p = fmpz_get_si(fmpq_numref(r));
                        ok = 1;
                    }
                }
                fmpq_clear(r);
            } else generic = 1;
            if (ok && (p > ALGF_MAX_POW || p < -ALGF_MAX_POW)) ok = 0;
            if (ok && p < 0) {
                if (nf_elem_is_zero(base, F->nf)) ok = 0;
                else { nf_elem_inv(inv, base, F->nf); nf_elem_swap(inv, base, F->nf); p = -p; }
            }
            if (ok) nf_elem_pow(out, base, (ulong)p, F->nf);
            nf_elem_clear(base, F->nf); nf_elem_clear(inv, F->nf);
            if (!generic) return ok;
        }
    }
    int k = algf_find(F, e);
    if (k < 0) return 0;
    nf_elem_set(out, F->aval + k, F->nf);
    return 1;
}

/* Monomial-basis coordinates of v: c = u * B^-1 (u its power-basis coords). */
static void algf_coords(AlgF* F, const nf_elem_t v, fmpq* c) {
    slong D = F->D;
    fmpq_poly_t pp; fmpq_poly_init(pp);
    nf_elem_get_fmpq_poly(pp, v, F->nf);
    fmpq* u = _fmpq_vec_init(D);
    for (slong j = 0; j < D; j++) fmpq_poly_get_coeff_fmpq(u + j, pp, j);
    fmpq_t t; fmpq_init(t);
    for (slong j = 0; j < D; j++) {
        fmpq_zero(c + j);
        for (slong i = 0; i < D; i++) {
            if (fmpq_is_zero(u + i)) continue;
            fmpq_mul(t, u + i, fmpq_mat_entry(F->binv, i, j));
            fmpq_add(c + j, c + j, t);
        }
    }
    fmpq_clear(t);
    _fmpq_vec_clear(u, D);
    fmpq_poly_clear(pp);
}

/* Append the terms c_j * basis_j (* xk) of v to `terms` (grown as needed). */
static void algf_push_terms(AlgF* F, const nf_elem_t v, const Expr* xk,
                            Expr*** terms, size_t* nt, size_t* cap) {
    slong D = F->D;
    fmpq* c = _fmpq_vec_init(D);
    algf_coords(F, v, c);
    for (slong j = 0; j < D; j++) {
        if (fmpq_is_zero(c + j)) continue;
        if (*nt == *cap) {
            *cap = *cap ? 2 * *cap : 16;
            *terms = realloc(*terms, sizeof(Expr*) * *cap);
        }
        Expr* f[3]; size_t nfac = 0;
        f[nfac++] = expr_from_fmpq(c + j);
        f[nfac++] = expr_copy(F->basis[j]);
        if (xk) f[nfac++] = expr_copy((Expr*)xk);
        (*terms)[(*nt)++] = expr_new_function(expr_new_symbol(SYM_Times), f, nfac);
    }
    _fmpq_vec_clear(c, D);
}

static Expr* algf_sum(Expr** terms, size_t nt) {
    if (nt == 0) return expr_new_integer(0);
    Expr* s = nt == 1 ? terms[0] : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
    return eval_and_free(s);
}

static Expr* algf_render(AlgF* F, const nf_elem_t v) {
    if (nf_elem_is_zero(v, F->nf)) return expr_new_integer(0);
    Expr** terms = NULL; size_t nt = 0, cap = 0;
    algf_push_terms(F, v, NULL, &terms, &nt, &cap);
    Expr* r = algf_sum(terms, nt);
    free(terms);
    return r;
}

static void algf_vec_free(AlgF* F, nf_elem_struct* v, long n) {
    for (long i = 0; i < n; i++) nf_elem_clear(v + i, F->nf);
    free(v);
}

/* Map n constant entries into K.  On success *outv holds n initialised
 * nf_elems (release with algf_vec_free). */
static int algf_load(AlgF* F, Expr* const* flat, long n, nf_elem_struct** outv) {
    for (long i = 0; i < n; i++) if (!algf_collect(F, flat[i])) return 0;
    if (!algf_build(F)) return 0;
    nf_elem_struct* v = malloc(sizeof(nf_elem_struct) * (size_t)(n > 0 ? n : 1));
    for (long i = 0; i < n; i++) nf_elem_init(v + i, F->nf);
    for (long i = 0; i < n; i++)
        if (!algf_value(F, flat[i], v + i)) { algf_vec_free(F, v, n); return 0; }
    *outv = v;
    return 1;
}

/* Gauss-Jordan in place on an r x cols row-major matrix, pivoting in the
 * first `pcols` columns only; pivots are scaled to 1.  pivcol[i] (if given)
 * receives the pivot column of row i.  *sign flips on each row swap and *det
 * accumulates the product of the pivots before scaling.  Returns the rank. */
static int algf_gauss_jordan(AlgF* F, nf_elem_struct* A, int r, int cols, int pcols,
                             int* pivcol, int* sign, nf_elem_struct* det) {
    nf_elem_t inv, fac, t;
    nf_elem_init(inv, F->nf); nf_elem_init(fac, F->nf); nf_elem_init(t, F->nf);
    int rank = 0;
    for (int col = 0; col < pcols && rank < r; col++) {
        int prow = -1;
        for (int i = rank; i < r; i++)
            if (!nf_elem_is_zero(A + (long)i * cols + col, F->nf)) { prow = i; break; }
        if (prow < 0) continue;
        if (prow != rank) {
            for (int j = 0; j < cols; j++)
                nf_elem_swap(A + (long)rank * cols + j, A + (long)prow * cols + j, F->nf);
            if (sign) *sign = -*sign;
        }
        nf_elem_struct* prw = A + (long)rank * cols;
        if (det) nf_elem_mul(det, det, prw + col, F->nf);
        nf_elem_inv(inv, prw + col, F->nf);
        for (int j = col; j < cols; j++)
            if (!nf_elem_is_zero(prw + j, F->nf)) nf_elem_mul(prw + j, prw + j, inv, F->nf);
        for (int i = 0; i < r; i++) {
            if (i == rank) continue;
            nf_elem_struct* row = A + (long)i * cols;
            if (nf_elem_is_zero(row + col, F->nf)) continue;
            nf_elem_set(fac, row + col, F->nf);
            for (int j = col; j < cols; j++) {
                if (nf_elem_is_zero(prw + j, F->nf)) continue;
                nf_elem_mul(t, fac, prw + j, F->nf);
                nf_elem_sub(row + j, row + j, t, F->nf);
            }
        }
        if (pivcol) pivcol[rank] = col;
        rank++;
    }
    nf_elem_clear(inv, F->nf); nf_elem_clear(fac, F->nf); nf_elem_clear(t, F->nf);
    return rank;
}

static Expr* algf_render_matrix(AlgF* F, nf_elem_struct* A, int r, int c, int stride, int off) {
    Expr** rows = malloc(sizeof(Expr*) * (size_t)r);
    for (int i = 0; i < r; i++) {
        Expr** el = malloc(sizeof(Expr*) * (size_t)c);
        for (int j = 0; j < c; j++) el[j] = algf_render(F, A + (long)i * stride + off + j);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)c);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)r);
    free(rows);
    return out;
}

int flint_qqbar_alg_mat_solve(Expr* const* A, int r, int c, Expr* const* B, int k,
                              Expr*** out) {
    if (!A || !B || r <= 0 || c <= 0 || k <= 0 || !out) return 0;
    int cols = c + k;
    long total = (long)r * cols;
    Expr** flat = malloc(sizeof(Expr*) * (size_t)total);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) flat[(long)i * cols + j] = A[(long)i * c + j];
        for (int j = 0; j < k; j++) flat[(long)i * cols + c + j] = B[(long)i * k + j];
    }
    AlgF F; memset(&F, 0, sizeof(F));
    nf_elem_struct* M = NULL;
    int ok = algf_load(&F, flat, total, &M);
    free(flat);
    if (!ok) { algf_clear(&F); return 0; }
    int* piv = malloc(sizeof(int) * (size_t)r);
    int rank = algf_gauss_jordan(&F, M, r, cols, c, piv, NULL, NULL);
    int status = 1;
    for (int i = rank; i < r && status == 1; i++)
        for (int j = 0; j < k; j++)
            if (!nf_elem_is_zero(M + (long)i * cols + c + j, F.nf)) { status = -1; break; }
    if (status == 1) {
        long ns = (long)c * k;
        Expr** sol = calloc((size_t)ns, sizeof(Expr*));
        for (int i = 0; i < rank; i++)
            for (int j = 0; j < k; j++)
                sol[(long)piv[i] * k + j] = algf_render(&F, M + (long)i * cols + c + j);
        for (long i = 0; i < ns; i++) if (!sol[i]) sol[i] = expr_new_integer(0);
        *out = sol;
    }
    free(piv);
    algf_vec_free(&F, M, total);
    algf_clear(&F);
    return status;
}

Expr* flint_qqbar_alg_mat_rref(Expr* const* flat, int r, int c) {
    if (!flat || r <= 0 || c <= 0) return NULL;
    AlgF F; memset(&F, 0, sizeof(F));
    nf_elem_struct* M = NULL;
    long total = (long)r * c;
    if (!algf_load(&F, flat, total, &M)) { algf_clear(&F); return NULL; }
    algf_gauss_jordan(&F, M, r, c, c, NULL, NULL, NULL);
    Expr* out = algf_render_matrix(&F, M, r, c, c, 0);
    algf_vec_free(&F, M, total);
    algf_clear(&F);
    return out;
}

Expr* flint_qqbar_alg_mat_det(Expr* const* flat, int n) {
    if (!flat || n <= 0) return NULL;
    AlgF F; memset(&F, 0, sizeof(F));
    nf_elem_struct* M = NULL;
    long total = (long)n * n;
    if (!algf_load(&F, flat, total, &M)) { algf_clear(&F); return NULL; }
    nf_elem_t det; nf_elem_init(det, F.nf); nf_elem_one(det, F.nf);
    int sign = 1;
    int rank = algf_gauss_jordan(&F, M, n, n, n, NULL, &sign, det);
    if (rank < n) nf_elem_zero(det, F.nf);
    else if (sign < 0) nf_elem_neg(det, det, F.nf);
    Expr* out = algf_render(&F, det);
    nf_elem_clear(det, F.nf);
    algf_vec_free(&F, M, total);
    algf_clear(&F);
    return out;
}

Expr* flint_qqbar_alg_mat_inverse(Expr* const* flat, int n, int* singular) {
    if (singular) *singular = 0;
    if (!flat || n <= 0) return NULL;
    AlgF F; memset(&F, 0, sizeof(F));
    nf_elem_struct* M0 = NULL;
    long total = (long)n * n;
    if (!algf_load(&F, flat, total, &M0)) { algf_clear(&F); return NULL; }
    int cols = 2 * n;
    long tot2 = (long)n * cols;
    nf_elem_struct* M = malloc(sizeof(nf_elem_struct) * (size_t)tot2);
    for (long i = 0; i < tot2; i++) { nf_elem_init(M + i, F.nf); nf_elem_zero(M + i, F.nf); }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) nf_elem_set(M + (long)i * cols + j, M0 + (long)i * n + j, F.nf);
        nf_elem_one(M + (long)i * cols + n + i, F.nf);
    }
    algf_vec_free(&F, M0, total);
    int rank = algf_gauss_jordan(&F, M, n, cols, n, NULL, NULL, NULL);
    Expr* out = NULL;
    if (rank < n) { if (singular) *singular = 1; }
    else out = algf_render_matrix(&F, M, n, n, cols, n);
    algf_vec_free(&F, M, tot2);
    algf_clear(&F);
    return out;
}

/* ---- an exact zero test over the field of a fixed set of expressions ---- */

struct FlintAlgField { AlgF F; };

FlintAlgField* flint_qqbar_alg_field_new(Expr* const* exprs, size_t n) {
    if (!exprs || n == 0) return NULL;
    FlintAlgField* A = calloc(1, sizeof(*A));
    int ok = 1;
    for (size_t i = 0; i < n && ok; i++) ok = algf_collect(&A->F, exprs[i]);
    if (ok) ok = algf_build(&A->F);
    if (!ok) { algf_clear(&A->F); free(A); return NULL; }
    return A;
}

int flint_qqbar_alg_field_is_zero(FlintAlgField* A, const Expr* e) {
    if (!A || !e) return -1;
    nf_elem_t v; nf_elem_init(v, A->F.nf);
    int r = algf_value(&A->F, e, v) ? (nf_elem_is_zero(v, A->F.nf) ? 1 : 0) : -1;
    nf_elem_clear(v, A->F.nf);
    return r;
}

void flint_qqbar_alg_field_free(FlintAlgField* A) {
    if (!A) return;
    algf_clear(&A->F);
    free(A);
}

#else /* !USE_FLINT */

int   flint_qqbar_is_constant_algebraic(const Expr* e) { (void)e; return 0; }
Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod m) { (void)e; (void)m; return NULL; }
Expr* flint_qqbar_reduce_coeffs(const Expr* e, QQBarMethod m) { (void)e; (void)m; return NULL; }
int   flint_qqbar_equal(const Expr* a, const Expr* b) { (void)a; (void)b; return -1; }
int   flint_qqbar_compare(const Expr* a, const Expr* b) { (void)a; (void)b; return -2; }
int   flint_qqbar_is_real(const Expr* e) { (void)e; return -1; }
int   flint_qqbar_part_sign(const Expr* e, int imag) { (void)e; (void)imag; return -2; }
Expr* flint_qqbar_algebraic_number(const Expr* g, const Expr* c) { (void)g; (void)c; return NULL; }
Expr* flint_qqbar_to_number_field(const Expr* a, const Expr* t) { (void)a; (void)t; return NULL; }
Expr* flint_qqbar_to_number_field_self(const Expr* x) { (void)x; return NULL; }
Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n, int s) { (void)as; (void)n; (void)s; return NULL; }
Expr* flint_qqbar_to_number_field_common_fn(Expr* (*get)(const void*, size_t, bool*),
                                            const void* ctx, size_t n, int s, bool o) { (void)get; (void)ctx; (void)n; (void)s; (void)o; return NULL; }
Expr* flint_qqbar_integral_basis(const Expr* a) { (void)a; return NULL; }
int   flint_qqbar_algebraic_integer_q(const Expr* x) { (void)x; return -1; }
int   flint_qqbar_algebraic_number_denominator(const Expr* x, Expr** out) { (void)x; (void)out; return -1; }
int   flint_qqbar_algebraic_number_norm(const Expr* a, const Expr* t, Expr** out) { (void)a; (void)t; (void)out; return -1; }
int   flint_qqbar_algebraic_number_trace(const Expr* a, const Expr* t, Expr** out) { (void)a; (void)t; (void)out; return -1; }
Expr* flint_qqbar_algnum_add(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }
Expr* flint_qqbar_algnum_mul(const Expr* a, const Expr* b) { (void)a; (void)b; return NULL; }
Expr* flint_qqbar_algnum_pow(const Expr* a, long p) { (void)a; (void)p; return NULL; }
Expr* flint_qqbar_algnum_add_rational(const Expr* a, const Expr* r) { (void)a; (void)r; return NULL; }
Expr* flint_qqbar_algnum_scale_rational(const Expr* a, const Expr* r) { (void)a; (void)r; return NULL; }
Expr* flint_qqbar_gen_minpoly_coeffs(const Expr* gen) { (void)gen; return NULL; }
Expr* flint_qqbar_nf_mat_rref(Expr** flat, int r, int c) {
    (void)flat; (void)r; (void)c; return NULL;
}

int   flint_qqbar_alg_mat_solve(Expr* const* A, int r, int c, Expr* const* B, int k, Expr*** out) {
    (void)A; (void)r; (void)c; (void)B; (void)k; (void)out; return 0;
}
Expr* flint_qqbar_alg_mat_rref(Expr* const* f, int r, int c) { (void)f; (void)r; (void)c; return NULL; }
Expr* flint_qqbar_alg_mat_det(Expr* const* f, int n) { (void)f; (void)n; return NULL; }
Expr* flint_qqbar_alg_mat_inverse(Expr* const* f, int n, int* s) { (void)f; (void)n; if (s) *s = 0; return NULL; }
struct FlintAlgField { int unused; };
FlintAlgField* flint_qqbar_alg_field_new(Expr* const* e, size_t n) { (void)e; (void)n; return NULL; }
int   flint_qqbar_alg_field_is_zero(FlintAlgField* A, const Expr* e) { (void)A; (void)e; return -1; }
void  flint_qqbar_alg_field_free(FlintAlgField* A) { (void)A; }

#endif /* USE_FLINT */
