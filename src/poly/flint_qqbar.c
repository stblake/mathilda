/* flint_qqbar.c — exact algebraic-number canonicalisation for RootReduce.
 *
 * See flint_qqbar.h for the contract. The engine is FLINT's `qqbar` (exact
 * real/complex algebraic numbers). A constant algebraic-number Expr is folded
 * bottom-up into a single qqbar, whose minimal polynomial, degree, rationality,
 * quadratic-radical form and root index then produce the WL-faithful canonical
 * Expr: a rational, a quadratic radical, or Root[Function[minpoly&], k].
 */

#include "flint_qqbar.h"

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
#include <flint/qqbar.h>

/* The maximal-order (Round 2 / Pohst-Zassenhaus) engine and its integral-basis
 * accessors, used by flint_qqbar_integral_basis below. numberfield_internal.h
 * guards its own FLINT includes and exposes nf_ok_basis / nf_ok_denom. */
#include "numberfield.h"
#include "numberfield_internal.h"

/* Intermediate qqbar degree above which we give up and fall back (identity /
 * the parametric engine). WL's degree-21 examples stay comfortably under. */
#define QQBAR_DEGREE_CAP 120
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
    /* both non-real: Re ascending */
    int c = qqbar_cmp_re(a, b);
    if (c != 0) return c;
    /* |Im| ascending */
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

/* 1-based WL index of x among the roots of its minimal polynomial. */
static slong wl_root_index(const qqbar_t x) {
    slong d = qqbar_degree(x);
    if (d <= 0) return 1;
    qqbar_ptr roots = _qqbar_vec_init(d);
    qqbar_roots_fmpz_poly(roots, QQBAR_POLY(x), QQBAR_ROOTS_IRREDUCIBLE);
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
            if (qqbar_degree(acc) > QQBAR_DEGREE_CAP) { ok = 0; break; }
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
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) ok = 0;
        qqbar_clear(base);
        return ok;
    }
    if (head_is(e, "Sqrt") && n == 1) {
        qqbar_t base; qqbar_init(base);
        int ok = to_qqbar(e->data.function.args[0], base);
        if (ok) qqbar_sqrt(out, base);
        if (ok && qqbar_degree(out) > QQBAR_DEGREE_CAP) ok = 0;
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
            if (qqbar_degree(acc) > QQBAR_DEGREE_CAP) { ok = 0; break; }
        }
        fmpq_clear(c);
        if (ok) qqbar_set(out, acc);
        qqbar_clear(acc); qqbar_clear(g);
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
    /* degree >= 3: Root[Function[minpoly&], k] (held). */
    Expr* fn = minpoly_slot_expr(x);
    slong k = wl_root_index(x);
    return expr_new_function(expr_new_symbol(SYM_Root),
               (Expr*[]){ fn, expr_new_integer((int64_t)k) }, 2);
}

/* ------------------------------------------------------------------ */
/*  Method -> "NumberField": common-number-field re-expression         */
/* ------------------------------------------------------------------ */

/* True if x is expressible in Q(alpha). */
static int in_field(const qqbar_t x, const qqbar_t alpha) {
    fmpq_poly_t f; fmpq_poly_init(f);
    int ok = qqbar_express_in_field(f, alpha, x, 100000, 0, 64);
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

/* Same-field binary op on two AlgebraicNumbers with expr_eq generators.
 * op 0 = add, 1 = mul. Returns a fresh combined value, or NULL to decline. */
static Expr* algnum_binop(const Expr* a, const Expr* b, int op) {
    if (!head_is(a, "AlgebraicNumber") || !head_is(b, "AlgebraicNumber")) return NULL;
    if (a->data.function.arg_count != 2 || b->data.function.arg_count != 2) return NULL;
    const Expr* g = a->data.function.args[0];
    if (!expr_eq(g, b->data.function.args[0])) return NULL;

    qqbar_t phi; qqbar_init(phi);
    if (!to_qqbar(g, phi)) { qqbar_clear(phi); return NULL; }
    slong n = qqbar_degree(phi);

    fmpq_poly_t M, pa, pb, r;
    fmpq_poly_init(M); fmpq_poly_init(pa); fmpq_poly_init(pb); fmpq_poly_init(r);
    fmpq_poly_set_fmpz_poly(M, QQBAR_POLY(phi));
    Expr* result = NULL;
    if (coeffs_to_fmpq_poly(a->data.function.args[1], pa) &&
        coeffs_to_fmpq_poly(b->data.function.args[1], pb)) {
        if (op == 0) { fmpq_poly_add(r, pa, pb); }
        else         { fmpq_poly_mul(r, pa, pb); fmpq_poly_rem(r, r, M); }
        result = poly_to_algnum(phi, r, n);
    }
    fmpq_poly_clear(M); fmpq_poly_clear(pa); fmpq_poly_clear(pb); fmpq_poly_clear(r);
    qqbar_clear(phi);
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

int flint_qqbar_equal(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -1;
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -1;
    if (to_qqbar(a, qa) && to_qqbar(b, qb)) r = qqbar_equal(qa, qb) ? 1 : 0;
    qqbar_clear(qa); qqbar_clear(qb);
    return r;
}

int flint_qqbar_compare(const Expr* a, const Expr* b) {
    if (!flint_qqbar_is_constant_algebraic(a) || !flint_qqbar_is_constant_algebraic(b))
        return -2;
    qqbar_t qa, qb; qqbar_init(qa); qqbar_init(qb);
    int r = -2;
    if (to_qqbar(a, qa) && to_qqbar(b, qb) && qqbar_is_real(qa) && qqbar_is_real(qb))
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

Expr* flint_qqbar_algebraic_number(const Expr* gen, const Expr* coeffs) {
    if (!gen || !coeffs || !head_is(coeffs, "List")) return NULL;

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
    if (qqbar_express_in_field(f, phi, av, 100000, 0, 64))
        result = poly_to_algnum(phi, f, n);   /* NULL stays NULL when a not in Q(theta) */
    fmpq_poly_clear(f);
    qqbar_clear(phi); qqbar_clear(av); qqbar_clear(tv);
    return result;
}

Expr* flint_qqbar_to_number_field_self(const Expr* x) {
    return flint_qqbar_to_number_field(x, x);
}

Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n,
                                         int smallest) {
    (void)smallest;                            /* Automatic and All coincide here */
    if (!as || n == 0) return NULL;

    qqbar_ptr vals = _qqbar_vec_init((slong)n);
    int ok = 1;
    for (size_t i = 0; i < n && ok; i++)
        if (!to_qqbar(as[i], vals + i)) ok = 0;

    qqbar_t alpha, cand, ca; qqbar_init(alpha); qqbar_init(cand); qqbar_init(ca);
    if (ok) qqbar_set(alpha, vals + 0);
    for (size_t i = 1; i < n && ok; i++) {
        int found = 0;
        for (slong c = 1; c <= 16 && !found; c++) {   /* primitive element a + c*b */
            qqbar_mul_si(ca, vals + i, c);
            qqbar_add(cand, alpha, ca);
            if (qqbar_degree(cand) <= QQBAR_DEGREE_CAP &&
                in_field(vals + i, cand) && in_field(alpha, cand)) {
                qqbar_set(alpha, cand); found = 1;
            }
        }
        if (!found) ok = 0;
    }
    qqbar_clear(cand); qqbar_clear(ca);

    Expr* result = NULL;
    if (ok) {
        fmpz_t lc; fmpz_init(lc);
        qqbar_t phi; qqbar_init(phi);
        algint_generator(alpha, phi, lc);
        fmpz_clear(lc);
        slong nn = qqbar_degree(phi);

        Expr** items = malloc(sizeof(Expr*) * n);
        int good = 1;
        for (size_t i = 0; i < n && good; i++) {
            fmpq_poly_t f; fmpq_poly_init(f);
            if (qqbar_express_in_field(f, phi, vals + i, 100000, 0, 64))
                items[i] = poly_to_algnum(phi, f, nn);
            else { items[i] = NULL; good = 0; }
            fmpq_poly_clear(f);
        }
        if (good) {
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
    for (slong i = 0; i < n; i++) {
        fmpq_poly_t row; fmpq_poly_init(row);
        for (slong j = 0; j < n; j++) {
            fmpq_set_fmpz_frac(coef, fmpz_mat_entry(H, n - 1 - i, n - 1 - j), Dz);
            fmpq_poly_set_coeff_fmpq(row, j, coef);
        }
        items[i] = poly_to_algnum(phi, row, n);     /* rational row -> Integer/Rational */
        fmpq_poly_clear(row);
    }
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
int flint_qqbar_algebraic_integer_q(const Expr* x) {
    if (!x) return 0;
    qqbar_t v; qqbar_init(v);
    if (!to_qqbar(x, v)) { qqbar_clear(v); return 0; }
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
            if (qqbar_express_in_field(f, phi, av, 100000, 0, 64)) {
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
            if (qqbar_express_in_field(f, phi, av, 100000, 0, 64)) {
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
    qqbar_t phi; qqbar_init(phi);
    if (!to_qqbar(g, phi)) { qqbar_clear(phi); return NULL; }
    slong n = qqbar_degree(phi);

    fmpq_poly_t M, base, acc, tmp;
    fmpq_poly_init(M); fmpq_poly_init(base); fmpq_poly_init(acc); fmpq_poly_init(tmp);
    fmpq_poly_set_fmpz_poly(M, QQBAR_POLY(phi));

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
            result = poly_to_algnum(phi, acc, n);
        }
    }
    fmpq_poly_clear(M); fmpq_poly_clear(base); fmpq_poly_clear(acc); fmpq_poly_clear(tmp);
    qqbar_clear(phi);
    return result;
}

/* a (AlgebraicNumber) combined with a rational scalar r: op 0 = a + r, 1 = a*r. */
static Expr* algnum_rational(const Expr* a, const Expr* r, int op) {
    if (!head_is(a, "AlgebraicNumber") || a->data.function.arg_count != 2) return NULL;
    fmpq_t rq; fmpq_init(rq);
    if (!fmpq_from_expr(r, rq)) { fmpq_clear(rq); return NULL; }
    const Expr* g = a->data.function.args[0];
    qqbar_t phi; qqbar_init(phi);
    if (!to_qqbar(g, phi)) { qqbar_clear(phi); fmpq_clear(rq); return NULL; }
    slong n = qqbar_degree(phi);
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
        result = poly_to_algnum(phi, p, n);
    }
    fmpq_poly_clear(p); qqbar_clear(phi); fmpq_clear(rq);
    return result;
}

Expr* flint_qqbar_algnum_add_rational(const Expr* a, const Expr* r) {
    return algnum_rational(a, r, 0);
}

Expr* flint_qqbar_algnum_scale_rational(const Expr* a, const Expr* r) {
    return algnum_rational(a, r, 1);
}

#else /* !USE_FLINT */

int   flint_qqbar_is_constant_algebraic(const Expr* e) { (void)e; return 0; }
Expr* flint_qqbar_canonical(const Expr* e, QQBarMethod m) { (void)e; (void)m; return NULL; }
Expr* flint_qqbar_reduce_coeffs(const Expr* e, QQBarMethod m) { (void)e; (void)m; return NULL; }
int   flint_qqbar_equal(const Expr* a, const Expr* b) { (void)a; (void)b; return -1; }
int   flint_qqbar_compare(const Expr* a, const Expr* b) { (void)a; (void)b; return -2; }
int   flint_qqbar_is_real(const Expr* e) { (void)e; return -1; }
Expr* flint_qqbar_algebraic_number(const Expr* g, const Expr* c) { (void)g; (void)c; return NULL; }
Expr* flint_qqbar_to_number_field(const Expr* a, const Expr* t) { (void)a; (void)t; return NULL; }
Expr* flint_qqbar_to_number_field_self(const Expr* x) { (void)x; return NULL; }
Expr* flint_qqbar_to_number_field_common(const Expr* const* as, size_t n, int s) { (void)as; (void)n; (void)s; return NULL; }
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

#endif /* USE_FLINT */
