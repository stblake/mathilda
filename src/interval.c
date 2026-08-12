/* Interval arithmetic — see interval.h for the design overview.
 *
 * Representation: Interval[{lo,hi}, ...] as an EXPR_FUNCTION with head
 * SYM_Interval, canonical form = pairs ordered lo<=hi, sorted by lo, and
 * pairwise non-overlapping and non-touching (closed-interval union).
 *
 * Exactness: endpoint operations are performed by the ordinary evaluator, so
 * exact operands yield exact/symbolic endpoints (Integer+Integer -> Integer,
 * Sin[2] stays Sin[2]); an inexact (Real/MPFR) result is then nudged one ULP
 * outward (iv_widen) so the enclosure stays rigorous. -Infinity is the standard
 * Times[-1, Infinity] form.
 */

#include "interval.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "numeric.h"
#include "sym_names.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* The minimum of Gamma (and LogGamma) on (0, inf), at x = 1.4616321449683623;
 * Gamma is decreasing below it and increasing above it. */
#define IV_GAMMA_MIN 1.4616321449683623

/* Outward-rounding direction for inexact endpoints. */
typedef enum { RND_DOWN = 0, RND_UP = 1 } IvRound;

/* ---------------------------------------------------------------------------
 * Small endpoint helpers
 * ------------------------------------------------------------------------- */

static Expr* iv_pos_inf(void) { return expr_new_symbol(SYM_Infinity); }

static Expr* iv_neg_inf(void) {
    Expr* a[2] = { expr_new_integer(-1), expr_new_symbol(SYM_Infinity) };
    return expr_new_function(expr_new_symbol(SYM_Times), a, 2);
}

static bool iv_is_pos_inf(Expr* e) { return is_infinity_sym(e); }
static bool iv_is_neg_inf(Expr* e) { return is_neg_infinity_form(e); }
static bool iv_is_inf(Expr* e) { return iv_is_pos_inf(e) || iv_is_neg_inf(e); }

static bool iv_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

static bool iv_is_exact_zero(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) == 0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr);
#endif
    { int64_t n, d; if (is_rational(e, &n, &d)) return n == 0; }
    return false;
}

/* Copy an infinite endpoint, preserving its sign. */
static Expr* iv_copy_inf(Expr* e) {
    return iv_is_neg_inf(e) ? iv_neg_inf() : iv_pos_inf();
}

/* Numeric approximation of an endpoint as a double, using N[] for symbolic
 * numerics (Pi, Sin[2], ...). Infinities map to +/-HUGE_VAL. */
static bool iv_to_double(Expr* e, double* out) {
    if (iv_is_pos_inf(e)) { *out = HUGE_VAL; return true; }
    if (iv_is_neg_inf(e)) { *out = -HUGE_VAL; return true; }
    Expr* call = expr_new_function(expr_new_symbol("N"),
                                   (Expr*[]){ expr_copy(e) }, 1);
    Expr* r = eval_and_free(call);
    bool ok = false;
    if (r) {
        if (r->type == EXPR_INTEGER) { *out = (double)r->data.integer; ok = true; }
        else if (r->type == EXPR_REAL) { *out = r->data.real; ok = true; }
        else if (r->type == EXPR_BIGINT) { *out = mpz_get_d(r->data.bigint); ok = true; }
#ifdef USE_MPFR
        else if (r->type == EXPR_MPFR) { *out = mpfr_get_d(r->data.mpfr, MPFR_RNDN); ok = true; }
#endif
        else { int64_t n, d; if (is_rational(r, &n, &d)) { *out = (double)n / (double)d; ok = true; } }
    }
    expr_free(r);
    return ok;
}

#ifdef USE_MPFR
/* Numericalize an endpoint into `out` (already mpfr_init2'd). Handles numeric
 * endpoints directly and symbolic-numerics (Pi, Sin[2], ...) via N[e, 40].
 * Returns false for genuinely non-numeric (a free variable). */
static bool iv_endpoint_to_mpfr(Expr* e, mpfr_t out) {
    mpfr_t im; mpfr_init2(im, mpfr_get_prec(out));
    bool ok = get_approx_mpfr(e, out, im, NULL) && mpfr_zero_p(im) && mpfr_number_p(out);
    if (!ok) {
        Expr* ne = eval_and_free(expr_new_function(expr_new_symbol("N"),
                       (Expr*[]){ expr_copy(e), expr_new_integer(40) }, 2));
        if (ne) {
            ok = get_approx_mpfr(ne, out, im, NULL) && mpfr_zero_p(im) && mpfr_number_p(out);
            expr_free(ne);
        }
    }
    mpfr_clear(im);
    return ok;
}
#endif

/* Compare two interval endpoints. Returns the sign of (a - b): -1, 0, or +1,
 * and sets *decidable to whether the comparison resolved.
 *
 * This comparison must be EXACT (bit-for-bit on the numeric values), NOT the
 * ~2^-46 relative-tolerance comparison the Less/Greater/Equal heads use for
 * machine reals: the outward-rounded min/max in the kernels rely on it to pick
 * the truly outermost of several endpoints that differ by only a ULP or two, so
 * a tolerant "equal" there would silently break the enclosure guarantee. */
int interval_endpoint_cmp(Expr* a, Expr* b, bool* decidable) {
    *decidable = true;
    bool a_ni = iv_is_neg_inf(a), a_pi = iv_is_pos_inf(a);
    bool b_ni = iv_is_neg_inf(b), b_pi = iv_is_pos_inf(b);
    if (a_ni && b_ni) return 0;
    if (a_pi && b_pi) return 0;
    if (a_ni) return -1;          /* -Inf < b */
    if (a_pi) return 1;           /* +Inf > b */
    if (b_ni) return 1;           /* a > -Inf */
    if (b_pi) return -1;          /* a < +Inf */
    if (expr_eq(a, b)) return 0;

#ifdef USE_MPFR
    {
        /* Exact comparison at high precision — 200 bits dwarfs a double's 53,
         * so endpoints separated by a ULP compare with their true sign. */
        mpfr_t ma, mb;
        mpfr_inits2((mpfr_prec_t)200, ma, mb, (mpfr_ptr)0);
        bool oka = iv_endpoint_to_mpfr(a, ma);
        bool okb = oka && iv_endpoint_to_mpfr(b, mb);
        int c = (oka && okb) ? mpfr_cmp(ma, mb) : 0;
        bool both = oka && okb;
        mpfr_clears(ma, mb, (mpfr_ptr)0);
        if (both) return c < 0 ? -1 : (c > 0 ? 1 : 0);
    }
#endif
    /* No MPFR (or a symbolic endpoint MPFR couldn't reach): exact double
     * comparison via N[]. This is still tolerance-free. */
    double da, db;
    if (iv_to_double(a, &da) && iv_to_double(b, &db)) {
        if (da < db) return -1;
        if (da > db) return 1;
        return 0;
    }
    *decidable = false;
    return 0;
}

/* Sign of an endpoint: 1, -1, 0, or -2 (unknown). */
static int iv_sign(Expr* e) {
    if (iv_is_pos_inf(e)) return 1;
    if (iv_is_neg_inf(e)) return -1;
    if (e->type == EXPR_INTEGER) return e->data.integer > 0 ? 1 : (e->data.integer < 0 ? -1 : 0);
    if (e->type == EXPR_REAL) return e->data.real > 0 ? 1 : (e->data.real < 0 ? -1 : 0);
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    { int64_t n, d; if (is_rational(e, &n, &d)) return n > 0 ? 1 : (n < 0 ? -1 : 0); }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_zero_p(e->data.mpfr)) return 0;
        return mpfr_sgn(e->data.mpfr) > 0 ? 1 : -1;
    }
#endif
    { double d; if (iv_to_double(e, &d)) return d > 0 ? 1 : (d < 0 ? -1 : 0); }
    return -2;
}

/* Nudge an inexact endpoint one ULP outward; exact values pass through. */
static Expr* iv_widen(Expr* v, IvRound dir) {
    if (!v) return v;
    if (v->type == EXPR_REAL) {
        double d = v->data.real;
        if (isfinite(d)) {
            double w = nextafter(d, dir == RND_DOWN ? -INFINITY : INFINITY);
            expr_free(v);
            return expr_new_real(w);
        }
        return v;
    }
#ifdef USE_MPFR
    if (v->type == EXPR_MPFR) {
        if (mpfr_number_p(v->data.mpfr)) {
            if (dir == RND_DOWN) mpfr_nextbelow(v->data.mpfr);
            else mpfr_nextabove(v->data.mpfr);
        }
        return v;
    }
#endif
    return v;
}

/* Evaluate head[a,b] and widen an inexact result outward. */
static Expr* iv_binop(const char* head, Expr* a, Expr* b, IvRound dir) {
    Expr* call = expr_new_function(expr_new_symbol(head),
                                   (Expr*[]){ expr_copy(a), expr_copy(b) }, 2);
    Expr* r = eval_and_free(call);
    if (iv_is_inexact(r)) r = iv_widen(r, dir);
    return r;
}

/* Evaluate head[x] and widen an inexact result outward. */
static Expr* iv_unary(const char* head, Expr* x, IvRound dir) {
    Expr* call = expr_new_function(expr_new_symbol(head),
                                   (Expr*[]){ expr_copy(x) }, 1);
    Expr* r = eval_and_free(call);
    if (iv_is_inexact(r)) r = iv_widen(r, dir);
    return r;
}

/* Owned copy of the smaller / larger of a,b, or NULL if undecidable. */
static Expr* iv_min2(Expr* a, Expr* b) {
    bool d; int c = interval_endpoint_cmp(a, b, &d);
    if (!d) return NULL;
    return expr_copy(c <= 0 ? a : b);
}
static Expr* iv_max2(Expr* a, Expr* b) {
    bool d; int c = interval_endpoint_cmp(a, b, &d);
    if (!d) return NULL;
    return expr_copy(c >= 0 ? a : b);
}

/* ---------------------------------------------------------------------------
 * Pair collector + canonicalization
 * ------------------------------------------------------------------------- */

typedef struct { Expr** los; Expr** his; size_t n, cap; } IvBuild;

static void ivb_init(IvBuild* b) { b->los = NULL; b->his = NULL; b->n = 0; b->cap = 0; }

static void ivb_push(IvBuild* b, Expr* lo, Expr* hi) {   /* adopts lo, hi */
    if (b->n == b->cap) {
        b->cap = b->cap ? b->cap * 2 : 4;
        b->los = realloc(b->los, b->cap * sizeof(Expr*));
        b->his = realloc(b->his, b->cap * sizeof(Expr*));
    }
    b->los[b->n] = lo;
    b->his[b->n] = hi;
    b->n++;
}

static void ivb_free(IvBuild* b) {
    for (size_t i = 0; i < b->n; i++) { expr_free(b->los[i]); expr_free(b->his[i]); }
    free(b->los); free(b->his);
    b->los = b->his = NULL; b->n = b->cap = 0;
}

/* Build Interval[{lo,hi},...] from ordered, disjoint pairs, adopting endpoints.
 * Frees the builder's arrays. */
static Expr* ivb_finish(IvBuild* b) {
    Expr** args = malloc(b->n * sizeof(Expr*));
    for (size_t i = 0; i < b->n; i++) {
        args[i] = expr_new_function(expr_new_symbol(SYM_List),
                                    (Expr*[]){ b->los[i], b->his[i] }, 2);
    }
    Expr* iv = expr_new_function(expr_new_symbol(SYM_Interval), args, b->n);
    free(args);
    free(b->los); free(b->his);
    b->los = b->his = NULL; b->n = b->cap = 0;
    return iv;
}

/* Canonicalize the given pairs (consuming ownership of every los[i]/his[i] and
 * the two arrays). Orders each pair, sorts by lower endpoint, and merges
 * overlapping/touching pairs. Always returns a fresh Interval. */
static Expr* iv_canonicalize_pairs(Expr** los, Expr** his, size_t n) {
    /* 1. order endpoints within each pair */
    for (size_t i = 0; i < n; i++) {
        bool d; int c = interval_endpoint_cmp(los[i], his[i], &d);
        if (d && c > 0) { Expr* t = los[i]; los[i] = his[i]; his[i] = t; }
    }
    /* 2. stable insertion sort by lower endpoint */
    for (size_t i = 1; i < n; i++) {
        for (size_t j = i; j > 0; j--) {
            bool d; int c = interval_endpoint_cmp(los[j - 1], los[j], &d);
            if (d && c > 0) {
                Expr* t = los[j - 1]; los[j - 1] = los[j]; los[j] = t;
                t = his[j - 1]; his[j - 1] = his[j]; his[j] = t;
            } else break;
        }
    }
    /* 3. merge overlapping/touching runs */
    IvBuild out; ivb_init(&out);
    size_t i = 0;
    while (i < n) {
        Expr* clo = los[i];
        Expr* chi = his[i];
        i++;
        while (i < n) {
            bool d; int c = interval_endpoint_cmp(los[i], chi, &d);   /* next.lo vs cur.hi */
            if (d && c <= 0) {
                /* overlap or touch: extend hi to max(chi, his[i]) */
                bool dm; int cm = interval_endpoint_cmp(chi, his[i], &dm);
                if (dm) {
                    if (cm < 0) { expr_free(chi); chi = his[i]; }
                    else { expr_free(his[i]); }
                    expr_free(los[i]);
                    i++;
                } else break;   /* undecidable extent: keep pairs separate */
            } else break;
        }
        ivb_push(&out, clo, chi);
    }
    free(los); free(his);
    return ivb_finish(&out);
}

/* ---------------------------------------------------------------------------
 * Predicate, constructors, accessors
 * ------------------------------------------------------------------------- */

bool is_interval(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Interval &&
           e->data.function.arg_count >= 1;
}

Expr* make_interval_pair(Expr* lo, Expr* hi) {
    Expr* pair = expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){ lo, hi }, 2);
    return expr_new_function(expr_new_symbol(SYM_Interval),
                             (Expr*[]){ pair }, 1);
}

size_t interval_pair_count(const Expr* iv) {
    return is_interval(iv) ? iv->data.function.arg_count : 0;
}

static Expr* iv_pair(const Expr* iv, size_t k) { return iv->data.function.args[k]; }

Expr* interval_pair_lo(const Expr* iv, size_t k) {
    Expr* p = iv_pair(iv, k);
    if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
        p->data.function.head->data.symbol.name == SYM_List &&
        p->data.function.arg_count >= 1)
        return p->data.function.args[0];
    return p;
}

Expr* interval_pair_hi(const Expr* iv, size_t k) {
    Expr* p = iv_pair(iv, k);
    if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
        p->data.function.head->data.symbol.name == SYM_List &&
        p->data.function.arg_count >= 2)
        return p->data.function.args[1];
    return p;
}

Expr* interval_from_scalar(const Expr* x) {
    return make_interval_pair(expr_copy((Expr*)x), expr_copy((Expr*)x));
}

/* ---------------------------------------------------------------------------
 * builtin_interval — canonicalizer
 * ------------------------------------------------------------------------- */

Expr* builtin_interval(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return NULL;

    Expr** los = malloc(n * sizeof(Expr*));
    Expr** his = malloc(n * sizeof(Expr*));

    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        bool is_list = a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
                       a->data.function.head->data.symbol.name == SYM_List;
        if (is_list) {
            if (a->data.function.arg_count != 2) {
                /* malformed pair: leave symbolic */
                for (size_t j = 0; j < i; j++) { expr_free(los[j]); expr_free(his[j]); }
                free(los); free(his);
                return NULL;
            }
            los[i] = expr_copy(a->data.function.args[0]);
            his[i] = expr_copy(a->data.function.args[1]);
        } else {
            los[i] = expr_copy(a);
            his[i] = expr_copy(a);
        }
    }

    Expr* canon = iv_canonicalize_pairs(los, his, n);
    if (expr_eq(canon, res)) { expr_free(canon); return NULL; }
    return canon;
}

/* ---------------------------------------------------------------------------
 * Arithmetic kernels
 * ------------------------------------------------------------------------- */

/* Extended-real endpoint sum with outward rounding. */
static Expr* iv_add_ext(Expr* x, Expr* y, IvRound dir) {
    if (iv_is_inf(x) || iv_is_inf(y)) {
        if (iv_is_inf(x)) return iv_copy_inf(x);
        return iv_copy_inf(y);
    }
    return iv_binop("Plus", x, y, dir);
}

Expr* interval_add(const Expr* A, const Expr* B) {
    if (!is_interval(A) || !is_interval(B)) return NULL;
    size_t na = interval_pair_count(A), nb = interval_pair_count(B);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        for (size_t j = 0; j < nb; j++) {
            Expr* lo = iv_add_ext(interval_pair_lo(A, i), interval_pair_lo(B, j), RND_DOWN);
            Expr* hi = iv_add_ext(interval_pair_hi(A, i), interval_pair_hi(B, j), RND_UP);
            if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
            ivb_push(&out, lo, hi);
        }
    }
    /* hand the collected pairs to canonicalization */
    size_t n = out.n;
    Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Negate a single endpoint (exact — negation is exact in floating point). */
static Expr* iv_neg_endpoint(Expr* e) {
    if (iv_is_pos_inf(e)) return iv_neg_inf();
    if (iv_is_neg_inf(e)) return iv_pos_inf();
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_new_integer(-1), expr_copy(e) }, 2));
}

Expr* interval_negate(const Expr* A) {
    if (!is_interval(A)) return NULL;
    size_t na = interval_pair_count(A);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* lo = iv_neg_endpoint(interval_pair_hi(A, i));
        Expr* hi = iv_neg_endpoint(interval_pair_lo(A, i));
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Extended-real endpoint product with outward rounding. */
static Expr* iv_mul_ext(Expr* x, Expr* y, IvRound dir) {
    if (iv_is_exact_zero(x) || iv_is_exact_zero(y)) return expr_new_integer(0);
    if (iv_is_inf(x) || iv_is_inf(y)) {
        int sx = iv_sign(x), sy = iv_sign(y);
        if (sx == -2 || sy == -2) return NULL;
        return (sx * sy > 0) ? iv_pos_inf() : iv_neg_inf();
    }
    return iv_binop("Times", x, y, dir);
}

/* Product of two pairs, writing the lo and hi endpoints to outlo and outhi. */
static bool iv_pair_mul(Expr* a1, Expr* a2, Expr* b1, Expr* b2,
                        Expr** outlo, Expr** outhi) {
    Expr* pa[4] = { a1, a1, a2, a2 };
    Expr* pb[4] = { b1, b2, b1, b2 };
    Expr* lod[4] = {0}; Expr* hiu[4] = {0};
    bool ok = true;
    for (int i = 0; i < 4 && ok; i++) {
        lod[i] = iv_mul_ext(pa[i], pb[i], RND_DOWN);
        hiu[i] = iv_mul_ext(pa[i], pb[i], RND_UP);
        if (!lod[i] || !hiu[i]) ok = false;
    }
    Expr* lo = NULL; Expr* hi = NULL;
    if (ok) {
        lo = expr_copy(lod[0]);
        hi = expr_copy(hiu[0]);
        for (int i = 1; i < 4 && lo && hi; i++) {
            Expr* m = iv_min2(lo, lod[i]); expr_free(lo); lo = m;
            Expr* x = iv_max2(hi, hiu[i]); expr_free(hi); hi = x;
        }
        if (!lo || !hi) ok = false;
    }
    for (int i = 0; i < 4; i++) { expr_free(lod[i]); expr_free(hiu[i]); }
    if (!ok) { expr_free(lo); expr_free(hi); return false; }
    *outlo = lo; *outhi = hi;
    return true;
}

Expr* interval_multiply(const Expr* A, const Expr* B) {
    if (!is_interval(A) || !is_interval(B)) return NULL;
    size_t na = interval_pair_count(A), nb = interval_pair_count(B);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        for (size_t j = 0; j < nb; j++) {
            Expr* lo; Expr* hi;
            if (!iv_pair_mul(interval_pair_lo(A, i), interval_pair_hi(A, i),
                             interval_pair_lo(B, j), interval_pair_hi(B, j),
                             &lo, &hi)) { ivb_free(&out); return NULL; }
            ivb_push(&out, lo, hi);
        }
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Reciprocal of a single endpoint (1/x = x^-1), rounding outward; 1/(+-Inf)=0. */
static Expr* iv_recip_endpoint(Expr* x, IvRound dir) {
    if (iv_is_inf(x)) return expr_new_integer(0);
    Expr* m1 = expr_new_integer(-1);
    Expr* r = iv_binop("Power", x, m1, dir);
    expr_free(m1);
    return r;
}

Expr* interval_reciprocal(const Expr* A) {
    if (!is_interval(A)) return NULL;
    size_t na = interval_pair_count(A);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(A, i);
        Expr* b = interval_pair_hi(A, i);
        int sa = iv_sign(a), sb = iv_sign(b);
        if (sa == -2 || sb == -2) { ivb_free(&out); return NULL; }

        if (sa > 0 || sb < 0) {
            /* no zero inside: reciprocal is [1/b, 1/a] */
            Expr* lo = iv_recip_endpoint(b, RND_DOWN);
            Expr* hi = iv_recip_endpoint(a, RND_UP);
            ivb_push(&out, lo, hi);
        } else if (sa == 0 && sb > 0) {
            ivb_push(&out, iv_recip_endpoint(b, RND_DOWN), iv_pos_inf());
        } else if (sa < 0 && sb == 0) {
            ivb_push(&out, iv_neg_inf(), iv_recip_endpoint(a, RND_UP));
        } else if (sa < 0 && sb > 0) {
            /* straddles 0: (-Inf, 1/a] U [1/b, +Inf) */
            ivb_push(&out, iv_neg_inf(), iv_recip_endpoint(a, RND_UP));
            ivb_push(&out, iv_recip_endpoint(b, RND_DOWN), iv_pos_inf());
        } else {
            /* [0,0]: 1/0 is unbounded */
            ivb_push(&out, iv_neg_inf(), iv_pos_inf());
        }
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* x^n for a single endpoint, rounding outward. n passed as an int. */
static Expr* iv_pow_endpoint(Expr* x, int64_t n, IvRound dir) {
    if (iv_is_inf(x)) {
        /* (+Inf)^n = +Inf ; (-Inf)^n = +Inf (even) / -Inf (odd), n>0 */
        if (iv_is_pos_inf(x)) return iv_pos_inf();
        return (n % 2 == 0) ? iv_pos_inf() : iv_neg_inf();
    }
    Expr* ne = expr_new_integer(n);
    Expr* r = iv_binop("Power", x, ne, dir);
    expr_free(ne);
    return r;
}

Expr* interval_power_int(const Expr* A, int64_t n) {
    if (!is_interval(A)) return NULL;
    if (n == 0) return make_interval_pair(expr_new_integer(1), expr_new_integer(1));
    if (n < 0) {
        Expr* rec = interval_reciprocal(A);
        if (!rec) return NULL;
        Expr* r = interval_power_int(rec, -n);
        expr_free(rec);
        return r;
    }
    bool odd = (n % 2 != 0);
    size_t na = interval_pair_count(A);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(A, i);
        Expr* b = interval_pair_hi(A, i);
        int sa = iv_sign(a), sb = iv_sign(b);
        if (sa == -2 || sb == -2) { ivb_free(&out); return NULL; }
        if (odd) {
            ivb_push(&out, iv_pow_endpoint(a, n, RND_DOWN), iv_pow_endpoint(b, n, RND_UP));
        } else if (sa >= 0) {
            ivb_push(&out, iv_pow_endpoint(a, n, RND_DOWN), iv_pow_endpoint(b, n, RND_UP));
        } else if (sb <= 0) {
            ivb_push(&out, iv_pow_endpoint(b, n, RND_DOWN), iv_pow_endpoint(a, n, RND_UP));
        } else {
            /* straddles 0, even power: [0, max(a^n, b^n)] */
            Expr* pa = iv_pow_endpoint(a, n, RND_UP);
            Expr* pb = iv_pow_endpoint(b, n, RND_UP);
            Expr* hi = iv_max2(pa, pb);
            expr_free(pa); expr_free(pb);
            if (!hi) { ivb_free(&out); return NULL; }
            ivb_push(&out, expr_new_integer(0), hi);
        }
    }
    size_t n2 = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n2);
}

Expr* interval_power_pos_exp(const Expr* A, const Expr* p) {
    if (!is_interval(A)) return NULL;
    if (iv_sign((Expr*)p) <= 0) return NULL;          /* need p > 0 */
    if (iv_sign(interval_pair_lo(A, 0)) < 0) return NULL; /* need A >= 0 */
    size_t na = interval_pair_count(A);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(A, i);
        Expr* b = interval_pair_hi(A, i);
        /* x^p is increasing for x >= 0, p > 0 */
        Expr* lo = iv_binop("Power", a, (Expr*)p, RND_DOWN);
        Expr* hi = iv_binop("Power", b, (Expr*)p, RND_UP);
        if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* ---------------------------------------------------------------------------
 * Unary function threading
 * ------------------------------------------------------------------------- */

/* True if an endpoint is not a real value -- used to reject a monotone thread
 * when a function is applied outside its real domain. An explicit Complex[re,im]
 * is caught directly; a symbolic result (ArcCos[5], Log[-2] = Log[2] + I Pi, ...)
 * that is really complex is caught by numericalizing and testing the imaginary
 * part. Real number kinds and +/-Infinity short-circuit as real. */
static bool iv_is_complex_result(Expr* e) {
    if (!e) return false;
    if (is_complex(e, NULL, NULL)) return true;
    if (iv_is_pos_inf(e) || iv_is_neg_inf(e)) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return false;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return false;
#endif
    { int64_t nn, dd; if (is_rational(e, &nn, &dd)) return false; }
    /* Symbolic: numericalize and inspect the imaginary part. */
    Expr* n = eval_and_free(expr_new_function(expr_new_symbol("N"),
                                              (Expr*[]){ expr_copy(e) }, 1));
    bool nonreal = false;
    Expr *re, *im;
    if (n && is_complex(n, &re, &im)) {
        double iv2;
        if (iv_to_double(im, &iv2)) nonreal = (iv2 > 1e-9 || iv2 < -1e-9);
        else nonreal = true;
    }
    expr_free(n);
    return nonreal;
}

Expr* interval_thread_monotone(const Expr* iv, IvMonotoneDir dir, const char* head) {
    if (!is_interval(iv)) return NULL;
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        Expr* lo; Expr* hi;
        if (dir == IV_INC) {
            lo = iv_unary(head, a, RND_DOWN);
            hi = iv_unary(head, b, RND_UP);
        } else {
            lo = iv_unary(head, b, RND_DOWN);
            hi = iv_unary(head, a, RND_UP);
        }
        /* Leaving the reals (e.g. Log of a negative endpoint) -> stay symbolic. */
        if (!lo || !hi || iv_is_complex_result(lo) || iv_is_complex_result(hi)) {
            expr_free(lo); expr_free(hi); ivb_free(&out); return NULL;
        }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Does some x = phase + period*k satisfy la <= x <= lb, with a small inclusion
 * bias so a critical point on the boundary counts (keeps the enclosure safe)? */
static bool iv_crit_in(double la, double lb, double phase, double period) {
    if (!isfinite(la) || !isfinite(lb)) return true;   /* infinite span hits every crit point */
    double eps = 1e-9 * (fabs(la) + fabs(lb) + 1.0);
    double klo = ceil((la - phase) / period - eps);
    double khi = floor((lb - phase) / period + eps);
    return klo <= khi;
}

typedef enum { IV_SIN, IV_COS } IvTrigKind;

static Expr* iv_trig(const Expr* iv, IvTrigKind kind) {
    const char* head = (kind == IV_SIN) ? "Sin" : "Cos";
    double crest_phase = (kind == IV_SIN) ? M_PI / 2.0 : 0.0;         /* f = +1 */
    double trough_phase = (kind == IV_SIN) ? 3.0 * M_PI / 2.0 : M_PI; /* f = -1 */
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        double la, lb;
        if (!iv_to_double(a, &la) || !iv_to_double(b, &lb)) { ivb_free(&out); return NULL; }

        bool crest = iv_crit_in(la, lb, crest_phase, 2.0 * M_PI);
        bool trough = iv_crit_in(la, lb, trough_phase, 2.0 * M_PI);

        Expr* hi;
        if (crest) hi = expr_new_integer(1);
        else {
            Expr* fa = iv_unary(head, a, RND_UP);
            Expr* fb = iv_unary(head, b, RND_UP);
            hi = (fa && fb) ? iv_max2(fa, fb) : NULL;
            expr_free(fa); expr_free(fb);
        }
        Expr* lo;
        if (trough) lo = expr_new_integer(-1);
        else {
            Expr* fa = iv_unary(head, a, RND_DOWN);
            Expr* fb = iv_unary(head, b, RND_DOWN);
            lo = (fa && fb) ? iv_min2(fa, fb) : NULL;
            expr_free(fa); expr_free(fb);
        }
        if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Tan (poles at pi/2 + k*pi, increasing between them) and Cot (poles at k*pi,
 * decreasing between them): a pole inside the pair makes the range unbounded,
 * otherwise the function is monotone across it. */
static Expr* iv_tan_cot(const Expr* iv, bool is_cot) {
    const char* head = is_cot ? "Cot" : "Tan";
    double pole_phase = is_cot ? 0.0 : (M_PI / 2.0);
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        double la, lb;
        if (!iv_to_double(a, &la) || !iv_to_double(b, &lb)) { ivb_free(&out); return NULL; }
        if (iv_crit_in(la, lb, pole_phase, M_PI)) {
            ivb_push(&out, iv_neg_inf(), iv_pos_inf());
        } else {
            /* Tan increases, Cot decreases, across a pole-free pair. */
            Expr* lo = iv_unary(head, is_cot ? b : a, RND_DOWN);
            Expr* hi = iv_unary(head, is_cot ? a : b, RND_UP);
            if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
            ivb_push(&out, lo, hi);
        }
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

static Expr* iv_abs(const Expr* iv) {
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        int sa = iv_sign(a), sb = iv_sign(b);
        if (sa == -2 || sb == -2) { ivb_free(&out); return NULL; }
        if (sa >= 0) {
            ivb_push(&out, expr_copy(a), expr_copy(b));
        } else if (sb <= 0) {
            ivb_push(&out, iv_neg_endpoint(b), iv_neg_endpoint(a));
        } else {
            Expr* na_ = iv_neg_endpoint(a);   /* |a| = -a since a < 0 */
            Expr* hi = iv_max2(na_, b);
            expr_free(na_);
            if (!hi) { ivb_free(&out); return NULL; }
            ivb_push(&out, expr_new_integer(0), hi);
        }
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Cosh: even, with its minimum 1 at 0 (like Abs but with Cosh applied). */
static Expr* iv_cosh(const Expr* iv) {
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        int sa = iv_sign(a), sb = iv_sign(b);
        if (sa == -2 || sb == -2) { ivb_free(&out); return NULL; }
        Expr* lo; Expr* hi;
        if (sa >= 0) {
            lo = iv_unary("Cosh", a, RND_DOWN); hi = iv_unary("Cosh", b, RND_UP);
        } else if (sb <= 0) {
            lo = iv_unary("Cosh", b, RND_DOWN); hi = iv_unary("Cosh", a, RND_UP);
        } else {
            lo = expr_new_integer(1);   /* Cosh[0] = 1, the exact minimum */
            Expr* ca = iv_unary("Cosh", a, RND_UP);
            Expr* cb = iv_unary("Cosh", b, RND_UP);
            hi = (ca && cb) ? iv_max2(ca, cb) : NULL;
            expr_free(ca); expr_free(cb);
        }
        if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Reciprocal of an owned base interval (frees base); NULL-safe. Used to build
 * Sec = 1/Cos, Csc = 1/Sin, Sech = 1/Cosh, Csch = 1/Sinh, Coth = 1/Tanh -- the
 * reciprocal kernel handles the pole (where the base straddles 0) as a disjoint
 * union, so it captures each function's singular behaviour for free. */
static Expr* iv_recip_of(Expr* base) {
    if (!base) return NULL;
    Expr* r = interval_reciprocal(base);
    expr_free(base);
    return r;
}

/* Inverse-reciprocal functions are inverse_base(1/x): ArcCot = ArcTan(1/x),
 * ArcSec = ArcCos(1/x), ArcCsc = ArcSin(1/x), ArcCoth = ArcTanh(1/x),
 * ArcSech = ArcCosh(1/x), ArcCsch = ArcSinh(1/x). Composing the (already
 * threaded) inverse base with the reciprocal kernel handles every subtlety for
 * free: reciprocal turns a straddling interval into a disjoint union (the
 * discontinuity of ArcCot / ArcCsch), and the inverse base's [-1,1] / [1,inf)
 * domain check sends out-of-domain intervals to symbolic. */
static Expr* iv_apply_of_recip(const Expr* iv, const char* inv_head) {
    Expr* rec = interval_reciprocal(iv);
    if (!rec) return NULL;
    Expr* r = interval_apply_function(inv_head, rec);
    expr_free(rec);
    return r;
}

/* Evaluate head[x] (or head[prefix, x] when has_prefix), widen inexact outward. */
static Expr* iv_fn_endpoint(const char* head, int64_t prefix, bool has_prefix,
                            Expr* x, IvRound dir) {
    Expr* call;
    if (has_prefix)
        call = expr_new_function(expr_new_symbol(head),
                   (Expr*[]){ expr_new_integer(prefix), expr_copy(x) }, 2);
    else
        call = expr_new_function(expr_new_symbol(head), (Expr*[]){ expr_copy(x) }, 1);
    Expr* r = eval_and_free(call);
    if (iv_is_inexact(r)) r = iv_widen(r, dir);
    return r;
}

/* Thread a function monotone in direction `dir`, but only when every pair's
 * numeric range lies strictly inside the open region (rlo, rhi) where that
 * monotonicity holds; else NULL (stays symbolic). This keeps a function threaded
 * on only part of its domain (Gamma, Zeta, ...) from ever producing a wrong
 * bound near an extremum or a pole. */
static Expr* iv_thread_region(const Expr* iv, IvMonotoneDir dir, const char* head,
                              int64_t prefix, bool has_prefix,
                              double rlo, double rhi) {
    size_t na = interval_pair_count(iv);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        Expr* a = interval_pair_lo(iv, i);
        Expr* b = interval_pair_hi(iv, i);
        double la, lb;
        if (!iv_to_double(a, &la) || !iv_to_double(b, &lb)) { ivb_free(&out); return NULL; }
        if (!(la > rlo && lb < rhi)) { ivb_free(&out); return NULL; }
        Expr* lo; Expr* hi;
        if (dir == IV_INC) {
            lo = iv_fn_endpoint(head, prefix, has_prefix, a, RND_DOWN);
            hi = iv_fn_endpoint(head, prefix, has_prefix, b, RND_UP);
        } else {
            lo = iv_fn_endpoint(head, prefix, has_prefix, b, RND_DOWN);
            hi = iv_fn_endpoint(head, prefix, has_prefix, a, RND_UP);
        }
        if (!lo || !hi || iv_is_complex_result(lo) || iv_is_complex_result(hi)) {
            expr_free(lo); expr_free(hi); ivb_free(&out); return NULL;
        }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Gamma / LogGamma: convex on (0, inf) with a minimum at IV_GAMMA_MIN; monotone
 * decreasing below it and increasing above. Thread on whichever branch the
 * interval lies safely within (a small margin around the minimum stays
 * symbolic); intervals reaching <= 0 or straddling the minimum stay symbolic. */
static Expr* iv_gamma_like(const Expr* iv, const char* head) {
    Expr* r = iv_thread_region(iv, IV_INC, head, 0, false, IV_GAMMA_MIN + 1e-6, INFINITY);
    if (r) return r;
    return iv_thread_region(iv, IV_DEC, head, 0, false, 0.0, IV_GAMMA_MIN - 1e-6);
}

Expr* interval_apply_function(const char* head, const Expr* iv) {
    if (!is_interval(iv)) return NULL;
    /* non-monotone (range analysis) */
    if (strcmp(head, "Sin") == 0)  return iv_trig(iv, IV_SIN);
    if (strcmp(head, "Cos") == 0)  return iv_trig(iv, IV_COS);
    if (strcmp(head, "Tan") == 0)  return iv_tan_cot(iv, false);
    if (strcmp(head, "Cot") == 0)  return iv_tan_cot(iv, true);
    if (strcmp(head, "Abs") == 0)  return iv_abs(iv);
    if (strcmp(head, "Cosh") == 0) return iv_cosh(iv);
    /* reciprocal trig / hyperbolic, built as 1/base (pole handled by reciprocal) */
    if (strcmp(head, "Sec") == 0)  return iv_recip_of(iv_trig(iv, IV_COS));
    if (strcmp(head, "Csc") == 0)  return iv_recip_of(iv_trig(iv, IV_SIN));
    if (strcmp(head, "Sech") == 0) return iv_recip_of(iv_cosh(iv));
    if (strcmp(head, "Csch") == 0) return iv_recip_of(interval_thread_monotone(iv, IV_INC, "Sinh"));
    if (strcmp(head, "Coth") == 0) return iv_recip_of(interval_thread_monotone(iv, IV_INC, "Tanh"));
    /* monotone increasing (domain-invalid endpoints -> complex result -> NULL) */
    if (strcmp(head, "Exp") == 0)      return interval_thread_monotone(iv, IV_INC, "Exp");
    if (strcmp(head, "Log") == 0)      return interval_thread_monotone(iv, IV_INC, "Log");
    if (strcmp(head, "Sqrt") == 0)     return interval_thread_monotone(iv, IV_INC, "Sqrt");
    if (strcmp(head, "Sinh") == 0)     return interval_thread_monotone(iv, IV_INC, "Sinh");
    if (strcmp(head, "Tanh") == 0)     return interval_thread_monotone(iv, IV_INC, "Tanh");
    if (strcmp(head, "ArcTan") == 0)   return interval_thread_monotone(iv, IV_INC, "ArcTan");
    if (strcmp(head, "ArcSin") == 0)   return interval_thread_monotone(iv, IV_INC, "ArcSin");
    if (strcmp(head, "ArcSinh") == 0)  return interval_thread_monotone(iv, IV_INC, "ArcSinh");
    if (strcmp(head, "ArcCosh") == 0)  return interval_thread_monotone(iv, IV_INC, "ArcCosh");
    if (strcmp(head, "ArcTanh") == 0)  return interval_thread_monotone(iv, IV_INC, "ArcTanh");
    /* monotone increasing on all of R */
    if (strcmp(head, "Erf") == 0)      return interval_thread_monotone(iv, IV_INC, "Erf");
    if (strcmp(head, "Sign") == 0)     return interval_thread_monotone(iv, IV_INC, "Sign");
    if (strcmp(head, "Floor") == 0)    return interval_thread_monotone(iv, IV_INC, "Floor");
    if (strcmp(head, "Ceiling") == 0)  return interval_thread_monotone(iv, IV_INC, "Ceiling");
    /* inverse-reciprocal = inverse_base(1/x) */
    if (strcmp(head, "ArcCot") == 0)   return iv_apply_of_recip(iv, "ArcTan");
    if (strcmp(head, "ArcSec") == 0)   return iv_apply_of_recip(iv, "ArcCos");
    if (strcmp(head, "ArcCsc") == 0)   return iv_apply_of_recip(iv, "ArcSin");
    if (strcmp(head, "ArcCoth") == 0)  return iv_apply_of_recip(iv, "ArcTanh");
    if (strcmp(head, "ArcSech") == 0)  return iv_apply_of_recip(iv, "ArcCosh");
    if (strcmp(head, "ArcCsch") == 0)  return iv_apply_of_recip(iv, "ArcSinh");
    /* monotone decreasing */
    if (strcmp(head, "ArcCos") == 0)   return interval_thread_monotone(iv, IV_DEC, "ArcCos");
    if (strcmp(head, "Erfc") == 0)     return interval_thread_monotone(iv, IV_DEC, "Erfc");
    /* monotone only on a sub-region (else symbolic) */
    if (strcmp(head, "Gamma") == 0)    return iv_gamma_like(iv, "Gamma");
    if (strcmp(head, "LogGamma") == 0) return iv_gamma_like(iv, "LogGamma");
    if (strcmp(head, "Zeta") == 0)     return iv_thread_region(iv, IV_DEC, "Zeta", 0, false, 1.0, INFINITY);
    return NULL;
}

Expr* interval_polygamma(int64_t n, const Expr* iv) {
    if (!is_interval(iv) || n < 0) return NULL;
    /* PolyGamma(n, .) is monotone on (0, inf): its derivative is PolyGamma(n+1,.)
     * with sign (-1)^n, so it increases for even n and decreases for odd n. */
    IvMonotoneDir dir = (n % 2 == 0) ? IV_INC : IV_DEC;
    return iv_thread_region(iv, dir, "PolyGamma", n, true, 0.0, INFINITY);
}

int interval_compare_intervals(const Expr* A, const Expr* B, bool* decidable) {
    *decidable = false;
    size_t na = interval_pair_count(A), nb = interval_pair_count(B);
    if (na == 0 || nb == 0) return 0;
    Expr* loA = interval_pair_lo(A, 0);
    Expr* hiA = interval_pair_hi(A, na - 1);
    Expr* loB = interval_pair_lo(B, 0);
    Expr* hiB = interval_pair_hi(B, nb - 1);
    bool d;
    int c = interval_endpoint_cmp(hiA, loB, &d);
    if (d && c < 0) { *decidable = true; return -1; }   /* A entirely below B */
    c = interval_endpoint_cmp(loA, hiB, &d);
    if (d && c > 0) { *decidable = true; return 1; }     /* A entirely above B */
    /* Equal only when both are the same single point. */
    if (na == 1 && nb == 1 &&
        expr_eq(loA, hiA) && expr_eq(loB, hiB) && expr_eq(loA, loB)) {
        *decidable = true; return 0;
    }
    return 0;
}

/* base^Interval for a positive scalar base: base^x is monotone in x, so the
 * extremes are at the endpoints. */
Expr* interval_scalar_pow_interval(const Expr* base, const Expr* B) {
    if (!is_interval(B)) return NULL;
    if (iv_sign((Expr*)base) <= 0) return NULL;   /* need base > 0 */
    size_t nb = interval_pair_count(B);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < nb; i++) {
        Expr* c = interval_pair_lo(B, i);
        Expr* d = interval_pair_hi(B, i);
        Expr* pc_lo = iv_binop("Power", (Expr*)base, c, RND_DOWN);
        Expr* pd_lo = iv_binop("Power", (Expr*)base, d, RND_DOWN);
        Expr* pc_hi = iv_binop("Power", (Expr*)base, c, RND_UP);
        Expr* pd_hi = iv_binop("Power", (Expr*)base, d, RND_UP);
        Expr* lo = (pc_lo && pd_lo) ? iv_min2(pc_lo, pd_lo) : NULL;
        Expr* hi = (pc_hi && pd_hi) ? iv_max2(pc_hi, pd_hi) : NULL;
        expr_free(pc_lo); expr_free(pd_lo); expr_free(pc_hi); expr_free(pd_hi);
        if (!lo || !hi) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
        ivb_push(&out, lo, hi);
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* Interval^Interval for a strictly-positive base interval: over each pair-pair
 * the extremes are at the four corners a^c, a^d, b^c, b^d. */
Expr* interval_power_interval(const Expr* A, const Expr* B) {
    if (!is_interval(A) || !is_interval(B)) return NULL;
    if (iv_sign(interval_pair_lo(A, 0)) <= 0) return NULL;   /* need A > 0 */
    size_t na = interval_pair_count(A), nb = interval_pair_count(B);
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < na; i++) {
        for (size_t j = 0; j < nb; j++) {
            Expr* a = interval_pair_lo(A, i); Expr* b = interval_pair_hi(A, i);
            Expr* c = interval_pair_lo(B, j); Expr* d = interval_pair_hi(B, j);
            Expr* corner_lo[4], *corner_hi[4];
            Expr* xa[4] = { a, a, b, b }; Expr* xe[4] = { c, d, c, d };
            bool ok = true;
            for (int k = 0; k < 4; k++) {
                corner_lo[k] = iv_binop("Power", xa[k], xe[k], RND_DOWN);
                corner_hi[k] = iv_binop("Power", xa[k], xe[k], RND_UP);
                if (!corner_lo[k] || !corner_hi[k]) ok = false;
            }
            Expr* lo = NULL; Expr* hi = NULL;
            if (ok) {
                lo = expr_copy(corner_lo[0]); hi = expr_copy(corner_hi[0]);
                for (int k = 1; k < 4 && lo && hi; k++) {
                    Expr* m = iv_min2(lo, corner_lo[k]); expr_free(lo); lo = m;
                    Expr* x = iv_max2(hi, corner_hi[k]); expr_free(hi); hi = x;
                }
                if (!lo || !hi) ok = false;
            }
            for (int k = 0; k < 4; k++) { expr_free(corner_lo[k]); expr_free(corner_hi[k]); }
            if (!ok) { expr_free(lo); expr_free(hi); ivb_free(&out); return NULL; }
            ivb_push(&out, lo, hi);
        }
    }
    size_t n = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, n);
}

/* ---------------------------------------------------------------------------
 * Comparison against a scalar
 * ------------------------------------------------------------------------- */

int interval_compare_scalar(const Expr* iv, Expr* c, bool* decidable) {
    *decidable = false;
    size_t n = interval_pair_count(iv);
    if (n == 0) return 0;
    Expr* lo = interval_pair_lo(iv, 0);
    Expr* hi = interval_pair_hi(iv, n - 1);
    bool d;
    int chi = interval_endpoint_cmp(hi, c, &d);
    if (d && chi < 0) { *decidable = true; return -1; }   /* entirely below c */
    int clo = interval_endpoint_cmp(lo, c, &d);
    if (d && clo > 0) { *decidable = true; return 1; }    /* entirely above c */
    return 0;
}

/* ---------------------------------------------------------------------------
 * Companion set operations: IntervalUnion / IntervalIntersection / MemberQ
 * ------------------------------------------------------------------------- */

/* Append every pair of `arg` (an interval, or a scalar as a point) to `out`.
 * Returns false if `arg` is neither. */
static bool iv_collect(Expr* arg, IvBuild* out) {
    if (is_interval(arg)) {
        size_t k = interval_pair_count(arg);
        for (size_t i = 0; i < k; i++)
            ivb_push(out, expr_copy(interval_pair_lo(arg, i)), expr_copy(interval_pair_hi(arg, i)));
        return true;
    }
    if (expr_is_numeric_like(arg) || iv_is_inf(arg)) {
        ivb_push(out, expr_copy(arg), expr_copy(arg));
        return true;
    }
    return false;
}

Expr* builtin_intervalunion(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    IvBuild out; ivb_init(&out);
    for (size_t i = 0; i < n; i++) {
        if (!iv_collect(res->data.function.args[i], &out)) { ivb_free(&out); return NULL; }
    }
    if (out.n == 0) return expr_new_function(expr_new_symbol(SYM_Interval), NULL, 0);
    size_t m = out.n; Expr** los = out.los; Expr** his = out.his;
    out.los = out.his = NULL; out.n = out.cap = 0;
    return iv_canonicalize_pairs(los, his, m);
}

/* Intersect the pair-list {los,his}[0..n) with one pair [c,d], appending the
 * overlaps to `dst`. */
static void iv_intersect_pair(Expr** los, Expr** his, size_t n, Expr* c, Expr* d, IvBuild* dst) {
    for (size_t i = 0; i < n; i++) {
        Expr* lo = iv_max2(los[i], c);    /* max of the two lower endpoints */
        Expr* hi = iv_min2(his[i], d);    /* min of the two upper endpoints */
        if (lo && hi) {
            bool dd; int c2 = interval_endpoint_cmp(lo, hi, &dd);
            if (dd && c2 <= 0) { ivb_push(dst, lo, hi); continue; }  /* non-empty overlap */
        }
        expr_free(lo); expr_free(hi);
    }
}

Expr* builtin_intervalintersection(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return NULL;
    /* Accumulator starts as the first interval's pairs. */
    IvBuild acc; ivb_init(&acc);
    if (!iv_collect(res->data.function.args[0], &acc)) { ivb_free(&acc); return NULL; }
    for (size_t i = 1; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        IvBuild nextiv; ivb_init(&nextiv);
        if (!iv_collect(arg, &nextiv)) { ivb_free(&nextiv); ivb_free(&acc); return NULL; }
        IvBuild result; ivb_init(&result);
        for (size_t j = 0; j < nextiv.n; j++)
            iv_intersect_pair(acc.los, acc.his, acc.n, nextiv.los[j], nextiv.his[j], &result);
        ivb_free(&nextiv);
        ivb_free(&acc);
        acc = result;
    }
    if (acc.n == 0) { ivb_free(&acc); return expr_new_function(expr_new_symbol(SYM_Interval), NULL, 0); }
    size_t m = acc.n; Expr** los = acc.los; Expr** his = acc.his;
    acc.los = acc.his = NULL; acc.n = acc.cap = 0;
    return iv_canonicalize_pairs(los, his, m);
}

Expr* builtin_intervalmemberq(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* iv = res->data.function.args[0];
    Expr* x  = res->data.function.args[1];
    if (!is_interval(iv)) return NULL;
    size_t n = interval_pair_count(iv);

    if (is_interval(x)) {
        /* subset test: every pair of x must sit inside some pair of iv */
        size_t m = interval_pair_count(x);
        for (size_t j = 0; j < m; j++) {
            Expr* xl = interval_pair_lo(x, j);
            Expr* xh = interval_pair_hi(x, j);
            bool inside = false, undecided = false;
            for (size_t i = 0; i < n; i++) {
                bool d1, d2;
                int c1 = interval_endpoint_cmp(interval_pair_lo(iv, i), xl, &d1);
                int c2 = interval_endpoint_cmp(xh, interval_pair_hi(iv, i), &d2);
                if (!d1 || !d2) { undecided = true; continue; }
                if (c1 <= 0 && c2 <= 0) { inside = true; break; }
            }
            if (!inside) {
                if (undecided) return NULL;         /* can't decide */
                return expr_new_symbol(SYM_False);  /* this piece escapes iv */
            }
        }
        return expr_new_symbol(SYM_True);
    }

    /* scalar membership: x in some [lo, hi] */
    bool undecided = false;
    for (size_t i = 0; i < n; i++) {
        bool d1, d2;
        int c1 = interval_endpoint_cmp(interval_pair_lo(iv, i), x, &d1);
        int c2 = interval_endpoint_cmp(x, interval_pair_hi(iv, i), &d2);
        if (!d1 || !d2) { undecided = true; continue; }
        if (c1 <= 0 && c2 <= 0) return expr_new_symbol(SYM_True);
    }
    if (undecided) return NULL;
    return expr_new_symbol(SYM_False);
}

/* ---------------------------------------------------------------------------
 * Registration
 * ------------------------------------------------------------------------- */

void interval_init(void) {
    symtab_add_builtin("Interval", builtin_interval);
    SymbolDef* def = symtab_get_def("Interval");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Interval",
        "Interval[{min, max}] represents the range of real values between min and "
        "max, inclusive. Interval[{a1,b1}, {a2,b2}, ...] is the union of the "
        "ranges. Arithmetic and elementary functions thread through intervals, "
        "producing rigorous enclosures; exact endpoints are kept exact.");

    symtab_add_builtin("IntervalUnion", builtin_intervalunion);
    def = symtab_get_def("IntervalUnion");
    if (def) def->attributes |= ATTR_PROTECTED | ATTR_ORDERLESS | ATTR_FLAT;
    symtab_set_docstring("IntervalUnion",
        "IntervalUnion[i1, i2, ...] gives the interval representing the union of "
        "the intervals ij.");

    symtab_add_builtin("IntervalIntersection", builtin_intervalintersection);
    def = symtab_get_def("IntervalIntersection");
    if (def) def->attributes |= ATTR_PROTECTED | ATTR_ORDERLESS | ATTR_FLAT;
    symtab_set_docstring("IntervalIntersection",
        "IntervalIntersection[i1, i2, ...] gives the interval representing the "
        "intersection of the intervals ij (Interval[] if they are disjoint).");

    symtab_add_builtin("IntervalMemberQ", builtin_intervalmemberq);
    def = symtab_get_def("IntervalMemberQ");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("IntervalMemberQ",
        "IntervalMemberQ[interval, x] gives True if x lies within interval, and "
        "False otherwise. IntervalMemberQ[interval, other] tests whether the "
        "interval other is wholly contained in interval.");
}
