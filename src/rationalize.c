/* PicoCAS — Rationalize implementation.
 *
 * Algorithm sketch
 * ----------------
 * The core of Rationalize is the classical "best rational approximation"
 * problem. Two flavours are needed:
 *
 *   (A) Rationalize[x, dx] — given a tolerance dx ≥ 0, find the rational
 *       p/q with the smallest q satisfying |p/q − x| ≤ dx.
 *
 *       This is solved by the Stern–Brocot / continued-fraction "simplest
 *       rational in interval" algorithm. We compute the continued fraction
 *       of x by repeatedly extracting integer parts of the lower and
 *       upper bound, descending until the integer parts diverge — at which
 *       point the smallest integer in the lower interval is the simplest
 *       rational. We use GMP throughout so denominators that exceed
 *       int64_t (achievable for very small dx, e.g. dx = 0) are handled.
 *
 *   (B) Rationalize[x] — no tolerance: walk the convergents of the
 *       continued-fraction expansion of x, returning the first p/q whose
 *       error satisfies |p/q − x| < 10^-4 / q^2 (the same threshold
 *       Mathematica uses). If none exist within the precision of x, the
 *       caller is told (returns false) and is expected to leave x alone.
 *
 * dx = 0 is treated as a thin tolerance derived from the ulp of x — this
 * lets `Rationalize[N[Pi], 0]` produce 245850922/78256779 rather than
 * the bit-exact (and uninteresting) 884279719003555/281474976710656.
 *
 * Threading
 * ---------
 * Rationalize maps over expression structure. The descent picks one of
 * two strategies per node:
 *   - Compound expressions that are themselves a NumericQ (e.g. Sqrt[2],
 *     Pi, Exp[Sqrt[2]]) are first numericalised via the existing N
 *     pipeline, then rationalised end-to-end. This is what makes
 *     `Rationalize[Exp[Sqrt[2]], 2^-12]` collapse to 218/53.
 *   - Otherwise the head and arguments are recursively rationalised
 *     and then re-evaluated, so e.g. `Rationalize[1.2 + 6.7 x]` becomes
 *     `6/5 + (67 x)/10` after Plus/Times re-canonicalise.
 *
 * In default mode (no dx) only inexact leaves (EXPR_REAL / EXPR_MPFR) get
 * converted. Exact numbers and symbolic atoms pass through. This matches
 * Mathematica: Rationalize[Pi] = Pi, Rationalize[3.14] = 157/50, ...
 *
 * Memory
 * ------
 * The core algorithm allocates / clears its own mpz_t scratch. The caller
 * owns the returned Expr from internal_rationalize_expr / builtin_*. The
 * builtin returns NULL when its input is malformed (caller retains the
 * res Expr); on a successful conversion it returns a fresh Expr and the
 * evaluator releases res itself — the builtin must not free res.
 */

#include "rationalize.h"

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "numeric.h"
#include "symtab.h"
#include "sym_names.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------------
 *  Numeric extraction
 *
 *  Pull a `double` out of any concrete numeric Expr we know how to handle.
 *  Returns true on success.
 * ---------------------------------------------------------------------- */
static bool extract_double(const Expr* e, double* out) {
    if (!e || !out) return false;
    int64_t n, d;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer;       return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);     return true;
        case EXPR_REAL:    *out = e->data.real;                  return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_rational((Expr*)e, &n, &d)) {
                *out = (double)n / (double)d;
                return true;
            }
            return false;
        default:
            return false;
    }
}

/* True iff `e` is an inexact numeric leaf (EXPR_REAL or EXPR_MPFR). */
static bool is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True iff `e` is one of the hold heads N already passes through unchanged. */
static bool is_hold_head_sym(const Expr* head) {
    if (!head || head->type != EXPR_SYMBOL) return false;
    const char* s = head->data.symbol;
    return s == SYM_Hold
        || s == SYM_HoldForm
        || s == SYM_HoldComplete
        || s == SYM_HoldPattern
        || s == SYM_Unevaluated;
}

/* True iff `e` would be considered a NumericQ — borrowed from core.c's
 * is_numeric_quantity but inlined here since that helper is file-static.
 * Recognises Integer/Real/BigInt, named constants, and NumericFunction
 * heads with all-numeric arguments. */
static bool is_numeric_quantity(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol;
        return name == SYM_Pi || name == SYM_E || name == SYM_I
            || name == SYM_EulerGamma || name == SYM_GoldenRatio
            || name == SYM_Catalan || name == SYM_Degree;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head_name = e->data.function.head->data.symbol;
        if (head_name == SYM_Complex || head_name == SYM_Rational) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!is_numeric_quantity(e->data.function.args[i])) return false;
            }
            return true;
        }
        SymbolDef* def = symtab_get_def(head_name);
        if (def && (def->attributes & ATTR_NUMERICFUNCTION)) {
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                if (!is_numeric_quantity(e->data.function.args[i])) return false;
            }
            return true;
        }
    }
    return false;
}

/* ------------------------------------------------------------------------
 *  Build an Expr from an mpz_t pair (n, d) representing n/d.
 *
 *  Promotes to BigInt when out of int64_t range; otherwise produces a
 *  plain Integer or Rational[n,d] via make_rational. Sign is normalised
 *  so the denominator is positive.
 * ---------------------------------------------------------------------- */
static Expr* expr_from_mpz_int(const mpz_t z) {
    if (mpz_fits_slong_p(z)) {
        long v = mpz_get_si(z);
        /* On 64-bit Unix, long == int64_t; on Windows it's 32. The fits
         * check above guarantees v is representable, so the cast is safe. */
        return expr_new_integer((int64_t)v);
    }
    return expr_new_bigint_from_mpz(z);
}

static Expr* expr_from_mpz_rational(const mpz_t n, const mpz_t d) {
    /* Normalise sign so d > 0. */
    mpz_t nn, dd, g;
    mpz_init_set(nn, n);
    mpz_init_set(dd, d);
    if (mpz_sgn(dd) < 0) {
        mpz_neg(nn, nn);
        mpz_neg(dd, dd);
    }
    mpz_init(g);
    mpz_gcd(g, nn, dd);
    if (mpz_cmp_ui(g, 1) > 0) {
        mpz_divexact(nn, nn, g);
        mpz_divexact(dd, dd, g);
    }
    Expr* result;
    if (mpz_cmp_ui(dd, 1) == 0) {
        result = expr_from_mpz_int(nn);
    } else {
        Expr* n_expr = expr_from_mpz_int(nn);
        Expr* d_expr = expr_from_mpz_int(dd);
        Expr* args[2] = { n_expr, d_expr };
        result = expr_new_function(expr_new_symbol("Rational"), args, 2);
    }
    mpz_clear(nn); mpz_clear(dd); mpz_clear(g);
    return result;
}

/* ------------------------------------------------------------------------
 *  Default-mode rationalisation: convergents with the c/q^2 bound.
 *
 *  Walk the continued fraction of |x|, compute successive convergents
 *  h_n / k_n, and return the first one with |h/k - x| < 10^-4 / k^2.
 *
 *  We stop the CF expansion at a generous depth; for double-precision
 *  inputs nothing useful comes out past ~50 partial quotients. If none
 *  satisfies the bound, return false so the caller can leave x alone.
 * ---------------------------------------------------------------------- */
static bool rationalize_default_double(double x, mpz_t out_n, mpz_t out_d) {
    if (!isfinite(x)) return false;

    const double C = 1.0e-4;
    const int    MAX_TERMS = 64;       /* plenty for double precision */
    const double TINY      = 1.0e-300; /* fractional-part underflow guard */

    int sign = (x < 0.0) ? -1 : 1;
    double ax = fabs(x);

    /* Past a certain denominator the floating-point noise of x dominates
     * the c/q^2 bound — every "convergent" produced from there on is
     * really an artefact of the input's last bits, not a true Diophantine
     * approximation. The natural cutoff is q_max = sqrt(C / eps_x), where
     * eps_x is the relative precision of x. For an IEEE double that's
     * about sqrt(1e-4 / 2.22e-16) ≈ 2.12e6.
     *
     * Without this cap, Rationalize[N[Pi]] would happily report a deep
     * convergent whose "matching" comes entirely from double-precision
     * round-off rather than any genuine smallness of the next CF
     * coefficient — Mathematica's behaviour is to leave N[Pi] unchanged. */
    const double k_max = 2.5e6;

    /* CF state. h0/k0 = current convergent; h1/k1 = previous. */
    mpz_t h0, h1, k0, k1, h, k;
    mpz_init_set_ui(h0, 1); mpz_init_set_ui(h1, 0);
    mpz_init_set_ui(k0, 0); mpz_init_set_ui(k1, 1);
    mpz_init(h); mpz_init(k);

    double xi = ax;
    bool found = false;

    for (int i = 0; i < MAX_TERMS; i++) {
        if (!isfinite(xi)) break;
        double a_double = floor(xi);
        if (a_double < 0 || a_double > 1.0e18) break;  /* sanity */
        unsigned long a = (unsigned long)a_double;

        /* h_n = a * h_{n-1} + h_{n-2}; same for k. */
        mpz_mul_ui(h, h0, a);
        mpz_add(h, h, h1);
        mpz_mul_ui(k, k0, a);
        mpz_add(k, k, k1);

        /* Bound check: |h/k - ax| < C / k^2, only meaningful while
         * k stays below the precision-noise threshold. */
        double k_d = mpz_get_d(k);
        double h_d = mpz_get_d(h);
        if (k_d > k_max) break;
        double err = fabs(h_d / k_d - ax);
        if (k_d > 0 && err < C / (k_d * k_d)) {
            found = true;
            break;
        }

        /* Advance: h1 ← h0, h0 ← h, similarly for k. */
        mpz_set(h1, h0); mpz_set(h0, h);
        mpz_set(k1, k0); mpz_set(k0, k);

        /* Continued-fraction step: xi ← 1 / (xi - a). */
        double frac = xi - a_double;
        if (frac <= TINY) {
            /* Exact CF termination — value is rational at h0/k0. Check
             * the bound on h0/k0. */
            double k0d = mpz_get_d(k0);
            double h0d = mpz_get_d(h0);
            double err0 = fabs(h0d / k0d - ax);
            if (k0d > 0 && err0 < C / (k0d * k0d)) {
                mpz_set(h, h0);
                mpz_set(k, k0);
                found = true;
            }
            break;
        }
        xi = 1.0 / frac;
    }

    if (found) {
        mpz_init_set(out_n, h);
        if (sign < 0) mpz_neg(out_n, out_n);
        mpz_init_set(out_d, k);
    }

    mpz_clear(h0); mpz_clear(h1); mpz_clear(k0); mpz_clear(k1);
    mpz_clear(h);  mpz_clear(k);
    return found;
}

/* ------------------------------------------------------------------------
 *  Tolerance-mode rationalisation: simplest rational in [lo, hi].
 *
 *  Classical recursive descent on the continued fraction of the interval.
 *  At each step we extract floor(lo) and floor(hi):
 *
 *    - If they differ, an integer lies in (lo, hi] — return it.
 *    - If lo is itself an integer and lies in [lo, hi] — return it.
 *    - Otherwise both have the same floor a; recurse on (1/(hi-a), 1/(lo-a))
 *      and combine: result = a + 1/sub.
 *
 *  Implemented iteratively with a stack of integer parts, then unwound in
 *  reverse to assemble the resulting fraction.
 *
 *  The caller passes lo, hi in any order and possibly negative; the helper
 *  normalises to lo < hi and assumes positive operands at the top level.
 * ---------------------------------------------------------------------- */
static void rationalize_interval(double lo, double hi, mpz_t out_n, mpz_t out_d) {
    if (lo > hi) { double t = lo; lo = hi; hi = t; }

    /* Sign: factor it out, work on positive interval. */
    int sign = 1;
    if (hi < 0) {
        sign = -1;
        double t = lo; lo = -hi; hi = -t;
    }

    /* If 0 ∈ [lo, hi], the simplest rational is 0. */
    if (lo <= 0.0) {
        mpz_init_set_ui(out_n, 0);
        mpz_init_set_ui(out_d, 1);
        return;
    }

    /* Stack of CF integer parts. 64 deep is comfortably more than enough
     * for any double-precision interval. */
    unsigned long a_stack[80];
    int           depth = 0;
    const int     MAX_DEPTH = 80;
    const double  HUGE_FLOOR = 1.0e18;  /* sanity ceiling */

    while (depth < MAX_DEPTH) {
        double f_lo = floor(lo);
        double f_hi = floor(hi);
        if (f_lo > HUGE_FLOOR || f_hi > HUGE_FLOOR) break;

        if (f_lo != f_hi || lo == f_lo) {
            /* Smallest integer ≥ lo, but strictly: any integer in [lo, hi].
             * If lo is itself integer, it's the answer; otherwise pick
             * f_lo + 1 (which equals f_hi when they differ). */
            unsigned long a = (lo == f_lo) ? (unsigned long)f_lo
                                            : (unsigned long)f_lo + 1UL;
            mpz_init_set_ui(out_n, a);
            mpz_init_set_ui(out_d, 1);
            /* Unwind the stack: at each level, replace (n, d) by
             * (a_i * n + d, n). */
            for (int i = depth - 1; i >= 0; i--) {
                mpz_t newn;
                mpz_init(newn);
                mpz_mul_ui(newn, out_n, a_stack[i]);
                mpz_add(newn, newn, out_d);
                mpz_swap(out_d, out_n);   /* d ← old n */
                mpz_swap(out_n, newn);    /* n ← a*old_n + old_d */
                mpz_clear(newn);
            }
            if (sign < 0) mpz_neg(out_n, out_n);
            return;
        }

        /* Same floor a — descend. lo' = 1 / (hi - a), hi' = 1 / (lo - a). */
        unsigned long a = (unsigned long)f_lo;
        a_stack[depth++] = a;
        double new_lo = 1.0 / (hi - f_lo);
        double new_hi = 1.0 / (lo - f_lo);
        lo = new_lo;
        hi = new_hi;
        /* If reciprocation overflowed (lo == f_lo) we'd be stuck; bail. */
        if (!isfinite(lo) || !isfinite(hi)) break;
    }

    /* Fallback: depth exceeded or pathology. Return 0/1 — should be
     * unreachable in practice for sane inputs. */
    mpz_init_set_ui(out_n, 0);
    mpz_init_set_ui(out_d, 1);
}

/* ------------------------------------------------------------------------
 *  Zero-tolerance rationalisation: walk convergents in GMP.
 *
 *  The interval-based algorithm above relies on doubles for its lo/hi
 *  bookkeeping, and at very tight tolerances (e.g. dx ≈ ½ ulp) the
 *  intermediate `hi - f_lo` losses precision long before the deep
 *  convergent that Mathematica picks (depth ~15 for π). For dx = 0 we
 *  walk continued-fraction convergents directly with mpz arithmetic,
 *  testing each against an ulp-derived tolerance.
 * ---------------------------------------------------------------------- */
static bool rationalize_zero_double(double x, mpz_t out_n, mpz_t out_d) {
    if (x == 0.0) {
        mpz_init_set_ui(out_n, 0);
        mpz_init_set_ui(out_d, 1);
        return true;
    }
    if (!isfinite(x)) return false;

    /* Tolerance: ½ ulp(x). Anything within this band rounds back to x in
     * IEEE-754, so it's the natural definition of "the same number". */
    double abx = fabs(x);
    double up = nextafter(abx, INFINITY);
    double dn = nextafter(abx, 0.0);
    double tol = 0.25 * (up - dn);            /* ≈ ½ ulp */
    if (tol <= 0.0) tol = abx * DBL_EPSILON * 0.5;

    int sign = (x < 0.0) ? -1 : 1;
    double xi = abx;

    mpz_t h0, h1, k0, k1, h, k;
    mpz_init_set_ui(h0, 1); mpz_init_set_ui(h1, 0);
    mpz_init_set_ui(k0, 0); mpz_init_set_ui(k1, 1);
    mpz_init(h); mpz_init(k);

    bool found = false;
    const int    MAX_TERMS = 80;
    const double TINY      = 1.0e-300;

    for (int i = 0; i < MAX_TERMS; i++) {
        if (!isfinite(xi)) break;
        double a_double = floor(xi);
        if (a_double < 0 || a_double > 1.0e18) break;
        unsigned long a = (unsigned long)a_double;

        mpz_mul_ui(h, h0, a);  mpz_add(h, h, h1);
        mpz_mul_ui(k, k0, a);  mpz_add(k, k, k1);

        double k_d = mpz_get_d(k);
        double h_d = mpz_get_d(h);
        double err = fabs(h_d / k_d - abx);
        if (err <= tol) {
            found = true;
            break;
        }

        mpz_set(h1, h0); mpz_set(h0, h);
        mpz_set(k1, k0); mpz_set(k0, k);

        double frac = xi - a_double;
        if (frac <= TINY) {
            /* Exact rational reached at h0/k0 — accept it. */
            mpz_set(h, h0); mpz_set(k, k0);
            found = true;
            break;
        }
        xi = 1.0 / frac;
    }

    if (found) {
        mpz_init_set(out_n, h);
        if (sign < 0) mpz_neg(out_n, out_n);
        mpz_init_set(out_d, k);
    }

    mpz_clear(h0); mpz_clear(h1); mpz_clear(k0); mpz_clear(k1);
    mpz_clear(h);  mpz_clear(k);
    return found;
}

/* ------------------------------------------------------------------------
 *  Public core: dispatch on mode.
 * ---------------------------------------------------------------------- */
bool internal_rationalize_double(double x, double dx, RationalizeMode mode,
                                 mpz_t out_n, mpz_t out_d) {
    if (!isfinite(x)) return false;

    if (mode == RATIONALIZE_DEFAULT) {
        return rationalize_default_double(x, out_n, out_d);
    }
    if (mode == RATIONALIZE_ZERO) {
        return rationalize_zero_double(x, out_n, out_d);
    }
    if (dx < 0.0) dx = -dx;
    rationalize_interval(x - dx, x + dx, out_n, out_d);
    return true;
}

/* ------------------------------------------------------------------------
 *  Expression-level threading.
 * ---------------------------------------------------------------------- */

/* Convert a numeric leaf to its rationalised form. On success returns a
 * fresh Expr (Integer or Rational[n,d]); on failure returns NULL meaning
 * "leave the leaf as-is". */
static Expr* rationalize_number_leaf(const Expr* e, double dx, RationalizeMode mode) {
    double xv;
    if (!extract_double(e, &xv)) return NULL;

    /* In default mode, exact numbers (Integer, BigInt, Rational) are
     * returned unchanged — only inexact (Real / MPFR) gets rewritten. */
    if (mode == RATIONALIZE_DEFAULT && !is_inexact_leaf(e)) {
        return expr_copy((Expr*)e);
    }

    mpz_t n, d;
    if (!internal_rationalize_double(xv, dx, mode, n, d)) {
        return NULL;  /* default mode failure: leaf passes through */
    }
    Expr* out = expr_from_mpz_rational(n, d);
    mpz_clear(n); mpz_clear(d);
    return out;
}

/* Rationalise a numeric quantity (NumericQ) end-to-end: numericalise to
 * a double then rationalise. Used in tolerance / zero modes when the
 * subexpression is provably numeric (e.g. Sqrt[2], Exp[Sqrt[2]], Pi). */
static Expr* rationalize_numeric_quantity(const Expr* e, double dx,
                                          RationalizeMode mode) {
    Expr* n_expr = numericalize(e, numeric_machine_spec());
    if (!n_expr) return NULL;

    /* If numericalisation produced a Complex with real components,
     * rationalise both parts. */
    Expr *re_part, *im_part;
    if (is_complex(n_expr, &re_part, &im_part)) {
        Expr* nre = rationalize_number_leaf(re_part, dx, mode);
        Expr* nim = rationalize_number_leaf(im_part, dx, mode);
        if (!nre) nre = expr_copy(re_part);
        if (!nim) nim = expr_copy(im_part);
        Expr* r = make_complex(nre, nim);
        expr_free(n_expr);
        /* make_complex returns Complex[re, im]; let the evaluator simplify
         * (it strips the wrapper when im == 0 etc.). */
        return eval_and_free(r);
    }

    /* Plain numeric leaf path. */
    Expr* leaf = rationalize_number_leaf(n_expr, dx, mode);
    if (leaf) {
        expr_free(n_expr);
        return leaf;
    }
    /* numericalize succeeded but rationalize_number_leaf rejected — only
     * happens in default mode if the value is not inexact, which it must
     * be after numericalize. Fall back to original. */
    expr_free(n_expr);
    return expr_copy((Expr*)e);
}

Expr* internal_rationalize_expr(const Expr* e, double dx, RationalizeMode mode) {
    if (!e) return NULL;

    /* Numeric leaves: handle directly. Note `is_inexact_leaf` catches REAL
     * and MPFR; integers / bigints flow through this path too because
     * extract_double will succeed and rationalize_number_leaf will return
     * a copy in default mode and the integer cast in non-default modes. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT
        || e->type == EXPR_REAL
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) {
        Expr* r = rationalize_number_leaf(e, dx, mode);
        return r ? r : expr_copy((Expr*)e);
    }

    /* Symbols (Pi, x, etc.) — pass through in default mode; in
     * tolerance / zero modes, named numeric constants are converted. */
    if (e->type == EXPR_SYMBOL) {
        if (mode != RATIONALIZE_DEFAULT && is_numeric_quantity(e)) {
            return rationalize_numeric_quantity(e, dx, mode);
        }
        return expr_copy((Expr*)e);
    }

    if (e->type == EXPR_STRING) return expr_copy((Expr*)e);

    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Rational[n, d] — already rational; just copy. */
    int64_t rn, rd;
    if (is_rational((Expr*)e, &rn, &rd)) {
        return expr_copy((Expr*)e);
    }

    /* Complex[re, im] — recurse componentwise. This handles the
     * `Rationalize[N[Pi + 33/211 E I], 0]` case directly since that
     * evaluates to Complex[N[Pi], N[33 E / 211]]. */
    Expr *re_arg, *im_arg;
    if (is_complex((Expr*)e, &re_arg, &im_arg)) {
        Expr* nre = internal_rationalize_expr(re_arg, dx, mode);
        Expr* nim = internal_rationalize_expr(im_arg, dx, mode);
        Expr* c = make_complex(nre, nim);
        return eval_and_free(c);
    }

    /* Hold-form heads pass through. */
    if (is_hold_head_sym(e->data.function.head)) {
        return expr_copy((Expr*)e);
    }

    /* In tolerance / zero modes, fully-numeric subexpressions get
     * numericalised and rationalised as a unit. This is what makes
     * Rationalize[Sqrt[2], 0.01] → 99/70 work. */
    if (mode != RATIONALIZE_DEFAULT && is_numeric_quantity(e)) {
        return rationalize_numeric_quantity(e, dx, mode);
    }

    /* General compound: rebuild with each part rationalised, then
     * re-evaluate so Plus/Times can re-canonicalise. */
    const size_t n = e->data.function.arg_count;
    Expr* new_head = internal_rationalize_expr(e->data.function.head, dx, mode);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) {
            expr_free(new_head);
            return expr_copy((Expr*)e);
        }
        for (size_t i = 0; i < n; i++) {
            new_args[i] = internal_rationalize_expr(e->data.function.args[i], dx, mode);
        }
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);
    return eval_and_free(rebuilt);
}

/* ------------------------------------------------------------------------
 *  Builtin entry point.
 * ---------------------------------------------------------------------- */
Expr* builtin_rationalize(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    const size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    Expr* x = res->data.function.args[0];

    /* Ownership: the evaluator frees `res` after a non-NULL return
     * (see eval.c, the builtin-dispatch site). We must not free it
     * ourselves. Returning NULL keeps `res` alive and unevaluated. */
    if (argc == 1) {
        return internal_rationalize_expr(x, 0.0, RATIONALIZE_DEFAULT);
    }

    Expr* dx_expr = res->data.function.args[1];
    double dx;
    if (!extract_double(dx_expr, &dx)) {
        /* Non-numeric tolerance — leave unevaluated. */
        return NULL;
    }
    if (dx < 0.0) {
        /* Negative tolerance is invalid in Mathematica; leave unevaluated. */
        return NULL;
    }

    RationalizeMode mode = (dx == 0.0) ? RATIONALIZE_ZERO : RATIONALIZE_TOLERANCE;
    return internal_rationalize_expr(x, dx, mode);
}

/* ------------------------------------------------------------------------
 *  Helpers exposed to other modules.
 *
 *  The exact-only symbolic functions (PolynomialGCD, PolynomialLCM,
 *  Together, Cancel, Apart, Factor, FactorSquareFree, Simplify, Series,
 *  Limit, …) cannot reason about EXPR_REAL / EXPR_MPFR coefficients —
 *  their algorithms assume rational arithmetic. The standard fix in a
 *  Mathematica-style CAS is "force-rationalize on the way in,
 *  numericalize on the way out": treat 1.5 as 3/2 internally so the
 *  algebraic algorithms work, then re-numericalize the result so the
 *  caller still observes an approximate-in / approximate-out
 *  signature.
 *
 *  These three helpers package that pattern.
 * ---------------------------------------------------------------------- */

bool internal_contains_inexact(const Expr* e) {
    if (!e) return false;
    if (is_inexact_leaf(e)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (internal_contains_inexact(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (internal_contains_inexact(e->data.function.args[i])) return true;
    }
    return false;
}

bool internal_args_contain_inexact(const Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (internal_contains_inexact(res->data.function.args[i])) return true;
    }
    return false;
}

/* Force-rationalise a single inexact leaf. We try the relaxed
 * convergent-bound algorithm first — that's what recovers `1/9.` as `1/9`
 * rather than its bit-exact representation `562949953421312/5066549580791809`.
 * Only when no convergent satisfies the c/q² bound (e.g. for an actual
 * irrational like N[Pi]) do we fall back to the bit-exact ½-ulp algorithm,
 * which always succeeds for finite inputs. The fallback ensures the
 * downstream exact-arithmetic code never sees a stray Real / MPFR. */
static Expr* force_rationalize_leaf(const Expr* e) {
    if (!is_inexact_leaf(e)) return expr_copy((Expr*)e);

    double xv;
    if (!extract_double(e, &xv)) return expr_copy((Expr*)e);

    mpz_t n, d;
    bool ok = internal_rationalize_double(xv, 0.0, RATIONALIZE_DEFAULT, n, d);
    if (!ok) {
        ok = internal_rationalize_double(xv, 0.0, RATIONALIZE_ZERO, n, d);
    }
    if (!ok) return expr_copy((Expr*)e);

    Expr* out = expr_from_mpz_rational(n, d);
    mpz_clear(n); mpz_clear(d);
    return out;
}

Expr* internal_force_rationalize(const Expr* e) {
    if (!e) return NULL;

    /* Inexact numeric leaves get the "lenient → bit-exact" treatment. */
    if (is_inexact_leaf(e)) {
        return force_rationalize_leaf(e);
    }

    /* Exact atoms pass through. Notably, named constants (Pi, E, ...) are
     * NOT rationalised — keeping them symbolic preserves the algebraic
     * power of the downstream pass and matches how a user would expect
     * Cancel[Pi/x + 1.5] to behave. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT
        || e->type == EXPR_SYMBOL || e->type == EXPR_STRING) {
        return expr_copy((Expr*)e);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    int64_t rn, rd;
    if (is_rational((Expr*)e, &rn, &rd)) return expr_copy((Expr*)e);

    /* Complex[re, im] — rationalise components, normalise via make_complex. */
    Expr *re_arg, *im_arg;
    if (is_complex((Expr*)e, &re_arg, &im_arg)) {
        Expr* nre = internal_force_rationalize(re_arg);
        Expr* nim = internal_force_rationalize(im_arg);
        Expr* c = make_complex(nre, nim);
        return eval_and_free(c);
    }

    /* Hold-forms pass through unchanged. */
    if (is_hold_head_sym(e->data.function.head)) {
        return expr_copy((Expr*)e);
    }

    /* General compound: recurse into head and args, rebuild, re-evaluate
     * so Plus/Times/Power can re-canonicalise. */
    const size_t n = e->data.function.arg_count;
    Expr* new_head = internal_force_rationalize(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) {
            expr_free(new_head);
            return expr_copy((Expr*)e);
        }
        for (size_t i = 0; i < n; i++) {
            new_args[i] = internal_force_rationalize(e->data.function.args[i]);
        }
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);
    return eval_and_free(rebuilt);
}

Expr* internal_rationalize_then_numericalize(Expr* res,
                                             Expr* (*core)(Expr*)) {
    if (!res || res->type != EXPR_FUNCTION || !core) return NULL;

    /* Build a clone of `res` with every argument force-rationalised.
     * The head is structurally copied (the symbol is the dispatch key
     * for the inner builtin call). */
    const size_t n = res->data.function.arg_count;
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = internal_force_rationalize(res->data.function.args[i]);
        }
    }
    Expr* clone = expr_new_function(expr_copy(res->data.function.head),
                                    new_args, n);
    free(new_args);

    /* Dispatch to the exact-arithmetic core on the rationalised clone.
     *
     * Ownership note: in this codebase the evaluator (eval.c) is the one
     * that frees a builtin's `res` after a non-NULL return — the builtin
     * itself merely produces a result. Since we are calling `core`
     * directly (bypassing the evaluator), the clone is *ours* either way:
     * we must free it after we extract or discard the result. We must
     * also leave `res` untouched so the outer evaluator can free it. */
    Expr* exact_result = core(clone);
    if (!exact_result) {
        expr_free(clone);
        return NULL;
    }

    /* Round-trip the exact result back to numeric so callers see
     * inexact-in / inexact-out. */
    Expr* numeric = numericalize(exact_result, numeric_machine_spec());
    expr_free(exact_result);
    expr_free(clone);
    return numeric;
}

/* ------------------------------------------------------------------------
 *  Registration.
 * ---------------------------------------------------------------------- */
void rationalize_init(void) {
    symtab_add_builtin("Rationalize", builtin_rationalize);
    symtab_get_def("Rationalize")->attributes |= ATTR_PROTECTED;
}
