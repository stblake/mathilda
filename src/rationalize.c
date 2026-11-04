/* Mathilda — Rationalize implementation.
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
    const char* s = head->data.symbol.name;
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
        const char* name = e->data.symbol.name;
        return name == SYM_Pi || name == SYM_E || name == SYM_I
            || name == SYM_EulerGamma || name == SYM_GoldenRatio
            || name == SYM_Catalan || name == SYM_Degree
            || name == SYM_GoldenAngle || name == SYM_Glaisher
            || name == SYM_Khinchin;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* head_name = e->data.function.head->data.symbol.name;
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
        result = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    }
    mpz_clear(nn); mpz_clear(dd); mpz_clear(g);
    return result;
}

/* ------------------------------------------------------------------------
 *  The exact core.
 *
 *  Every mode works on the EXACT value of the number, in GMP rationals.
 *  A double or an MPFR value is a dyadic rational m 2^e, so nothing is lost
 *  by converting it, and nothing below ever rounds.  The earlier core ran
 *  its continued fractions in doubles, and that broke at both ends of the
 *  range and at every precision above machine:
 *
 *    - Rationalize[1.*^-300, 0] was 0 -- the fractional-part underflow
 *      guard (frac <= 1e-300) accepted the zeroth convergent 0/1 -- so
 *      Rationalize[1.*^300 + 1.*^-300 x, 0] dropped the x term;
 *    - Rationalize[1.*^20, 0] and Rationalize[1.*^300, 0] stayed inexact
 *      (a partial quotient above 1e18 aborted the walk);
 *    - Rationalize[N[Pi, 30], 0] was accurate to 1e-17 only (the MPFR value
 *      was read through a double), and Rationalize[N[Pi, 30], 10^-20] did
 *      not evaluate at all (a tolerance with a bignum denominator was not
 *      read);
 *    - Rationalize[10^30 + 1, 1/2] went through a double as well.
 *
 *  The precision of a leaf is its own: 53 bits for a machine real, the MPFR
 *  precision otherwise.
 * ---------------------------------------------------------------------- */

/* Simplest rational (smallest denominator) in the closed interval
 * [lo, hi], 0 < lo <= hi, by the continued-fraction descent: take the
 * integer part a of lo; stop at lo itself if it is an integer, or at a + 1
 * if that lies in the interval; otherwise recurse on (1/(hi-a), 1/(lo-a)).
 * The convergent recurrence accumulates the answer.  `lo` and `hi` are
 * clobbered. */
static void simplest_rational_between(mpq_t lo, mpq_t hi, mpz_t num, mpz_t den) {
    mpz_t a, h0, h1, k0, k1, t;
    mpz_inits(a, h0, h1, k0, k1, t, NULL);
    mpz_set_ui(h0, 0); mpz_set_ui(h1, 1);       /* h_{-2}, h_{-1} */
    mpz_set_ui(k0, 1); mpz_set_ui(k1, 0);       /* k_{-2}, k_{-1} */
    mpq_t q;
    mpq_init(q);
    for (;;) {
        mpz_fdiv_q(a, mpq_numref(lo), mpq_denref(lo));
        bool done = false;
        if (mpz_cmp_ui(mpq_denref(lo), 1) == 0) {
            done = true;                                   /* lo is an integer */
        } else {
            mpq_set_z(q, a);
            mpz_add_ui(mpq_numref(q), mpq_numref(q), 1);   /* q = a + 1 */
            if (mpq_cmp(q, hi) <= 0) { mpz_add_ui(a, a, 1); done = true; }
        }
        /* h = a h1 + h0, k = a k1 + k0 */
        mpz_mul(t, a, h1); mpz_add(t, t, h0); mpz_set(h0, h1); mpz_set(h1, t);
        mpz_mul(t, a, k1); mpz_add(t, t, k0); mpz_set(k0, k1); mpz_set(k1, t);
        if (done) break;
        mpq_set_z(q, a);
        mpq_sub(lo, lo, q);                        /* lo - a > 0 */
        mpq_sub(hi, hi, q);                        /* hi - a >= lo - a */
        mpq_inv(q, hi);
        mpq_inv(hi, lo);
        mpq_set(lo, q);
    }
    mpz_set(num, h1);
    mpz_set(den, k1);
    mpq_clear(q);
    mpz_clears(a, h0, h1, k0, k1, t, NULL);
}

/* The simplest rational with magnitude in [|x| - wlo, |x| + whi] (wlo, whi
 * >= 0), with the sign of x, into the initialised (num, den).  An interval
 * holding 0 answers 0. */
static void rz_simplest_around2(const mpq_t x, const mpq_t wlo, const mpq_t whi,
                                mpz_t num, mpz_t den) {
    mpq_t lo, hi;
    mpq_inits(lo, hi, NULL);
    mpq_abs(lo, x);
    mpq_set(hi, lo);
    mpq_sub(lo, lo, wlo);
    mpq_add(hi, hi, whi);
    if (mpq_sgn(lo) <= 0) {
        mpz_set_ui(num, 0);
        mpz_set_ui(den, 1);
    } else {
        simplest_rational_between(lo, hi, num, den);
        if (mpq_sgn(x) < 0) mpz_neg(num, num);
    }
    mpq_clears(lo, hi, NULL);
}

static void rz_simplest_around(const mpq_t x, const mpq_t w, mpz_t num, mpz_t den) {
    rz_simplest_around2(x, w, w, num, den);
}

/* log2 of a positive mpz, to double accuracy, at any magnitude. */
static double rz_log2_z(const mpz_t z) {
    long e;
    double d = mpz_get_d_2exp(&e, z);
    return log2(d) + (double)e;
}

/* Rationalize[x] with no tolerance: the first continued-fraction
 * convergent h/k of |x| with |h/k - x| < 10^-4 / k^2, walked exactly.
 * The walk stops once k passes the noise threshold
 * k_max = sqrt(10^-4 / (eps |x|)), eps = 2^(1 - prec): past it a
 * convergent's "match" is an artefact of the last bits of x, not a genuine
 * Diophantine approximation (so N[Pi] is left alone, as in the reference).
 * The integer convergent (k = 1) is always eligible: Rationalize[1.*^12] is
 * 1000000000000 even though k_max < 1 there.  A zero convergent never
 * stands for a nonzero x (Rationalize[0.000001] is not 0).  Returns false
 * when no convergent qualifies. */
static bool rz_default(const mpq_t x, long prec, mpz_t num, mpz_t den) {
    if (mpq_sgn(x) == 0) { mpz_set_ui(num, 0); mpz_set_ui(den, 1); return true; }
    mpz_t A, B, a, r, h0, h1, k0, k1, t, u;
    mpz_inits(A, B, a, r, h0, h1, k0, k1, t, u, NULL);
    mpz_abs(A, mpq_numref(x));
    mpz_set(B, mpq_denref(x));
    mpz_t XA, XB;                                /* |x| = XA / XB, kept */
    mpz_init_set(XA, A);
    mpz_init_set(XB, B);
    const double log2_kmax = 0.5 * (log2(1.0e-4) - (double)(1 - prec)
                                    - (rz_log2_z(XA) - rz_log2_z(XB)));
    mpz_set_ui(h0, 0); mpz_set_ui(h1, 1);        /* h_{-2}, h_{-1} */
    mpz_set_ui(k0, 1); mpz_set_ui(k1, 0);        /* k_{-2}, k_{-1} */
    bool found = false;
    while (mpz_sgn(B) != 0) {
        mpz_fdiv_qr(a, r, A, B);                 /* A/B = a + r/B */
        mpz_mul(t, a, h1); mpz_add(t, t, h0); mpz_set(h0, h1); mpz_set(h1, t);
        mpz_mul(t, a, k1); mpz_add(t, t, k0); mpz_set(k0, k1); mpz_set(k1, t);
        if (mpz_cmp_ui(k1, 1) > 0 && rz_log2_z(k1) > log2_kmax) break;
        if (mpz_sgn(h1) > 0) {
            /* |h/k - XA/XB| < 1/(10^4 k^2)  <=>  |h XB - k XA| k 10^4 < XB */
            mpz_mul(t, h1, XB);
            mpz_mul(u, k1, XA);
            mpz_sub(t, t, u);
            mpz_abs(t, t);
            mpz_mul(t, t, k1);
            mpz_mul_ui(t, t, 10000);
            if (mpz_cmp(t, XB) < 0) { found = true; break; }
        }
        mpz_set(A, B);
        mpz_set(B, r);
    }
    if (found) {
        mpz_set(num, h1);
        if (mpq_sgn(x) < 0) mpz_neg(num, num);
        mpz_set(den, k1);
    }
    mpz_clears(A, B, a, r, h0, h1, k0, k1, t, u, XA, XB, NULL);
    return found;
}

/* One number, exactly: its value `x`, and for an inexact one its precision
 * in bits and the exponents of the gaps to its neighbours: the next larger
 * magnitude is 2^ulp_e away, the next smaller 2^ulp_lo_e (half as far when
 * |x| is a power of two, except among the subnormals, whose spacing is
 * fixed).  Exact numbers (Integer, BigInt, Rational with any size of parts)
 * report prec = 0.  False for anything that is not a finite number. */
static bool rz_number2(const Expr* e, mpq_t x, long* prec, long* ulp_e, long* ulp_lo_e) {
    long dummy;
    if (!ulp_lo_e) ulp_lo_e = &dummy;
    *prec = 0;
    *ulp_e = 0;
    *ulp_lo_e = 0;
    switch (e->type) {
        case EXPR_INTEGER: mpq_set_si(x, (long)e->data.integer, 1UL); return true;
        case EXPR_BIGINT:  mpq_set_z(x, e->data.bigint);              return true;
        case EXPR_REAL: {
            double v = e->data.real;
            if (!isfinite(v)) return false;
            mpq_set_d(x, v);                       /* exact */
            int E = 0;
            if (v != 0.0) (void)frexp(v, &E);      /* |v| = f 2^E, f in [1/2, 1) */
            const long umin = DBL_MIN_EXP - DBL_MANT_DIG;   /* subnormal spacing */
            long ue = (long)E - DBL_MANT_DIG;
            if (ue < umin) ue = umin;
            long le = ue;
            if (v != 0.0 && fabs(v) == ldexp(1.0, E - 1) && ue - 1 >= umin) le = ue - 1;
            *prec = DBL_MANT_DIG;
            *ulp_e = ue;
            *ulp_lo_e = le;
            return true;
        }
#ifdef USE_MPFR
        case EXPR_MPFR: {
            mpfr_srcptr v = e->data.mpfr;
            if (!mpfr_number_p(v)) return false;
            *prec = (long)mpfr_get_prec(v);
            if (mpfr_zero_p(v)) { mpq_set_ui(x, 0, 1); return true; }
            mpz_t m;
            mpz_init(m);
            mpfr_exp_t ex = mpfr_get_z_2exp(m, v);   /* v = m 2^ex, m has prec bits */
            mpz_abs(m, m);
            *ulp_lo_e = (mpz_popcount(m) == 1) ? (long)ex - 1 : (long)ex;
            if (mpfr_sgn(v) < 0) mpz_neg(m, m);
            mpq_set_z(x, m);
            if (ex >= 0) mpq_mul_2exp(x, x, (mp_bitcnt_t)ex);
            else         mpq_div_2exp(x, x, (mp_bitcnt_t)(-ex));
            mpz_clear(m);
            *ulp_e = (long)ex;
            return true;
        }
#endif
        case EXPR_FUNCTION:
            if (is_rational_like(e)) {
                mpz_t n, d;
                mpz_inits(n, d, NULL);
                expr_to_mpz(e->data.function.args[0], n);
                expr_to_mpz(e->data.function.args[1], d);
                bool ok = mpz_sgn(d) != 0;
                if (ok) {
                    mpz_set(mpq_numref(x), n);
                    mpz_set(mpq_denref(x), d);
                    mpq_canonicalize(x);
                }
                mpz_clears(n, d, NULL);
                return ok;
            }
            return false;
        default:
            return false;
    }
}

static bool rz_number(const Expr* e, mpq_t x, long* prec, long* ulp_e) {
    return rz_number2(e, x, prec, ulp_e, NULL);
}

/* The rationalised value of the number `x` (precision `prec`, ulp 2^ulp_e;
 * prec == 0 for an exact one) into the initialised (num, den), in `mode`
 * with tolerance `dx` (TOLERANCE only).  False only in DEFAULT mode, when
 * no convergent qualifies and the caller keeps the leaf. */
static bool rz_core(const mpq_t x, long prec, long ulp_e, long ulp_lo_e,
                    RationalizeMode mode, const mpq_t dx, mpz_t num, mpz_t den) {
    if (mode == RATIONALIZE_TOLERANCE) {
        rz_simplest_around(x, dx, num, den);
        return true;
    }
    /* An exact number is already rational; so is an inexact one whose ulp
     * is at least 1 (every such value is an integer, held exactly). */
    if (prec == 0 || ulp_e >= 0 || mpq_sgn(x) == 0) {
        mpz_set(num, mpq_numref(x));
        mpz_set(den, mpq_denref(x));
        return true;
    }
    if (mode == RATIONALIZE_DEFAULT) return rz_default(x, prec, num, den);
    /* ZERO: the simplest rational that rounds back to x at its own precision
     * (so it is "the same number"): strictly within half the gap to each
     * neighbour.  Open at both ends, since a value exactly halfway may round
     * to the neighbour (ties to even: 1.5 2^-1074 is 2^-1073), and the gap
     * below a power of two is half the gap above.  Shrinking each width by
     * 2^-64 of itself makes the closed search interval an open one. */
    mpq_t wlo, whi, t;
    mpq_inits(wlo, whi, t, NULL);
    mpq_set_ui(whi, 1, 1);
    mpq_div_2exp(whi, whi, (mp_bitcnt_t)(1 - ulp_e));      /* 2^(ulp_e - 1) */
    mpq_set_ui(wlo, 1, 1);
    mpq_div_2exp(wlo, wlo, (mp_bitcnt_t)(1 - ulp_lo_e));   /* 2^(ulp_lo_e - 1) */
    mpq_div_2exp(t, whi, 64); mpq_sub(whi, whi, t);
    mpq_div_2exp(t, wlo, 64); mpq_sub(wlo, wlo, t);
    rz_simplest_around2(x, wlo, whi, num, den);
    mpq_clears(wlo, whi, t, NULL);
    return true;
}

/* ------------------------------------------------------------------------
 *  Public core: one double.
 * ---------------------------------------------------------------------- */
bool internal_rationalize_double(double x, double dx, RationalizeMode mode,
                                 mpz_t out_n, mpz_t out_d) {
    if (!isfinite(x)) return false;
    if (mode == RATIONALIZE_TOLERANCE && !isfinite(dx)) return false;
    Expr leaf;
    memset(&leaf, 0, sizeof leaf);
    leaf.type = EXPR_REAL;
    leaf.data.real = x;
    mpq_t q, w;
    mpq_inits(q, w, NULL);
    long prec, ulp_e, ulp_lo_e;
    (void)rz_number2(&leaf, q, &prec, &ulp_e, &ulp_lo_e);
    mpq_set_d(w, fabs(dx));
    mpz_init(out_n);
    mpz_init(out_d);
    bool ok = rz_core(q, prec, ulp_e, ulp_lo_e, mode, w, out_n, out_d);
    if (!ok) { mpz_clear(out_n); mpz_clear(out_d); }
    mpq_clears(q, w, NULL);
    return ok;
}

/* ------------------------------------------------------------------------
 *  Expression-level threading.
 * ---------------------------------------------------------------------- */

/* Convert a numeric leaf to its rationalised form. On success returns a
 * fresh Expr (Integer or Rational[n,d]); on failure returns NULL meaning
 * "leave the leaf as-is".  In DEFAULT mode an exact number is returned
 * unchanged; in ZERO mode too (it is its own exact value). */
static Expr* rationalize_number_leaf(const Expr* e, const mpq_t dx, RationalizeMode mode) {
    if (mode != RATIONALIZE_TOLERANCE && !is_inexact_leaf(e)) {
        return expr_copy((Expr*)e);
    }
    mpq_t x;
    mpq_init(x);
    long prec, ulp_e, ulp_lo_e;
    if (!rz_number2(e, x, &prec, &ulp_e, &ulp_lo_e)) { mpq_clear(x); return NULL; }
    mpz_t n, d;
    mpz_inits(n, d, NULL);
    Expr* out = NULL;
    if (rz_core(x, prec, ulp_e, ulp_lo_e, mode, dx, n, d)) out = expr_from_mpz_rational(n, d);
    mpz_clears(n, d, NULL);
    mpq_clear(x);
    return out;
}

/* Rationalise a numeric quantity (NumericQ) end-to-end: numericalise it,
 * then rationalise the number.  Used in tolerance / zero modes when the
 * subexpression is provably numeric (e.g. Sqrt[2], Exp[Sqrt[2]], Pi).
 * Machine precision unless the tolerance is finer than a double resolves:
 * Rationalize[Sqrt[2], 10^-30] numericalises at enough bits that the
 * answer is within 10^-30 of Sqrt[2] itself, not of N[Sqrt[2]]. */
static Expr* rationalize_numeric_quantity(const Expr* e, const mpq_t dx,
                                          RationalizeMode mode) {
    Expr* n_expr = numericalize(e, numeric_machine_spec());
    if (!n_expr) return NULL;
#ifdef USE_MPFR
    if (mode == RATIONALIZE_TOLERANCE && mpq_sgn(dx) > 0) {
        double mag = 1.0;                          /* |x| at machine precision */
        Expr *re0, *im0;
        double v;
        if (is_complex(n_expr, &re0, &im0)) {
            double a = 0, b = 0;
            if (extract_double(re0, &a) && extract_double(im0, &b)) mag = fmax(fabs(a), fabs(b));
        } else if (extract_double(n_expr, &v)) {
            mag = fabs(v);
        }
        if (!(mag >= 1.0) || !isfinite(mag)) mag = 1.0;
        double log2_dx = rz_log2_z(mpq_numref(dx)) - rz_log2_z(mpq_denref(dx));
        double need = log2(mag) - log2_dx + 16.0;  /* bits to resolve dx, with margin */
        if (need > (double)DBL_MANT_DIG && need < 1.0e7) {
            NumericSpec spec;
            spec.mode = NUMERIC_MODE_MPFR;
            spec.bits = (long)ceil(need);
            Expr* hp = numericalize(e, spec);
            if (hp) { expr_free(n_expr); n_expr = hp; }
        }
    }
#endif

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
    /* Not a finite number (an overflowed or Indeterminate value): keep
     * the original. */
    expr_free(n_expr);
    return expr_copy((Expr*)e);
}

static Expr* rationalize_walk(const Expr* e, const mpq_t dx, RationalizeMode mode) {
    if (!e) return NULL;

    /* Numeric leaves: handle directly. */
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

    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    /* Rational[n, d] (either part a bignum) — a number in its own right:
     * unchanged unless a tolerance asks for a simpler one. */
    if (is_rational_like(e)) {
        Expr* r = rationalize_number_leaf(e, dx, mode);
        return r ? r : expr_copy((Expr*)e);
    }

    /* Complex[re, im] — recurse componentwise. This handles the
     * `Rationalize[N[Pi + 33/211 E I], 0]` case directly since that
     * evaluates to Complex[N[Pi], N[33 E / 211]]. */
    Expr *re_arg, *im_arg;
    if (is_complex((Expr*)e, &re_arg, &im_arg)) {
        Expr* nre = rationalize_walk(re_arg, dx, mode);
        Expr* nim = rationalize_walk(im_arg, dx, mode);
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
    Expr* new_head = rationalize_walk(e->data.function.head, dx, mode);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) {
            expr_free(new_head);
            return expr_copy((Expr*)e);
        }
        for (size_t i = 0; i < n; i++) {
            new_args[i] = rationalize_walk(e->data.function.args[i], dx, mode);
        }
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);
    return eval_and_free(rebuilt);
}

Expr* internal_rationalize_expr(const Expr* e, double dx, RationalizeMode mode) {
    if (!e) return NULL;
    mpq_t w;
    mpq_init(w);
    if (isfinite(dx)) mpq_set_d(w, fabs(dx));
    Expr* r = rationalize_walk(e, w, mode);
    mpq_clear(w);
    return r;
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

    /* The tolerance, exactly: 10^-20 (a bignum denominator) or N[10^-40, 50]
     * are read as given, not through a double. */
    mpq_t dx;
    mpq_init(dx);
    long prec, ulp_e;
    if (!rz_number(res->data.function.args[1], dx, &prec, &ulp_e) || mpq_sgn(dx) < 0) {
        /* Non-numeric or negative tolerance — leave unevaluated. */
        mpq_clear(dx);
        return NULL;
    }
    RationalizeMode mode = (mpq_sgn(dx) == 0) ? RATIONALIZE_ZERO : RATIONALIZE_TOLERANCE;
    Expr* r = rationalize_walk(x, dx, mode);
    mpq_clear(dx);
    return r;
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

/* Force-rationalise a single inexact leaf.
 *
 * Common to all modes: we always try `RATIONALIZE_DEFAULT` first (the
 * relaxed convergent-bound algorithm) so a clean rational like `1.5`
 * recovers as `3/2` rather than its bit-exact double form, and `1/9.`
 * recovers as `1/9`.  The mode parameter only affects the *fallback*
 * for cases where no convergent meets the c/q² bound (typically honest
 * transcendentals like N[Pi]):
 *
 *   bits == 0           -- legacy default.  Fall back to RATIONALIZE_ZERO
 *                          (½-ulp at the input's double precision).  The
 *                          resulting rational round-trips exactly through
 *                          double.
 *
 *   1 ≤ bits < 52       -- sub-machine precision (only seen when the
 *                          caller has an MPFR leaf with explicit low
 *                          precision).  Fall back to RATIONALIZE_TOLERANCE
 *                          with dx = 2^(-bits): the smallest-denominator
 *                          rational in [x - 2^(-bits), x + 2^(-bits)].
 *                          This is no tighter than the requested precision,
 *                          giving short denominators when the caller
 *                          accepts a coarse approximation.
 *
 *   bits ≥ 52           -- at or above machine precision.  TOLERANCE with
 *                          dx ≤ ulp(x) loses precision in the interval
 *                          arithmetic, so we use
 *                          the same RATIONALIZE_ZERO fallback as the
 *                          legacy path.  The output is bit-exact at
 *                          double precision; if the caller wants more
 *                          precision than that, they will recover it
 *                          when numericalising the eventual result at
 *                          MPFR. */
#ifdef USE_MPFR
/* An MPFR leaf made exact at `bits` of precision (its own, if lower): the
 * simplest rational within a quarter unit in the last place, so it rounds
 * back to the same value, and short when the value is (2.`30 -> 2,
 * 0.1`30 -> 1/10).  Going through a double, as the machine path does, threw
 * away every digit past the 16th -- Factor[N[Pi, 50] x + 1] came back as
 * 245850922/78256779 (issue #124). */
static Expr* force_rationalize_mpfr(const Expr* e, long bits) {
    mpfr_srcptr x = e->data.mpfr;
    if (!mpfr_number_p(x)) return expr_copy((Expr*)e);
    if (mpfr_zero_p(x)) return expr_new_integer(0);
    long p = (long)mpfr_get_prec(x);
    long s = (bits > 0 && bits < p) ? p - bits : 0;   /* coarser than its own */
    mpz_t m, num, den;
    mpz_inits(m, num, den, NULL);
    mpfr_exp_t ex = mpfr_get_z_2exp(m, x);           /* x = m 2^ex exactly */
    int sign = mpz_sgn(m);
    mpz_abs(m, m);
    /* [lo, hi] = (4m -+ 2^s) 2^(ex-2) */
    mpq_t lo, hi;
    mpq_inits(lo, hi, NULL);
    mpz_mul_2exp(m, m, 2);
    mpz_t w;
    mpz_init_set_ui(w, 1);
    mpz_mul_2exp(w, w, (mp_bitcnt_t)s);
    mpz_sub(mpq_numref(lo), m, w);
    mpz_add(mpq_numref(hi), m, w);
    mpz_clear(w);
    long sh = (long)ex - 2;
    mpz_set_ui(mpq_denref(lo), 1);
    mpz_set_ui(mpq_denref(hi), 1);
    if (sh >= 0) {
        mpz_mul_2exp(mpq_numref(lo), mpq_numref(lo), (mp_bitcnt_t)sh);
        mpz_mul_2exp(mpq_numref(hi), mpq_numref(hi), (mp_bitcnt_t)sh);
    } else {
        mpz_mul_2exp(mpq_denref(lo), mpq_denref(lo), (mp_bitcnt_t)(-sh));
        mpz_mul_2exp(mpq_denref(hi), mpq_denref(hi), (mp_bitcnt_t)(-sh));
    }
    mpq_canonicalize(lo);
    mpq_canonicalize(hi);
    Expr* out;
    if (mpq_sgn(lo) <= 0) {
        out = expr_new_integer(0);                   /* the interval holds 0 */
    } else {
        simplest_rational_between(lo, hi, num, den);
        if (sign < 0) mpz_neg(num, num);
        out = expr_from_mpz_rational(num, den);
    }
    mpq_clears(lo, hi, NULL);
    mpz_clears(m, num, den, NULL);
    return out;
}
#endif

static Expr* force_rationalize_leaf(const Expr* e, long bits) {
    if (!is_inexact_leaf(e)) return expr_copy((Expr*)e);
#ifdef USE_MPFR
    /* An arbitrary-precision leaf asked for more than a double holds. */
    if (e->type == EXPR_MPFR && bits > 53) return force_rationalize_mpfr(e, bits);
#endif

    double xv;
    if (!extract_double(e, &xv)) return expr_copy((Expr*)e);

    mpz_t n, d;
    bool ok = internal_rationalize_double(xv, 0.0, RATIONALIZE_DEFAULT, n, d);
    if (!ok) {
        if (bits > 0 && bits < 52) {
            double dx = ldexp(1.0, -bits);
            ok = internal_rationalize_double(xv, dx, RATIONALIZE_TOLERANCE, n, d);
        } else {
            ok = internal_rationalize_double(xv, 0.0, RATIONALIZE_ZERO, n, d);
        }
    }
    if (!ok && isfinite(xv) && fabs(xv) >= 9007199254740992.0) {
        /* |x| >= 2^53: every such double is an integer, so it converts
         * exactly.  Returning the Real unchanged sent the round-trip
         * callers (Factor, Together, ...) back into themselves forever. */
        mpz_init(n); mpz_init_set_ui(d, 1);
        mpz_set_d(n, xv);
        ok = true;
    }
    if (!ok) return expr_copy((Expr*)e);

    Expr* out = expr_from_mpz_rational(n, d);
    mpz_clear(n); mpz_clear(d);
    return out;
}

Expr* internal_force_rationalize_bits(const Expr* e, long bits) {
    if (!e) return NULL;

    /* Inexact numeric leaves get the "lenient → bit-exact" treatment
     * (bits == 0) or the precision-respecting interval search (bits > 0). */
    if (is_inexact_leaf(e)) {
        return force_rationalize_leaf(e, bits);
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
        Expr* nre = internal_force_rationalize_bits(re_arg, bits);
        Expr* nim = internal_force_rationalize_bits(im_arg, bits);
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
    Expr* new_head = internal_force_rationalize_bits(e->data.function.head, bits);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) {
            expr_free(new_head);
            return expr_copy((Expr*)e);
        }
        for (size_t i = 0; i < n; i++) {
            new_args[i] = internal_force_rationalize_bits(e->data.function.args[i], bits);
        }
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);
    return eval_and_free(rebuilt);
}

/* SetPrecision[x, Infinity] / SetAccuracy[x, Infinity]: every approximate
 * number in e replaced by its exact value (a machine or MPFR number IS a
 * dyadic rational m 2^e, so this rounds nothing): 0.1 ->
 * 3602879701896397/36028797018963968.  Exact numbers and symbolic constants
 * (Pi stays Pi) are untouched, held forms are not entered (as SetPrecision's
 * other forms), and the rebuilt expression is re-evaluated so Plus / Times /
 * Complex re-canonicalise. */
Expr* internal_exact_value(const Expr* e) {
    if (!e) return NULL;
    if (is_inexact_leaf(e)) {
        mpq_t x;
        mpq_init(x);
        long prec, ulp_e;
        Expr* out;
        if (rz_number(e, x, &prec, &ulp_e)) out = expr_from_mpz_rational(mpq_numref(x), mpq_denref(x));
        else out = expr_copy((Expr*)e);            /* not finite: unchanged */
        mpq_clear(x);
        return out;
    }
    if (e->type != EXPR_FUNCTION || is_rational_like(e) ||
        is_hold_head_sym(e->data.function.head))
        return expr_copy((Expr*)e);
    const size_t n = e->data.function.arg_count;
    Expr* new_head = internal_exact_value(e->data.function.head);
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) { expr_free(new_head); return expr_copy((Expr*)e); }
        for (size_t i = 0; i < n; i++) new_args[i] = internal_exact_value(e->data.function.args[i]);
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);
    return eval_and_free(rebuilt);
}

Expr* internal_force_rationalize(const Expr* e) {
    return internal_force_rationalize_bits(e, 0);
}

/* ------------------------------------------------------------------------
 *  Opaque kernels in the round trip.
 *
 *  The final numericalise of the round trip must turn the *coefficients* of
 *  the exact answer back into machine numbers -- and nothing else.  A kernel
 *  the input carried verbatim, c[1] or x[-2 + g], is an indeterminate to the
 *  exact core, and numericalising it gave c[1.0] and x[-1.0]: a different
 *  symbol, so Factor[0.4 a + c[1] + d] silently changed what it was
 *  factoring.  So each such kernel K of the input is recorded with its
 *  rationalised form K' (the form the core saw and hands back); in the exact
 *  result every K' is parked behind a placeholder symbol across the
 *  numericalise, and then replaced by the original K.  This also gives back
 *  f[0.5] where the input had it, rather than the f[1/2] the core saw.
 *
 *  A kernel is a compound with a symbolic head other than the arithmetic
 *  and structural heads the walk descends through, and which is not itself a
 *  numeric quantity: Sin[2] or Sqrt[2] still numericalise to 0.909..,
 *  1.414.., as the coefficients they are.
 * ---------------------------------------------------------------------- */
typedef struct {
    Expr** orig;    /* the input's kernel, verbatim                       */
    Expr** exact;   /* its rationalised form, as the core sees it         */
    Expr** ph;      /* placeholder symbol parked in its place             */
    size_t n, cap;
} RnKernels;

/* Heads the kernel walk descends through rather than treats as kernels. */
static bool rn_is_structural_head(const char* h) {
    return h == SYM_Plus || h == SYM_Times || h == SYM_Power || h == SYM_List
        || h == SYM_Complex || h == SYM_Rational
        || h == SYM_Rule || h == SYM_RuleDelayed
        || h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
        || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual
        || h == SYM_And || h == SYM_Or || h == SYM_Not;
}

/* Does `e` carry an inexact number as a *coefficient* -- outside every
 * kernel?  The walk descends the same structural heads as the kernel
 * collector; a numeric quantity (Sin[2.], 1.5) is a coefficient and is
 * scanned whole, while any other compound (f[1.], x[0.5]) is an
 * indeterminate whose inner numbers the exact core never computes with. */
static bool rn_coefficient_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type != EXPR_FUNCTION) return internal_contains_inexact(e);
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || is_hold_head_sym(h))
        return internal_contains_inexact(e);
    if (rn_is_structural_head(h->data.symbol.name)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (rn_coefficient_inexact(e->data.function.args[i])) return true;
        return false;
    }
    return is_numeric_quantity(e) && internal_contains_inexact(e);
}

bool internal_args_contain_inexact_coefficient(const Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        if (rn_coefficient_inexact(res->data.function.args[i])) return true;
    return false;
}

static void rn_collect_kernels(const Expr* e, RnKernels* k, long bits) {
    if (!e || e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || is_hold_head_sym(h)) return;
    if (rn_is_structural_head(h->data.symbol.name)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            rn_collect_kernels(e->data.function.args[i], k, bits);
        return;
    }
    if (is_numeric_quantity(e)) return;
    for (size_t i = 0; i < k->n; i++)
        if (expr_eq((Expr*)e, k->orig[i])) return;
    Expr* ex = internal_force_rationalize_bits(e, bits);
    /* Two distinct inputs that rationalise alike (f[0.5], f[1/2]) cannot
     * both be restored from one exact form: drop the ambiguous key. */
    for (size_t i = 0; i < k->n; i++) {
        if (expr_eq(ex, k->exact[i])) {
            expr_free(ex);
            expr_free(k->orig[i]); expr_free(k->exact[i]);
            k->orig[i] = k->exact[i] = NULL;
            return;
        }
    }
    if (k->n == k->cap) {
        k->cap = k->cap ? 2 * k->cap : 8;
        k->orig  = realloc(k->orig,  sizeof(Expr*) * k->cap);
        k->exact = realloc(k->exact, sizeof(Expr*) * k->cap);
        k->ph    = realloc(k->ph,    sizeof(Expr*) * k->cap);
    }
    char buf[48];
    snprintf(buf, sizeof buf, "RationalizeKernel$%zu", k->n);
    k->orig[k->n]  = expr_copy((Expr*)e);
    k->exact[k->n] = ex;
    k->ph[k->n]    = expr_new_symbol(buf);
    k->n++;
}

static void rn_kernels_free(RnKernels* k) {
    for (size_t i = 0; i < k->n; i++) {
        expr_free(k->orig[i]); expr_free(k->exact[i]); expr_free(k->ph[i]);
    }
    free(k->orig); free(k->exact); free(k->ph);
}

/* Does `e` mention one of the placeholder names?  Then parking would be
 * ambiguous, and the caller keeps the plain round trip. */
static bool rn_mentions_placeholder(const Expr* e, const RnKernels* k) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < k->n; i++)
            if (e->data.symbol.name == k->ph[i]->data.symbol.name) return true;
        return false;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (rn_mentions_placeholder(e->data.function.head, k)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (rn_mentions_placeholder(e->data.function.args[i], k)) return true;
    return false;
}

/* Replace (outermost first) each `from[i]` occurring in `e` by `to[i]`.
 * Structural only; the caller evaluates. */
static Expr* rn_swap(const Expr* e, Expr** from, Expr** to, size_t n, bool* hit) {
    for (size_t i = 0; i < n; i++) {
        if (from[i] && expr_eq((Expr*)e, from[i])) { *hit = true; return expr_copy(to[i]); }
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t m = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++)
        args[i] = rn_swap(e->data.function.args[i], from, to, n, hit);
    Expr* r = expr_new_function(rn_swap(e->data.function.head, from, to, n, hit), args, m);
    free(args);
    return r;
}

/* Numericalise an exact round-trip result, keeping the input's kernels as
 * they were written (see the block comment above).  `bits` is the round
 * trip's working precision: 0 for machine, else the MPFR precision. */
static Expr* rn_numericalize_keeping_kernels(const Expr* exact_result, const Expr* res,
                                             long bits) {
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (bits > 0) { spec.mode = NUMERIC_MODE_MPFR; spec.bits = bits; }
#endif
    RnKernels k = {0};
    for (size_t i = 0; i < res->data.function.arg_count; i++)
        rn_collect_kernels(res->data.function.args[i], &k, bits);

    Expr* out = NULL;
    if (k.n > 0 && !rn_mentions_placeholder(exact_result, &k)) {
        bool hit = false;
        Expr* parked = rn_swap(exact_result, k.exact, k.ph, k.n, &hit);
        if (hit) {
            Expr* num = numericalize(parked, spec);
            bool back = false;
            out = eval_and_free(rn_swap(num, k.ph, k.orig, k.n, &back));
            expr_free(num);
        }
        expr_free(parked);
    }
    rn_kernels_free(&k);
    return out ? out : numericalize((Expr*)exact_result, spec);
}

/* Smallest precision, in bits, over the inexact leaves of `e`'s arguments
 * when EVERY one is arbitrary-precision; 0 when any is a machine number
 * (or there is none).  Machine precision is the floor of the contagion, so
 * one machine leaf makes the whole round trip a machine one (issue #124:
 * Factor[N[Pi, 50] x + 1] came back at machine precision). */
static long rn_round_trip_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_REAL) return -1;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
#endif
    if (e->type != EXPR_FUNCTION) return 0;
    long best = rn_round_trip_bits(e->data.function.head);
    if (best < 0) return -1;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        long b = rn_round_trip_bits(e->data.function.args[i]);
        if (b < 0) return -1;
        if (b > 0 && (best == 0 || b < best)) best = b;
    }
    return best;
}

Expr* internal_rationalize_then_numericalize_post(Expr* res, Expr* (*core)(Expr*),
                                                  Expr* (*post)(Expr*)) {
    if (!res || res->type != EXPR_FUNCTION || !core) return NULL;

    /* Working precision of the round trip: the input's own, when it is
     * entirely arbitrary-precision, else machine (bits == 0, the legacy
     * convergent-bound rationalisation). */
    long bits = 0;
    for (size_t i = 0; i < res->data.function.arg_count && bits >= 0; i++) {
        long b = rn_round_trip_bits(res->data.function.args[i]);
        if (b < 0) bits = -1;
        else if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    if (bits <= 53) bits = 0;

    /* Build a clone of `res` with every argument force-rationalised.
     * The head is structurally copied (the symbol is the dispatch key
     * for the inner builtin call). */
    const size_t n = res->data.function.arg_count;
    Expr** new_args = NULL;
    if (n > 0) {
        new_args = (Expr**)malloc(sizeof(Expr*) * n);
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = internal_force_rationalize_bits(res->data.function.args[i], bits);
        }
    }
    Expr* clone = expr_new_function(expr_copy(res->data.function.head),
                                    new_args, n);
    free(new_args);

    /* A leaf that could not be made exact would bring `core` straight back
     * here with the same input: decline instead of recursing without end. */
    if (internal_args_contain_inexact(clone)) {
        expr_free(clone);
        return NULL;
    }

    /* Dispatch to the exact-arithmetic core on the rationalised clone.
     *
     * Ownership note: in this codebase the evaluator (eval.c) is the one
     * that frees a builtin's `res` after a non-NULL return — the builtin
     * itself merely produces a result. Since we are calling `core`
     * directly (bypassing the evaluator), the clone is *ours* either way:
     * we must free it after we extract or discard the result. We must
     * also leave `res` untouched so the outer evaluator can free it. */
    Expr* exact_result = core(clone);
    expr_free(clone);
    if (!exact_result) return NULL;
    if (post) exact_result = post(exact_result);

    /* Round-trip the exact result back to numeric so callers see
     * inexact-in / inexact-out. */
    Expr* numeric = rn_numericalize_keeping_kernels(exact_result, res, bits);
    expr_free(exact_result);
    return numeric;
}

Expr* internal_rationalize_then_numericalize(Expr* res,
                                             Expr* (*core)(Expr*)) {
    return internal_rationalize_then_numericalize_post(res, core, NULL);
}

/* ------------------------------------------------------------------------
 *  Registration.
 * ---------------------------------------------------------------------- */
void rationalize_init(void) {
    symtab_add_builtin("Rationalize", builtin_rationalize);
    symtab_get_def("Rationalize")->attributes |= ATTR_PROTECTED;
}
