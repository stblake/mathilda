/*
 * Mathilda Comparisons
 *
 * Implements structural comparison functions similar to the Wolfram Language.
 */

#include "ndarray.h"
#include "comparisons.h"
#include "symtab.h"
#include "sym_names.h"
#include "attr.h"
#include "expr.h"
#include "arithmetic.h"
#include "eval.h"
#include "numeric.h"
#include "internal.h"
#include "interval.h"
#include "zero_test.h"
#include "assoc_struct.h"   /* Equal/Unequal between associations */
#include "equal_tol.h"      /* Internal`$EqualTolerance / $SameQTolerance */
#include "pack.h"           /* pack_unpack: SameQ on a packed list */
#include <stdio.h>
#include <gmp.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

/* The smallest positive double (2^-1074); DBL_TRUE_MIN is C11. */
#define DBL_TRUE_MIN_OR_DENORM 4.9406564584124654e-324

static bool get_numeric_value(Expr* e, double* val, int64_t* num, int64_t* den, bool* is_exact) {
    if (e->type == EXPR_INTEGER) {
        *val = (double)e->data.integer;
        *num = e->data.integer;
        *den = 1;
        *is_exact = true;
        return true;
    } else if (e->type == EXPR_REAL) {
        *val = e->data.real;
        *is_exact = false;
        return true;
    } else if (is_rational(e, num, den)) {
        *val = (double)(*num) / (*den);
        *is_exact = true;
        return true;
    }

    /* A BigInt or a Rational with a BigInt part, met by a machine number (two
     * exact operands are compared exactly before this is reached). Read it
     * correctly rounded: past DBL_MAX that is +-inf, which orders right
     * against every finite double and is never tolerantly Equal to one; an
     * exact value that underflows keeps its sign and a magnitude below
     * DBL_MIN, so it orders below every normal machine number (and 0. ==
     * 10^-400 stays False, as 0. == 10^-300 is). 1.5 < 10^400 and
     * 1.*^308 < 10^309 were left unevaluated: N of these gives an extended-
     * range MPFR, which the fallback below declines. *num and *den are not set:
     * the value is used only on the inexact (tolerant) path. */
    {
        double d;
        if ((e->type == EXPR_BIGINT || is_rational_like(e)) && numeric_machine_value(e, &d)) {
            if (d == 0.0 || numeric_machine_subnormal(d)) {
                int sg = expr_numeric_sign(e);
                if (sg != 0 && d == 0.0) d = sg > 0 ? DBL_TRUE_MIN_OR_DENORM : -DBL_TRUE_MIN_OR_DENORM;
            }
            *val = d;
            *num = 0;
            *den = 1;
            *is_exact = false;
            return true;
        }
    }

    /* Fallback: drive the expression through N[] (machine precision).
     * Captures NumericQ-style values like Sqrt[3], Pi, or 1 + Sqrt[2]/3
     * that aren't raw numeric literals but have a definite real value.
     * If numericalize can't reduce the expression to a real number, we
     * leave it as unevaluated for the caller. The result is treated as
     * inexact so comparisons use the floating-point tolerance branch.
     * Only a numeric quantity is numericalised: N maps into the arguments
     * of an inert call, so ff[5] == x evaluated ff[5.] -- firing a
     * definition the reference never reaches (FeynCalc's GA[5.] guard
     * printed GA::noint inside every `... == x` test). */
    if (!expr_is_numeric_quantity(e)) return false;
    Expr* approx = numericalize(e, numeric_machine_spec());
    if (!approx) return false;
    bool ok = false;
    if (approx->type == EXPR_INTEGER) {
        *val = (double)approx->data.integer;
        *num = approx->data.integer;
        *den = 1;
        *is_exact = false;
        ok = true;
    } else if (approx->type == EXPR_REAL) {
        *val = approx->data.real;
        *is_exact = false;
        ok = true;
    } else if (is_rational(approx, num, den)) {
        *val = (double)(*num) / (*den);
        *is_exact = false;
        ok = true;
    }
    expr_free(approx);
    return ok;
}


/* An integer-like Expr into an already-initialised mpz (expr_to_mpz
 * initialises its target, so it goes through a temporary). */
static void cmp_set_mpz(mpz_t z, const Expr* e) {
    mpz_t t;
    expr_to_mpz(e, t);
    mpz_set(z, t);
    mpz_clear(t);
}

/* An exact rational (Integer, BigInt, or Rational with either kind of part)
 * read into `q`. */
static bool cmp_exact_rational(const Expr* e, mpq_t q) {
    if (expr_is_integer_like((Expr*)e)) {
        cmp_set_mpz(mpq_numref(q), e);
        mpz_set_ui(mpq_denref(q), 1);
        return true;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Rational) return false;
    const Expr* n = e->data.function.args[0];
    const Expr* d = e->data.function.args[1];
    if (!expr_is_integer_like((Expr*)n) || !expr_is_integer_like((Expr*)d)) return false;
    cmp_set_mpz(mpq_numref(q), n);
    cmp_set_mpz(mpq_denref(q), d);
    mpq_canonicalize(q);
    return true;
}

/* +1 for Infinity / DirectedInfinity[d] with d > 0, -1 for -Infinity /
 * DirectedInfinity[d] with d < 0 (Times[-1, Infinity] is how -Infinity
 * evaluates), 0 for anything else -- ComplexInfinity, I Infinity and
 * DirectedInfinity[] are not on the real line. */
static int real_infinity_sign(const Expr* e) {
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Infinity ? 1 : 0;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr* const* A = e->data.function.args;
    const Expr* dir = NULL;
    if (h == SYM_DirectedInfinity && n == 1) dir = A[0];
    else if (h == SYM_Times && n == 2 && A[1]->type == EXPR_SYMBOL &&
             A[1]->data.symbol.name == SYM_Infinity) dir = A[0];
    if (!dir) return 0;
    if (dir->type == EXPR_INTEGER) return (dir->data.integer > 0) - (dir->data.integer < 0);
    if (dir->type == EXPR_REAL) return (dir->data.real > 0) - (dir->data.real < 0);
    int64_t p, q;
    if (is_rational(dir, &p, &q)) return (p > 0) - (p < 0);
    return 0;
}

/* True iff `e` holds no approximate number anywhere: an exact numeric
 * quantity such as Pi + 10^-20, Exp[1 + 10^-20] or Pi/3^700. */
static bool cmp_is_exact_quantity(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return false;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return false;
#endif
    if (e->type != EXPR_FUNCTION) return true;
    for (size_t i = 0; i < e->data.function.arg_count; ++i)
        if (!cmp_is_exact_quantity(e->data.function.args[i])) return false;
    return cmp_is_exact_quantity(e->data.function.head);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------------
 *  Tolerant comparison of approximate numbers
 *
 *  Two approximate numbers are Equal when they agree to within the
 *  uncertainty of the less accurate of them, ignoring its last few binary
 *  digits. Concretely, with every inexact operand x assigned
 *
 *      precision  p(x) = its mantissa bits     (53 for a machine real,
 *                                               0 for an approximate zero)
 *      accuracy   a(x) = p(x) - log2|x|        (an approximate zero: its
 *                                               bit precision; machine 0.:
 *                                               53 + 1022 bits)
 *
 *  and an exact operand (Integer, Rational, Pi, Sqrt[2], ...) infinite
 *  accuracy, the operand g of LEAST accuracy governs:
 *
 *      x == y   iff   |x - y| <= 2^(k(g) - a(g)),
 *      k(g) = K * min(p(g), 53) / 53        (bits ignored)
 *
 *  with K = Internal`$EqualTolerance * log2 10 (7 by default; equal_tol.h)
 *
 *  so a number of at least machine precision has its last 7 bits (~2.107
 *  decimal digits) ignored, and a lower-precision one proportionally fewer
 *  -- an approximate zero none, it compares by accuracy alone. For two
 *  machine reals of similar size this is the familiar relative 2^-46.
 *  Complex operands compare part by part, each part with its own precision.
 *
 *  The ordering heads are consistent with it: inside the tolerance x < y is
 *  False and x <= y True; outside it the sign of x - y decides.
 *
 *  Before this rule existed an arbitrary-precision operand was read as a
 *  double and compared with the machine tolerance, so N[1 - 10^-20, 30] == 1
 *  was True and N[1 - 10^-20, 30] < 1 False. Two machine reals (or a machine
 *  real and an exact number) never come here: compare_numeric keeps its
 *  double path for them.
 * ---------------------------------------------------------------------- */

enum { CMP_PART_EXACT = 0, CMP_PART_MACHINE = 1, CMP_PART_MPFR = 2 };

/* One real component of an operand: its value at the working precision, and
 * how much it is worth -- exact, or inexact with `bits` of precision. */
typedef struct {
    mpfr_t v;
    int kind;
    long bits;
} CmpPart;

/* Bits of a part's accuracy, a(x) above. Only called on inexact parts. */
static double cmp_part_accuracy_bits(const CmpPart* p) {
    if (mpfr_zero_p(p->v))
        return p->kind == CMP_PART_MACHINE ? (double)(DBL_MANT_DIG - DBL_MIN_EXP + 1)
                                           : (double)p->bits;
    long e;
    double m = mpfr_get_d_2exp(&e, p->v, MPFR_RNDN);
    return (double)p->bits - ((double)e + log2(fabs(m)));
}

/* Load one real atom into `p`: an exact rational of either size, a machine
 * real, or an MPFR. `force_exact` marks a value that stands for an exact
 * quantity (a numericalized Pi, Sqrt[2], ...). False for anything else. */
static bool cmp_load_atom(const Expr* e, CmpPart* p, bool force_exact) {
    if (e->type == EXPR_REAL) {
        if (!isfinite(e->data.real)) return false;
        mpfr_set_d(p->v, e->data.real, MPFR_RNDN);
        p->kind = CMP_PART_MACHINE; p->bits = DBL_MANT_DIG;
    } else if (e->type == EXPR_MPFR) {
        if (!mpfr_number_p(e->data.mpfr)) return false;
        mpfr_set(p->v, e->data.mpfr, MPFR_RNDN);
        /* its precision, not the guard-padded mantissa width (numeric.h) */
        p->kind = CMP_PART_MPFR;
        p->bits = numeric_width_precision_bits((long)mpfr_get_prec(e->data.mpfr));
    } else {
        mpq_t q;
        mpq_init(q);
        bool ok = cmp_exact_rational(e, q);
        if (ok) mpfr_set_q(p->v, q, MPFR_RNDN);
        mpq_clear(q);
        if (!ok) return false;
        p->kind = CMP_PART_EXACT; p->bits = 0;
    }
    if (force_exact) { p->kind = CMP_PART_EXACT; p->bits = 0; }
    return true;
}

/* A real number or Complex[re, im] of real atoms into (re, im). */
static bool cmp_load_number(const Expr* e, CmpPart* re, CmpPart* im, bool force_exact) {
    Expr *r, *i;
    if (e->type == EXPR_FUNCTION && is_complex((Expr*)e, &r, &i))
        return cmp_load_atom(r, re, force_exact) && cmp_load_atom(i, im, force_exact);
    if (!cmp_load_atom(e, re, force_exact)) return false;
    mpfr_set_zero(im->v, +1);
    im->kind = CMP_PART_EXACT; im->bits = 0;
    return true;
}

/* An operand into (re, im) at `wbits`: a number directly, a closed-form
 * numeric quantity (Pi, 1 + Sqrt[2], E^I) through N at the working
 * precision, keeping its exactness. */
static bool cmp_load_operand(const Expr* e, long wbits, CmpPart* re, CmpPart* im) {
    if (cmp_load_number(e, re, im, false)) return true;
    if (!expr_is_numeric_quantity(e)) return false;
    NumericSpec s; s.mode = NUMERIC_MODE_MPFR; s.bits = wbits;
    Expr* n = numericalize(e, s);
    if (!n) return false;
    bool ok = cmp_load_number(n, re, im, true);
    expr_free(n);
    return ok;
}

/* Compare two parts. Returns false when both are exact (no tolerance to
 * apply: the caller's exact paths own that), else sets *out to 0 when they
 * are equal within the governing tolerance, otherwise to the sign of x - y. */
static bool cmp_parts(const CmpPart* x, const CmpPart* y, long wbits, int* out) {
    bool ix = x->kind != CMP_PART_EXACT, iy = y->kind != CMP_PART_EXACT;
    if (!ix && !iy) return false;
    const CmpPart* g;
    double acc;
    if (ix && iy) {
        double ax = cmp_part_accuracy_bits(x), ay = cmp_part_accuracy_bits(y);
        g = ax <= ay ? x : y;
        acc = ax <= ay ? ax : ay;
    } else {
        g = ix ? x : y;
        acc = cmp_part_accuracy_bits(g);
    }
    double prec = mpfr_zero_p(g->v) ? 0.0 : (double)g->bits;
    double ignored = mth_equal_tolerance_bits()
                   * (prec < DBL_MANT_DIG ? prec : DBL_MANT_DIG) / DBL_MANT_DIG;

    mpfr_t d, tol;
    mpfr_init2(d, (mpfr_prec_t)wbits);
    mpfr_init2(tol, 64);
    mpfr_sub(d, x->v, y->v, MPFR_RNDN);
    int sign = mpfr_sgn(d);
    mpfr_abs(d, d, MPFR_RNDN);
    mpfr_set_d(tol, ignored - acc, MPFR_RNDN);
    mpfr_ui_pow(tol, 2, tol, MPFR_RNDN);
    *out = (sign == 0 || mpfr_cmp(d, tol) <= 0) ? 0 : (sign > 0 ? 1 : -1);
    mpfr_clears(d, tol, (mpfr_ptr)0);
    return true;
}

/* The tolerant comparison of a and b when at least one carries an MPFR
 * value. Returns false (caller keeps its own path) when the rule does not
 * apply: no MPFR operand, an operand that is not a number, or a pair of
 * exact parts. On success *out is
 *   want_order:  -1 / 0 / +1   (both operands real; a complex one declines)
 *   otherwise:    0 equal, 1 unequal   (complex parts compared one by one). */
static bool compare_approximate(Expr* a, Expr* b, bool want_order, int* out) {
    if (!numeric_any_mpfr(a, b)) return false;
    long wbits = numeric_combined_bits(a, b, 0) + 32;
    CmpPart ar, ai, br, bi;
    mpfr_inits2((mpfr_prec_t)wbits, ar.v, ai.v, br.v, bi.v, (mpfr_ptr)0);
    bool ok = cmp_load_operand(a, wbits, &ar, &ai) && cmp_load_operand(b, wbits, &br, &bi);
    if (ok) {
        int cr = 0, ci = 0;
        bool real_pair = mpfr_zero_p(ai.v) && ai.kind == CMP_PART_EXACT
                      && mpfr_zero_p(bi.v) && bi.kind == CMP_PART_EXACT;
        if (want_order) {
            ok = real_pair && cmp_parts(&ar, &br, wbits, &cr);
            if (ok) *out = cr;
        } else {
            bool dr = cmp_parts(&ar, &br, wbits, &cr);
            bool di = real_pair || cmp_parts(&ai, &bi, wbits, &ci);
            /* A pair of exact parts (Complex[0, x] against I y: both real
             * parts 0) carries no tolerance. Equal working values settle it;
             * unequal ones decline, as they may be two roundings of
             * the same exact value. */
            if (!dr && !mpfr_equal_p(ar.v, br.v)) ok = false;
            if (!di && !mpfr_equal_p(ai.v, bi.v)) ok = false;
            if (ok) *out = (cr == 0 && ci == 0) ? 0 : 1;
        }
    }
    mpfr_clears(ar.v, ai.v, br.v, bi.v, (mpfr_ptr)0);
    return ok;
}


/* `e` at `bits` into the real `v`, refusing a value that came back as a
 * machine real (some head ignored the requested precision, so the value
 * carries 53 bits whatever was asked) or with a non-zero imaginary part. */
static bool cmp_real_at(const Expr* e, long bits, mpfr_t v) {
    NumericSpec s; s.mode = NUMERIC_MODE_MPFR; s.bits = bits;
    Expr* n = numericalize(e, s);
    if (!n) return false;
    bool ok = false;
    if (n->type == EXPR_MPFR && mpfr_number_p(n->data.mpfr)) {
        mpfr_set(v, n->data.mpfr, MPFR_RNDN);
        ok = true;
    } else if (n->type == EXPR_INTEGER || n->type == EXPR_BIGINT ||
               is_rational_like(n)) {
        mpq_t q;
        mpq_init(q);
        ok = cmp_exact_rational(n, q);
        if (ok) mpfr_set_q(v, q, MPFR_RNDN);
        mpq_clear(q);
    }
    expr_free(n);
    return ok;
}

/* True iff `e` contains no Plus anywhere: a product, power or elementary
 * function of exact constants (Exp[-Exp[10]], Pi^2/7, Sin[355]).  Such a value
 * involves no subtraction, so its numerical value at any MPFR precision is
 * correct to that RELATIVE precision and its sign is the sign of that value --
 * there is no cancellation that a low precision could get wrong. */
static bool cmp_plus_free(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && h->data.symbol.name == SYM_Plus) return false;
    if (!cmp_plus_free(h)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; ++i)
        if (!cmp_plus_free(e->data.function.args[i])) return false;
    return true;
}

/* The structural magnitude of `e` at `bits` into `out`: the size its value
 * may have cancelled down from.
 *   Plus           the sum of its terms' magnitudes
 *   Times          the product of its factors'
 *   Power[b, p]    magnitude(b)^p, p a positive rational (Sqrt[c] included)
 *   Plus-free      |e| itself (no subtraction, full relative precision)
 * Anything else that contains a Plus (1/(a - b), Log[1 + a - b], ...) has no
 * such bound, and the function declines (false): the caller then has no MPFR
 * order and keeps its confirmed-by-zero-test path.
 *
 * The resolution threshold of cmp_exact_order is measured against this.
 * Measured against |value| instead, the value of a cancelling sum
 * (PolyLog[2, 1/3] - 2006534954/5479143815, 6.9*10^-21, whose machine value is
 * -5.6*10^-17) or the rounding noise of an exact zero -- Sqrt[5 + 2 Sqrt[6]] -
 * Sqrt[2] - Sqrt[3], or 169 (-1 + Sqrt[16 - 6 Sqrt[7]] + Sqrt[11 - 4 Sqrt[7]])
 * from Reduce -- is "resolved" against itself. */
static bool cmp_l1_at(const Expr* e, long bits, mpfr_t out) {
    if (cmp_plus_free(e)) {
        mpfr_t t;
        mpfr_init2(t, (mpfr_prec_t)bits);
        bool ok = cmp_real_at(e, bits, t);
        if (ok) mpfr_abs(out, t, MPFR_RNDN);
        mpfr_clear(t);
        return ok;
    }
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t n = e->data.function.arg_count;
    Expr* const* A = e->data.function.args;
    if (h == SYM_Plus || h == SYM_Times) {
        mpfr_t t;
        mpfr_init2(t, (mpfr_prec_t)bits);
        if (h == SYM_Plus) mpfr_set_zero(out, 1); else mpfr_set_ui(out, 1, MPFR_RNDN);
        bool ok = true;
        for (size_t i = 0; ok && i < n; ++i) {
            ok = cmp_l1_at(A[i], bits, t);
            if (ok) {
                if (h == SYM_Plus) mpfr_add(out, out, t, MPFR_RNDU);
                else               mpfr_mul(out, out, t, MPFR_RNDU);
            }
        }
        mpfr_clear(t);
        return ok;
    }
    if (h == SYM_Power && n == 2) {
        mpq_t q;
        mpq_init(q);
        bool ok = cmp_exact_rational(A[1], q) && mpq_sgn(q) > 0 && cmp_l1_at(A[0], bits, out);
        if (ok) {
            mpfr_t p;
            mpfr_init2(p, (mpfr_prec_t)bits);
            mpfr_set_q(p, q, MPFR_RNDN);
            mpfr_pow(out, out, p, MPFR_RNDU);
            mpfr_clear(p);
        }
        mpq_clear(q);
        return ok;
    }
    return false;
}

/* Outcome of cmp_exact_order. */
enum { CMP_ORDER_NA = 0, CMP_ORDER_RESOLVED = 1, CMP_ORDER_UNRESOLVED = 2 };

/* The order of two EXACT real quantities.
 *
 * First from their exact difference: a - b, evaluated, is frequently a
 * rational (Pi + 10^-400 - Pi is 10^-400) or a Plus-free constant
 * (Pi + Exp[-Exp[10]] - Pi is Exp[-Exp[10]]), and then its sign is the order
 * at once -- exactly, or from a single MPFR read that cannot cancel
 * (cmp_plus_free).
 *
 * Otherwise at increasing precision.  A double cannot order them: N[Pi +
 * 10^-20] and N[Pi] are the same double, so the machine comparison called
 * them equal (Pi + 10^-20 == Pi was True), and Pi/3^700 has no double at all
 * (it is below 2^-1074), so Pi/3^700 > 0 stayed unevaluated.  At p bits each
 * side is good to a relative 2^-p, so a difference above 2^(16 - p) of the
 * larger side is resolved and its sign is the order.  Every magnitude stays in
 * MPFR, so the exponent range is no limit.
 *
 * Returns CMP_ORDER_RESOLVED with *out set; CMP_ORDER_UNRESOLVED when the
 * difference is below the top rung's resolution -- NOT "equal": that verdict
 * made Pi + 10^-400 == Pi True and Sqrt[2] + 10^-400 > Sqrt[2] False, and the
 * caller leaves such a pair to the exact zero test; CMP_ORDER_NA when a side
 * does not evaluate to a real number at the requested precision (the caller
 * keeps its own rule). */
static int cmp_exact_order(const Expr* a, const Expr* b, int* out) {
    {
        Expr* d = eval_and_free(
            internal_subtract((Expr*[]){ expr_copy((Expr*)a), expr_copy((Expr*)b) }, 2));
        int verdict = CMP_ORDER_NA;
        if (d) {
            mpq_t q;
            mpq_init(q);
            if (cmp_exact_rational(d, q)) {
                int sg = mpq_sgn(q);
                *out = (sg > 0) - (sg < 0);
                verdict = CMP_ORDER_RESOLVED;
            } else if (cmp_plus_free(d) && expr_is_numeric_quantity(d)) {
                mpfr_t v;
                mpfr_init2(v, 128);
                if (cmp_real_at(d, 128, v) && !mpfr_zero_p(v)) {
                    *out = mpfr_sgn(v) > 0 ? 1 : -1;
                    verdict = CMP_ORDER_RESOLVED;
                }
                mpfr_clear(v);
            }
            mpq_clear(q);
            expr_free(d);
        }
        if (verdict == CMP_ORDER_RESOLVED) return verdict;
    }

    static const long rungs[] = { 128, 1024 };
    int verdict = CMP_ORDER_NA;
    mpfr_t x, y, d, t;
    mpfr_inits2(64, x, y, d, t, (mpfr_ptr)0);
    for (size_t r = 0; r < sizeof rungs / sizeof rungs[0]; ++r) {
        long bits = rungs[r];
        mpfr_set_prec(x, (mpfr_prec_t)bits);
        mpfr_set_prec(y, (mpfr_prec_t)bits);
        mpfr_set_prec(d, (mpfr_prec_t)bits);
        mpfr_set_prec(t, (mpfr_prec_t)bits);
        if (!cmp_real_at(a, bits, x) || !cmp_real_at(b, bits, y)) { verdict = CMP_ORDER_NA; break; }
        verdict = CMP_ORDER_UNRESOLVED;
        mpfr_sub(d, x, y, MPFR_RNDN);
        if (mpfr_zero_p(d)) continue;
        /* t = 2^(16 - bits) max(L1(a), L1(b)), the structural magnitudes */
        {
            mpfr_t la, lb;
            mpfr_inits2((mpfr_prec_t)bits, la, lb, (mpfr_ptr)0);
            bool lok = cmp_l1_at(a, bits, la) && cmp_l1_at(b, bits, lb);
            if (lok) mpfr_max(t, la, lb, MPFR_RNDU);
            mpfr_clears(la, lb, (mpfr_ptr)0);
            if (!lok) { verdict = CMP_ORDER_NA; break; }
        }
        mpfr_mul_2si(t, t, 16 - bits, MPFR_RNDN);
        if (mpfr_cmpabs(d, t) > 0) {
            *out = mpfr_sgn(d) > 0 ? 1 : -1;
            verdict = CMP_ORDER_RESOLVED;
            break;
        }
    }
    mpfr_clears(x, y, d, t, (mpfr_ptr)0);
    return verdict;
}
#endif

/* How a compare_numeric order was reached -- what it may be trusted for.
 *   CMP_CERT_EXACT     exact arithmetic (GMP), or an order of two exact
 *                      quantities resolved in MPFR (cmp_exact_order) or read
 *                      off their exact difference: final.
 *   CMP_CERT_TOLERANT  an approximate operand, decided by the tolerance of
 *                      Internal`$EqualTolerance: final (the difference of two
 *                      approximate numbers is itself approximate, so no exact
 *                      test can overrule it).
 *   CMP_CERT_MACHINE   exact operands read only as doubles: a nonzero order
 *                      may be a cancellation residue (Sqrt[5 + 2 Sqrt[6]] -
 *                      Sqrt[2] - Sqrt[3] reads ~1e-16), so the callers confirm a
 *                      nearly-equal or unequal verdict with the exact zero test. */
enum { CMP_CERT_EXACT = 0, CMP_CERT_TOLERANT = 1, CMP_CERT_MACHINE = 2 };

static int compare_numeric_ex(Expr* a, Expr* b, bool* can_compare, int* cert);

static int compare_numeric(Expr* a, Expr* b, bool* can_compare) {
    int cert;
    return compare_numeric_ex(a, b, can_compare, &cert);
}

static int compare_numeric_ex(Expr* a, Expr* b, bool* can_compare, int* cert) {
    double va, vb;
    int64_t na, da, nb, db;
    bool exact_a, exact_b;

    *can_compare = false;
    *cert = CMP_CERT_EXACT;

    /* Two machine integers: compare in place.  Every guard and the GMP route
     * below give exactly this order for them, at the cost of two mpz
     * allocations per i == j / j > i test in a Table body. */
    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
        *can_compare = true;
        return (a->data.integer > b->data.integer) - (a->data.integer < b->data.integer);
    }

    /* The real infinities order against each other and against every real
     * number, as in the reference: 0 < Infinity, -Infinity < Pi, Infinity <=
     * Infinity. A complex or symbolic other side stays undecided. */
    {
        int ia = real_infinity_sign(a), ib = real_infinity_sign(b);
        if (ia || ib) {
            if (ia && ib) { *can_compare = true; return (ia > ib) - (ia < ib); }
            /* any real number is finite, so the infinity side wins exactly */
            if (!get_numeric_value(ia ? b : a, &va, &na, &da, &exact_a)) return 0;
            *can_compare = true;
            return ia ? ia : -ib;
        }
    }

    /* Interval vs scalar: decide only when the interval lies entirely on one
     * side of the scalar (disjoint), so Interval[{5,8}] > Pi -> True. An
     * interval straddling the scalar, or two intervals, stays symbolic. */
    if (is_interval(a) || is_interval(b)) {
        *cert = CMP_CERT_TOLERANT;
        if (is_interval(a) && is_interval(b)) return interval_compare_intervals(a, b, can_compare);
        if (is_interval(a)) return interval_compare_scalar(a, b, can_compare);
        bool d; int s = interval_compare_scalar(b, a, &d);
        *can_compare = d;
        return -s;
    }

    /* Exact integer-like comparison via GMP. Avoids precision loss when
     * either operand is a BigInt whose magnitude exceeds 2^53 (the limit
     * for distinguishing adjacent values through double). Without this,
     * Less[10^30, 10^30 + 1] coerces both sides to the same double and
     * answers False. */
    if (expr_is_integer_like(a) && expr_is_integer_like(b)) {
        mpz_t ma, mb;
        mpz_init(ma); mpz_init(mb);
        expr_to_mpz(a, ma);
        expr_to_mpz(b, mb);
        int cmp = mpz_cmp(ma, mb);
        mpz_clear(ma); mpz_clear(mb);
        *can_compare = true;
        if (cmp < 0) return -1;
        if (cmp > 0) return 1;
        return 0;
    }

    /* Two exact rationals compare exactly. The int64 n/d read below cannot
     * hold a BigInt part, so 1 + 10^-20 used to reach the tolerant double
     * comparison and test Equal to 1. */
    if ((a->type == EXPR_FUNCTION || b->type == EXPR_FUNCTION) &&
        (a->type == EXPR_FUNCTION || expr_is_integer_like(a)) &&
        (b->type == EXPR_FUNCTION || expr_is_integer_like(b))) {
        mpq_t qa, qb;
        mpq_init(qa); mpq_init(qb);
        bool ok = cmp_exact_rational(a, qa) && cmp_exact_rational(b, qb);
        int cmp = ok ? mpq_cmp(qa, qb) : 0;
        mpq_clear(qa); mpq_clear(qb);
        if (ok) {
            *can_compare = true;
            return (cmp > 0) - (cmp < 0);
        }
    }

#ifdef USE_MPFR
    {
        int cmp;
        if (compare_approximate(a, b, true, &cmp)) {
            *can_compare = true;
            *cert = CMP_CERT_TOLERANT;
            return cmp;
        }
    }
#endif

    bool got_a = get_numeric_value(a, &va, &na, &da, &exact_a);
    bool got_b = got_a && get_numeric_value(b, &vb, &nb, &db, &exact_b);

#ifdef USE_MPFR
    /* Two exact quantities, at least one beyond an exact rational (Pi + 10^-20,
     * Pi/3^700).  Their doubles decide only a clear order: both in the normal
     * range (or an exact 0) and apart by more than 2^-40 relative.  Anything
     * closer -- or a value the double range cannot hold, which N returns as an
     * extended-range MPFR and get_numeric_value declines -- is ordered at
     * increasing precision (cmp_exact_order). */
    if (cmp_is_exact_quantity(a) && cmp_is_exact_quantity(b) &&
        expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b)) {
        bool clear = false;
        /* Doubles are trusted only for operands that cannot cancel: a sum's
         * machine value carries no relative accuracy (PolyLog[2, 1/3] -
         * 2006534954/5479143815 is +6.9*10^-21 and reads -5.6*10^-17, so it
         * was "clearly" below 0). */
        if (got_b && cmp_plus_free(a) && cmp_plus_free(b) && isfinite(va) && isfinite(vb) &&
            (va == 0.0 || fabs(va) >= DBL_MIN) && (vb == 0.0 || fabs(vb) >= DBL_MIN)) {
            double mx = fabs(va) > fabs(vb) ? fabs(va) : fabs(vb);
            clear = fabs(va - vb) > mx * 0x1p-40;
        }
        int cmp = 0;
        int how = clear ? CMP_ORDER_NA : cmp_exact_order(a, b, &cmp);
        if (how == CMP_ORDER_RESOLVED) {
            *can_compare = true;
            return cmp;
        }
        /* Closer than the top MPFR rung resolves: no numeric order to give.
         * The caller's exact zero test decides equality; the ordering heads
         * stay unevaluated rather than read a 1024-bit tie as "equal". */
        if (how == CMP_ORDER_UNRESOLVED) return 0;
    }
#endif

    if (!got_b) return 0;

    *can_compare = true;
    *cert = (exact_a && exact_b) ? CMP_CERT_EXACT
          : (cmp_is_exact_quantity(a) && cmp_is_exact_quantity(b)) ? CMP_CERT_MACHINE
          : CMP_CERT_TOLERANT;
    if (exact_a && exact_b) {
        long double val_a = (long double)na / da;
        long double val_b = (long double)nb / db;
        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
        return 0;
    } else {
        /* the last Internal`$EqualTolerance digits ignored: 2^-46 relative
         * by default (equal_tol.h) */
        return mth_machine_tolerant_cmp(va, vb, mth_machine_equal_rel());
    }
}

/* A literal datum whose comparison verdict is final: a number atom, a string,
 * or an exact rational of ANY size.  A rational with a BigInt part counts:
 * compare_numeric decides it exactly in GMP, and the exact zero test it would
 * otherwise be re-checked by reads 2^2000/3^1300 at machine precision, where
 * a denominator past the double range makes it 0 -- so Equal[2^2000/3^1300, 0]
 * was True (and every rational beyond 1e308 in the denominator tested equal
 * to 0). */
static bool is_raw_data(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_STRING) return true;
    return is_rational_like(e);
}

/*
 * IEEE 754 / ISO 60559 makes a NaN *unordered* with respect to every value,
 * itself included: Equal, Less, LessEqual, Greater and GreaterEqual all yield
 * False, and Unequal yields True.  Indeterminate is Mathilda's NaN — 0/0,
 * Infinity - Infinity and 0 ComplexInfinity all normalise to it — so the
 * comparison heads give it exactly that treatment.  The consequence that
 * matters most is that Equal must not reach its structural-identity shortcut
 * for such operands, since expr_eq would otherwise report
 * Indeterminate == Indeterminate as True.
 *
 * An EXPR_REAL carrying a machine NaN is the same object one level down (see
 * the isnan branch in expr_eq), so it counts here too.
 *
 * SameQ and UnsameQ are deliberately exempt: they test structural identity of
 * expressions, not mathematical equality of values — the same distinction IEEE
 * draws between its comparison predicates and totalOrder, which does order a
 * NaN with itself.  Indeterminate === Indeterminate therefore stays True.
 */
static bool is_nan_operand(Expr* e) {
    if (!e) return false;
    if (is_indeterminate_sym(e)) return true;
    return e->type == EXPR_REAL && isnan(e->data.real);
}

/* True if any argument of the comparison `res` is unordered in the IEEE sense.
 * Every argument of a chained comparison takes part in at least one adjacent
 * pair, so a single NaN anywhere decides the whole chain. */
static bool has_nan_operand(Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (is_nan_operand(res->data.function.args[i])) return true;
    }
    return false;
}

/* ------------------------------------------------------------------------
 *  SameQ on approximate numbers
 *
 *  SameQ is structural identity, except that two approximate numbers are the
 *  same when they differ in at most their last k binary digits, k =
 *  Internal`$SameQTolerance * log2 10 (1 bit by default), at the LOWER of the
 *  two precisions -- 1. === 1. + 2^-52, 1`30 === 1`20 and 1. === 1`30 are all
 *  True, 1. === 1. + 2^-51 is False.  The rule applies at every leaf of the
 *  two expressions ({1.} === {1. + 2^-52} and f[1.] === f[1. + 2^-52] are
 *  True).  expr_eq itself stays exact: hashing, the matcher and the dedup
 *  tables depend on that, which is why MemberQ, Count and patterns do not
 *  carry the tolerance (nor do they in the reference).
 * ---------------------------------------------------------------------- */

/* True for an approximate atom (a finite machine real or MPFR); *bits gets
 * its precision. */
static bool sameq_inexact_atom(const Expr* e, double* bits) {
    if (e->type == EXPR_REAL) {
        if (!isfinite(e->data.real)) return false;
        *bits = DBL_MANT_DIG;
        return true;
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (!mpfr_number_p(e->data.mpfr)) return false;
        *bits = (double)numeric_width_precision_bits((long)mpfr_get_prec(e->data.mpfr));
        return true;
    }
#endif
    return false;
}

/* Two approximate atoms the same within k ignored bits at the lower
 * precision. */
static bool sameq_inexact_pair(const Expr* a, const Expr* b, double k) {
    double pa, pb;
    if (!sameq_inexact_atom(a, &pa) || !sameq_inexact_atom(b, &pb)) return false;
    double p = pa < pb ? pa : pb;
    if (a->type == EXPR_REAL && b->type == EXPR_REAL) {
        double x = a->data.real, y = b->data.real;
        if (x == y) return true;
        double mx = fmax(fabs(x), fabs(y));
        double tol = mx * exp2(k - (double)DBL_MANT_DIG);
        /* a subnormal's last digits are absolute: one ulp is 2^-1074 */
        double floor_tol = exp2(k - 1.0) * (DBL_MIN * DBL_EPSILON);   /* 2^-1074 */
        return fabs(x - y) <= (tol > floor_tol ? tol : floor_tol);
    }
#ifdef USE_MPFR
    long wbits = (long)(pa > pb ? pa : pb) + 64;
    mpfr_t x, y, d, t;
    mpfr_inits2((mpfr_prec_t)wbits, x, y, d, t, (mpfr_ptr)0);
    if (a->type == EXPR_REAL) mpfr_set_d(x, a->data.real, MPFR_RNDN); else mpfr_set(x, a->data.mpfr, MPFR_RNDN);
    if (b->type == EXPR_REAL) mpfr_set_d(y, b->data.real, MPFR_RNDN); else mpfr_set(y, b->data.mpfr, MPFR_RNDN);
    mpfr_sub(d, x, y, MPFR_RNDN);
    mpfr_abs(d, d, MPFR_RNDN);
    if (mpfr_cmpabs(x, y) >= 0) mpfr_abs(t, x, MPFR_RNDN); else mpfr_abs(t, y, MPFR_RNDN);
    /* t = max |.| * 2^(k - p) */
    mpfr_t sc; mpfr_init2(sc, 64);
    mpfr_set_d(sc, k - p, MPFR_RNDN);
    mpfr_ui_pow(sc, 2, sc, MPFR_RNDN);
    mpfr_mul(t, t, sc, MPFR_RNDN);
    bool same = mpfr_cmp(d, t) <= 0;
    mpfr_clear(sc);
    mpfr_clears(x, y, d, t, (mpfr_ptr)0);
    return same;
#else
    (void)p;
    return false;
#endif
}

/* SameQ of a and b with `k` ignored bits on approximate leaves. */
static bool sameq_tolerant(const Expr* a, const Expr* b, double k) {
    if (expr_eq(a, b)) return true;
    if (k <= 0.0) return false;
    double pa, pb;
    if (sameq_inexact_atom(a, &pa) && sameq_inexact_atom(b, &pb))
        return sameq_inexact_pair(a, b, k);
    /* a packed list against anything: compare as the List it is */
    if (is_packed_list(a) || is_packed_list(b)) {
        Expr* ua = is_packed_list(a) ? pack_unpack(a) : NULL;
        Expr* ub = is_packed_list(b) ? pack_unpack(b) : NULL;
        bool same = (ua || !is_packed_list(a)) && (ub || !is_packed_list(b))
                 && sameq_tolerant(ua ? ua : a, ub ? ub : b, k);
        if (ua) expr_free(ua);
        if (ub) expr_free(ub);
        return same;
    }
    if (a->type != EXPR_FUNCTION || b->type != EXPR_FUNCTION) return false;
    size_t n = a->data.function.arg_count;
    if (n != b->data.function.arg_count) return false;
    if (!sameq_tolerant(a->data.function.head, b->data.function.head, k)) return false;
    for (size_t i = 0; i < n; i++)
        if (!sameq_tolerant(a->data.function.args[i], b->data.function.args[i], k)) return false;
    return true;
}

/*
 * builtin_sameq: Implements SameQ[x, y, ...].
 * SameQ returns True if all its arguments are identical, and False otherwise.
 * Identical means the same structure and atoms, approximate numbers within
 * Internal`$SameQTolerance (see above).
 */
Expr* builtin_sameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }

    // SameQ[] and SameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }

    double k = mth_sameq_tolerance_bits();
    Expr* first = res->data.function.args[0];
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        if (!sameq_tolerant(first, res->data.function.args[i], k)) {
            return expr_new_symbol(SYM_False);
        }
    }

    return expr_new_symbol(SYM_True);
}

/*
 * builtin_unsameq: Implements UnsameQ[e1, e2, ...].
 * UnsameQ returns True if no two of the ei are identical (in SameQ's sense),
 * and False otherwise.
 */
Expr* builtin_unsameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }

    // UnsameQ[] and UnsameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }

    double k = mth_sameq_tolerance_bits();
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            if (sameq_tolerant(res->data.function.args[i], res->data.function.args[j], k)) {
                return expr_new_symbol(SYM_False);
            }
        }
    }

    return expr_new_symbol(SYM_True);
}

/* Exact zero-test on (a - b) for two numeric-quantity operands.  Mirrors the
 * cancellation guard already used inside builtin_equal: build a - b, evaluate,
 * and consult the exact zero-test procedure.  Returns ZERO_TEST_TRUE (the two
 * are equal), ZERO_TEST_FALSE (provably unequal) or ZERO_TEST_UNKNOWN.
 *
 * This is what lets the comparison heads decide constants that carry no real
 * ordering — I == 0, 2 I Pi != 0, ½(-1-I Sqrt[3]) != 1 — which compare_numeric
 * leaves undecided.  Callers gate it on expr_is_numeric_quantity for BOTH
 * operands, so a free symbol (x == 0) never reaches here and stays symbolic. */
static ZeroTestResult numeric_pair_zero_test(Expr* a, Expr* b) {
    Expr* diff = eval_and_free(
        internal_subtract((Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
    ZeroTestResult zt = diff ? zero_test_decide(diff) : ZERO_TEST_UNKNOWN;
    if (diff) expr_free(diff);
    return zt;
}

/* Equality of two associations, as Mathematica 15 decides it:
 *
 *     <|a -> 1, b -> 2|> == <|b -> 2, a -> 1|>    False  (order matters)
 *     <|a -> 1|> == <|a -> 1.|>                   True   (values by Equal)
 *     <|1 -> a|> == <|1. -> a|>                   False  (keys by SameQ)
 *     <|a :> 1|> == <|a -> 1|>                    False  (Rule vs RuleDelayed)
 *     <|a -> x|> == <|a -> y|>                    stays unevaluated
 *
 * The entries must agree pairwise, in order, in their Rule/RuleDelayed head and
 * (structurally) in their key; the values are then compared with Equal itself.
 * Returns ASSOC_EQ_TRUE / ASSOC_EQ_FALSE / ASSOC_EQ_UNKNOWN, or ASSOC_EQ_NA
 * when a and b are not both well-formed associations (the caller's ordinary
 * rules then apply; <|a -> 1|> == 1 stays unevaluated). */
enum { ASSOC_EQ_NA = -2, ASSOC_EQ_UNKNOWN = -1, ASSOC_EQ_FALSE = 0, ASSOC_EQ_TRUE = 1 };

static int assoc_equal_verdict(Expr* a, Expr* b) {
    if (!assoc_is_wellformed(a) || !assoc_is_wellformed(b)) return ASSOC_EQ_NA;
    size_t n = a->data.function.arg_count;
    if (n != b->data.function.arg_count) return ASSOC_EQ_FALSE;
    for (size_t i = 0; i < n; i++) {
        Expr* ea = a->data.function.args[i];
        Expr* eb = b->data.function.args[i];
        if (ea->data.function.head->data.symbol.name != eb->data.function.head->data.symbol.name ||
            !expr_eq(ea->data.function.args[0], eb->data.function.args[0]))
            return ASSOC_EQ_FALSE;
    }
    int verdict = ASSOC_EQ_TRUE;
    for (size_t i = 0; i < n; i++) {
        Expr* va = a->data.function.args[i]->data.function.args[1];
        Expr* vb = b->data.function.args[i]->data.function.args[1];
        if (expr_eq(va, vb)) continue;
        Expr* pair[2] = { expr_copy(va), expr_copy(vb) };
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Equal), pair, 2));
        bool t = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
        bool f = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_False;
        if (r) expr_free(r);
        if (f) return ASSOC_EQ_FALSE;
        if (!t) verdict = ASSOC_EQ_UNKNOWN;
    }
    return verdict;
}

/* Equal on one pair (a, b), borrowing both operands: builds a transient
 * Equal[a, b] that shares a and b, runs builtin_equal on it, and detaches them
 * again before freeing the wrapper. ASSOC_EQ_TRUE / _FALSE / _UNKNOWN. */
static int equal_pair_verdict(Expr* a, Expr* b) {
    Expr* pair[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_Equal), pair, 2);
    Expr* r = builtin_equal(call);
    call->data.function.args[0] = NULL;
    call->data.function.args[1] = NULL;
    expr_free(call);
    int v = ASSOC_EQ_UNKNOWN;
    if (r && r->type == EXPR_SYMBOL) {
        if (r->data.symbol.name == SYM_True) v = ASSOC_EQ_TRUE;
        else if (r->data.symbol.name == SYM_False) v = ASSOC_EQ_FALSE;
    }
    if (r) expr_free(r);
    return v;
}

static bool is_list_head(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List;
}

/* Equality of two lists, as Mathematica decides it: element by element with
 * Equal itself, so the machine-real tolerance and every other scalar rule
 * apply at each position, recursively through nested lists.
 *
 *     {0.1 + 0.2, 1.} == {0.3, 1.}     True   (tolerance, elementwise)
 *     {1, 2} == {1, 2, 3}              False  (lengths differ)
 *     {x, 1} == {y, 2}                 False  (one position decidably False)
 *     {x} == {y}                       stays unevaluated
 *
 * Returns ASSOC_EQ_NA unless both a and b are Lists (a structurally identical
 * pair never gets here: the caller's expr_eq shortcut answers it first). */
static int list_equal_verdict(Expr* a, Expr* b) {
    if (!is_list_head(a) || !is_list_head(b)) return ASSOC_EQ_NA;
    size_t n = a->data.function.arg_count;
    if (n != b->data.function.arg_count) return ASSOC_EQ_FALSE;
    int verdict = ASSOC_EQ_TRUE;
    for (size_t i = 0; i < n; i++) {
        Expr* ea = a->data.function.args[i];
        Expr* eb = b->data.function.args[i];
        if (expr_eq(ea, eb)) continue;
        int v = equal_pair_verdict(ea, eb);
        if (v == ASSOC_EQ_FALSE) return ASSOC_EQ_FALSE;
        if (v != ASSOC_EQ_TRUE) verdict = ASSOC_EQ_UNKNOWN;
    }
    return verdict;
}

/*
 * builtin_equal: Implements Equal[lhs, rhs, ...].
 * Equal returns True if its arguments are identical (SameQ) or numerically equal (2 == 2.0).
 * Otherwise it returns NULL (unevaluated).
 */
/* Decide through the Lists; when that does not settle, leave the call as the
 * user wrote it (NDArray[...] == x stays NDArray[...] == x). */
static Expr* equal_via_lists(Expr* res) {
    Expr* r = ndarray_delist_and_reeval(res);
    if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL
        && expr_eq(r->data.function.head, res->data.function.head)) {
        expr_free(r);
        return NULL;
    }
    return r;
}

/* (Kept below equal_via_lists: check_packed_aware.py attributes an
 * is_ndarray() call to the preceding Expr* definition.)
 * A visible NDArray[...] operand: compare through its nested List, as the
 * packed surface already is (list_equal_verdict), so NDArray[{3., 2.}] ==
 * {3., 2.} is True rather than unevaluated -- the representations must agree. */
static bool has_visible_ndarray_operand(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        const Expr* a = res->data.function.args[i];
        if (is_ndarray(a) && !is_packed_list(a)) return true;
    }
    return false;
}

/* Equal on a pair with an arbitrary-precision operand, by the tolerance of
 * compare_approximate (complex parts one by one). True when it decides,
 * with *unequal set; the verdict is final -- an exact zero test of the
 * difference must not override it, since the difference of two approximate
 * numbers is itself approximate. */
static bool approximate_equal_verdict(Expr* a, Expr* b, bool* unequal) {
#ifdef USE_MPFR
    int v;
    if (compare_approximate(a, b, false, &v)) { *unequal = v != 0; return true; }
#else
    (void)a; (void)b; (void)unequal;
#endif
    return false;
}

/* An approximate number: a machine real, an MPFR, or Complex[] with one. */
static bool is_inexact_number(Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr *re, *im;
    return e->type == EXPR_FUNCTION && is_complex(e, &re, &im)
        && (is_inexact_number(re) || is_inexact_number(im));
}

/* The tolerance makes approximate equality intransitive: in a chain
 * x1 == x2 == x3 whose neighbours each agree, x1 and x3 may still differ by
 * more than the tolerance, and then the chain is False -- every pair must
 * be equal. Exact equality is transitive, so only pairs of approximate
 * numbers are re-checked. True when some non-adjacent pair of the n values
 * v[0], v[stride], v[2 stride], ... is unequal (stride 2 walks the values of
 * an Inequality[v0, Equal, v1, Equal, v2] chain). */
static bool chain_has_unequal_inexact_pair(Expr* const* v, size_t n, size_t stride) {
    for (size_t i = 0; i + 2 < n; i++) {
        Expr* a = v[i * stride];
        if (!is_inexact_number(a)) continue;
        for (size_t j = i + 2; j < n; j++) {
            Expr* b = v[j * stride];
            if (!is_inexact_number(b)) continue;
            bool unequal;
            if (approximate_equal_verdict(a, b, &unequal)) {
                if (unequal) return true;
                continue;
            }
            bool can_compare = false;
            int cmp = compare_numeric(a, b, &can_compare);
            if (can_compare && cmp != 0) return true;
        }
    }
    return false;
}

/* The equality verdict of ONE pair (a, b): ASSOC_EQ_TRUE, ASSOC_EQ_FALSE or
 * ASSOC_EQ_UNKNOWN.  Equal and Unequal both read every pair through here, so
 * a pair Equal decides is decided the same way by Unequal and they are
 * complements by construction.  (They used to carry two copies of this logic,
 * and only Equal's ran the exact zero test, so 22/7 == Pi and 22/7 != Pi were
 * both True when that test misfired.)
 *
 *   - a NaN operand (Indeterminate, a machine NaN) is IEEE-unordered: unequal
 *     to everything, itself included -- checked before structural identity;
 *   - structurally identical: equal;
 *   - two associations / two lists: entry- / element-wise;
 *   - an approximate operand: by the tolerance (final);
 *   - otherwise by compare_numeric.  An order it reached EXACTLY (GMP, or two
 *     exact quantities resolved in MPFR or by their exact difference) or by
 *     the tolerance is final.  One read off machine doubles of EXACT operands
 *     may be a cancellation residue, so "unequal" is confirmed by the exact
 *     zero test, which overrides it only when it proves the difference zero
 *     (Sqrt[5 + 2 Sqrt[6]] == Sqrt[2] + Sqrt[3]);
 *   - a pair compare_numeric cannot order (complex constants, or exact
 *     quantities closer than MPFR resolves): two raw data are unequal; two
 *     closed-form numeric quantities are decided by the exact zero test,
 *     UNKNOWN when it cannot. */
static int equal_pair_decide(Expr* a, Expr* b) {
    if (is_nan_operand(a) || is_nan_operand(b)) return ASSOC_EQ_FALSE;
    if (expr_eq(a, b)) return ASSOC_EQ_TRUE;
    int av;
    if ((av = assoc_equal_verdict(a, b)) != ASSOC_EQ_NA) return av;
    if ((av = list_equal_verdict(a, b)) != ASSOC_EQ_NA) return av;
    bool unequal = false;
    if (approximate_equal_verdict(a, b, &unequal))
        return unequal ? ASSOC_EQ_FALSE : ASSOC_EQ_TRUE;

    bool can_compare = false;
    int cert = CMP_CERT_EXACT;
    int cmp = compare_numeric_ex(a, b, &can_compare, &cert);
    bool raw = is_raw_data(a) && is_raw_data(b);
    if (can_compare) {
        if (cmp == 0) return ASSOC_EQ_TRUE;
        if (raw || cert != CMP_CERT_MACHINE) return ASSOC_EQ_FALSE;
        Expr* diff = eval_and_free(
            internal_subtract((Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
        ZeroTestResult zt = diff ? zero_test_decide(diff) : ZERO_TEST_UNKNOWN;
        if (diff) expr_free(diff);
        return zt == ZERO_TEST_TRUE ? ASSOC_EQ_TRUE : ASSOC_EQ_FALSE;
    }
    if (raw) return ASSOC_EQ_FALSE;
    if (expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b)) {
        ZeroTestResult zt = numeric_pair_zero_test(a, b);
        if (zt == ZERO_TEST_TRUE) return ASSOC_EQ_TRUE;
        if (zt == ZERO_TEST_FALSE) return ASSOC_EQ_FALSE;
    }
    return ASSOC_EQ_UNKNOWN;
}

Expr* builtin_equal(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }

    // Equal[] and Equal[x] return True.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    if (has_visible_ndarray_operand(res)) return equal_via_lists(res);

    /* IEEE unordered: a NaN operand makes every pair it takes part in False. */
    if (has_nan_operand(res)) {
        return expr_new_symbol(SYM_False);
    }

    bool all_equal = true;
    for (size_t i = 0; i + 1 < res->data.function.arg_count; i++) {
        int v = equal_pair_decide(res->data.function.args[i], res->data.function.args[i + 1]);
        if (v == ASSOC_EQ_FALSE) return expr_new_symbol(SYM_False);
        if (v != ASSOC_EQ_TRUE) all_equal = false;
    }

    if (all_equal) {
        if (res->data.function.arg_count > 2
            && chain_has_unequal_inexact_pair(res->data.function.args,
                                              res->data.function.arg_count, 1))
            return expr_new_symbol(SYM_False);
        return expr_new_symbol(SYM_True);
    }

    return NULL;
}

/*
 * builtin_unequal: Implements Unequal[e1, e2, ...].
 * Unequal returns False if any two arguments are equal, and True if every
 * pair is decidably unequal -- each pair by equal_pair_decide, the verdict
 * Equal uses, so the two heads cannot disagree on a pair.
 */
Expr* builtin_unequal(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }

    // Unequal[] and Unequal[x] return True.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    if (has_visible_ndarray_operand(res)) return equal_via_lists(res);

    bool all_definitely_unequal = true;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            int v = equal_pair_decide(res->data.function.args[i], res->data.function.args[j]);
            if (v == ASSOC_EQ_TRUE) return expr_new_symbol(SYM_False);
            if (v != ASSOC_EQ_FALSE) all_definitely_unequal = false;
        }
    }

    if (all_definitely_unequal) {
        return expr_new_symbol(SYM_True);
    }

    return NULL;
}

/* True iff the numeric quantity `e` has a real machine value (not Complex,
 * not Indeterminate): the ordering heads compare only real numbers. */
static bool cmp_is_real_constant(Expr* e) {
    double v; int64_t n, d; bool ex;
    return get_numeric_value(e, &v, &n, &d, &ex) && !isnan(v);
}

/* The sign of a real EXACT numeric constant (no free symbol, no approximate
 * number), as the comparison heads decide e vs 0: an order resolved exactly is
 * final, one read off doubles is confirmed by the exact zero test.  True with
 * *sgn in {-1, 0, 1} when decided; false when e is not such a constant, is not
 * real at 128 bits, or the order cannot be settled.  This is the fallback
 * Sign / Abs use when the algebraic (FLINT) certificate declines --
 * Sign[-1 + Zeta[3]] stayed unevaluated. */
bool comparisons_real_constant_sign(Expr* e, int* sgn) {
    if (!e || !expr_is_numeric_quantity(e) || !cmp_is_exact_quantity(e)) return false;
#ifdef USE_MPFR
    {
        mpfr_t v;
        mpfr_init2(v, 128);
        bool real = cmp_real_at(e, 128, v);
        mpfr_clear(v);
        if (!real) return false;
    }
#else
    if (!cmp_is_real_constant(e)) return false;
#endif
    Expr* zero = expr_new_integer(0);
    bool can = false;
    int cert = CMP_CERT_EXACT;
    int cmp = compare_numeric_ex(e, zero, &can, &cert);
    expr_free(zero);
    if (can && cert != CMP_CERT_MACHINE) { *sgn = cmp; return true; }
    ZeroTestResult zt = zero_test_decide(e);
    if (zt == ZERO_TEST_TRUE) { *sgn = 0; return true; }
    if (zt == ZERO_TEST_FALSE && can && cmp != 0) { *sgn = cmp; return true; }
    return false;
}

/* Inequality helpers */
static Expr* evaluate_inequality(Expr* res, int expected_cmp_1, int expected_cmp_2) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count < 2) return expr_new_symbol(SYM_True);

    /* IEEE unordered: an ordering predicate against a NaN is always False. */
    if (has_nan_operand(res)) return expr_new_symbol(SYM_False);

    for (size_t i = 0; i < res->data.function.arg_count - 1; i++) {
        Expr* a = res->data.function.args[i];
        Expr* b = res->data.function.args[i+1];
        
        bool can_compare;
        int cert = CMP_CERT_EXACT;
        int cmp = compare_numeric_ex(a, b, &can_compare, &cert);

        if (!can_compare) {
            /* Two exact quantities closer than MPFR resolves have no numeric
             * order; when the exact zero test proves them EQUAL, the order is
             * equality (Sqrt[5 + 2 Sqrt[6]] <= Sqrt[2] + Sqrt[3]).  Otherwise
             * leave the comparison unevaluated. */
            if (cmp_is_exact_quantity(a) && cmp_is_exact_quantity(b) &&
                expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b) &&
                cmp_is_real_constant(a) && cmp_is_real_constant(b) &&
                numeric_pair_zero_test(a, b) == ZERO_TEST_TRUE)
                cmp = 0;
            else
                return NULL; // Unevaluated
        }

        /* An order read off machine values of EXACT operands can be a
         * cancellation residue: E^w - 3 w at w = -ProductLog[-1/3] is exactly
         * 0 but evaluates to ~1e-35, which made E^w - 3 w > 0 True.  As in
         * Equal, a nearly-equal exact pair is confirmed by the exact zero test
         * before an order is committed to -- but only an order read off
         * doubles (CMP_CERT_MACHINE).  An order resolved in MPFR or from the
         * exact difference is final: the zero test overruling it made
         * Pi < 80143857/25510582 and Pi > 80143857/25510582 both False. */
        if (cmp != 0 && cert == CMP_CERT_MACHINE && (!is_raw_data(a) || !is_raw_data(b))) {
            double va, vb; int64_t n1, d1, n2, d2; bool ea, eb;
            /* "Nearly equal" is RELATIVE (two small but distinct numbers, a
             * 1e-13 residual against 10^-12, are not), except against an exact
             * 0, where only the size of the other side can say it. */
            bool zero_side = (a->type == EXPR_INTEGER && a->data.integer == 0)
                          || (b->type == EXPR_INTEGER && b->data.integer == 0);
            if (get_numeric_value(a, &va, &n1, &d1, &ea) && get_numeric_value(b, &vb, &n2, &d2, &eb)
                && (zero_side ? fabs(va - vb) <= 1e-10
                              : fabs(va - vb) <= 1e-12 * fmax(fabs(va), fabs(vb)))) {
                if (numeric_pair_zero_test(a, b) == ZERO_TEST_TRUE) cmp = 0;
            }
        }
        
        if (cmp != expected_cmp_1 && cmp != expected_cmp_2) {
            return expr_new_symbol(SYM_False);
        }
    }
    return expr_new_symbol(SYM_True);
}

Expr* builtin_less(Expr* res) {
    return evaluate_inequality(res, -1, -1);
}

Expr* builtin_greater(Expr* res) {
    return evaluate_inequality(res, 1, 1);
}

Expr* builtin_lessequal(Expr* res) {
    return evaluate_inequality(res, -1, 0);
}

Expr* builtin_greaterequal(Expr* res) {
    return evaluate_inequality(res, 1, 0);
}

/* Map an operator-symbol head (Less/LessEqual/Greater/GreaterEqual/Equal/Unequal)
 * to a pairwise decision. Returns:
 *    1  if `a OP b` is definitely True,
 *    0  if it is definitely False,
 *   -1  if the comparison is not decidable (operand non-numeric or the
 *       op is not a recognised chain head).
 * Equal/Unequal defer to their builtins, so x == x is True even when x has
 * no numeric value and {0.1 + 0.2} == {0.3} applies the list tolerance. */
/* Run the Equal / Unequal builtin on the pair (a, b) without consuming
 * either: 1 for True, 0 for False, -1 when it stays unevaluated. Routing a
 * chained `==` / `!=` through the builtin keeps a chain's pairwise verdict
 * identical to the binary head's -- the list tolerance, association rules and
 * exact zero test included -- instead of a weaker private copy of it. */
static int decide_via_builtin(Expr* (*fn)(Expr*), const char* head, Expr* a, Expr* b) {
    Expr* pair[2] = { expr_copy(a), expr_copy(b) };
    Expr* call = expr_new_function(expr_new_symbol(head), pair, 2);
    Expr* r = fn(call);
    expr_free(call);
    int d = -1;
    if (r && r->type == EXPR_SYMBOL) {
        if (r->data.symbol.name == SYM_True) d = 1;
        else if (r->data.symbol.name == SYM_False) d = 0;
    }
    if (r) expr_free(r);
    return d;
}

static int decide_pair(const char* op, Expr* a, Expr* b) {
    if (op == SYM_Equal)   return decide_via_builtin(builtin_equal, SYM_Equal, a, b);
    if (op == SYM_Unequal) return decide_via_builtin(builtin_unequal, SYM_Unequal, a, b);
    /* IEEE unordered: every ordering head is False against a NaN operand. */
    if (is_nan_operand(a) || is_nan_operand(b)) {
        if (op == SYM_Less || op == SYM_LessEqual
         || op == SYM_Greater || op == SYM_GreaterEqual) return 0;
        return -1;
    }
    if (op == SYM_Less || op == SYM_LessEqual
     || op == SYM_Greater || op == SYM_GreaterEqual) {
        bool can; int cmp = compare_numeric(a, b, &can);
        if (!can) return -1;
        if (op == SYM_Less)         return cmp <  0 ? 1 : 0;
        if (op == SYM_LessEqual)    return cmp <= 0 ? 1 : 0;
        if (op == SYM_Greater)      return cmp >  0 ? 1 : 0;
        /* GreaterEqual */          return cmp >= 0 ? 1 : 0;
    }
    return -1;
}

/* Inequality[v0, op0, v1, op1, v2, ...]: Mathematica's variadic chained
 * comparison head, produced by the parser for `a < b <= c == d > e` and
 * similar chains. Semantics:
 *   - True iff every adjacent pair `v_i op_i v_{i+1}` is True.
 *   - False if any adjacent pair is decidably False.
 *   - Otherwise, drop pairs that we can prove True, and either collapse
 *     to a single binary head (when one pair remains) or rebuild a
 *     smaller Inequality with the surviving pairs.
 *
 * The argument layout is strict: 2n+1 args, even slots are values, odd
 * slots are bare operator symbols. Anything else returns NULL so the
 * caller's input is preserved unevaluated. */
Expr* builtin_inequality(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_True);
    if (n == 1) return expr_new_symbol(SYM_True);  /* one value, no pair */
    if ((n & 1u) == 0) return NULL;              /* must be 2k+1 */

    /* Validate operator slots upfront — every odd index must be one of the
     * five chain heads or Unequal (`a < b != c`). */
    for (size_t i = 1; i < n; i += 2) {
        Expr* op = res->data.function.args[i];
        if (op->type != EXPR_SYMBOL) return NULL;
        const char* s = op->data.symbol.name;
        if (s != SYM_Less && s != SYM_LessEqual
         && s != SYM_Greater && s != SYM_GreaterEqual
         && s != SYM_Equal && s != SYM_Unequal) return NULL;
    }

    /* Decide each pair; any False short-circuits the whole chain. We also
     * remember which pairs were "undecided" so we can rebuild a residual
     * Inequality if needed. */
    size_t npairs = (n - 1) / 2;
    bool* undecided = (bool*)calloc(npairs, sizeof(bool));
    if (!undecided) return NULL;
    bool any_undecided = false;
    for (size_t k = 0; k < npairs; k++) {
        Expr* a  = res->data.function.args[2*k];
        Expr* op = res->data.function.args[2*k + 1];
        Expr* b  = res->data.function.args[2*k + 2];
        int d = decide_pair(op->data.symbol.name, a, b);
        if (d == 0) {                    /* definitely False */
            free(undecided);
            return expr_new_symbol(SYM_False);
        }
        if (d < 0) {
            undecided[k] = true;
            any_undecided = true;
        }
    }
    if (!any_undecided) {
        free(undecided);
        /* Every maximal run of Equal links (v_s == ... == v_e) must hold
         * pairwise, not just between neighbours: see
         * chain_has_unequal_inexact_pair. */
        for (size_t k = 0; k < npairs; ) {
            size_t e = k;
            while (e < npairs && res->data.function.args[2*e + 1]->data.symbol.name == SYM_Equal) e++;
            if (e - k >= 2 && chain_has_unequal_inexact_pair(res->data.function.args + 2*k, e - k + 1, 2))
                return expr_new_symbol(SYM_False);
            k = e > k ? e : k + 1;
        }
        return expr_new_symbol(SYM_True);
    }

    /* Build a residual Inequality from the undecided pairs. Each surviving
     * pair contributes (lhs, op, rhs). When the previous emitted RHS is
     * structurally equal to the next pair's LHS (the common case when
     * undecided pairs are adjacent, but also when a dropped pair was
     * `b == b`), we omit the duplicate value so the chain stays well-
     * formed (2n+1 args). */
    size_t cap = 2 * npairs + 1;
    Expr** out_args = (Expr**)malloc(sizeof(Expr*) * cap);
    if (!out_args) { free(undecided); return NULL; }
    size_t out_n = 0;
    Expr* last_emitted_val = NULL;
    for (size_t k = 0; k < npairs; k++) {
        if (!undecided[k]) continue;
        Expr* a  = res->data.function.args[2*k];
        Expr* op = res->data.function.args[2*k + 1];
        Expr* b  = res->data.function.args[2*k + 2];
        if (!last_emitted_val || !expr_eq(last_emitted_val, a)) {
            out_args[out_n++] = expr_copy(a);
        }
        out_args[out_n++] = expr_copy(op);
        out_args[out_n++] = expr_copy(b);
        last_emitted_val = b;
    }
    free(undecided);

    if (out_n == 3) {
        /* Exactly one undecided pair survived — return it as the binary
         * comparison so downstream consumers see the familiar shape. */
        Expr* binary_args[2] = { out_args[0], out_args[2] };
        const char* head_sym = out_args[1]->data.symbol.name;
        expr_free(out_args[1]);
        Expr* result = expr_new_function(expr_new_symbol(head_sym), binary_args, 2);
        free(out_args);
        return result;
    }
    Expr* residual = expr_new_function(expr_new_symbol(SYM_Inequality), out_args, out_n);
    free(out_args);
    return residual;
}

/*
 * comparisons_init: Registers comparison-related builtins in the symbol table.
 */
void comparisons_init(void) {
    equal_tol_init();
    symtab_add_builtin("SameQ", builtin_sameq);
    symtab_add_builtin("UnsameQ", builtin_unsameq);
    symtab_add_builtin("Equal", builtin_equal);
    symtab_add_builtin("Unequal", builtin_unequal);
    symtab_add_builtin("Less", builtin_less);
    symtab_add_builtin("Greater", builtin_greater);
    symtab_add_builtin("LessEqual", builtin_lessequal);
    symtab_add_builtin("GreaterEqual", builtin_greaterequal);
    symtab_add_builtin("Inequality", builtin_inequality);
    set_attributes("Inequality", ATTR_PROTECTED);
}
