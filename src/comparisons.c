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
#include <stdio.h>
#include <gmp.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

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
 *      k(g) = 7 * min(p(g), 53) / 53        (bits ignored)
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
    double ignored = 7.0 * (prec < DBL_MANT_DIG ? prec : DBL_MANT_DIG) / DBL_MANT_DIG;

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
#endif

static int compare_numeric(Expr* a, Expr* b, bool* can_compare) {
    double va, vb;
    int64_t na, da, nb, db;
    bool exact_a, exact_b;

    *can_compare = false;

    /* The real infinities order against each other and against every real
     * number, as in the reference: 0 < Infinity, -Infinity < Pi, Infinity <=
     * Infinity. A complex or symbolic other side stays undecided. */
    {
        int ia = real_infinity_sign(a), ib = real_infinity_sign(b);
        if (ia || ib) {
            if (ia && ib) { *can_compare = true; return (ia > ib) - (ia < ib); }
            if (!get_numeric_value(ia ? b : a, &va, &na, &da, &exact_a)) return 0;
            *can_compare = true;
            return ia ? ia : -ib;
        }
    }

    /* Interval vs scalar: decide only when the interval lies entirely on one
     * side of the scalar (disjoint), so Interval[{5,8}] > Pi -> True. An
     * interval straddling the scalar, or two intervals, stays symbolic. */
    if (is_interval(a) || is_interval(b)) {
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
            return cmp;
        }
    }
#endif

    if (!get_numeric_value(a, &va, &na, &da, &exact_a)) return 0;
    if (!get_numeric_value(b, &vb, &nb, &db, &exact_b)) return 0;

    *can_compare = true;
    if (exact_a && exact_b) {
        long double val_a = (long double)na / da;
        long double val_b = (long double)nb / db;
        if (val_a < val_b) return -1;
        if (val_a > val_b) return 1;
        return 0;
    } else {
        double diff = va - vb;
        if (diff < 0) diff = -diff;
        double max_val = (va < 0 ? -va : va);
        if ((vb < 0 ? -vb : vb) > max_val) max_val = (vb < 0 ? -vb : vb);
        
        // 2^(-46) = 1.4210854715202004e-14
        if (diff <= max_val * 1.4210854715202004e-14 || diff == 0.0) {
            return 0;
        }
        if (va < vb) return -1;
        return 1;
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

/*
 * builtin_sameq: Implements SameQ[x, y, ...].
 * SameQ returns True if all its arguments are identical, and False otherwise.
 * Identical means they have the same structure and atoms.
 */
Expr* builtin_sameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // SameQ[] and SameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    Expr* first = res->data.function.args[0];
    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        if (!expr_eq(first, res->data.function.args[i])) {
            return expr_new_symbol(SYM_False);
        }
    }
    
    return expr_new_symbol(SYM_True);
}

/*
 * builtin_unsameq: Implements UnsameQ[e1, e2, ...].
 * UnsameQ returns True if no two of the ei are identical, and False otherwise.
 */
Expr* builtin_unsameq(Expr* res) {
    if (res->type != EXPR_FUNCTION) {
        return NULL;
    }
    
    // UnsameQ[] and UnsameQ[x] return True by convention.
    if (res->data.function.arg_count < 2) {
        return expr_new_symbol(SYM_True);
    }
    
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        for (size_t j = i + 1; j < res->data.function.arg_count; j++) {
            if (expr_eq(res->data.function.args[i], res->data.function.args[j])) {
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
    for (size_t i = 0; i < res->data.function.arg_count - 1; i++) {
        Expr* a = res->data.function.args[i];
        Expr* b = res->data.function.args[i+1];

        bool equal = false;
        bool definitely_unequal = false;

        int av;
        if (expr_eq(a, b)) {
            equal = true;
        } else if ((av = assoc_equal_verdict(a, b)) != ASSOC_EQ_NA
                   || (av = list_equal_verdict(a, b)) != ASSOC_EQ_NA) {
            if (av == ASSOC_EQ_FALSE) return expr_new_symbol(SYM_False);
            if (av == ASSOC_EQ_TRUE) equal = true;
            else all_equal = false;
            continue;
        } else if (approximate_equal_verdict(a, b, &definitely_unequal)) {
            if (definitely_unequal) return expr_new_symbol(SYM_False);
            equal = true;
        } else {
            bool can_compare = false;
            int cmp = compare_numeric(a, b, &can_compare);
            if (can_compare) {
                if (cmp == 0) equal = true;
                else definitely_unequal = true;
            }
        }

        if (!equal) {
            if (definitely_unequal) {
                /* compare_numeric found a nonzero difference.  On EXACT operands
                 * (not raw machine floats) catastrophic cancellation over an
                 * algebraic tower — e.g. Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2],
                 * which is exactly 0 — collapses to a tiny-but-nonzero residual
                 * at machine precision, so the numeric verdict is unreliable.
                 * Confirm with the exact zero test before committing to
                 * inequality; only override to equal when it PROVES the
                 * difference zero. */
                if (!is_raw_data(a) || !is_raw_data(b)) {
                    Expr* diff = eval_and_free(
                        internal_subtract((Expr*[]){ expr_copy(a), expr_copy(b) }, 2));
                    ZeroTestResult zt = diff ? zero_test_decide(diff) : ZERO_TEST_UNKNOWN;
                    if (diff) expr_free(diff);
                    if (zt == ZERO_TEST_TRUE) { equal = true; }
                    else return expr_new_symbol(SYM_False);
                } else {
                    return expr_new_symbol(SYM_False);
                }
            } else if (is_raw_data(a) && is_raw_data(b)) {
                return expr_new_symbol(SYM_False);
            } else if (expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b)) {
                /* compare_numeric could not order the pair, yet both are
                 * closed-form numeric constants (e.g. I, 2 I Pi, a complex
                 * radical).  The exact zero-test decides equality soundly. */
                ZeroTestResult zt = numeric_pair_zero_test(a, b);
                if (zt == ZERO_TEST_TRUE) equal = true;
                else if (zt == ZERO_TEST_FALSE) return expr_new_symbol(SYM_False);
                /* ZERO_TEST_UNKNOWN: leave undecided (falls through to NULL). */
            }
            if (!equal) all_equal = false;
        }
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
 * Unequal returns False if any two arguments are numerically equal,
 * and True if all pairs are determined to be unequal.
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
            Expr* a = res->data.function.args[i];
            Expr* b = res->data.function.args[j];
            
            bool equal = false;
            bool definitely_unequal = false;

            if (is_nan_operand(a) || is_nan_operand(b)) {
                /* IEEE unordered: NaN != anything, itself included. Checked
                 * ahead of expr_eq, which reports two Indeterminates equal. */
                definitely_unequal = true;
            } else if (expr_eq(a, b)) {
                // Structural identity
                equal = true;
            } else if ((assoc_is_wellformed(a) && assoc_is_wellformed(b))
                       || (is_list_head(a) && is_list_head(b))) {
                /* Entry-/element-wise, as Equal decides it (assoc_equal_verdict,
                 * list_equal_verdict): {1, 2} != {1, 3} is True and
                 * {0.1 + 0.2} != {0.3} False. */
                int av = is_list_head(a) ? list_equal_verdict(a, b)
                                         : assoc_equal_verdict(a, b);
                if (av == ASSOC_EQ_TRUE) equal = true;
                else if (av == ASSOC_EQ_FALSE) definitely_unequal = true;
            } else if (approximate_equal_verdict(a, b, &definitely_unequal)) {
                equal = !definitely_unequal;
            } else {
                bool can_compare;
                int cmp = compare_numeric(a, b, &can_compare);
                if (can_compare) {
                    if (cmp == 0) {
                        equal = true;
                    } else {
                        definitely_unequal = true;
                    }
                } else if (is_raw_data(a) && is_raw_data(b)) {
                    // Different raw data types or values (and not equal numerically)
                    definitely_unequal = true;
                } else if (expr_is_numeric_quantity(a) && expr_is_numeric_quantity(b)) {
                    /* Unorderable but closed-form numeric constants (I, 2 I Pi,
                     * a complex radical): the exact zero-test decides. */
                    ZeroTestResult zt = numeric_pair_zero_test(a, b);
                    if (zt == ZERO_TEST_TRUE) equal = true;
                    else if (zt == ZERO_TEST_FALSE) definitely_unequal = true;
                }
            }
            
            if (equal) {
                return expr_new_symbol(SYM_False);
            }
            if (!definitely_unequal) {
                all_definitely_unequal = false;
            }
        }
    }
    
    if (all_definitely_unequal) {
        return expr_new_symbol(SYM_True);
    }
    
    return NULL;
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
        int cmp = compare_numeric(a, b, &can_compare);
        
        if (!can_compare) {
            return NULL; // Unevaluated
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
