#include "list_common.h"
#include "internal.h"

#include <math.h>       /* isnan, isinf -- C99, no feature-test macro */

bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Overflow;
}

bool is_listq(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

bool is_minus_infinity(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times &&
        e->data.function.arg_count == 2) {
        Expr* a1 = e->data.function.args[0];
        Expr* a2 = e->data.function.args[1];
        if (a1->type == EXPR_INTEGER && a1->data.integer == -1 && is_infinity(a2)) return true;
        if (a2->type == EXPR_INTEGER && a2->data.integer == -1 && is_infinity(a1)) return true;
    }
    return false;
}

Expr* make_minus_infinity(void) {
    Expr* args[2] = { expr_new_integer(-1), expr_new_symbol(SYM_Infinity) };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

bool is_real_numeric(Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) return true;
    if (is_rational(e, NULL, NULL)) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    Expr* re, *im;
    if (is_complex(e, &re, &im)) {
        if (im->type == EXPR_INTEGER && im->data.integer == 0) return true;
        if (im->type == EXPR_REAL && im->data.real == 0.0) return true;
        if (im->type == EXPR_BIGINT && mpz_sgn(im->data.bigint) == 0) return true;
#ifdef USE_MPFR
        if (im->type == EXPR_MPFR && mpfr_zero_p(im->data.mpfr)) return true;
#endif
    }
    return false;
}

/* ------------------------------------------------------------------------- *
 * Exact numeric ordering. See the block comment in list_common.h for why
 * none of expr_compare / is_real_numeric / expr_numeric_sign can serve here.
 * ------------------------------------------------------------------------- */

/* True for a real number the helpers below can order exactly.
 *
 * Deliberately NOT list_common.h's is_real_numeric: that routes rationals
 * through is_rational (arithmetic.c:105-117), which requires int64 numerator
 * AND denominator, so an ordinary exact input like 1/10^25 -- a
 * Rational[1, BigInt] -- is rejected and the whole call declines. is_rational_like
 * (arithmetic.c:125) is the bigint-aware predicate. */
bool list_real_number_q(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_REAL) return !isnan(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return !mpfr_nan_p(e->data.mpfr);
#endif
    return is_rational_like(e);
}

/* Sign of a real number: -1, 0, +1. Sets *ok false for anything this file
 * no sign can be read from, which callers turn into an unevaluated result.
 *
 * The *ok flag is the whole point. expr_numeric_sign (arithmetic.c) returns a
 * bare 0 both for "genuinely zero" and for "I do not recognise this", and it
 * recognises neither MPFR nor a bigint-component Rational -- so using it here
 * would silently report a tie for two distances it simply could not read. */
int list_numeric_sign(Expr* e, bool* ok) {
    *ok = true;
    if (!e) { *ok = false; return 0; }
    if (e->type == EXPR_INTEGER)
        return (e->data.integer > 0) - (e->data.integer < 0);
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL) {
        if (isnan(e->data.real)) { *ok = false; return 0; }
        return (e->data.real > 0.0) - (e->data.real < 0.0);
    }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        if (mpfr_nan_p(e->data.mpfr)) { *ok = false; return 0; }
        return mpfr_sgn(e->data.mpfr);
    }
#endif
    /* Rational[num, den], bigint components included. The denominator is
     * conventionally positive, but both signs are read so a hand-built
     * Rational with a negative denominator still orders correctly. */
    if (is_rational_like(e) && e->type == EXPR_FUNCTION) {
        bool ok_n = true, ok_d = true;
        int sn = list_numeric_sign(e->data.function.args[0], &ok_n);
        int sd = list_numeric_sign(e->data.function.args[1], &ok_d);
        if (!ok_n || !ok_d) { *ok = false; return 0; }
        return sn * sd;
    }
    *ok = false;
    return 0;
}

/* Lift an exactly-representable real into a canonical mpq. Handles Integer,
 * BigInt, bigint-component Rational, and -- the point of this helper -- Real,
 * via mpq_set_d, which is exact: a double IS a rational, m * 2^e.
 *
 * Shaped after mod_quot_expr_to_mpq (core.c:2569), which is static there and
 * has no Real case. */
static bool lc_to_mpq(const Expr* e, mpq_t out) {
    mpq_init(out);
    if (e->type == EXPR_INTEGER) { mpq_set_si(out, (long)e->data.integer, 1UL); return true; }
    if (e->type == EXPR_BIGINT)  { mpq_set_z(out, e->data.bigint); return true; }
    if (e->type == EXPR_REAL) {
        if (isnan(e->data.real) || isinf(e->data.real)) { mpq_clear(out); return false; }
        mpq_set_d(out, e->data.real);          /* exact */
        return true;
    }
    if (is_rational_like(e) && e->type == EXPR_FUNCTION) {
        mpz_t num, den;
        mpz_init(num); mpz_init(den);
        expr_to_mpz(e->data.function.args[0], num);
        expr_to_mpz(e->data.function.args[1], den);
        if (mpz_sgn(den) == 0) { mpz_clears(num, den, NULL); mpq_clear(out); return false; }
        mpq_set_num(out, num);
        mpq_set_den(out, den);
        mpq_canonicalize(out);
        mpz_clears(num, den, NULL);
        return true;
    }
    mpq_clear(out);
    return false;
}

/* Order two distances by NUMERIC VALUE: -1, 0, +1. Sets *ok false when the
 * comparison cannot be decided.
 *
 * NOT expr_compare, which is a canonical total order and is wrong here in both
 * directions. It breaks a value tie between different ExprTypes on the type
 * enum (sort.c:376), so Nearest[{0, 2.0}, 1] would drop the 2.0 -- distances 1
 * and 1.0 are equal in value but Integer sorts before Real -- defeating the
 * all-ties contract this file exists to honour. And for atoms that are not both
 * integer-like it falls back to comparing get_numeric_value() doubles
 * (sort.c:372-377), so two exact rationals differing below double resolution
 * compare equal and a strictly-farther element is reported as tied.
 *
 * Subtracting instead is exact where it must be and inexact only where the
 * input already was: 1 - 1.0 is 0.0 (a real tie, as Mathematica has it), while
 * 1/3 - (1/3 + 1/10^18) is the exact Rational[-1, 10^18]. */
int list_numeric_cmp(Expr* a, Expr* b, bool* ok) {
    *ok = true;
    /* Same-type fast paths for the two common cases, avoiding an allocation and
     * two evaluator passes per comparison. Both are exact for their type. */
    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER)
        return (a->data.integer > b->data.integer) - (a->data.integer < b->data.integer);
    if (a->type == EXPR_REAL && b->type == EXPR_REAL) {
        if (isnan(a->data.real) || isnan(b->data.real)) { *ok = false; return 0; }
        return (a->data.real > b->data.real) - (a->data.real < b->data.real);
    }
    /* MIXED EXACT / INEXACT: compare as exact rationals, never by subtracting.
     *
     * Subtraction here loses the exact operand's precision, because
     * internal_subtract widens the pair to a double. That is not merely
     * imprecise, it makes the comparator INTRANSITIVE, and an intransitive
     * comparator feeding a sort produces order-dependent output. With plain
     * int64 values:
     *
     *     a = 2^60, b = 2.0^60, c = 2^60 + 1
     *     a - b -> 0.0   so a == b
     *     c - b -> 0.0   so c == b
     *     a - c -> -1    so a <  c
     *
     * A double is exactly a rational, so lifting both sides with mpq and
     * comparing is exact and total. It also matches Mathematica on the two
     * cases that pin the semantics: 1 and 1.0 still tie (1.0 lifts to exactly
     * 1), while 0.1 and 1/10 no longer do -- and Mathematica's
     * Nearest[{0.9, 11/10}, 1] is {0.9}, not both. */
    bool a_exact = (a->type != EXPR_REAL);
    bool b_exact = (b->type != EXPR_REAL);
    if (a_exact != b_exact) {
        mpq_t qa, qb;
        if (lc_to_mpq(a, qa)) {
            if (lc_to_mpq(b, qb)) {
                int r = mpq_cmp(qa, qb);
                mpq_clear(qa); mpq_clear(qb);
                return (r > 0) - (r < 0);
            }
            mpq_clear(qa);
        }
        /* Not liftable -- an MPFR operand, or a non-finite Real. Fall through
         * to the subtraction path, which is what handled these before. */
    }

    Expr* sub_args[2] = { expr_copy(a), expr_copy(b) };
    Expr* diff = eval_and_free(internal_subtract(sub_args, 2));
    int s = list_numeric_sign(diff, ok);
    expr_free(diff);
    return s;
}
