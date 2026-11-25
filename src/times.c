#include "times.h"
#include "radcanon.h"
#include "attr.h"       /* call_head_attributes — a cleared Orderless */
#include "sort.h"       /* expr_orderless_sort — memoized Orderless canonicalisation */
#include "arithmetic.h"
#include "complex.h"
#include "interval.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "trig_canon.h"
#include "series.h"
#include "ndarray.h"
#include "checked_int.h"
#include "flint_qqbar.h"
#include "expand.h"       /* expr_expand: O4 exponent sums */
#include "symtab.h"       /* SymbolDef: the opaque-factor fast path */
#include "print.h"        /* expr_to_string_fullform: fast-path self-check */
#include "core.h"         /* numeric_certified_real_sign: c Infinity */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>    /* isfinite, fabs */
#include <float.h>   /* DBL_MIN */
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>


static bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Overflow;
}


/* True for positive numeric expressions: positive int/bigint, positive real,
 * or Rational[n, d] with positive numerator (denominators are conventionally
 * positive in Mathilda). Used to guard the radical-fusion rewrite
 * a^q * b^(-q) -> (a/b)^q, which is only valid on the principal branch when
 * both bases are strictly positive. */
static bool is_positive_numeric_expr(Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) > 0;
    if (e->type == EXPR_REAL)    return e->data.real > 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return n > 0;
    return false;
}

/* Sound (conservative) positivity oracle for the radical-fusion guard: returns
 * true only when e is provably strictly positive.  Extends
 * is_positive_numeric_expr to the positive real constants Pi, E and EulerGamma
 * and to Times/Power/Plus built from provably-positive operands, so that e.g.
 * Sqrt[Pi/2] Sqrt[2/Pi] fuses to Sqrt[1] = 1 (both radicands are positive, and
 * their product is a positive numeric).  Any operand it cannot prove positive
 * makes the whole expression fall through to false — the guard is never
 * loosened on the principal-branch validity of a^q b^q = (a b)^q. */
static bool is_positive_constant_expr(Expr* e) {
    if (!e) return false;
    if (is_positive_numeric_expr(e)) return true;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        return s == SYM_Pi || s == SYM_E || s == SYM_EulerGamma;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        Expr** args = e->data.function.args;
        size_t n = e->data.function.arg_count;
        if (h == SYM_Times || h == SYM_Plus) {
            /* product / sum of provably-positive operands is positive */
            for (size_t i = 0; i < n; i++)
                if (!is_positive_constant_expr(args[i])) return false;
            return n > 0;
        }
        if (h == SYM_Power && n == 2) {
            /* positive base raised to a real exponent stays positive */
            Expr* ex = args[1];
            bool real_exp = (ex->type == EXPR_INTEGER || ex->type == EXPR_REAL
                             || ex->type == EXPR_BIGINT || is_rational(ex, NULL, NULL));
#ifdef USE_MPFR
            real_exp = real_exp || ex->type == EXPR_MPFR;
#endif
            return real_exp && is_positive_constant_expr(args[0]);
        }
    }
    return false;
}

/* Sign of a real number atom: -1, 0 or 1 (0 also for a non-number). */
static int real_atom_sign(const Expr* c) {
    if (c->type == EXPR_INTEGER) return (c->data.integer > 0) - (c->data.integer < 0);
    if (c->type == EXPR_BIGINT) return mpz_sgn(c->data.bigint);
    if (c->type == EXPR_REAL) return (c->data.real > 0.0) - (c->data.real < 0.0);
    if (c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
        c->data.function.head->data.symbol.name == SYM_Rational && c->data.function.arg_count == 2)
        return real_atom_sign(c->data.function.args[0]);
    return 0;
}

static bool is_complex_atom(const Expr* c) {
    return c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
           c->data.function.head->data.symbol.name == SYM_Complex && c->data.function.arg_count == 2;
}

/* A term of a sum reads as negative, as the printer writes it with a minus:
 * a negative real leading coefficient, or a complex one whose real part is
 * negative, or zero with a negative imaginary part (-I t, (-1 + I) t). */
static bool term_looks_negative(const Expr* t) {
    if (t->type != EXPR_FUNCTION || t->data.function.head->type != EXPR_SYMBOL ||
        t->data.function.head->data.symbol.name != SYM_Times || t->data.function.arg_count < 2)
        return false;
    const Expr* c = t->data.function.args[0];
    if (is_complex_atom(c)) {
        int re = real_atom_sign(c->data.function.args[0]);
        return re < 0 || (re == 0 && real_atom_sign(c->data.function.args[1]) < 0);
    }
    return real_atom_sign(c) < 0;
}

/* The exponent sign a merged power is written without, as in the reference:
 *  - a product with a negative real leading coefficient (a, not -a), or a
 *    negative imaginary one (I t, not -I t); any other complex coefficient
 *    never flips: (2/3)^((-1 + I) t) and (2/3)^((1 - I) t) stay;
 *  - a sum whose LAST term does not read as negative (a - b, not b - a;
 *    1 - a, not -1 + a; -1 - Sqrt[3], not 1 + Sqrt[3]; 1 - I t, not
 *    -1 + I t);
 *  - a complex number with a negative imaginary part (I, not -I; -1 + I, not
 *    1 - I). */
static bool exponent_looks_negative(const Expr* e) {
    if (!e) return false;
    if (is_complex_atom(e)) return real_atom_sign(e->data.function.args[1]) < 0;
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count < 2) return false;
    const Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    if (h->data.symbol.name == SYM_Plus)
        return !term_looks_negative(e->data.function.args[e->data.function.arg_count - 1]);
    if (h->data.symbol.name != SYM_Times) return false;
    const Expr* c = e->data.function.args[0];
    if (is_complex_atom(c))
        return real_atom_sign(c->data.function.args[0]) == 0 &&
               real_atom_sign(c->data.function.args[1]) < 0;
    return real_atom_sign(c) < 0;
}

/* Does e contain a machine or arbitrary-precision real?  The reference leaves
 * 2^(2.5 t) 3^(-2.5 t) and 2^(1.5 + a) 3^(-1.5 - a) unmerged. */
static bool contains_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_inexact(e->data.function.args[i])) return true;
    return false;
}

/* A numeric base is positive when that is provable, or when its machine
 * value is a real clear of zero (Sinh[1], Log[2]): a numeric quantity has one
 * value, so its sign is decided by evaluating it. */
static bool numeric_base_positive(Expr* b) {
    if (is_positive_constant_expr(b)) return true;
    Expr* v = numericalize(b, numeric_machine_spec());
    bool pos = v && v->type == EXPR_REAL && v->data.real > 1e-12;
    if (v) expr_free(v);
    return pos;
}

/* ComplexInfinity absorbs every other factor of a product -- a number of
 * any phase, a symbol, any finite expression -- since multiplying an
 * infinity of undetermined direction by a nonzero quantity leaves it an
 * infinity of undetermined direction: I ComplexInfinity, x ComplexInfinity
 * and ComplexInfinity Log[5] are all ComplexInfinity (audit Times#3,
 * Divide#4, Coth#2, Csch#2, Tanh#1, Log#5, Pochhammer#2). Not when a factor
 * is an explicit zero (that is 0 ComplexInfinity, Indeterminate, decided by
 * the caller), nor beside the structured factors that carry their own
 * product rules (Interval, SeriesData, Overflow[]). */
static bool times_cinf_absorbs(Expr** args, size_t n) {
    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (is_complex_infinity_sym(a) || is_infinity_sym(a)) continue;
        if ((a->type == EXPR_INTEGER && a->data.integer == 0) ||
            (a->type == EXPR_REAL && a->data.real == 0.0))
            return false;
#ifdef USE_MPFR
        if (a->type == EXPR_MPFR && mpfr_zero_p(a->data.mpfr)) return false;
#endif
        if (is_interval(a) || is_overflow(a)) return false;
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_SeriesData)
            return false;
    }
    return true;
}

/* An exact or inexact zero part of a Complex. */
static bool times_is_zero_part(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer == 0) ||
           (e->type == EXPR_REAL && e->data.real == 0.0)
#ifdef USE_MPFR
           || (e->type == EXPR_MPFR && mpfr_zero_p(e->data.mpfr))
#endif
           ;
}

/* A directed infinity's direction is a unit number: Times[c, Infinity] with
 * a complex c is rewritten to Times[c/Abs[c], Infinity] (2 I Infinity is
 * I Infinity, so I Infinity + I Infinity is I Infinity). A direction on the
 * imaginary axis becomes the exact +-I whatever the precision of c, as a
 * real direction becomes the exact +-1. NULL when c is already a unit (the
 * canonical form) or not a complex number. */
static Expr* times_normalize_direction(Expr* c) {
    Expr *re, *im;
    if (!is_complex(c, &re, &im)) return NULL;
    if (times_is_zero_part(re)) {
        int s = expr_numeric_sign(im);
#ifdef USE_MPFR
        if (im->type == EXPR_MPFR) s = mpfr_sgn(im->data.mpfr);
#endif
        if (s == 0) return NULL;
        if (re->type == EXPR_INTEGER && im->type == EXPR_INTEGER &&
            (im->data.integer == 1 || im->data.integer == -1))
            return NULL;                                  /* already +-I */
        return make_complex(expr_new_integer(0), expr_new_integer(s));
    }
    /* Off the axes: c / |c|, unless |c| is already 1. */
    Expr* a2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(re), expr_new_integer(2) }, 2),
                   expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(im), expr_new_integer(2) }, 2) }, 2));
    bool unit = false;
    if (a2->type == EXPR_INTEGER) unit = a2->data.integer == 1;
    else if (a2->type == EXPR_REAL) unit = fabs(a2->data.real - 1.0) < 1e-12;
#ifdef USE_MPFR
    else if (a2->type == EXPR_MPFR) {
        double d = mpfr_get_d(a2->data.mpfr, MPFR_RNDN);
        unit = fabs(d - 1.0) < 1e-12;
    }
#endif
    else if (a2->type != EXPR_FUNCTION) unit = true;      /* undecidable: leave it */
    expr_free(a2);
    if (unit) return NULL;
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Sign), (Expr*[]){ expr_copy(c) }, 1));
}

static Expr* multiply_numbers(Expr* a, Expr* b) {
    if (is_overflow(a) || is_overflow(b)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
#ifdef USE_MPFR
    if (numeric_any_mpfr(a, b) && !is_complex(a, NULL, NULL) && !is_complex(b, NULL, NULL)) {
        /* An exact 0 times an arbitrary-precision number is the exact 0
         * (0 N[Pi, 20] is 0): the product is zero whatever the factor's
         * digits, so it has no precision to report. This is what keeps the
         * zero part of Complex[0, y] exact through a product, Complex[0, 1]
         * N[Pi, 20] being Complex[0, 3.14...`20]. A machine Real keeps the
         * machine rule (0 2.5 is 0.). A non-finite MPFR is not a number. */
        if (((a->type == EXPR_INTEGER && a->data.integer == 0 && b->type == EXPR_MPFR
              && mpfr_number_p(b->data.mpfr))
             || (b->type == EXPR_INTEGER && b->data.integer == 0 && a->type == EXPR_MPFR
                 && mpfr_number_p(a->data.mpfr))))
            return expr_new_integer(0);
        /* Likewise a machine 0. times a finite MPFR -- in practice a value
         * the machine/MPFR boundary policy keeps outside the double range,
         * since contagion lowers every in-range MPFR first -- is the machine
         * 0.: Exp[1000.] 0. was an MPFR zero, 0.e-15. */
        if ((a->type == EXPR_REAL && a->data.real == 0.0 && b->type == EXPR_MPFR
             && mpfr_number_p(b->data.mpfr))
            || (b->type == EXPR_REAL && b->data.real == 0.0 && a->type == EXPR_MPFR
                && mpfr_number_p(a->data.mpfr)))
            return expr_new_real(0.0);
        Expr* r = numeric_mpfr_mul(a, b, 0);
        if (r) return r;
    }
#endif
    if (a->type == EXPR_REAL || b->type == EXPR_REAL) {
        /* Correctly rounded operand values: an exact BigInt or Rational of
         * any width (mpz_get_d truncated, and a Rational[BigInt, BigInt]
         * read as 0). Anything else is not this branch's to fold. */
        double va, vb;
        if (!numeric_machine_value(a, &va) || !numeric_machine_value(b, &vb))
            return NULL;
        /* Overflow (10.^308 * 10), an exact operand past DBL_MAX, or an
         * underflow -- a zero or subnormal product of nonzero operands
         * (1.*^-200 * 1.*^-200 silently gave 0.) -- is recomputed at machine
         * precision in MPFR, whose exponent range is far wider: the policy
         * Power/Exp already follow (Exp[1000.]). */
        double r = va * vb;
        if (!isfinite(r) || numeric_machine_subnormal(r) || (r == 0.0 && va != 0.0 && vb != 0.0)) {
            Expr* p = numeric_machine_binop_promote(a, b, true);
            if (p) return p;
        }
        return expr_new_real(r);
    }
    /* Mixed BigInt fast path: at least one side is BigInt, and *both*
     * sides fit the (Integer | BigInt | Rational[int64, int64]) shape.
     * Rationals whose components have overflowed int64 fail the
     * is_rational check and would corrupt expr_to_mpz below — those
     * fall through to the generic GMP block at the bottom. */
    if ((a->type == EXPR_BIGINT || b->type == EXPR_BIGINT) &&
        (expr_is_integer_like(a) || is_rational(a, NULL, NULL)) &&
        (expr_is_integer_like(b) || is_rational(b, NULL, NULL))) {
        int64_t n1 = 1, d1 = 1, n2 = 1, d2 = 1;
        bool a_is_rat = is_rational(a, &n1, &d1);
        bool b_is_rat = is_rational(b, &n2, &d2);

        /* Borrow a BigInt operand's limbs rather than copying them (a
         * Rational operand's int64 numerator goes through scratch). */
        mpz_t as, bs, r;
        mpz_srcptr av, bv;
        if (a_is_rat) { mth_mpz_init_set_i64(as, n1); av = as; }
        else av = expr_mpz_borrow(a, as);
        if (b_is_rat) { mth_mpz_init_set_i64(bs, n2); bv = bs; }
        else bv = expr_mpz_borrow(b, bs);

        mpz_init(r);
        mpz_mul(r, av, bv);
        if (a_is_rat) mpz_clear(as); else expr_mpz_release(a, as);
        if (b_is_rat) mpz_clear(bs); else expr_mpz_release(b, bs);

        // Now handle the denominators!
        int64_t den = d1 * d2;
        if (den == 1)
            return expr_bigint_normalize(expr_new_bigint_take(r));
        
        // Compute GCD of r and den
        mpz_t m_den, m_gcd;
        mpz_inits(m_den, m_gcd, NULL);
        mth_mpz_set_i64(m_den, den);
        mpz_gcd(m_gcd, r, m_den);
        
        mpz_divexact(r, r, m_gcd);
        mpz_divexact(m_den, m_den, m_gcd);
        
        int64_t final_den = 0;
        (void)expr_mpz_get_int64(m_den, &final_den);   /* divides den: fits */
        mpz_clears(m_den, m_gcd, NULL);
        
        if (final_den == 1) {
            Expr* res = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clear(r);
            return res;
        }
        
        // Return Times[num, Power[den, -1]] if we can't form a Rational (which only holds int64_t)
        // Wait, if final_den is small, we can form Rational ? No, Rational holds 2 args. If num is BigInt, we can't use Rational.
        // Or can we? Rational usually holds integers.
        // Actually, if we return Times[BigInt, Power[final_den, -1]], it'll infinite loop.
        // Because multiply_numbers will be called again!
        // We MUST return something that won't trigger multiply_numbers.
        // Wait! In Mathilda, Rational CAN hold BigInts?! Let's allow it:
        Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clear(r);
        Expr* r_den = expr_new_integer(final_den);
        Expr* r_args[2] = { r_num, r_den };
        return expr_new_function(expr_new_symbol(SYM_Rational), r_args, 2);
    }

    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
        __int128_t res = (__int128_t)a->data.integer * b->data.integer;
        if (res > INT64_MAX || res < INT64_MIN) {
            mpz_t av, bv, r;
            expr_to_mpz(a, av);
            expr_to_mpz(b, bv);
            mpz_init(r);
            mpz_mul(r, av, bv);
            mpz_clear(av); mpz_clear(bv);
            Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clear(r);
            return result;
        }
        return expr_new_integer((int64_t)res);
    }
    int64_t n1, d1, n2, d2;
    if (is_rational(a, &n1, &d1) && is_rational(b, &n2, &d2)) {
        __int128_t num = (__int128_t)n1 * n2;
        __int128_t den = (__int128_t)d1 * d2;
        if (num > INT64_MAX || num < INT64_MIN || den > INT64_MAX || den < INT64_MIN) {
            /* Promote to BigInt rational to keep symbolic computations
             * exact when intermediates exceed 64 bits. */
            mpz_t mn1, md1, mn2, md2, mnum, mden, g;
            mth_mpz_init_set_i64(mn1, n1); mth_mpz_init_set_i64(md1, d1);
            mth_mpz_init_set_i64(mn2, n2); mth_mpz_init_set_i64(md2, d2);
            mpz_inits(mnum, mden, g, NULL);
            mpz_mul(mnum, mn1, mn2);
            mpz_mul(mden, md1, md2);
            if (mpz_sgn(mden) < 0) { mpz_neg(mden, mden); mpz_neg(mnum, mnum); }
            mpz_gcd(g, mnum, mden);
            if (mpz_sgn(g) != 0) { mpz_divexact(mnum, mnum, g); mpz_divexact(mden, mden, g); }
            Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(mnum));
            Expr* r_den = expr_bigint_normalize(expr_new_bigint_from_mpz(mden));
            mpz_clears(mn1, md1, mn2, md2, mnum, mden, g, NULL);
            if (r_den->type == EXPR_INTEGER && r_den->data.integer == 1) {
                expr_free(r_den);
                return r_num;
            }
            Expr* r_args[2] = { r_num, r_den };
            return expr_new_function(expr_new_symbol(SYM_Rational), r_args, 2);
        }
        return make_rational((int64_t)num, (int64_t)den);
    }

    /* Generic GMP rational fallback: handle any combination of
     * Integer / BigInt / Rational[Integer-or-BigInt, Integer-or-BigInt].
     * Without this, a Rational with BigInt components (which can arise
     * during polynomial arithmetic, resultant computation, etc.) cannot
     * be multiplied with another such value via the int64 fast paths
     * above, and we would return NULL — which the caller in
     * builtin_times then dereferences. */
    {
        bool a_ok = false, b_ok = false;
        mpz_t an, ad, bn, bd;
        if (expr_is_integer_like(a)) {
            expr_to_mpz(a, an);
            mpz_init_set_si(ad, 1);
            a_ok = true;
        } else if (a->type == EXPR_FUNCTION &&
                   a->data.function.head->type == EXPR_SYMBOL &&
                   a->data.function.head->data.symbol.name == SYM_Rational &&
                   a->data.function.arg_count == 2 &&
                   expr_is_integer_like(a->data.function.args[0]) &&
                   expr_is_integer_like(a->data.function.args[1])) {
            expr_to_mpz(a->data.function.args[0], an);
            expr_to_mpz(a->data.function.args[1], ad);
            a_ok = true;
        }
        if (a_ok) {
            if (expr_is_integer_like(b)) {
                expr_to_mpz(b, bn);
                mpz_init_set_si(bd, 1);
                b_ok = true;
            } else if (b->type == EXPR_FUNCTION &&
                       b->data.function.head->type == EXPR_SYMBOL &&
                       b->data.function.head->data.symbol.name == SYM_Rational &&
                       b->data.function.arg_count == 2 &&
                       expr_is_integer_like(b->data.function.args[0]) &&
                       expr_is_integer_like(b->data.function.args[1])) {
                expr_to_mpz(b->data.function.args[0], bn);
                expr_to_mpz(b->data.function.args[1], bd);
                b_ok = true;
            }
        }
        if (a_ok && b_ok) {
            mpz_t mnum, mden, g;
            mpz_inits(mnum, mden, g, NULL);
            mpz_mul(mnum, an, bn);
            mpz_mul(mden, ad, bd);
            if (mpz_sgn(mden) < 0) { mpz_neg(mden, mden); mpz_neg(mnum, mnum); }
            mpz_gcd(g, mnum, mden);
            if (mpz_sgn(g) != 0) {
                mpz_divexact(mnum, mnum, g);
                mpz_divexact(mden, mden, g);
            }
            Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(mnum));
            Expr* r_den = expr_bigint_normalize(expr_new_bigint_from_mpz(mden));
            mpz_clears(an, ad, bn, bd, mnum, mden, g, NULL);
            if (r_den->type == EXPR_INTEGER && r_den->data.integer == 1) {
                expr_free(r_den);
                return r_num;
            }
            Expr* r_args[2] = { r_num, r_den };
            return expr_new_function(expr_new_symbol(SYM_Rational), r_args, 2);
        }
        if (a_ok) { mpz_clear(an); mpz_clear(ad); }
        if (b_ok) { mpz_clear(bn); mpz_clear(bd); }
    }

    return NULL;
}

Expr* make_times(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* BasePower (the (base, exponent) pair) is shared with trig_canon.h, which
 * mutates the array we build here. */

/* ---- AlgebraicNumber field multiplication ------------------------------- */

static int times_is_algnum(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_AlgebraicNumber &&
           e->data.function.arg_count == 2;
}

static int times_is_rational_number(const Expr* e) {
    return e && (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
                 (e->type == EXPR_FUNCTION && e->data.function.head &&
                  e->data.function.head->type == EXPR_SYMBOL &&
                  e->data.function.head->data.symbol.name == SYM_Rational));
}

/* Combine same-generator AlgebraicNumber factors of a product via number-field
 * multiplication; when there is exactly one generator, fold a rational
 * coefficient into it.  Returns the rebuilt (evaluated) result, or NULL. */
static Expr* algnum_times_combine(Expr** args, size_t n) {
    int any = 0;
    for (size_t i = 0; i < n; i++) if (times_is_algnum(args[i])) { any = 1; break; }
    if (!any) return NULL;

    Expr**       gens   = malloc(sizeof(Expr*) * n);
    Expr**       accs   = malloc(sizeof(Expr*) * n);
    const Expr** rats   = malloc(sizeof(Expr*) * n);
    const Expr** others = malloc(sizeof(Expr*) * n);
    size_t ng = 0, nrat = 0, no = 0;
    int changed = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* a = args[i];
        if (times_is_algnum(a)) {
            size_t gi; int found = 0;
            for (gi = 0; gi < ng; gi++)
                if (expr_eq(gens[gi], a->data.function.args[0])) { found = 1; break; }
            if (!found) { gens[ng] = a->data.function.args[0]; accs[ng] = expr_copy(a); ng++; }
            else {
                Expr* s = flint_qqbar_algnum_mul(accs[gi], a);
                if (s) { expr_free(accs[gi]); accs[gi] = s; changed = 1; }
                else others[no++] = a;
            }
        } else if (times_is_rational_number(a)) {
            rats[nrat++] = a;
        } else {
            others[no++] = a;
        }
    }

    int scaled = 0;
    if (ng == 1 && nrat > 0 && times_is_algnum(accs[0])) {
        Expr* rprod;
        if (nrat == 1) rprod = expr_copy((Expr*)rats[0]);
        else {
            Expr** rc = malloc(sizeof(Expr*) * nrat);
            for (size_t k = 0; k < nrat; k++) rc[k] = expr_copy((Expr*)rats[k]);
            rprod = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), rc, nrat));
            free(rc);
        }
        Expr* s = flint_qqbar_algnum_scale_rational(accs[0], rprod);
        expr_free(rprod);
        if (s) { expr_free(accs[0]); accs[0] = s; scaled = 1; changed = 1; }
    }

    if (!changed) {
        for (size_t gi = 0; gi < ng; gi++) expr_free(accs[gi]);
        free(gens); free(accs); free(rats); free(others);
        return NULL;
    }

    size_t total = ng + no + (scaled ? 0 : nrat);
    Expr** terms = malloc(sizeof(Expr*) * (total ? total : 1));
    size_t t = 0;
    for (size_t gi = 0; gi < ng; gi++)  terms[t++] = accs[gi];
    for (size_t k = 0; k < no; k++)     terms[t++] = expr_copy((Expr*)others[k]);
    if (!scaled)
        for (size_t k = 0; k < nrat; k++) terms[t++] = expr_copy((Expr*)rats[k]);
    free(gens); free(accs); free(rats); free(others);

    Expr* result;
    if (t == 0)      result = expr_new_integer(1);
    else if (t == 1) result = terms[0];
    else             result = expr_new_function(expr_new_symbol(SYM_Times), terms, t);
    free(terms);
    return eval_and_free(result);
}

/* >0 while the radical-fusion pass probes whether a fused Power would split
 * back apart (see there). Single-threaded like the evaluator itself. */
static int g_radical_fusion_probe = 0;

/* ---- O3: a power of a product splits off its positive numeric factors ----
 *
 * In a product, the reference writes Power[Times[c, rest], e] (e not a
 * number) as c^e rest^e when c is a positive number or numeric quantity:
 * y (2 x)^b is 2^b x^b y, Pi^a (Pi x)^b is Pi^b x^(a + b), and a negative
 * number contributes its absolute value, (-2 x)^b x^a being 2^b (-x)^b x^a.
 * A base with no symbolic factor ((2 Pi)^b) splits only when another factor
 * of the product shares one of its bases, so the powers can combine:
 * Pi^a (2 Pi)^b is 2^b Pi^(a + b), while y (2 Pi)^b stays. A lone (2 x)^b is
 * not a product and stays too. */
static bool o3_is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return is_rational((Expr*)e, NULL, NULL) || is_complex((Expr*)e, NULL, NULL);
}

static int o3_real_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return (e->data.integer > 0) - (e->data.integer < 0);
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (e->type == EXPR_REAL) return (e->data.real > 0) - (e->data.real < 0);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_sgn(e->data.mpfr);
#endif
    /* Rational[n, d] with int64 or bignum parts (d > 0): the sign of n */
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational && e->data.function.arg_count == 2)
        return o3_real_sign(e->data.function.args[0]);
    return 0;
}

static const Expr* o3_base_of(const Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Power && e->data.function.arg_count == 2)
        return e->data.function.args[0];
    return e;
}

static bool o3_is_power_of_times(const Expr* a) {
    if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL ||
        a->data.function.head->data.symbol.name != SYM_Power || a->data.function.arg_count != 2)
        return false;
    const Expr* b = a->data.function.args[0];
    if (b->type != EXPR_FUNCTION || b->data.function.head->type != EXPR_SYMBOL ||
        b->data.function.head->data.symbol.name != SYM_Times) return false;
    return !o3_is_number(a->data.function.args[1]);
}

/* The split of factor args[k] = Power[Times[...], e], or NULL when it stays. */
static Expr* o3_split_one(Expr* const* args, size_t n, size_t k) {
    Expr* pw = args[k];
    Expr* B = pw->data.function.args[0];
    Expr* e = pw->data.function.args[1];
    size_t m = B->data.function.arg_count;
    bool symbolic = false;
    for (size_t j = 0; j < m; j++)
        if (!expr_is_numeric_quantity(B->data.function.args[j])) { symbolic = true; break; }
    if (!symbolic) {
        /* (B): only when a factor's base is another factor's base. */
        bool shared = false;
        for (size_t j = 0; j < m && !shared; j++) {
            const Expr* fb = o3_base_of(B->data.function.args[j]);
            for (size_t i = 0; i < n && !shared; i++)
                if (i != k && expr_eq((Expr*)o3_base_of(args[i]), (Expr*)fb)) shared = true;
        }
        if (!shared) return NULL;
    }
    Expr** out = malloc(sizeof(Expr*) * (m + 2));
    Expr** rest = malloc(sizeof(Expr*) * (m + 1));
    size_t no = 0, nr = 0;
    for (size_t j = 0; j < m; j++) {
        Expr* f = B->data.function.args[j];
        int s = o3_real_sign(f);
        if (s > 0 || (s == 0 && !o3_is_number(f) && expr_is_numeric_quantity(f) && numeric_base_positive(f))) {
            int64_t p, q;
            if (is_rational(f, &p, &q) && p == 1)      /* (1/q)^e is q^-e */
                out[no++] = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_new_integer(q),
                               expr_new_function(expr_new_symbol(SYM_Times),
                                   (Expr*[]){ expr_new_integer(-1), expr_copy(e) }, 2) }, 2);
            else
                out[no++] = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(f), expr_copy(e) }, 2);
        } else if (symbolic && is_complex(f, NULL, NULL) && o3_real_sign(f->data.function.args[0]) == 0 &&
                   f->data.function.args[0]->type == EXPR_INTEGER &&
                   !(f->data.function.args[1]->type == EXPR_INTEGER &&
                     (f->data.function.args[1]->data.integer == 1 || f->data.function.args[1]->data.integer == -1))) {
            /* 2 I in (2 I x)^b: 2 out, I (or -I) stays with x */
            Expr* im = f->data.function.args[1];
            int sg = o3_real_sign(im);
            Expr* absim = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                              (Expr*[]){ expr_new_integer(sg), expr_copy(im) }, 2));
            out[no++] = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ absim, expr_copy(e) }, 2);
            rest[nr++] = make_complex(expr_new_integer(0), expr_new_integer(sg));
        } else if (s < 0 && symbolic && !(f->type == EXPR_INTEGER && f->data.integer == -1)) {
            /* -2 in (-2 x)^b: 2 out, -1 stays with x */
            Expr* absf = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                             (Expr*[]){ expr_new_integer(-1), expr_copy(f) }, 2));
            out[no++] = expr_new_function(expr_new_symbol(SYM_Power),
                (Expr*[]){ absf, expr_copy(e) }, 2);
            rest[nr++] = expr_new_integer(-1);
        } else {
            rest[nr++] = expr_copy(f);
        }
    }
    if (no == 0) {
        for (size_t j = 0; j < nr; j++) expr_free(rest[j]);
        free(out); free(rest);
        return NULL;
    }
    if (nr == 1) {
        out[no++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ rest[0], expr_copy(e) }, 2);
    } else if (nr > 1) {
        Expr* rt = expr_new_function(expr_new_symbol(SYM_Times), rest, nr);
        out[no++] = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ rt, expr_copy(e) }, 2);
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), out, no);
    free(out); free(rest);
    return r;
}

/* ---- O4: numeric powers over a coprime basis -----------------------------
 *
 * The reference writes a product of numeric powers whose bases share factors
 * over a common coprime basis when an exponent is symbolic: 4^a 6^b is
 * 2^(2 a + b) 3^b, 8 2^a 4^b is 2^(3 + a + 2 b), 4^a/2 is 2^(-1 + 2 a),
 * 12^a 2^-a is 6^a, and (2/3)^a (9/4)^b is (2/3)^(a - 2 b) (the equal-exponent
 * merge below recombines bases whose exponents end up equal or opposite).
 * The basis comes from gcds alone (no integer factoring), so huge bases cost
 * only gcds: FeynCalc's FactorList2 check multiplies m^-7 (p/m^4)^a
 * (r/m^3)^b, which must collapse to p^a r^b.  Bases that share nothing stay
 * as they are (3 4^a, 2^a 3^b).
 *
 * A negative rational base takes part through its magnitude: (-n)^e is
 * (-1)^e n^e for every complex e (arg(-1) + arg(n) = Pi lies in (-Pi, Pi]),
 * so the sign rides along as (-1)^e while n^e joins the basis.  That sign
 * goes back onto a basis element whose exponent is still exactly e, else it
 * stays a separate (-1)^e: (-6)^a/3 is (-2)^a 3^(-1 + a), (-6)^a/4 is
 * (-3)^a 2^(-2 + a), (-2)^a/2 is (-1)^a 2^(-1 + a), and (-2)^a (-8)^b is
 * (-1)^(a + b) 2^(a + 3 b).  Two negative powers under one exponent cannot
 * merge their magnitudes ((-2)^a (-3)^a is not 6^a), so they are kept apart
 * when deciding whether anything is shared: (-2)^a (-3)^a stays, while
 * (-2)^a (-4)^a is (-1)^(2 a) 2^(3 a).  A coefficient equal to the base is
 * the plain same-base merge, (-3)^a (-3) being (-3)^(1 + a). */
typedef struct { mpz_t v; Expr* e; int src; } O4Atom;  /* v^e, v > 1; src: the
                                    negative-base arg it came from, or -1 */

static bool o4_rat_parts(const Expr* x, mpz_t num, mpz_t den) {
    if (expr_is_integer_like(x)) {
        mpz_t t; expr_to_mpz(x, t); mpz_set(num, t); mpz_clear(t);
        mpz_set_ui(den, 1);
        return mpz_sgn(num) != 0;
    }
    if (x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL &&
        x->data.function.head->data.symbol.name == SYM_Rational && x->data.function.arg_count == 2 &&
        expr_is_integer_like(x->data.function.args[0]) && expr_is_integer_like(x->data.function.args[1])) {
        mpz_t t;
        expr_to_mpz(x->data.function.args[0], t); mpz_set(num, t); mpz_clear(t);
        expr_to_mpz(x->data.function.args[1], t); mpz_set(den, t); mpz_clear(t);
        return true;
    }
    return false;
}

/* Is a a candidate power: a rational base other than 0 and +-1, an exponent
 * that is not an integer?  *sym when the exponent is not numeric at all, *neg
 * when the base is negative. */
static bool o4_power(const Expr* a, bool* sym, bool* neg) {
    if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL ||
        a->data.function.head->data.symbol.name != SYM_Power || a->data.function.arg_count != 2) return false;
    const Expr* b = a->data.function.args[0];
    const Expr* e = a->data.function.args[1];
    if (expr_is_integer_like(e)) return false;
    mpz_t n, d; mpz_init(n); mpz_init(d);
    bool ok = o4_rat_parts(b, n, d) &&
              !(mpz_cmpabs_ui(n, 1) == 0 && mpz_cmp_ui(d, 1) == 0);
    bool negative = ok && mpz_sgn(n) < 0;
    mpz_clear(n); mpz_clear(d);
    if (!ok) return false;
    *sym = !expr_is_numeric_like(e);
    *neg = negative;
    return true;
}

static void o4_add(O4Atom** at, size_t* na, size_t* cap, const mpz_t v, Expr* e, int src) {
    if (mpz_cmp_ui(v, 1) == 0) { expr_free(e); return; }
    if (*na == *cap) { *cap = *cap ? 2 * *cap : 8; *at = realloc(*at, *cap * sizeof(O4Atom)); }
    mpz_init_set((*at)[*na].v, v);
    (*at)[*na].e = e;
    (*at)[*na].src = src;
    (*na)++;
}

/* The same-base cases O4 leaves to the grouping pass, or merges itself.
 *
 * Sets *defer when two negative powers share a base, which grouping merges
 * ((-2)^a (-2)^b -> (-2)^(a + b)) before O4 should split them.
 *
 * Returns the merged product, evaluated, when the rational coefficient (all
 * number factors multiplied together, as they are before anything else)
 * equals the base of a negative power: (-3)^a (-3) is (-3)^(1 + a).  The
 * reference meets the factors in canonical order (numeric bases increasing),
 * so the coefficient merges only with the FIRST power whose base shares a
 * prime with it; when that is another power, the split comes first.  So (-4)^(2 + n) (-4) 2^(2 x) is
 * (-1)^(3 + n) 2^(6 + 2 n + 2 x) (merged, then split), while
 * (-8)^a (-1/2)^(2 a) (-1/2) is -(-1)^(3 a) 2^(-1 + a) (split as written),
 * and (-1/2)^a (-1/2) 4 has the coefficient -2, so no merge at all. */
static Expr* o4_negative_same_base(Expr* const* args, size_t n, bool* defer) {
    *defer = false;
    for (size_t i = 0; i < n; i++) {
        bool s, ng;
        if (!o4_power(args[i], &s, &ng) || !ng || !s) continue;
        for (size_t j = i + 1; j < n; j++) {
            bool s2, ng2;
            if (o4_power(args[j], &s2, &ng2) && ng2 &&
                expr_eq(args[j]->data.function.args[0], args[i]->data.function.args[0]))
                *defer = true;
        }
    }
    if (*defer) return NULL;

    /* the coefficient */
    mpq_t q, t; mpq_init(q); mpq_init(t);
    mpq_set_ui(q, 1, 1);
    mpz_t num, den; mpz_init(num); mpz_init(den);
    size_t ncoef = 0;
    for (size_t k = 0; k < n; k++) {
        bool s, ng;
        if (o4_power(args[k], &s, &ng) || !o4_rat_parts(args[k], num, den)) continue;
        mpq_set_num(t, num); mpq_set_den(t, den); mpq_canonicalize(t);
        mpq_mul(q, q, t);
        ncoef++;
    }
    Expr* r = NULL;
    if (ncoef > 0 && mpq_sgn(q) < 0) {
        /* |num q| * den q: the primes of the coefficient */
        mpz_t cq, bq, g;
        mpz_init(cq); mpz_init(bq); mpz_init(g);
        mpz_mul(cq, mpq_numref(q), mpq_denref(q));
        mpz_abs(cq, cq);
        /* the interacting power with the smallest base (the reference's
         * canonical order puts numeric bases in increasing order) */
        size_t first = n;
        mpq_t fb; mpq_init(fb);
        for (size_t i = 0; i < n; i++) {
            bool s, ng;
            if (!o4_power(args[i], &s, &ng) || !s) continue;
            o4_rat_parts(args[i]->data.function.args[0], num, den);
            mpz_mul(bq, num, den);
            mpz_abs(bq, bq);
            mpz_gcd(g, bq, cq);
            if (mpz_cmp_ui(g, 1) <= 0) continue;
            mpq_set_num(t, num); mpq_set_den(t, den); mpq_canonicalize(t);
            if (first == n || mpq_cmp(t, fb) < 0) { first = i; mpq_set(fb, t); }
        }
        if (first < n && mpq_equal(fb, q)) {
            size_t i = first;
            Expr** fa = malloc(n * sizeof(Expr*));
            size_t nn = 0;
            for (size_t k = 0; k < n; k++) {
                bool s3, ng3;
                if (k == i)
                    fa[nn++] = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ expr_copy(args[i]->data.function.args[0]),
                                   expr_new_function(expr_new_symbol(SYM_Plus),
                                       (Expr*[]){ expr_new_integer(1),
                                                  expr_copy(args[i]->data.function.args[1]) }, 2) }, 2);
                else if (o4_power(args[k], &s3, &ng3) || !o4_rat_parts(args[k], num, den))
                    fa[nn++] = expr_copy(args[k]);
            }
            Expr* p = nn == 1 ? fa[0] : expr_new_function(expr_new_symbol(SYM_Times), fa, nn);
            free(fa);
            r = eval_and_free(p);
        }
        mpq_clear(fb);
        mpz_clear(cq); mpz_clear(bq); mpz_clear(g);
    }
    mpq_clear(q); mpq_clear(t); mpz_clear(num); mpz_clear(den);
    return r;
}

static Expr* o4_coprime_basis(Expr* const* args, size_t n) {
    size_t nsym = 0;
    size_t nneg = 0;
    for (size_t i = 0; i < n; i++) {
        bool s, ng;
        if (o4_power(args[i], &s, &ng)) { if (s) nsym++; if (ng) nneg++; }
    }
    if (nsym == 0) return NULL;
    static int off = -1;                 /* MATHILDA_NO_O4: A/B kill switch */
    if (off < 0) off = getenv("MATHILDA_NO_O4") != NULL;
    if (off) return NULL;
    if (nneg) {
        bool defer;
        Expr* merged = o4_negative_same_base(args, n, &defer);
        if (merged) return merged;
        if (defer) return NULL;
    }

    /* atoms: v^e from every candidate power and from the rational coefficient */
    O4Atom* at = NULL; size_t na = 0, cap = 0;
    bool* used = calloc(n, sizeof(bool));
    int sign = 1;
    mpz_t num, den; mpz_init(num); mpz_init(den);
    for (size_t i = 0; i < n; i++) {
        bool s, ng;
        if (o4_power(args[i], &s, &ng)) {
            o4_rat_parts(args[i]->data.function.args[0], num, den);
            if (ng) mpz_neg(num, num);           /* (-n)^e: n^e here, (-1)^e below */
            int src = ng ? (int)i : -1;
            Expr* e = args[i]->data.function.args[1];
            o4_add(&at, &na, &cap, num, expr_copy(e), src);
            o4_add(&at, &na, &cap, den, eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(e) }, 2)), src);
            used[i] = true;
        } else if (o4_rat_parts(args[i], num, den)) {
            if (mpz_sgn(num) < 0) { sign = -sign; mpz_neg(num, num); }
            o4_add(&at, &na, &cap, num, expr_new_integer(1), -1);
            o4_add(&at, &na, &cap, den, expr_new_integer(-1), -1);
            used[i] = true;
        }
    }
    /* Merging, in the reference's order:
     *  1. powers of one base add their exponents (2 2^a is 2^(1 + a),
     *     2^a 2^a is 2^(2 a)); a negative power's magnitude joins a positive
     *     power only to cancel it ((-1/2)^n 2^n is (-1)^n), and two negative
     *     powers never meet here (one base was left to the grouping pass);
     *  2. values under one exponent multiply together (2^a 8^a is 16^a);
     *     atoms of two different negative powers never do: (-2)^a (-3)^a is
     *     not 6^a times any sign;
     *  3. the coprime basis is taken of the values that remain, so 6^x 4^x 9
     *     is 8^x 3^(2 + x) (24 and 9 refine to {3, 8}) and 2^a 8^a 4^b is
     *     4^(2 a + b), as in the reference.
     * Anything to do: a merge in step 1, or two entries after step 2 whose
     * values share a factor.  Either needs a symbolic exponent: 17 Sqrt[17]
     * would otherwise be rewritten 17^(3/2), which evaluates straight back
     * to 17 Sqrt[17] -- an endless loop through Times. */
    bool share = false;
    mpz_t g; mpz_init(g);
    for (size_t i = 0; i < na; i++) {                        /* step 1 */
        if (!at[i].e) continue;
        for (size_t j = i + 1; j < na; j++) {
            if (!at[j].e || mpz_cmp(at[i].v, at[j].v) != 0) continue;
            /* never two negative sources (each keeps its own sign); a
             * negative source's atom meets a positive one only when they
             * cancel exactly: (-1/2)^n 2^n is (-1)^n, while (-2)^a 2^a is
             * left to merge as (-4)^a */
            if (at[i].src >= 0 && at[j].src >= 0) continue;
            Expr* sum = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                            (Expr*[]){ expr_copy(at[i].e), expr_copy(at[j].e) }, 2));
            bool cancels = sum->type == EXPR_INTEGER && sum->data.integer == 0;
            if ((at[i].src >= 0 || at[j].src >= 0) && !cancels) { expr_free(sum); continue; }
            if (!expr_is_numeric_like(at[i].e) || !expr_is_numeric_like(at[j].e)) share = true;
            expr_free(at[i].e);
            expr_free(at[j].e);
            at[i].e = sum;
            if (at[j].src >= 0) at[i].src = at[j].src;
            at[j].e = NULL;
        }
    }
    size_t nm = 0;                                           /* step 2 */
    mpz_t* mv = malloc((na ? na : 1) * sizeof(mpz_t));
    Expr** me = malloc((na ? na : 1) * sizeof(Expr*));
    bool* ms = malloc((na ? na : 1) * sizeof(bool));
    int* msrc = malloc((na ? na : 1) * sizeof(int));
    for (size_t i = 0; i < na; i++) {
        if (!at[i].e || (at[i].e->type == EXPR_INTEGER && at[i].e->data.integer == 0)) continue;
        size_t k = 0;
        while (k < nm && !(expr_eq(me[k], at[i].e) &&
                           (at[i].src < 0 || msrc[k] < 0 || msrc[k] == at[i].src))) k++;
        if (k == nm) {
            mpz_init_set(mv[nm], at[i].v); me[nm] = at[i].e;
            ms[nm] = !expr_is_numeric_like(at[i].e); msrc[nm] = at[i].src; nm++;
        } else {
            mpz_mul(mv[k], mv[k], at[i].v);
            if (at[i].src >= 0) msrc[k] = at[i].src;
        }
    }
    for (size_t i = 0; i < nm && !share; i++)
        for (size_t j = i + 1; j < nm && !share; j++) {
            if (!ms[i] && !ms[j]) continue;
            mpz_gcd(g, mv[i], mv[j]);
            if (mpz_cmp_ui(g, 1) > 0) share = true;
        }
    Expr* out = NULL;
    if (share) {
        /* step 3: coprime basis of the merged values by gcd refinement */
        size_t nb = 0, bcap = nm * 2 + 4;
        mpz_t* B = malloc(bcap * sizeof(mpz_t));
        for (size_t i = 0; i < nm; i++) {
            bool dup = mpz_cmp_ui(mv[i], 1) == 0;
            for (size_t k = 0; k < nb && !dup; k++) dup = mpz_cmp(B[k], mv[i]) == 0;
            if (!dup) mpz_init_set(B[nb++], mv[i]);
        }
        bool again = true;
        while (again) {
            again = false;
            for (size_t i = 0; i < nb && !again; i++)
                for (size_t j = i + 1; j < nb && !again; j++) {
                    mpz_gcd(g, B[i], B[j]);
                    if (mpz_cmp_ui(g, 1) <= 0) continue;
                    mpz_t x, y;
                    mpz_init(x); mpz_init(y);
                    mpz_divexact(x, B[i], g);
                    mpz_divexact(y, B[j], g);
                    mpz_set(B[i], g);
                    mpz_set(B[j], x);
                    if (nb == bcap) { bcap *= 2; B = realloc(B, bcap * sizeof(mpz_t)); }
                    mpz_init_set(B[nb++], y);
                    mpz_clear(x); mpz_clear(y);
                    /* drop 1s and duplicates */
                    size_t w = 0;
                    for (size_t r = 0; r < nb; r++) {
                        bool drop = mpz_cmp_ui(B[r], 1) == 0;
                        for (size_t q = 0; q < w && !drop; q++) drop = mpz_cmp(B[q], B[r]) == 0;
                        if (drop) { mpz_clear(B[r]); continue; }
                        if (w != r) { mpz_init_set(B[w], B[r]); mpz_clear(B[r]); }
                        w++;
                    }
                    nb = w;
                    again = true;
                }
        }
        /* exponent of every basis element */
        Expr** X = malloc((nb ? nb : 1) * sizeof(Expr*));
        /* (exponents are expanded below: -7 + 3 (1 - 4 e) + 4 (1 + 3 e) is 0) */
        for (size_t k = 0; k < nb; k++) X[k] = expr_new_integer(0);
        bool ok = true;
        mpz_t r; mpz_init(r);
        for (size_t i = 0; i < nm && ok; i++) {
            mpz_set(r, mv[i]);
            for (size_t k = 0; k < nb; k++) {
                long c = 0;
                while (mpz_divisible_p(r, B[k])) { mpz_divexact(r, r, B[k]); c++; }
                if (c)
                    X[k] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                        (Expr*[]){ X[k], eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer((int64_t)c), expr_copy(me[i]) }, 2)) }, 2));
            }
            if (mpz_cmp_ui(r, 1) != 0) ok = false;
        }
        mpz_clear(r);
        for (size_t k = 0; k < nb; k++) { Expr* x = expr_expand(X[k]); expr_free(X[k]); X[k] = x; }
        if (ok) {
            size_t cap2 = 2 * n + nb + 1, nn = 0;
            Expr** fa = malloc(cap2 * sizeof(Expr*));
            if (sign < 0) fa[nn++] = expr_new_integer(-1);
            /* each negative power's sign: onto a basis element of its own
             * value whose exponent is still its own, else into one
             * (-1)^(e1 + e2 + ...).  Its own value: (-6)^a/3 is
             * (-2)^a 3^(-1 + a), (-6)^(2 a) 4^(2 a) 9 is (-8)^(2 a) 3^(2 + 2 a)
             * (6 and 4 merged under one exponent), but in
             * (-9)^(1 + a) (-12)^a 4 the 4^(1 + a) comes from the 4, so the
             * signs stay as (-1)^(1 + 2 a).  One power, not several: the
             * grouping pass would fold a -1 coefficient into a merge of (-1)
             * powers ((-1)^(1 + 3 a)), where the reference keeps -(-1)^(3 a). */
            bool* signed_b = calloc(nb ? nb : 1, sizeof(bool));
            Expr* sgn_e = NULL;
            for (size_t i = 0; i < n; i++) {
                bool s, ng;
                if (!used[i] || !o4_power(args[i], &s, &ng) || !ng) continue;
                Expr* e = args[i]->data.function.args[1];
                Expr* ee = expr_expand(e);
                size_t k = 0;
                for (; k < nb; k++) {
                    if (signed_b[k] || !expr_eq(X[k], ee)) continue;
                    bool own = false;
                    for (size_t m = 0; m < nm && !own; m++)
                        own = msrc[m] == (int)i && mpz_divisible_p(mv[m], B[k]);
                    if (own) break;
                }
                expr_free(ee);
                if (k < nb) signed_b[k] = true;
                else sgn_e = sgn_e ? eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                         (Expr*[]){ sgn_e, expr_copy(e) }, 2))
                                   : expr_copy(e);
            }
            if (sgn_e) fa[nn++] = expr_new_function(expr_new_symbol(SYM_Power),
                                     (Expr*[]){ expr_new_integer(-1), sgn_e }, 2);
            for (size_t k = 0; k < nb; k++) {
                if (X[k]->type == EXPR_INTEGER && X[k]->data.integer == 0) continue;
                if (signed_b[k]) mpz_neg(B[k], B[k]);
                fa[nn++] = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_bigint_normalize(expr_new_bigint_from_mpz(B[k])), expr_copy(X[k]) }, 2);
            }
            free(signed_b);
            for (size_t i = 0; i < n; i++) if (!used[i]) fa[nn++] = expr_copy(args[i]);
            out = nn == 0 ? expr_new_integer(1)
                          : expr_new_function(expr_new_symbol(SYM_Times), fa, nn);
            free(fa);
        }
        for (size_t k = 0; k < nb; k++) { mpz_clear(B[k]); expr_free(X[k]); }
        free(B); free(X);
    }
    for (size_t i = 0; i < nm; i++) mpz_clear(mv[i]);
    free(mv); free(me); free(ms); free(msrc);
    mpz_clear(g); mpz_clear(num); mpz_clear(den);
    for (size_t i = 0; i < na; i++) { mpz_clear(at[i].v); if (at[i].e) expr_free(at[i].e); }
    free(at); free(used);
    return out ? eval_and_free(out) : NULL;
}

/* ---- O5: exact rational radicals in canonical form ------------------------
 *
 * The number factors of a product and its powers r^e of rationals with
 * rational, non-integer exponents are rewritten together into the one
 * canonical form of radcanon.c: -1/3 Sqrt[3/2] is -1/Sqrt[6], 2 Sqrt[3/2]
 * is Sqrt[6], Sqrt[6]/3 is Sqrt[2/3], 6^(1/3)/2 is 3^(1/3)/2^(2/3),
 * 1000000007 Sqrt[3/1000000007] is Sqrt[3000000021] (gcds only, no
 * factoring), (3 + 6 I)/Sqrt[3] is (1 + 2 I) Sqrt[3].  Every other factor
 * passes through untouched.
 *
 * Taking part: every power of a positive rational (other than 1) to a
 * Rational exponent, and the exact numbers (Integer, Rational, and Complex
 * with exact parts).  A negative base takes part only when it is the only
 * one and its exponent lies in (0, 1) -- the shape Power leaves -- so its
 * sign can be placed (radcanon.c); otherwise negative powers pass through.
 * (-1)^e takes part the same way, as a sign with no magnitude:
 * (-1)^(1/3) 2^(1/3) is (-2)^(1/3), while 2 (-1)^(2/3) 2^(1/3) stays.  An inexact number declines the whole pass (the product
 * is numericised anyway).
 *
 * Returns the evaluated canonical product, or NULL when the product is
 * already canonical -- which it is after one rewrite, so the evaluation of
 * the result meets NULL here and the pass cannot cycle. */
static bool o5_exact_number(const Expr* x, mpq_t re, mpq_t im) {
    mpz_t n, d;
    mpz_init(n); mpz_init(d);
    bool ok = false;
    if (o4_rat_parts(x, n, d) || (expr_is_integer_like(x))) {
        if (expr_is_integer_like(x)) { mpz_t t; expr_to_mpz(x, t); mpz_set(n, t); mpz_clear(t); mpz_set_ui(d, 1); }
        mpq_set_num(re, n); mpq_set_den(re, d); mpq_canonicalize(re);
        mpq_set_ui(im, 0, 1);
        ok = true;
    } else if (x->type == EXPR_FUNCTION && x->data.function.head->type == EXPR_SYMBOL &&
               x->data.function.head->data.symbol.name == SYM_Complex &&
               x->data.function.arg_count == 2) {
        const Expr* a = x->data.function.args[0];
        const Expr* b = x->data.function.args[1];
        mpz_t n2, d2;
        mpz_init(n2); mpz_init(d2);
        bool oa = expr_is_integer_like(a) ? (expr_to_mpz(a, n), mpz_set_ui(d, 1), true) : o4_rat_parts(a, n, d);
        bool ob = expr_is_integer_like(b) ? (expr_to_mpz(b, n2), mpz_set_ui(d2, 1), true) : o4_rat_parts(b, n2, d2);
        if (oa && ob) {
            mpq_set_num(re, n); mpq_set_den(re, d); mpq_canonicalize(re);
            mpq_set_num(im, n2); mpq_set_den(im, d2); mpq_canonicalize(im);
            ok = true;
        }
        mpz_clear(n2); mpz_clear(d2);
    }
    mpz_clear(n); mpz_clear(d);
    return ok;
}

/* A power r^e of a rational r other than 0 and +-1 to a Rational exponent:
 * fills num, den (den > 0) and e. */
static bool o5_radical(const Expr* a, mpz_t num, mpz_t den, mpq_t e) {
    if (a->type != EXPR_FUNCTION || a->data.function.head->type != EXPR_SYMBOL ||
        a->data.function.head->data.symbol.name != SYM_Power || a->data.function.arg_count != 2)
        return false;
    const Expr* x = a->data.function.args[1];
    if (x->type != EXPR_FUNCTION || x->data.function.head->type != EXPR_SYMBOL ||
        x->data.function.head->data.symbol.name != SYM_Rational || x->data.function.arg_count != 2)
        return false;
    mpz_t en, ed;
    mpz_init(en); mpz_init(ed);
    bool ok = o4_rat_parts(x, en, ed) && mpz_cmp_ui(ed, 1) > 0;
    if (ok) { mpq_set_num(e, en); mpq_set_den(e, ed); mpq_canonicalize(e); }
    mpz_clear(en); mpz_clear(ed);
    if (!ok) return false;
    const Expr* b = a->data.function.args[0];
    if (expr_is_integer_like(b)) {
        mpz_t t; expr_to_mpz(b, t); mpz_set(num, t); mpz_clear(t);
        mpz_set_ui(den, 1);
    } else if (!o4_rat_parts(b, num, den)) {
        return false;
    }
    /* not 0 or 1; -1 takes part as a pure sign carrier */
    return mpz_sgn(num) != 0 && !(mpz_cmp_ui(num, 1) == 0 && mpz_cmp_ui(den, 1) == 0);
}

static Expr* o5_rational_radicals(Expr* const* args, size_t n) {
    static int off = -1;                 /* MATHILDA_NO_O5: A/B kill switch */
    if (off < 0) off = getenv("MATHILDA_NO_O5") != NULL;
    if (off) return NULL;
    /* cheap scan: anything to do at all? */
    size_t nrad = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = args[i];
        if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.head->data.symbol.name == SYM_Power && a->data.function.arg_count == 2 &&
            a->data.function.args[1]->type == EXPR_FUNCTION &&
            a->data.function.args[1]->data.function.head->type == EXPR_SYMBOL &&
            a->data.function.args[1]->data.function.head->data.symbol.name == SYM_Rational)
            nrad++;
    }
    if (nrad == 0) return NULL;

    mpq_t cre, cim, xre, xim, t1, t2, e;
    mpq_inits(cre, cim, xre, xim, t1, t2, e, NULL);
    mpq_set_ui(cre, 1, 1);
    mpz_t num, den;
    mpz_init(num); mpz_init(den);
    int* role = calloc(n, sizeof(int));  /* 0 pass, 1 number, 2 radical, 3 negative radical */
    size_t nnum = 0, npos = 0, nneg = 0, neg_i = 0;
    bool ok = true;
    for (size_t i = 0; i < n && ok; i++) {
        const Expr* a = args[i];
        if (a->type == EXPR_REAL
#ifdef USE_MPFR
            || a->type == EXPR_MPFR
#endif
            ) { ok = false; break; }
        if (o5_exact_number(a, xre, xim)) {
            /* (cre + i cim)(xre + i xim) */
            mpq_mul(t1, cre, xre); mpq_mul(t2, cim, xim); mpq_sub(t1, t1, t2);
            mpq_mul(t2, cre, xim); mpq_mul(xim, cim, xre); mpq_add(cim, t2, xim);
            mpq_set(cre, t1);
            role[i] = 1; nnum++;
        } else if (is_complex((Expr*)a, NULL, NULL)) {
            ok = false;                  /* an inexact complex number */
        } else if (o5_radical(a, num, den, e)) {
            if (mpz_sgn(num) > 0) { role[i] = 2; npos++; }
            else { role[i] = 3; nneg++; neg_i = i; }
        }
    }
    if (ok && mpq_sgn(cre) == 0 && mpq_sgn(cim) == 0) ok = false;   /* 0: Times owns it */
    /* Negative powers.  One with exponent in (0, 1) takes part as the sign
     * carrier.  Several take part only when their magnitudes share a factor
     * with some other factor of the product: then every sign splits off,
     * (-n)^e = (-1)^e n^e, and the signs combine, (-1)^a (-1)^b =
     * (-1)^(a + b), both exact on the principal branch; the combined sign
     * (-1)^k (-1)^E, k = floor(sum), E in [0, 1), leaves at most one
     * carrier.  So (-3)^(1/3) (-1)^(2/3) 3^(2/3) is -3 and (-2)^(1/3)
     * (-4)^(1/3) is 2 (-1)^(2/3), while (-2)^(1/3) (-3)^(1/3), with
     * nothing shared, stays -- as in the reference. */
    bool multi_neg = false;
    if (ok && nneg > 0) {
        bool take = false;
        if (nneg == 1) {
            o5_radical(args[neg_i], num, den, e);
            take = mpq_sgn(e) > 0 && mpz_cmp(mpq_numref(e), mpq_denref(e)) < 0;
        } else {
            /* does any negative power's magnitude share a factor with another
             * factor's magnitude (or with the coefficient's content)? */
            mpz_t g, m1, m2;
            mpz_init(g); mpz_init(m1); mpz_init(m2);
            mpz_t cc;
            mpz_init(cc);
            mpz_mul(cc, mpq_numref(cre), mpq_denref(cre));
            if (mpq_sgn(cim) != 0) { mpz_mul(cc, cc, mpq_numref(cim)); mpz_mul(cc, cc, mpq_denref(cim)); }
            mpz_abs(cc, cc);
            for (size_t i = 0; i < n && !take; i++) {
                if (role[i] != 3) continue;
                o5_radical(args[i], num, den, e);
                mpz_mul(m1, num, den);
                mpz_abs(m1, m1);
                if (mpz_cmp_ui(m1, 1) == 0) continue;
                mpz_gcd(g, m1, cc);
                if (mpz_cmp_ui(g, 1) > 0) take = true;
                for (size_t j = 0; j < n && !take; j++) {
                    if (j == i || (role[j] != 2 && role[j] != 3)) continue;
                    o5_radical(args[j], num, den, e);
                    mpz_mul(m2, num, den);
                    mpz_abs(m2, m2);
                    mpz_gcd(g, m1, m2);
                    if (mpz_cmp_ui(g, 1) > 0) take = true;
                }
            }
            mpz_clear(g); mpz_clear(m1); mpz_clear(m2); mpz_clear(cc);
            multi_neg = take;
        }
        if (!take) {
            for (size_t i = 0; i < n; i++) if (role[i] == 3) role[i] = 0;
            nneg = 0;
        }
    }
    size_t nrad_in = npos + nneg;
    /* Fast certificate (int64 only).  The grouping pass below merges equal
     * bases and folds an integer power into the coefficient, in this same
     * call; so first do just that here, and when what remains is already
     * canonical, step aside and let the ordinary pass finish -- a rewrite
     * would evaluate the product a second time, which is what made
     * number-field arithmetic (3 Sqrt[2] x 5 Sqrt[2] y, products of
     * monomials over Q(Sqrt[2])) twice as slow.  "Canonical" here: no
     * negative power takes part, every remaining exponent lies in (-1, 1)
     * with pairwise distinct magnitudes (nothing merges or carries), a
     * non-integer rational base has a positive exponent (Power's
     * orientation), the bases are pairwise coprime, and a prime they share
     * with the coefficient has the same sign there (nothing moves).
     * Anything else -- or any overflow --
     * takes the full canonicalisation. */
    if (ok && nneg == 0 && npos <= 8 && mpq_sgn(cim) == 0 &&
        mpz_fits_slong_p(mpq_numref(cre)) && mpz_fits_slong_p(mpq_denref(cre))) {
        int64_t cn = (int64_t)mpz_get_si(mpq_numref(cre)), cd = (int64_t)mpz_get_si(mpq_denref(cre));
        int64_t bn_[8], bd_[8], en_[8], ed_[8];
        size_t m = 0;
        bool canon = true;
        for (size_t i = 0; i < n && canon; i++) {
            if (role[i] != 2) continue;
            const Expr* a = args[i];
            int64_t bn, bd, en, ed;
            if (!is_rational(a->data.function.args[0], &bn, &bd) ||
                !is_rational(a->data.function.args[1], &en, &ed) || bn <= 0 ||
                en == INT64_MIN || ed > (INT64_C(1) << 30) || (en < 0 ? -en : en) > (INT64_C(1) << 30)) {
                canon = false; break;
            }
            size_t k = 0;
            while (k < m && !(bn_[k] == bn && bd_[k] == bd)) k++;
            if (k == m) { bn_[m] = bn; bd_[m] = bd; en_[m] = en; ed_[m] = ed; m++; }
            else {                                   /* same base: add exponents */
                int64_t nn = en_[k] * ed + en * ed_[k], dd = ed_[k] * ed;
                int64_t g = gcd(nn < 0 ? -nn : nn, dd);
                if (g > 1) { nn /= g; dd /= g; }
                if (dd > (INT64_C(1) << 30) || nn > (INT64_C(1) << 30) || nn < -(INT64_C(1) << 30)) { canon = false; break; }
                en_[k] = nn; ed_[k] = dd;
            }
        }
        /* integer exponents fold into the coefficient */
        size_t w = 0;
        for (size_t k = 0; k < m && canon; k++) {
            if (ed_[k] == 1) {
                int64_t e = en_[k];
                int64_t pn = bn_[k], pd = bd_[k];
                if (e < 0) { int64_t t = pn; pn = pd; pd = t; e = -e; }
                for (int64_t r = 0; r < e && canon; r++) {
                    if (ci_mul_i64(cn, pn, &cn) || ci_mul_i64(cd, pd, &cd)) canon = false;
                }
                if (canon) {
                    int64_t g = gcd(cn < 0 ? -cn : cn, cd);
                    if (g > 1) { cn /= g; cd /= g; }
                }
                continue;
            }
            if (en_[k] == 0) continue;
            bn_[w] = bn_[k]; bd_[w] = bd_[k]; en_[w] = en_[k]; ed_[w] = ed_[k]; w++;
        }
        uint64_t cc = (uint64_t)(cn < 0 ? -cn : cn), cdd = (uint64_t)cd;
        for (size_t k = 0; k < w && canon; k++) {
            int64_t ae = en_[k] < 0 ? -en_[k] : en_[k];
            if (ae >= ed_[k] || (bd_[k] > 1 && en_[k] < 0)) { canon = false; break; }
            /* a prime shared with the coefficient stays put when its integer
             * exponent there has the radical's sign (trunc(k + f) = k):
             * 2 Sqrt[2], 3 Sqrt[3/2] stay; Sqrt[2]/2, 2 Sqrt[3/2] move */
            bool pos = en_[k] > 0;
            if ((cc > 1 && ((!pos && gcd((int64_t)cc, bn_[k]) != 1) || (pos && gcd((int64_t)cc, bd_[k]) != 1))) ||
                (cdd > 1 && ((pos && gcd((int64_t)cdd, bn_[k]) != 1) || (!pos && gcd((int64_t)cdd, bd_[k]) != 1)))) {
                canon = false; break;
            }
            for (size_t j = 0; j < k && canon; j++) {
                int64_t aj = en_[j] < 0 ? -en_[j] : en_[j];
                if (aj == ae && ed_[j] == ed_[k]) canon = false;
                else if (gcd(bn_[j], bn_[k]) != 1 || gcd(bn_[j], bd_[k]) != 1 ||
                         gcd(bd_[j], bn_[k]) != 1 || gcd(bd_[j], bd_[k]) != 1)
                    canon = false;
            }
        }
        if (canon) ok = false;                       /* the grouping pass finishes it */
    }
    /* a lone radical with a unit coefficient is Power's canonical output */
    if (ok && nrad_in == 1) {
        bool unit;
        if (mpq_sgn(cim) == 0) unit = mpz_cmpabs_ui(mpq_numref(cre), 1) == 0 && mpz_cmp_ui(mpq_denref(cre), 1) == 0;
        else {
            mpz_t g, l;
            mpz_init(g); mpz_init(l);
            mpz_gcd(g, mpq_numref(cre), mpq_numref(cim));
            mpz_lcm(l, mpq_denref(cre), mpq_denref(cim));
            unit = mpz_cmp_ui(g, 1) == 0 && mpz_cmp_ui(l, 1) == 0;
            mpz_clear(g); mpz_clear(l);
        }
        if (unit) ok = false;
    }
    if (!ok || nrad_in == 0) {
        free(role);
        mpz_clear(num); mpz_clear(den);
        mpq_clears(cre, cim, xre, xim, t1, t2, e, NULL);
        return NULL;
    }

    RcAtoms at;
    rc_atoms_init(&at);
    mpq_t neg_e, ne;
    mpq_init(neg_e); mpq_init(ne);
    for (size_t i = 0; i < n; i++) {
        if (role[i] != 2 && role[i] != 3) continue;
        o5_radical(args[i], num, den, e);
        if (role[i] == 3) { mpz_neg(num, num); mpq_add(neg_e, neg_e, e); }
        rc_atoms_add(&at, num, e);
        mpq_neg(ne, e);
        rc_atoms_add(&at, den, ne);
    }
    bool has_carrier = nneg > 0;
    if (multi_neg) {
        /* (-1)^sum = (-1)^k (-1)^E: the integer part flips the coefficient */
        mpz_t k;
        mpz_init(k);
        mpz_fdiv_q(k, mpq_numref(neg_e), mpq_denref(neg_e));
        if (mpz_odd_p(k)) { mpq_neg(cre, cre); mpq_neg(cim, cim); }
        mpq_set_z(ne, k);
        mpq_sub(neg_e, neg_e, ne);
        mpz_clear(k);
        has_carrier = mpq_sgn(neg_e) != 0;
    }
    Expr** F = NULL;
    size_t nf = rc_canonical(&at, cre, cim, has_carrier ? (const mpq_t*)&neg_e : NULL, &F);
    rc_atoms_clear(&at);
    mpq_clear(neg_e); mpq_clear(ne);

    /* unchanged?  at most one number, equal to the new one, and the same
     * radicals as a multiset */
    bool same = nnum <= 1;
    /* F holds the number (first, when not 1) and Powers only; a Rational with
     * a bignum part is a number too, which o3_is_number's int64 test misses */
    size_t fnum = (nf > 0 && !(F[0]->type == EXPR_FUNCTION &&
                               F[0]->data.function.head->type == EXPR_SYMBOL &&
                               F[0]->data.function.head->data.symbol.name == SYM_Power)) ? 1 : 0;
    if (same) {
        const Expr* in_num = NULL;
        for (size_t i = 0; i < n; i++) if (role[i] == 1) in_num = args[i];
        if (fnum) same = in_num && expr_eq((Expr*)in_num, F[0]);
        else same = in_num == NULL;
    }
    if (same) same = nf - fnum == nrad_in;
    if (same) {
        bool* used = calloc(n, sizeof(bool));
        for (size_t k = fnum; k < nf && same; k++) {
            size_t i = 0;
            for (; i < n; i++)
                if ((role[i] == 2 || role[i] == 3) && !used[i] && expr_eq(args[i], F[k])) break;
            if (i == n) same = false;
            else used[i] = true;
        }
        free(used);
    }
    Expr* out = NULL;
    if (same) {
        for (size_t k = 0; k < nf; k++) expr_free(F[k]);
    } else {
        Expr** fa = malloc(sizeof(Expr*) * (nf + n + 1));
        size_t m = 0;
        for (size_t k = 0; k < nf; k++) fa[m++] = F[k];
        for (size_t i = 0; i < n; i++) if (role[i] == 0) fa[m++] = expr_copy(args[i]);
        out = m == 0 ? expr_new_integer(1)
            : m == 1 ? fa[0]
            : expr_new_function(expr_new_symbol(SYM_Times), fa, m);
        free(fa);
        out = eval_and_free(out);
    }
    free(F);
    free(role);
    mpz_clear(num); mpz_clear(den);
    mpq_clears(cre, cim, xre, xim, t1, t2, e, NULL);
    return out;
}

/* The product with every splittable power replaced by its split, evaluated;
 * NULL when nothing splits. */
static Expr* o3_split_powers(Expr* const* args, size_t n) {
    bool any = false;
    for (size_t i = 0; i < n && !any; i++) any = o3_is_power_of_times(args[i]);
    if (!any) return NULL;
    static int off = -1;                 /* MATHILDA_NO_O3: A/B kill switch */
    if (off < 0) off = getenv("MATHILDA_NO_O3") != NULL;
    if (off) return NULL;
    /* The split factors are spliced into ONE flat product: wrapped in an
     * inner Times they would evaluate on their own first, and the numeric-base
     * merge would rebuild (2 Pi)^b from 2^b Pi^b before Pi^b met Pi. */
    size_t cap = n + 8, nn = 0;
    Expr** na = malloc(sizeof(Expr*) * cap);
    bool changed = false;
    for (size_t i = 0; i < n; i++) {
        Expr* s = o3_is_power_of_times(args[i]) ? o3_split_one(args, n, i) : NULL;
        size_t add = s ? s->data.function.arg_count : 1;
        if (nn + add > cap) { cap = 2 * (nn + add); na = realloc(na, sizeof(Expr*) * cap); }
        if (s) {
            changed = true;
            for (size_t j = 0; j < s->data.function.arg_count; j++) {
                na[nn++] = s->data.function.args[j];
                s->data.function.args[j] = NULL;
            }
            s->data.function.arg_count = 0;
            expr_free(s);
        } else {
            na[nn++] = expr_copy(args[i]);
        }
    }
    if (!changed) {
        for (size_t i = 0; i < nn; i++) expr_free(na[i]);
        free(na);
        return NULL;
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_Times), na, nn);
    free(na);
    return eval_and_free(r);
}

/* An opaque factor for the fast path below: a non-System symbol, or a call
 * whose head is a non-System symbol with no builtin, no value and no
 * NumericFunction attribute -- g[x[1], x[2]], f[a], c$12.  None of the
 * collector's special passes (numbers, complex/interval/series/NDArray
 * contagion, Infinity, radical fusion, numeric-base merging, trig ratios,
 * AlgebraicNumber, O3/O4 power splitting) can act on such a factor with
 * exponent 1: all of them key on System heads, Power, or numeric bases. */
static bool times_opaque_factor(Expr* a) {
    Expr* h;
    if (a->type == EXPR_SYMBOL) h = a;
    else if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL)
        h = a->data.function.head;
    else return false;
    SymbolDef* d = h->data.symbol.def;
    if (!d) { d = symtab_get_def(h->data.symbol.name); h->data.symbol.def = d; }
    return d && !d->is_system && !d->builtin_func && !d->own_values &&
           !(d->attributes & ATTR_NUMERICFUNCTION);
}

static int g_times_fast_off = 0;   /* differential self-check: run the full path */

/* Product of machine-integer coefficients and DISTINCT opaque factors -- every
 * term of a polynomial in kernels such as the issue #112 Dirac trace's
 * -4 g[x[1], x[14]] ... g[x[13], x[2]].  The full collector costs ~2 us here,
 * almost all of it passes that cannot fire; the answer is simply the folded
 * coefficient followed by the canonically sorted factors.  Any other shape (a
 * repeated factor, a zero, an overflowing coefficient, a Power, a number that
 * is not a machine integer) returns NULL for the full path.
 * MATHILDA_TIMES_FAST_CHECK=1 runs both paths and aborts on a mismatch. */
static Expr* times_opaque_fast(Expr* res) {
    size_t n = res->data.function.arg_count;
    int64_t coef = 1;
    size_t nf = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        if (a->type == EXPR_INTEGER) {
            if (a->data.integer == 0 || ci_mul_i64(coef, a->data.integer, &coef)) return NULL;
        } else if (times_opaque_factor(a)) {
            nf++;
        } else {
            return NULL;
        }
    }
    if (nf == 0 || !(call_head_attributes(res) & ATTR_ORDERLESS)) return NULL;

    Expr* stack_buf[16];
    Expr** f = (nf + 1 <= 16) ? stack_buf : malloc(sizeof(Expr*) * (nf + 1));
    size_t lead = (coef != 1) ? 1 : 0;      /* the coefficient sorts first */
    if (lead) f[0] = expr_new_integer(coef);
    Expr** fac = f + lead;
    size_t j0 = 0;
    for (size_t i = 0; i < n; i++)
        if (res->data.function.args[i]->type != EXPR_INTEGER)
            fac[j0++] = expr_copy(res->data.function.args[i]);
    size_t k = lead + nf;
    expr_orderless_sort(fac, nf);
    for (size_t i = 1; i < nf; i++) {
        if (expr_eq(fac[i - 1], fac[i])) {          /* x x -> x^2: full path */
            for (size_t j = 0; j < k; j++) expr_free(f[j]);
            if (f != stack_buf) free(f);
            return NULL;
        }
    }
    Expr* out = (k == 1) ? f[0] : expr_new_function(expr_new_symbol(SYM_Times), f, k);
    if (f != stack_buf) free(f);
    /* Certify a fresh product as a fixed point (as builtin_plus does) when a
     * second evaluation could only re-derive it: every factor is evaluated and
     * no rule can fire on Times. */
    if (k > 1 && symtab_up_value_count == 0) {
        SymbolDef* td = symtab_get_def(SYM_Times);
        bool ok = td && !td->down_values && !td->sub_values && !td->own_values;
        for (size_t i = 0; ok && i < k; i++)
            ok = eval_node_reusable(out->data.function.args[i]);
        if (ok) eval_mark_fixed_point(out);
    }

    static int check = -1;
    if (check < 0) check = getenv("MATHILDA_TIMES_FAST_CHECK") ? 1 : 0;
    if (check) {
        g_times_fast_off++;               /* the whole check runs the full path */
        Expr* slow = builtin_times(res);
        Expr* slow_ev = slow ? evaluate(slow) : NULL;
        Expr* fast_ev = evaluate(out);
        g_times_fast_off--;
        if (!slow_ev || !expr_eq(slow_ev, fast_ev)) {
            char* a = expr_to_string_fullform(res);
            fprintf(stderr, "MATHILDA_TIMES_FAST_CHECK mismatch on %s\n", a);
            abort();
        }
        expr_free(slow); expr_free(slow_ev); expr_free(fast_ev);
    }
    return out;
}

/* A factor of a product the collector has just assembled that a second pass of
 * the collector would hand back unchanged: an exact number in the leading
 * coefficient slot, a symbol without a value, or an evaluated fixed point whose
 * head none of the special passes keys on.  Deliberately conservative -- the
 * shapes excluded here are exactly the ones the order-dependent passes work on:
 *   - inexact numbers (contagion), Complex, Interval, the infinities and I;
 *   - List / Association / SparseArray (Listable threading), Times (Flat),
 *     SeriesData, DirectedInfinity, Rational and Complex outside the
 *     coefficient slot;
 *   - a Power of a numeric quantity to a non-integer exponent, and any Power
 *     whose base is a number (radical canonicalisation, radical fusion,
 *     numeric-base merging, the O3/O4 splits). */
static bool times_certifiable_factor(Expr* a, bool lead) {
    switch (a->type) {
        case EXPR_INTEGER:
            return lead && a->data.integer != 0 && a->data.integer != 1;
        case EXPR_BIGINT:
            return lead;
        case EXPR_SYMBOL: {
            const char* s = a->data.symbol.name;
            if (s == SYM_I || s == SYM_Infinity || s == SYM_ComplexInfinity ||
                s == SYM_Indeterminate)
                return false;
            return eval_node_reusable(a);
        }
        case EXPR_FUNCTION:
            break;
        default:
            return false;
    }
    Expr* h = a->data.function.head;
    if (h->type != EXPR_SYMBOL) return false;
    const char* hn = h->data.symbol.name;
    if (hn == SYM_Rational) {
        /* the folded coefficient: canonical by construction */
        return lead && a->data.function.arg_count == 2 &&
               expr_is_integer_like(a->data.function.args[0]) &&
               expr_is_integer_like(a->data.function.args[1]);
    }
    if (hn == SYM_Times || hn == SYM_List || hn == SYM_Association ||
        hn == SYM_SparseArray || hn == SYM_SeriesData || hn == SYM_Complex ||
        hn == SYM_Interval || hn == SYM_DirectedInfinity || hn == SYM_Overflow)
        return false;
    if (hn == SYM_Power) {
        if (a->data.function.arg_count != 2) return false;
        Expr* b = a->data.function.args[0];
        Expr* e = a->data.function.args[1];
        if (b->type != EXPR_SYMBOL && b->type != EXPR_FUNCTION) return false;
        if (b->type == EXPR_FUNCTION && b->data.function.head->type == EXPR_SYMBOL &&
            (b->data.function.head->data.symbol.name == SYM_Rational ||
             b->data.function.head->data.symbol.name == SYM_Complex))
            return false;
        if (e->type != EXPR_INTEGER && expr_is_numeric_quantity(b)) return false;
    }
    return eval_node_reusable(a);
}

/* Stamp the collector's freshly built product as a fixed point when a second
 * pass could only re-derive it (see times_certifiable_factor), as builtin_plus
 * does for its sums.  Without this the evaluator's fixed-point loop re-runs the
 * whole collector on every product it returns, purely to find it unchanged --
 * the dominant per-entry cost of an exact matrix built in a Table body
 * (-Cos[t] Sin[t]^k, Sin[Pi k/2]/(Pi k), 1/(3 Pi)).  Requires an Orderless
 * Times with no rules of its own and no UpValues anywhere.
 * MATHILDA_TIMES_FAST_CHECK=1 re-runs the collector on each certified product
 * and aborts if it is not reproduced. */
static void times_certify_result(Expr* input, Expr* out, bool orderless) {
    if (!orderless || symtab_up_value_count != 0 || out->type != EXPR_FUNCTION) return;
    size_t k = out->data.function.arg_count;
    if (k < 2) return;
    SymbolDef* td = out->data.function.head->data.symbol.def;   /* cached on the node */
    if (!td) td = symtab_get_def(SYM_Times);
    if (!td || td->down_values || td->sub_values || td->own_values) return;
    /* Every factor but the folded coefficient must be one of the input
     * factors, untouched.  A factor the passes built -- a merged group, a
     * re-emitted trig ratio, a rewritten radical -- may meet another factor
     * only on a second pass: (a^-1)^(1/3) (a^-1)^(2/3) a collects to
     * a^-1 a, which only the next pass cancels. */
    size_t nin = input->data.function.arg_count;
    for (size_t i = 0; i < k; i++) {
        Expr* f = out->data.function.args[i];
        if (!times_certifiable_factor(f, i == 0)) return;
        if (f->type == EXPR_INTEGER || f->type == EXPR_BIGINT || is_rational(f, NULL, NULL) ||
            (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL &&
             f->data.function.head->data.symbol.name == SYM_Rational))
            continue;                     /* the coefficient (lead slot only) */
        bool original = false;
        for (size_t j = 0; j < nin && !original; j++)
            original = input->data.function.args[j] == f;
        if (!original) return;
    }

    static int check = -1;
    static int checking = 0;
    if (check < 0) check = getenv("MATHILDA_TIMES_FAST_CHECK") ? 1 : 0;
    if (check && !checking) {
        checking++;
        Expr* again = builtin_times(out);
        checking--;
        if (again && !expr_eq(again, out)) {
            char* a = expr_to_string_fullform(out);
            char* b = expr_to_string_fullform(again);
            char* c = expr_to_string_fullform(input);
            fprintf(stderr, "MATHILDA_TIMES_FAST_CHECK certified %s (from %s) re-collects to %s\n",
                    a, c, b);
            abort();
        }
        expr_free(again);
    }
    /* The folded coefficient is usually a fresh Rational: stamp it too (it
     * is GROUND), so the product can carry the INERT bit and survive the
     * next iterator binding instead of being re-collected. */
    Expr* lead = out->data.function.args[0];
    if (lead->type == EXPR_FUNCTION && !eval_node_reusable(lead))
        eval_mark_fixed_point(lead);
    eval_mark_fixed_point(out);
}

Expr* builtin_times(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_integer(1);
    if (n == 1) return expr_copy(res->data.function.args[0]);

    /* Machine-integer fast path: Times of all int64 args is just their product,
     * with no symbolic factors to collect. Skips the NDArray / SeriesData /
     * inexact-contagion / (base,exp) grouping scans below -- the hot case in
     * integer arithmetic (Case B's coefficient products, integer loop bodies).
     * A leading 0 short-circuits. On overflow we fall through to the generic
     * path, which promotes to bigint; correctness is unchanged. */
    {
        bool all_int = true;
        for (size_t i = 0; i < n; i++)
            if (res->data.function.args[i]->type != EXPR_INTEGER) { all_int = false; break; }
        if (all_int) {
            int64_t acc = 1;
            bool overflow = false;
            for (size_t i = 0; i < n; i++) {
                int64_t v = res->data.function.args[i]->data.integer;
                if (v == 0) { acc = 0; overflow = false; break; }
                if (ci_mul_i64(acc, v, &acc)) { overflow = true; break; }
            }
            if (!overflow) return expr_new_integer(acc);
        }
    }

    if (!g_times_fast_off) {
        Expr* fast = times_opaque_fast(res);
        if (fast) return fast;
    }

    /* AlgebraicNumber field multiplication: combine same-generator factors. Placed
     * after the int64 fast path so hot integer products never pay for the scan. */
    {
        Expr* combined = algnum_times_combine(res->data.function.args, n);
        if (combined) return combined;
    }

    /* Underflow[] factors (numeric.c): 3. Underflow[] = 0., Underflow[]^k
     * collapses, with Overflow[] it is Indeterminate. NULL unless one is here. */
    {
        Expr* u = numeric_underflow_times(res);
        if (u) return u;
        /* Overflow[] factors: numeric ones absorbed, symbolic ones kept
         * (Overflow[] x), an exact 0 gives 0. NULL once reduced, after
         * which Overflow[] is an ordinary factor below. */
        Expr* o = numeric_overflow_times(res);
        if (o) return o;
    }

    /* Fused special-case guard (mirrors builtin_plus). The NDArray, SeriesData,
     * inexact-contagion and Infinity/Indeterminate pre-scans below each traverse
     * the whole factor list, and for an ordinary symbolic or rational product --
     * the dominant case, e.g. every `c[k] x` term a symbolic Table builds -- all
     * find nothing. One detection pass replaces those traversals (and the size-n
     * contagion malloc): if no factor could trip any pre-scan we jump straight to
     * the (base, exponent) collector. The predicate is a provably-correct superset
     * of "some pass would act" -- NDArray needs is_ndarray, series needs
     * is_series_data, contagion mutates only when arg_is_inexact holds, and the
     * infinity block acts only on an Infinity/ComplexInfinity/Indeterminate factor
     * -- so skipping the passes when it is false is exactly equivalent. */
    {
        bool needs_prescan = false;
        for (size_t i = 0; i < n; i++) {
            Expr* a = res->data.function.args[i];
            if (is_ndarray(a) || is_series_data(a) || arg_is_inexact(a) ||
                is_indeterminate_sym(a) || is_complex_infinity_sym(a) || is_infinity_sym(a)) {
                needs_prescan = true;
                break;
            }
        }
        if (!needs_prescan) goto times_grouping;
    }

    /* NDArray fast path: same-shape NDArray operands multiply elementwise over
     * raw buffers, with numpy-style broadcasting of numeric scalars (3 * NDArray;
     * this is also how -NDArray and NDArray - NDArray reduce). If the array
     * operands disagree in shape, warn (NDArray::shape) and leave the product
     * unevaluated, like Dot::dotsh. A NULL (a symbolic operand, or no NDArray at
     * all) falls through to the generic path, which treats any NDArrays as
     * opaque non-numeric factors. */
    {
        /* Cheap guard: only walk the NDArray machinery (three arg scans) when
         * an NDArray operand is actually present. Pure-scalar Times — the
         * dominant case — skips it with a single pass. */
        bool any_nd = false;
        for (size_t i = 0; i < n; i++)
            if (is_ndarray(res->data.function.args[i])) { any_nd = true; break; }
        if (any_nd) {
            bool degenerate = false;
            Expr* fast = ndarray_elementwise(res->data.function.args, n, false, &degenerate);
            if (fast) return fast;
            /* A lane or scalar outside the machine range: the element path
             * keeps it at machine precision in MPFR (numeric-tower review C3). */
            if (degenerate) {
                Expr* r = ndarray_delist_and_reeval(res);
                if (r) return r;
            }
            Expr* retry = ndarray_int64_delist_retry(res);
            if (retry) return retry;
            if (ndarray_warn_shape_mismatch(res->data.function.args, n, "multiplied"))
                return NULL;
            /* NDArray combined with a symbolic factor (c * NDArray): purely
             * numeric, so it can't be multiplied elementwise. An invisible
             * packed List threads instead -- see the matching note in
             * src/plus.c -- and only a visible NDArray[...] declines. */
            Expr* symretry = ndarray_symbolic_delist_retry(res);
            if (symretry) return symretry;
            Expr* cretry = ndarray_numeric_constant_retry(res);   /* Pi NDArray */
            if (cretry) return cretry;
            ndarray_warn_symbolic(res->data.function.args, n, "multiplied");
        }
    }

    /* SeriesData arithmetic: if any factor is a power series, fold the whole
     * Times through the series algebra. Runs before inexact contagion, which
     * would otherwise apply N to the whole SeriesData (turning its integer
     * nmin/nmax/den into reals and defeating recognition); mixed-precision
     * coefficients still combine correctly because so_mul evaluates each
     * product. NULL (incompatible operands) falls through to the generic
     * collector, leaving a symbolic Times. */
    for (size_t i = 0; i < n; i++) {
        if (is_series_data(res->data.function.args[i])) {
            Expr* r = series_combine_times(res->data.function.args, n);
            if (r) return r;
            break;
        }
    }

    /* Inexact contagion: if any factor is an inexact Real/MPFR, numericalize
     * exact numeric parts in-place (Pi -> 3.14159, Sqrt[2] -> 1.41421, ...).
     * This is what turns `1. Pi` into `3.14159` instead of leaving it as a
     * frozen `1. Pi` Times. */
    {
        enum { CONTAGION_SMALL = 8 };
        Expr* numed_stack[CONTAGION_SMALL];
        Expr** numed = (n <= CONTAGION_SMALL) ? numed_stack
                                              : malloc(sizeof(Expr*) * n);
        if (numeric_contagion_args(res->data.function.args, n, numed)) {
            for (size_t i = 0; i < n; i++) {
                expr_free(res->data.function.args[i]);
                res->data.function.args[i] = numed[i];
            }
        }
        if (n > CONTAGION_SMALL) free(numed);
    }

    /* Infinity / Indeterminate preprocessing.
     *
     * Mathematica semantics for Times:
     *   Indeterminate * anything            -> Indeterminate
     *   0 * Infinity, 0 * ComplexInfinity   -> Indeterminate (with message)
     *   c * Infinity   (c numeric, c > 0)   -> Infinity
     *   c * Infinity   (c numeric, c < 0)   -> -Infinity (Times[-1, Infinity])
     *   c * ComplexInfinity (c != 0)        -> ComplexInfinity
     *   Infinity * Infinity                 -> Infinity
     *   Infinity * ComplexInfinity          -> ComplexInfinity
     *
     * If the product contains symbolic (non-numeric, non-Infinity) factors we
     * cannot decide the sign and fall through to the normal symbolic handler.
     */
    {
        bool has_indet = false;
        size_t inf_count = 0;
        size_t cinf_count = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* arg = res->data.function.args[i];
            if (is_indeterminate_sym(arg))      { has_indet = true; }
            else if (is_complex_infinity_sym(arg)) { cinf_count++; }
            else if (is_infinity_sym(arg))      { inf_count++; }
        }
        if (has_indet) return expr_new_symbol(SYM_Indeterminate);
        if (cinf_count > 0 && times_cinf_absorbs(res->data.function.args, n))
            return expr_new_symbol(SYM_ComplexInfinity);
        if (inf_count == 1 && cinf_count == 0 && n == 2) {
            Expr* c = is_infinity_sym(res->data.function.args[0]) ? res->data.function.args[1]
                                                                 : res->data.function.args[0];
            Expr* u = times_normalize_direction(c);
            if (u) return expr_new_function(expr_new_symbol(SYM_Times),
                                            (Expr*[]){ u, expr_new_symbol(SYM_Infinity) }, 2);
        }
        if (inf_count > 0 || cinf_count > 0) {
            Expr* coeff = expr_new_integer(1);
            bool has_symbolic = false;
            for (size_t i = 0; i < n; i++) {
                Expr* arg = res->data.function.args[i];
                if (is_infinity_sym(arg) || is_complex_infinity_sym(arg)) continue;
                if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL ||
                    arg->type == EXPR_BIGINT || is_rational(arg, NULL, NULL)) {
                    Expr* nn = multiply_numbers(coeff, arg);
                    expr_free(coeff);
                    coeff = nn;
                    if (!coeff) { coeff = expr_new_integer(1); has_symbolic = true; break; }
                    continue;
                }
                has_symbolic = true;
                break;
            }
            if (!has_symbolic) {
                int sign = expr_numeric_sign(coeff);
                bool is_zero = (coeff->type == EXPR_INTEGER && coeff->data.integer == 0) ||
                               (coeff->type == EXPR_REAL && coeff->data.real == 0.0);
                expr_free(coeff);
                if (is_zero) {
                    const char* what = (cinf_count > 0) ? "ComplexInfinity" : "Infinity";
                    mth_message_gated(g_arith_warnings_muted, "Infinity", "indet",
                            "Indeterminate expression 0 %s encountered.", what);
                    return expr_new_symbol(SYM_Indeterminate);
                }
                if (cinf_count > 0) return expr_new_symbol(SYM_ComplexInfinity);
                /* inf_count > 0, coeff != 0 */
                if (sign > 0) return expr_new_symbol(SYM_Infinity);
                return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
            }
            expr_free(coeff);
            /* A symbolic factor does not save 0 * Infinity: the collector below
             * would let the 0 absorb everything, ComplexInfinity included, so
             * 0 C[2] ComplexInfinity came out 0 -- an IVP fit at a basis
             * singularity then read 0 == 1 and reported "no solution". */
            for (size_t i = 0; i < n; i++) {
                Expr* arg = res->data.function.args[i];
                if ((arg->type == EXPR_INTEGER && arg->data.integer == 0) ||
                    (arg->type == EXPR_REAL && arg->data.real == 0.0)) {
                    const char* what = (cinf_count > 0) ? "ComplexInfinity" : "Infinity";
                    mth_message_gated(g_arith_warnings_muted, "Infinity", "indet",
                            "Indeterminate expression 0 %s encountered.", what);
                    return expr_new_symbol(SYM_Indeterminate);
                }
            }
            /* A real numeric-quantity coefficient decides the sign as a
             * number does: Infinity Pi and Infinity/Log[10] are Infinity,
             * -Infinity Sqrt[2] is -Infinity (the sign is certified, never
             * read off a rounded value). A non-real one stays. */
            if (cinf_count == 0 && inf_count == 1) {
                size_t m = 0;
                bool all_num = true;
                for (size_t i = 0; i < n && all_num; i++) {
                    Expr* arg = res->data.function.args[i];
                    if (is_infinity_sym(arg)) continue;
                    if (!expr_is_numeric_quantity(arg)) all_num = false;
                    m++;
                }
                if (all_num && m > 0) {
                    Expr** fa = malloc(sizeof(Expr*) * m);
                    size_t k = 0;
                    for (size_t i = 0; i < n; i++)
                        if (!is_infinity_sym(res->data.function.args[i]))
                            fa[k++] = expr_copy(res->data.function.args[i]);
                    Expr* prod = m == 1 ? fa[0]
                               : expr_new_function(expr_new_symbol(SYM_Times), fa, m);
                    free(fa);
                    int sg = 0;
                    int kind = numeric_certified_real_sign(prod, &sg);
                    expr_free(prod);
                    if (kind == 1 && sg > 0) return expr_new_symbol(SYM_Infinity);
                    if (kind == 1 && sg < 0)
                        return expr_new_function(expr_new_symbol(SYM_Times),
                            (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
                }
            }
            /* fall through: keep symbolic, let the rest run */
        }
    }

    times_grouping:;   /* fused special-case guard jumps here for an ordinary product */
    {
        Expr* split = o3_split_powers(res->data.function.args, n);
        if (split) return split;
        Expr* basis = o4_coprime_basis(res->data.function.args, n);
        if (basis) return basis;
        Expr* rad = o5_rational_radicals(res->data.function.args, n);
        if (rad) return rad;
    }
    Expr* num_prod = expr_new_integer(1);
    Expr* complex_val = NULL;
    Expr* interval_val = NULL;   /* running product of Interval factors */

    /* Collector buffers on the stack for the common small product (mirrors
     * builtin_plus): avoids three malloc/free pairs per Times in tight numeric
     * loops. n <= TIMES_SMALL_N bounds ht_cap <= TIMES_SMALL_HT below. */
    enum { TIMES_SMALL_N = 8, TIMES_SMALL_HT = 16 };
    /* Times is Orderless by default, but a user may ClearAttributes it. Then,
     * as in the reference, numbers still fold into one leading coefficient
     * but like bases combine only with their immediate predecessor
     * (b a b stays; b a^2 a is b a^3) and the factors are not sorted. One
     * attribute load; the default path is unchanged. */
    const bool orderless = (call_head_attributes(res) & ATTR_ORDERLESS) != 0;
    bool heap_bufs = (n > TIMES_SMALL_N);
    BasePower groups_stack[TIMES_SMALL_N];
    BasePower* groups = heap_bufs ? malloc(sizeof(BasePower) * n) : groups_stack;
    size_t group_count = 0;

    /* Open-addressing hash table mapping base-hash -> group index, so that
     * collecting like bases (Power exponents summed) is O(1) expected per
     * factor rather than a linear scan of every existing group. Without it a
     * product of many distinct factors (e.g. Product of n symbolic terms) is
     * O(n*groups) = O(n^2). Mirrors builtin_plus. Encounter order is
     * unchanged, so the order-dependent radical passes below are unaffected. */
    size_t ht_cap = 8;
    while (ht_cap < n * 2) ht_cap <<= 1;
    size_t ht_mask = ht_cap - 1;
    int64_t  slot_group_stack[TIMES_SMALL_HT];
    uint64_t slot_hash_stack[TIMES_SMALL_HT];
    int64_t* slot_group = heap_bufs ? malloc(sizeof(int64_t) * ht_cap) : slot_group_stack;
    uint64_t* slot_hash = heap_bufs ? malloc(sizeof(uint64_t) * ht_cap) : slot_hash_stack;
    for (size_t s = 0; s < ht_cap; s++) slot_group[s] = -1;

    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (expr_is_numeric_like(arg) && !is_complex(arg, NULL, NULL) && !is_interval(arg)) {
            Expr* next = multiply_numbers(num_prod, arg);
            if (next == NULL) {
                /* Could not fold this numeric factor into num_prod;
                 * stash it as its own group so we don't lose it (and so
                 * we don't propagate NULL to the type checks below). */
                Expr* nb = expr_copy(arg);
                uint64_t hb = expr_hash(nb);
                size_t slot = (size_t)hb & ht_mask;
                while (slot_group[slot] != -1) slot = (slot + 1) & ht_mask;
                slot_group[slot] = (int64_t)group_count;
                slot_hash[slot] = hb;
                groups[group_count].base = nb;
                groups[group_count].exponent = expr_new_integer(1);
                group_count++;
            } else {
                expr_free(num_prod);
                num_prod = next;
            }
        } else if (is_complex(arg, NULL, NULL) || (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_I)) {
            Expr* c_arg;
            if (arg->type == EXPR_SYMBOL) {
                Expr* z0 = expr_new_integer(0);
                Expr* z1 = expr_new_integer(1);
                c_arg = make_complex(z0, z1);
            } else {
                c_arg = expr_copy(arg);
            }
            if (!complex_val) complex_val = c_arg;
            else {
                Expr *re1, *im1, *re2, *im2;
                is_complex(complex_val, &re1, &im1); is_complex(c_arg, &re2, &im2);
                Expr* re = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                    expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(re1), expr_copy(re2)}, 2),
                    expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_new_integer(-1), expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(im1), expr_copy(im2)}, 2)}, 2)
                }, 2));
                Expr* im = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                    expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(re1), expr_copy(im2)}, 2),
                    expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(re2), expr_copy(im1)}, 2)
                }, 2));
                expr_free(complex_val); expr_free(c_arg);
                complex_val = make_complex(re, im);
            }
        } else if (is_interval(arg)) {
            /* Accumulate Interval factors: interval_val *= arg. If a product is
             * undecidable (symbolic endpoints), keep the factor as its own group
             * so nothing is lost. */
            if (!interval_val) {
                interval_val = expr_copy(arg);
            } else {
                Expr* prod = interval_multiply(interval_val, arg);
                if (prod) { expr_free(interval_val); interval_val = prod; }
                else {
                    uint64_t hb = expr_hash(arg);
                    size_t slot = (size_t)hb & ht_mask;
                    while (slot_group[slot] != -1) slot = (slot + 1) & ht_mask;
                    slot_group[slot] = (int64_t)group_count;
                    slot_hash[slot] = hb;
                    groups[group_count].base = expr_copy(arg);
                    groups[group_count].exponent = expr_new_integer(1);
                    group_count++;
                }
            }
        } else {
            Expr* base = arg; Expr* exponent;
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol.name == SYM_Power && arg->data.function.arg_count == 2) {
                base = arg->data.function.args[0]; exponent = expr_copy(arg->data.function.args[1]);
            } else { exponent = expr_new_integer(1); }
            
            uint64_t hb = 0;
            size_t slot = 0;
            int found = -1;
            if (orderless) {
                hb = expr_hash(base);
                slot = (size_t)hb & ht_mask;
                while (slot_group[slot] != -1) {
                    if (slot_hash[slot] == hb && expr_eq(groups[slot_group[slot]].base, base)) {
                        found = (int)slot_group[slot];
                        break;
                    }
                    slot = (slot + 1) & ht_mask;
                }
            } else if (group_count > 0 && expr_eq(groups[group_count - 1].base, base)) {
                found = (int)(group_count - 1);    /* adjacent like base only */
            }
            if (found != -1) {
                Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){groups[found].exponent, exponent}, 2));
                groups[found].exponent = new_exp;
            } else if (!orderless) {
                groups[group_count].base = expr_copy(base);
                groups[group_count].exponent = exponent;
                group_count++;
            } else {
                /* `slot` is the empty slot the probe stopped on. */
                slot_group[slot] = (int64_t)group_count;
                slot_hash[slot] = hb;
                groups[group_count].base = expr_copy(base);
                groups[group_count].exponent = (arg == base) ? exponent : exponent;
                group_count++;
            }
        }
    }
    if (heap_bufs) free(slot_group);
    if (heap_bufs) free(slot_hash);

    /* A factor equal to the base of another power joins it, as in the
     * reference: m^2 (m^2)^(-1 + D/2) is (m^2)^(D/2), Sqrt[m^2] m^2 is
     * (m^2)^(3/2).  The grouping above keys m^2 as (m, 2) and (m^2)^a as
     * (m^2, a), so the two never met.  Only exact equality merges
     * (m^4 Sqrt[m^2] stays), and only bases that are themselves powers are
     * examined, so ordinary products pay one scan. */
    if (group_count >= 2) {
        bool any_power_base = false;
        for (size_t gi = 0; gi < group_count && !any_power_base; gi++) {
            Expr* gb = groups[gi].base;
            any_power_base = gb->type == EXPR_FUNCTION && gb->data.function.head->type == EXPR_SYMBOL
                && gb->data.function.head->data.symbol.name == SYM_Power
                && gb->data.function.arg_count == 2;
        }
        if (any_power_base) {
            for (size_t j = 0; j < group_count; j++) {
                Expr* bj = groups[j].base;
                if (!(bj->type == EXPR_FUNCTION && bj->data.function.head->type == EXPR_SYMBOL
                      && bj->data.function.head->data.symbol.name == SYM_Power
                      && bj->data.function.arg_count == 2)) continue;
                for (size_t i = 0; i < group_count; i++) {
                    if (i == j) continue;
                    if (!expr_eq(groups[i].base, bj->data.function.args[0])
                        || !expr_eq(groups[i].exponent, bj->data.function.args[1])) continue;
                    /* groups[i] is the whole factor bj, to the first power */
                    Expr* ne = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),
                                   (Expr*[]){ groups[j].exponent, expr_new_integer(1) }, 2));
                    groups[j].exponent = ne;
                    expr_free(groups[i].base);
                    expr_free(groups[i].exponent);
                    for (size_t k = i + 1; k < group_count; k++) groups[k - 1] = groups[k];
                    group_count--;
                    if (i < j) j--;
                    break;
                }
            }
        }
    }

    if (num_prod->type == EXPR_INTEGER && num_prod->data.integer == 0) {
        if (complex_val) expr_free(complex_val);
        if (interval_val) expr_free(interval_val);   /* 0 * Interval -> 0 */
        for(size_t j=0; j<group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
        if (heap_bufs) free(groups);
        return num_prod;
    }
    /* An inexact zero absorbs symbolic factors, as in the reference:
     * 0. d is 0., 0. (d y) is 0.  (A factor that is an infinity keeps the
     * product, which the Indeterminate handling elsewhere decides.) */
    if (num_prod->type == EXPR_REAL && num_prod->data.real == 0.0
        && !complex_val && !interval_val && group_count > 0) {
        bool has_inf = false;
        for (size_t j = 0; j < group_count && !has_inf; j++) {
            Expr* gb = groups[j].base;
            has_inf = (gb->type == EXPR_SYMBOL
                       && (gb->data.symbol.name == SYM_Infinity
                           || gb->data.symbol.name == SYM_ComplexInfinity
                           || gb->data.symbol.name == SYM_Indeterminate))
                   || (gb->type == EXPR_FUNCTION && gb->data.function.head->type == EXPR_SYMBOL
                       && gb->data.function.head->data.symbol.name == SYM_DirectedInfinity);
        }
        if (!has_inf) {
            for (size_t j = 0; j < group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
            if (heap_bufs) free(groups);
            return num_prod;
        }
    }

    /* (The coefficient/radical exchange that used to live here --
     * Sqrt[2]/2 -> 1/Sqrt[2], 8/Sqrt[2] -> 4 Sqrt[2] -- is the O5 pass
     * above, which handles every rational radical, not only integer bases.) */

    /* Generalized radical fusion: collapse Power[a, e_i] * Power[b, e_j]
     * into a single Power when a, b are positive numerics and one
     * exponent is an integer multiple (positive or negative) of the
     * other. Reduce both exponents to lowest terms; if the denominators
     * match (q_i == q_j) and one numerator divides the other, the
     * integer ratio k is well defined. We then collapse to
     *   Power[a_pri * a_sec^k, e_pri]
     * where pri is whichever group has the smaller-magnitude numerator
     * (so the residual Power displays the shortest exponent). Subsumes:
     *   k = -1: Sqrt[6]/Sqrt[2]  -> Sqrt[3]            (basic fusion)
     *           Sqrt[2]/Sqrt[3]  -> 1/Sqrt[3/2]
     *   k = +1: Sqrt[2]*Sqrt[3]  -> Sqrt[6]            (same-exp collapse)
     *           2^(1/3)*3^(1/3)  -> 6^(1/3)
     *   k = -2: 12^(1/3) * 2^(-2/3)   -> 3^(1/3)
     *   k = +2: 2^(1/3) * 8^(2/3)     -> 128^(1/3) -> 4 * 2^(1/3)
     * Restricted to positive bases for principal-branch correctness;
     * q > 1 so we don't compete with ordinary arithmetic on integer
     * exponents (which the base-grouping pass has already merged). */
    for (size_t i = 0; i < group_count; i++) {
        if (!is_positive_constant_expr(groups[i].base)) continue;
        int64_t pi_n, pi_d;
        if (!is_rational(groups[i].exponent, &pi_n, &pi_d)) continue;
        if (pi_d == 1) continue;
        if (pi_n == 0) continue;
        /* INT64_MIN numerators: INT64_MIN % -1 traps (SIGFPE) and llabs
         * overflows below. This fusion is optional -- skip them. */
        if (pi_n == INT64_MIN) continue;

        for (size_t j = i + 1; j < group_count; j++) {
            if (!is_positive_constant_expr(groups[j].base)) continue;
            /* two rational bases are O5's (radcanon.c): already canonical */
            {
                mpz_t a1, a2;
                mpz_init(a1); mpz_init(a2);
                bool both = (expr_is_integer_like(groups[i].base) || o4_rat_parts(groups[i].base, a1, a2)) &&
                            (expr_is_integer_like(groups[j].base) || o4_rat_parts(groups[j].base, a1, a2));
                mpz_clear(a1); mpz_clear(a2);
                if (both) continue;
            }
            int64_t pj_n, pj_d;
            if (!is_rational(groups[j].exponent, &pj_n, &pj_d)) continue;
            if (pj_d != pi_d) continue;
            if (pj_n == 0 || pj_n == INT64_MIN) continue;

            /* Determine integer ratio. Default: i is primary (k = pj/pi). */
            bool i_primary;
            int64_t k;
            if (pj_n % pi_n == 0) {
                i_primary = true;
                k = pj_n / pi_n;
            } else if (pi_n % pj_n == 0) {
                i_primary = false;
                k = pi_n / pj_n;
            } else {
                continue;
            }

            /* Tie-break when |numerators| are equal: prefer the
             * positive-exponent primary so the printed Power has a
             * positive exponent (matches existing fusion convention). */
            if (llabs(pi_n) == llabs(pj_n)) {
                if (pi_n > 0)        i_primary = true;
                else if (pj_n > 0)   i_primary = false;
                k = i_primary ? (pj_n / pi_n) : (pi_n / pj_n);
            }

            /* Don't fuse coprime positive-integer bases when |k| > 1.
             * Fusion in that case produces a strictly larger combined
             * base (a * b^|k| with no shared prime factor to cancel),
             * inverting the canonical split form. For example
             * Power[2, 1/3] * Power[3, 2/3] (k=2, gcd(2,3)=1) should
             * stay split as Mathematica's canonical 2^(1/3) 3^(2/3),
             * not collapse to 18^(1/3). The k = ±1 case is exempt --
             * that's the basic "same exponent" combination
             * (Sqrt[2]*Sqrt[3] -> Sqrt[6], 2^(1/3)*3^(1/3) -> 6^(1/3))
             * which adds information by sharing the exponent. Cases
             * with gcd > 1 (e.g. 12^(1/3) * 2^(-2/3) -> 3^(1/3)) still
             * fuse because the combined base reduces via the common
             * prime factor. */
            if (k != 1 && k != -1
                && groups[i].base->type == EXPR_INTEGER
                && groups[j].base->type == EXPR_INTEGER
                && groups[i].base->data.integer > 0
                && groups[j].base->data.integer > 0
                && gcd(groups[i].base->data.integer,
                       groups[j].base->data.integer) == 1) {
                continue;
            }
            size_t pri = i_primary ? i : j;
            size_t sec = i_primary ? j : i;

            /* new_base = pri.base * sec.base^k. */
            Expr* sec_pow_k = (k == 1)
                ? expr_copy(groups[sec].base)
                : eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(groups[sec].base), expr_new_integer(k) }, 2));
            Expr* new_base = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(groups[pri].base), sec_pow_k }, 2));

            /* The fused base must remain a positive numeric we can wrap
             * in one Power. If evaluation produced anything else
             * (shouldn't on positive numerics, but be safe) abandon the
             * fusion and leave the pair intact. */
            if (!is_positive_numeric_expr(new_base)) {
                expr_free(new_base);
                continue;
            }
            /* A fusion that Power would immediately undo is a cycle, not a
             * simplification: Power[b, e] splits a rational base carrying a
             * perfect power (b = 50/32487^2 -> 50^(1/3) / 32487^(2/3)) back into
             * the very factors this pass fused, and Times fuses them again --
             * 2^(1/3) 5^(2/3) 32487^(-2/3) recursed to $RecursionLimit. Fuse
             * only when the fused Power stays one Power (or a number). */
            if (k != 1 && k != -1) {
                /* Inside the probe the |k| > 1 fusion is off (the probe's own
                 * Times must not re-enter this test, which is the cycle). */
                if (g_radical_fusion_probe) { expr_free(new_base); continue; }
                g_radical_fusion_probe++;
                Expr* probe = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ expr_copy(new_base), expr_copy(groups[pri].exponent) }, 2));
                g_radical_fusion_probe--;
                bool splits = probe && probe->type == EXPR_FUNCTION &&
                              probe->data.function.head->type == EXPR_SYMBOL &&
                              probe->data.function.head->data.symbol.name == SYM_Times;
                expr_free(probe);
                if (splits) { expr_free(new_base); continue; }
            }
            Expr* new_exp = expr_copy(groups[pri].exponent);

            /* Replace at index i (smaller of the two), drop j. */
            expr_free(groups[i].base);
            expr_free(groups[i].exponent);
            expr_free(groups[j].base);
            expr_free(groups[j].exponent);
            groups[i].base = new_base;
            groups[i].exponent = new_exp;
            for (size_t s = j; s + 1 < group_count; s++) groups[s] = groups[s + 1];
            group_count--;
            /* Restart -- the fused base may pair with an earlier group.
             * Each fusion strictly decreases group_count so the overall
             * loop terminates in O(group_count) restarts. */
            i = (size_t)-1;
            break;
        }
    }

    /* Numeric bases under one symbolic exponent merge, as in the reference:
     *   2^a 3^a -> 6^a,   3^-a 5^a -> (5/3)^a,   2^a 3^a 5^-a -> (6/5)^a,
     *   (-2)^a 3^a -> (-6)^a,   Sinh[1]^-a Cosh[1]^a -> (Cosh[1]/Sinh[1])^a.
     * (b1 b2)^e = b1^e b2^e holds whenever at most one of the bases is not a
     * positive real, so a merged group admits at most one base this cannot
     * prove positive ((-2)^a (-3)^a stays).  Only a positive base may be
     * inverted to flip its exponent's sign ((1/z)^a is not z^-a for z < 0),
     * so the group keeps the exponent of its unsigned base, if any; a group
     * whose exponents all agree keeps that one (3^-a 5^-a -> 15^-a); a mixed
     * group takes the sign without an explicit negative factor (a, not -a).
     * Rational exponents are the radical fusion's business above; any other
     * numeric exponent merges the same way (2^Sqrt[2] 3^-Sqrt[2] is
     * (2/3)^Sqrt[2], 2^Pi 3^Pi is 6^Pi); bases that
     * are not numeric quantities (x^a y^a, f[1]^a) never merge.  A power of a
     * product stays a power, so the merge cannot cycle. */
    for (size_t i = 0; i < group_count; i++) {
        Expr* ei = groups[i].exponent;
        if (expr_is_integer_like(ei) || contains_inexact(ei) ||
            (ei->type == EXPR_FUNCTION && ei->data.function.head->type == EXPR_SYMBOL &&
             ei->data.function.head->data.symbol.name == SYM_Rational) ||
            !expr_is_numeric_quantity(groups[i].base))
            continue;
        Expr* neg_ei = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(-1), expr_copy(ei) }, 2));
        /* Members: group indices, the first being i; flip[k] when the
         * member's exponent is -ei. */
        size_t* mem = malloc(sizeof(size_t) * group_count);
        bool* flip = malloc(sizeof(bool) * group_count);
        size_t nmem = 0;
        int unsigned_k = -1;
        mem[nmem] = i; flip[nmem] = false;
        if (!numeric_base_positive(groups[i].base)) unsigned_k = 0;
        nmem++;
        for (size_t j = i + 1; j < group_count; j++) {
            if (!expr_is_numeric_quantity(groups[j].base)) continue;
            bool same = expr_eq(groups[j].exponent, ei);
            if (!same && !expr_eq(groups[j].exponent, neg_ei)) continue;
            if (!numeric_base_positive(groups[j].base)) {
                if (unsigned_k >= 0) continue;
                unsigned_k = (int)nmem;
            }
            mem[nmem] = j; flip[nmem] = !same;
            nmem++;
        }
        if (nmem > 1) {
            /* Target sign: the unsigned member's, else the non-negative-looking one. */
            bool mixed = false;
            for (size_t k = 1; k < nmem; k++) mixed = mixed || flip[k];
            bool target_flip;
            if (unsigned_k >= 0) target_flip = flip[unsigned_k];
            else target_flip = mixed && exponent_looks_negative(ei);
            Expr** fac = malloc(sizeof(Expr*) * nmem);
            for (size_t k = 0; k < nmem; k++) {
                Expr* b = expr_copy(groups[mem[k]].base);
                fac[k] = flip[k] == target_flip ? b
                    : expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ b, expr_new_integer(-1) }, 2);
            }
            Expr* merged = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), fac, nmem));
            free(fac);
            expr_free(groups[i].base);
            groups[i].base = merged;
            if (target_flip) {           /* ei aliases groups[i].exponent */
                expr_free(groups[i].exponent);
                groups[i].exponent = expr_copy(neg_ei);
            }
            /* Drop the other members, highest index first. */
            for (size_t k = nmem; k-- > 1; ) {
                size_t j = mem[k];
                expr_free(groups[j].base);
                expr_free(groups[j].exponent);
                for (size_t s2 = j; s2 + 1 < group_count; s2++) groups[s2] = groups[s2 + 1];
                group_count--;
            }
        }
        free(mem);
        free(flip);
        expr_free(neg_ei);
    }

    /* Trig / hyperbolic ratio canonicalization. After base grouping, look for
     * factors of the form Sin/Cos/Tan/Cot/Sec/Csc[arg] (or hyperbolic
     * counterparts) sharing an arg and re-emit the shortest form:
     *   Sin[x]/Cos[x] -> Tan[x],   1/Cos[x] -> Sec[x],   1/Tanh[x] -> Coth[x],
     *   Sin[x]*Csc[x] -> 1,        Cos[x]*Tan[x] -> Sin[x], etc. */
    trig_canon_groups(groups, &group_count);

    /* (Sqrt-coefficient absorption, c Power[r, +-1/2] -> sign (p/q) Sqrt[a/b],
     * is subsumed by the O5 pass above.) */

    /* Fold the numeric coefficient into the interval product (num_prod becomes 1
     * once absorbed). Done before the Complex fold so it also reads a clean
     * num_prod. If the product is undecidable, num_prod and interval_val stay as
     * separate factors. */
    if (interval_val && !(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) {
        Expr* np_iv = interval_from_scalar(num_prod);
        Expr* prod = interval_multiply(np_iv, interval_val);
        expr_free(np_iv);
        if (prod) {
            expr_free(interval_val); interval_val = prod;
            expr_free(num_prod); num_prod = expr_new_integer(1);
        }
    }

    if (complex_val && !(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) {
        Expr *re, *im; is_complex(complex_val, &re, &im);
        Expr* nr = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(num_prod), expr_copy(re)}, 2));
        Expr* ni = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){expr_copy(num_prod), expr_copy(im)}, 2));
        expr_free(complex_val); complex_val = make_complex(nr, ni);
        expr_free(num_prod); num_prod = expr_new_integer(1);
    }

    size_t final_count = 0;
    if (!(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) final_count++;
    if (complex_val) final_count++;
    if (interval_val) final_count++;
    for (size_t i = 0; i < group_count; i++) {
        if (!(groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 0)) final_count++;
    }

    if (final_count == 0) {
        expr_free(num_prod); if (complex_val) expr_free(complex_val);
        if (interval_val) expr_free(interval_val);
        for(size_t j=0; j<group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
        if (heap_bufs) free(groups);
        return expr_new_integer(1);
    }

    /* final_count <= 2 + group_count <= 2 + n; when !heap_bufs, n <= TIMES_SMALL_N,
     * so the result vector fits a stack buffer (expr_new_function memcpys it).
     * Saves a malloc/free per small Times in tight numeric loops. */
    Expr* final_args_stack[TIMES_SMALL_N + 3];
    Expr** final_args = heap_bufs ? malloc(sizeof(Expr*) * final_count)
                                  : final_args_stack;
    size_t idx = 0;
    if (!(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) final_args[idx++] = num_prod;
    else expr_free(num_prod);
    if (complex_val) final_args[idx++] = complex_val;
    if (interval_val) final_args[idx++] = interval_val;
    for (size_t i = 0; i < group_count; i++) {
        if (groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 0) {
            expr_free(groups[i].base); expr_free(groups[i].exponent); continue;
        }
        if (groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 1) {
            final_args[idx++] = groups[i].base; expr_free(groups[i].exponent);
        } else {
            /* A group that is still exactly one input factor b^e -- the same
             * base and exponent nodes, untouched by every pass above -- is
             * that factor.  When the factor is an evaluated fixed point,
             * re-evaluating a rebuilt Power[b, e] could only re-derive it:
             * hand back the factor instead (a full Power dispatch saved per
             * power in every product, e.g. -Cos[t] Sin[t]^k). */
            Expr* same = NULL;
            for (size_t a = 0; a < n && !same; a++) {
                Expr* f = res->data.function.args[a];
                if (f->type == EXPR_FUNCTION && f->data.function.arg_count == 2 &&
                    f->data.function.args[0] == groups[i].base &&
                    f->data.function.args[1] == groups[i].exponent &&
                    f->data.function.head->type == EXPR_SYMBOL &&
                    f->data.function.head->data.symbol.name == SYM_Power &&
                    eval_node_reusable(f))
                    same = f;
            }
            if (same) {
                expr_free(groups[i].base); expr_free(groups[i].exponent);
                final_args[idx++] = expr_copy(same);
            } else {
                final_args[idx++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){groups[i].base, groups[i].exponent}, 2));
            }
        }
    }
    if (heap_bufs) free(groups);
    if (idx == 1) {
        Expr* res_final = final_args[0];
        if (heap_bufs) free(final_args);
        return res_final;
    }
    /* Canonically order the collapsed factors here (group count is small) so
     * the evaluator can skip the generic ORDERLESS sort of the raw input; see
     * the ORDERLESS block in eval.c. The leading numeric coefficient still
     * sorts first (numbers precede symbols/powers in expr_compare). */
    if (idx >= 2 && orderless) expr_orderless_sort(final_args, idx);   /* memoized symbol sets */
    /* A product of exactly a coefficient -1 or -1/n and a sum distributes the
     * sign, as in the reference: -(-1 + x) is 1 - x, Times[-1/2, x + y] is
     * (-x - y)/2. Only these two shapes: Times[-3/2, x + y], Times[-1., x + y],
     * Times[-1, x + y, z] and -(x + y)^2 stay as they are. */
    if (idx == 2 && final_args[1]->type == EXPR_FUNCTION &&
        final_args[1]->data.function.head->type == EXPR_SYMBOL &&
        final_args[1]->data.function.head->data.symbol.name == SYM_Plus) {
        Expr* c = final_args[0];
        bool minus_one = c->type == EXPR_INTEGER && c->data.integer == -1;
        bool minus_unit_frac = c->type == EXPR_FUNCTION && c->data.function.head->type == EXPR_SYMBOL &&
            c->data.function.head->data.symbol.name == SYM_Rational && c->data.function.arg_count == 2 &&
            c->data.function.args[0]->type == EXPR_INTEGER && c->data.function.args[0]->data.integer == -1;
        if (minus_one || minus_unit_frac) {
            Expr* sum = final_args[1];
            size_t m = sum->data.function.arg_count;
            Expr** terms = malloc(sizeof(Expr*) * (m ? m : 1));
            for (size_t i = 0; i < m; i++)
                terms[i] = expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(sum->data.function.args[i]) }, 2);
            Expr* neg = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), terms, m));
            free(terms);
            expr_free(sum);
            Expr* out;
            if (minus_one) {
                expr_free(c);
                out = neg;
            } else {
                Expr* inv = make_rational(1, c->data.function.args[1]->type == EXPR_INTEGER
                                             ? c->data.function.args[1]->data.integer : 1);
                bool big_den = c->data.function.args[1]->type != EXPR_INTEGER;
                if (big_den) {
                    /* -1/n with a bignum n: 1/n is Rational[1, n]. */
                    expr_free(inv);
                    inv = expr_new_function(expr_new_symbol(SYM_Rational),
                        (Expr*[]){ expr_new_integer(1), expr_copy(c->data.function.args[1]) }, 2);
                }
                expr_free(c);
                out = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ inv, neg }, 2));
            }
            if (heap_bufs) free(final_args);
            return out;
        }
    }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Times), final_args, idx);
    if (heap_bufs) free(final_args);
    times_certify_result(res, result, orderless);
    return result;
}
