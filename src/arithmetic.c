
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include "checked_int.h"
#include "numeric.h"       /* numeric_contagion_args: one-kind Complex */
#include "ndarray.h"
#include "linkedlist.h"       /* ndarray_divide: one-rounding buffer quotient (K1) */
#include "mth_machine_pow.h" /* mth_cdiv_machine: the one machine complex quotient */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

int g_arith_warnings_muted = 0;

/* Count trailing zeros of a nonzero word. */
static inline int arith_ctz64(uint64_t v) {
#if defined(__GNUC__)
    return __builtin_ctzll(v);
#else
    int k = 0;
    while (!(v & 1)) { v >>= 1; k++; }
    return k;
#endif
}

int64_t gcd(int64_t a, int64_t b) {
    /* Magnitudes as unsigned words: no llabs overflow at INT64_MIN. */
    uint64_t x = a < 0 ? (uint64_t)0 - (uint64_t)a : (uint64_t)a;
    uint64_t y = b < 0 ? (uint64_t)0 - (uint64_t)b : (uint64_t)b;
    /* Operands that fit in 32 bits -- every Rational a Table of small exact
     * entries builds -- run the Euclidean loop on 32-bit words: a 64-bit
     * divide costs several times a 32-bit one on x86-64, and this gcd sits
     * under every make_rational. */
    if (x <= UINT32_MAX && y <= UINT32_MAX) {
        uint32_t u = (uint32_t)x, v = (uint32_t)y;
        while (v) { uint32_t t = u % v; u = v; v = t; }
        return (int64_t)u;
    }
    /* Wider operands (a factorial denominator, k^j / j! past 12!) take the
     * binary gcd: shifts and subtractions only. The Euclidean loop spent ~40
     * cycles per 64-bit division, ~1.5k cycles per make_rational. */
    if (x == 0) return (int64_t)y;
    if (y == 0) return (int64_t)x;
    int sh = arith_ctz64(x | y);
    x >>= arith_ctz64(x);
    do {
        y >>= arith_ctz64(y);
        if (x > y) { uint64_t t = x; x = y; y = t; }
        y -= x;
    } while (y);
    return (int64_t)(x << sh);
}

/* lcm(|a|, |b|), or 0 when either is 0. SATURATES at INT64_MAX when the lcm
 * does not fit in int64 (including an INT64_MIN operand, whose magnitude has
 * no int64 form): the old (a/g)*b silently wrapped to a wrong, possibly
 * negative, value. Saturation is stable under chaining -- lcm(INT64_MAX, q)
 * is again >= INT64_MAX -- so a fold like solverad's stays saturated, and its
 * degree guard rejects it. A caller that needs the exact value treats
 * INT64_MAX as "use the GMP LCM" (rat.c). */
int64_t lcm(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    int64_t aa, bb, r;
    if (ci_abs_i64(a, &aa) || ci_abs_i64(b, &bb)) return INT64_MAX;
    if (ci_mul_i64(aa / gcd(aa, bb), bb, &r)) return INT64_MAX;
    return r;
}

Expr* make_rational(int64_t n, int64_t d) {
    if (d == 0) return NULL; // Error
    if (n == 0) return expr_new_integer(0);

    /* INT64_MIN has no positive int64 counterpart: llabs() in gcd() and the
     * sign flip below overflow, and n / -1 raises a hardware divide trap
     * (SIGFPE) -- Mathar M7801 crashed the kernel this way.  Exact via GMP. */
    if (n == INT64_MIN || d == INT64_MIN) {
        mpz_t zn, zd;
        mpz_init(zn); mpz_init(zd);
        mpz_import(zn, 1, 1, sizeof(uint64_t), 0, 0, &(uint64_t){ n < 0 ? (uint64_t)0 - (uint64_t)n : (uint64_t)n });
        if (n < 0) mpz_neg(zn, zn);
        mpz_import(zd, 1, 1, sizeof(uint64_t), 0, 0, &(uint64_t){ d < 0 ? (uint64_t)0 - (uint64_t)d : (uint64_t)d });
        if (d < 0) mpz_neg(zd, zd);
        Expr* r = make_rational_mpz(zn, zd);
        mpz_clear(zn); mpz_clear(zd);
        return r;
    }

    int64_t common = gcd(n, d);
    n /= common;
    d /= common;

    if (d < 0) {
        n = -n;
        d = -d;
    }

    if (d == 1) return expr_new_integer(n);

    Expr* args[2];
    args[0] = expr_new_integer(n);
    args[1] = expr_new_integer(d);
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

/* make_rational(-n, d) without the int64 negation: -INT64_MIN wraps back to
 * INT64_MIN, so the printer's "flip a negative exponent" step re-produced the
 * same negative exponent forever (x^(-2^63/3) printed 1/(1/(1/(...)))). */
Expr* make_rational_neg(int64_t n, int64_t d) {
    if (n != INT64_MIN) return make_rational(-n, d);
    mpz_t zn, zd;
    mpz_init(zn); mpz_init(zd);
    mpz_set_ui(zn, 1);
    mpz_mul_2exp(zn, zn, 63);   /* -INT64_MIN = 2^63 */
    mpz_import(zd, 1, 1, sizeof(uint64_t), 0, 0,
               &(uint64_t){ d < 0 ? (uint64_t)0 - (uint64_t)d : (uint64_t)d });
    if (d < 0) mpz_neg(zd, zd);
    Expr* r = make_rational_mpz(zn, zd);
    mpz_clear(zn); mpz_clear(zd);
    return r;
}

/* Bignum analogue of make_rational: build a canonical Rational (or Integer /
 * BigInt when it reduces) from GMP numerator/denominator.  Reduces by gcd and
 * forces the denominator positive, matching the invariant every other rational
 * producer keeps.  Components are demoted to machine Integer via
 * expr_bigint_normalize when they fit.  Returns NULL only for a zero
 * denominator. */
Expr* make_rational_mpz(const mpz_t n, const mpz_t d) {
    if (mpz_sgn(d) == 0) return NULL;
    mpz_t nn, dd, g;
    mpz_init_set(nn, n);
    mpz_init_set(dd, d);
    if (mpz_sgn(dd) < 0) { mpz_neg(nn, nn); mpz_neg(dd, dd); }
    mpz_init(g);
    mpz_gcd(g, nn, dd);
    if (mpz_cmp_ui(g, 1) > 0) {
        mpz_divexact(nn, nn, g);
        mpz_divexact(dd, dd, g);
    }
    Expr* result;
    if (mpz_cmp_ui(dd, 1) == 0) {
        result = expr_bigint_normalize(expr_new_bigint_from_mpz(nn));
    } else {
        Expr* args[2];
        args[0] = expr_bigint_normalize(expr_new_bigint_from_mpz(nn));
        args[1] = expr_bigint_normalize(expr_new_bigint_from_mpz(dd));
        result = expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
    }
    mpz_clear(nn); mpz_clear(dd); mpz_clear(g);
    return result;
}

Expr* builtin_rational(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_expr = res->data.function.args[0];
    Expr* d_expr = res->data.function.args[1];
    
    if (n_expr->type == EXPR_INTEGER && d_expr->type == EXPR_INTEGER) {
        int64_t n = n_expr->data.integer;
        int64_t d = d_expr->data.integer;
        if (d == 0) {
            arith_warn("Power::infy: Infinite expression 1/0 encountered.\n");
            if (n == 0) {
                arith_warn("Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        
        /* Already canonical (d > 1, coprime): the common case -- every
         * evaluation of a Rational an arithmetic builtin produced lands here.
         * Decline without building and freeing the reduced copy. */
        if (d > 1 && n != 0 && n != INT64_MIN && gcd(n, d) == 1) return NULL;
        Expr* r = make_rational(n, d);
        if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL && r->data.function.head->data.symbol.name == SYM_Rational) {
            Expr* rn = r->data.function.args[0];
            Expr* rd = r->data.function.args[1];
            if (rn->type == EXPR_INTEGER && rd->type == EXPR_INTEGER && rn->data.integer == n && rd->data.integer == d) {
                // No simplification happened
                expr_free(r);
                return NULL;
            }
        }
        return r;
    }

    /* A BigInt part (e.g. Rational[2^70, 2^70], Rational[6, -2^70], or an
     * unnormalised BigInt that fits int64): canonicalise through GMP -- gcd
     * reduction, positive denominator, demotion. Decline only when the input
     * is already exactly the canonical tree (same values AND representation),
     * so an unnormalised part is still rewritten. */
    if (expr_is_integer_like(n_expr) && expr_is_integer_like(d_expr)) {
        mpz_t zn, zd;
        expr_to_mpz(n_expr, zn);
        expr_to_mpz(d_expr, zd);
        if (mpz_sgn(zd) == 0) {
            bool zero_num = mpz_sgn(zn) == 0;
            mpz_clear(zn); mpz_clear(zd);
            arith_warn("Power::infy: Infinite expression 1/0 encountered.\n");
            if (zero_num) {
                arith_warn("Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        Expr* r = make_rational_mpz(zn, zd);
        mpz_clear(zn); mpz_clear(zd);
        if (r && r->type == EXPR_FUNCTION && r->data.function.arg_count == 2
            && r->data.function.args[0]->type == n_expr->type
            && r->data.function.args[1]->type == d_expr->type
            && expr_eq(r->data.function.args[0], n_expr)
            && expr_eq(r->data.function.args[1], d_expr)) {
            expr_free(r);
            return NULL;
        }
        return r;
    }
    return NULL;
}

bool is_rational(const Expr* e, int64_t* n, int64_t* d) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        if (n) *n = e->data.integer;
        if (d) *d = 1;
        return true;
    }
    /* Function nodes built during evaluator transitions can transiently
     * carry a NULL head pointer; guard before dereferencing.  Without
     * this check, a Function with a NULL head reaching this fast-path
     * crashes the dispatcher (observed on Linux under nested-radical
     * Simplify; macOS heap layout makes the deref usually land in
     * mapped memory). */
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) {
        if (e->data.function.arg_count == 2 &&
            e->data.function.args[0] &&
            e->data.function.args[1] &&
            e->data.function.args[0]->type == EXPR_INTEGER &&
            e->data.function.args[1]->type == EXPR_INTEGER) {
            if (n) *n = e->data.function.args[0]->data.integer;
            if (d) *d = e->data.function.args[1]->data.integer;
            return true;
        }
    }
    return false;
}

/* Bigint-aware version: matches Integer, BigInt, or Rational[X, Y] where
 * both components are integer-like (Integer or BigInt).  Used in numeric
 * predicates so the Times/Plus folders recognise rationals whose
 * components have overflowed int64 — without this, a BigInt times a
 * Rational[1, BigInt] is left unsimplified and can defeat exact
 * polynomial division. */
bool is_rational_like(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        return true;
    }
    return false;
}

bool is_complex(Expr* e, Expr** re, Expr** im) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.arg_count == 2) {
        if (re) *re = e->data.function.args[0];
        if (im) *im = e->data.function.args[1];
        return true;
    }
    return false;
}

Expr* make_complex(Expr* re, Expr* im) {
    Expr* args[2] = { re, im };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

Expr* builtin_subtract(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    Expr* minus_one = expr_new_integer(-1);
    Expr* mb_args[2] = { minus_one, expr_copy(b) };
    Expr* minus_b = expr_new_function(expr_new_symbol(SYM_Times), mb_args, 2);

    Expr* p_args[2] = { expr_copy(a), minus_b };
    return expr_new_function(expr_new_symbol(SYM_Plus), p_args, 2);
}

/* The numeric kind of a Complex part, for builtin_complex's one-kind rule. */
enum { CPK_EXACT, CPK_MACHINE, CPK_MPFR, CPK_MPFR_ZERO, CPK_OTHER };
static int complex_part_kind(const Expr* e) {
    if (e->type == EXPR_REAL) return CPK_MACHINE;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) ? CPK_MPFR_ZERO : CPK_MPFR;
#endif
    if (is_rational_like(e)) return CPK_EXACT;
    return CPK_OTHER;
}

Expr* builtin_complex(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* re = res->data.function.args[0];
    Expr* im = res->data.function.args[1];

    if (im->type == EXPR_INTEGER && im->data.integer == 0) {
        return expr_copy(re);
    }

    /* One number, one kind: a machine Real or an exact number meeting a part
     * of another kind takes the contagion Plus / Times use
     * (numeric_contagion_args). Complex[1.0, 2] is Complex[1.0, 2.0] (it was
     * not SameQ to it), Complex[1.5`30, 2.] is machine (it printed 30 digits
     * but reported MachinePrecision), Complex[1/2, 1.5`30] is 30 digits. With
     * a machine 0. imaginary part an exact real part becomes its machine value
     * (Complex[2^70, 0.] was an exact bignum where Complex[2, 0.] was a Real).
     * Left alone:
     *  - an exact-zero part: canonical as it stands (Complex[0, 1.5], N[I, 30]);
     *  - an MPFR zero part: it carries ACCURACY, not precision, and lowering
     *    its width would change it (0``120 became 0.e-20);
     *  - two MPFR parts of different widths, and Complex[x_MPFR, 0.]: the
     *    special-function kernels build these deliberately (guard widths),
     *    and they were never mixed-kind. */
    int kre = complex_part_kind(re), kim = complex_part_kind(im);
    bool re_zero = re->type == EXPR_INTEGER && re->data.integer == 0;
    if (im->type == EXPR_REAL && im->data.real == 0.0) {
        if (kre == CPK_EXACT) {
            double v;
            if (numeric_machine_value(re, &v) && isfinite(v)) return expr_new_real(v);
        }
        return expr_copy(re);
    }
    if (!re_zero && kre != CPK_OTHER && kim != CPK_OTHER && kre != kim
        && kre != CPK_MPFR_ZERO && kim != CPK_MPFR_ZERO
        && (kre == CPK_MACHINE || kim == CPK_MACHINE || kre == CPK_EXACT || kim == CPK_EXACT)) {
        Expr* parts[2] = { re, im };
        Expr* out[2];
        if (numeric_contagion_args(parts, 2, out))
            return expr_new_function(expr_new_symbol(SYM_Complex), out, 2);
    }

    return NULL;
}

/* A Divide operand as a machine complex value: a Real, an Integer a double
 * holds exactly, or a Complex of machine Reals (*cplx). */
static bool divide_machine_complex_operand(const Expr* e, double* re, double* im, bool* cplx) {
    *im = 0.0; *cplx = false;
    if (e->type == EXPR_REAL) { *re = e->data.real; return isfinite(*re); }
    if (e->type == EXPR_INTEGER && e->data.integer >= -(INT64_C(1) << 53) &&
        e->data.integer <= (INT64_C(1) << 53)) {
        *re = (double)e->data.integer; return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Complex &&
        e->data.function.args[0]->type == EXPR_REAL &&
        e->data.function.args[1]->type == EXPR_REAL) {
        *re = e->data.function.args[0]->data.real;
        *im = e->data.function.args[1]->data.real;
        *cplx = true;
        return isfinite(*re) && isfinite(*im);
    }
    return false;
}

Expr* builtin_divide(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];

    /* A LinkedList operand threads as the List it stands for, one rounding per
     * element (linkedlist_ops.c: a raw quotient loop for machine reals, else
     * Divide per element), exactly as a packed operand's ndarray_divide. */
    if (num->type == EXPR_LINKEDLIST || den->type == EXPR_LINKEDLIST)
        return ll_thread_listable(res);

    /* A machine array is not a scalar, and the Real branch below cannot be
     * allowed to think it is.  That branch fires when EITHER operand is Real
     * and then reads the OTHER through a chain of scalar type tests, falling
     * back to 0.0 when none matches -- so with an array numerator and a Real
     * denominator it computed 0.0 / 2.0 and answered `0.` for the whole call.
     *
     * Harmless while the transparency gate materialised every packed argument
     * before this function ran; a live wrong answer the moment Divide was
     * marked packed-aware on 2026-08-02.  The rewrite at the bottom
     * (Times[num, Power[den, -1]]) is what handles an array, and it handles it
     * well -- both of those heads have buffer paths -- so the fix is to reach
     * it, which is also what makes marking the head aware worth anything.
     * Found by the packed-vs-plain differential in tests/test_packed_list.c. */
    if (num->type == EXPR_NDARRAY || den->type == EXPR_NDARRAY) {
        /* Real machine operands: one rounding per element, as the scalar
         * branch below (K1: the rewrite rounds 1/den, then the product). */
        bool degenerate = false;
        Expr* q = ndarray_divide(num, den, &degenerate);
        if (q) return q;
        if (degenerate) {
            Expr* r = ndarray_delist_and_reeval(res);
            if (r) return r;
        }
        goto rewrite;
    }

    /* The machine-Real fast path may only claim operands it can actually read
     * as a machine real (Real / Integer / BigInt).  A Rational, MPFR, or Complex
     * sibling of a machine real would otherwise be misread as 0 -- e.g.
     * Divide[2.5, 1/3] -> ComplexInfinity, Divide[N[1/4,30], 2.0] -> 0.,
     * Divide[2.0, 1+I] -> ComplexInfinity.  For those, fall through to the exact
     * Times[num, Power[den, -1]] rewrite below (Power/Times handle rational,
     * arbitrary-precision, and complex operands correctly). */
    bool num_mreal = (num->type == EXPR_REAL || num->type == EXPR_INTEGER || num->type == EXPR_BIGINT);
    bool den_mreal = (den->type == EXPR_REAL || den->type == EXPR_INTEGER || den->type == EXPR_BIGINT);
    if ((num->type == EXPR_REAL || den->type == EXPR_REAL) && num_mreal && den_mreal) {
        double vnum = (num->type == EXPR_REAL) ? num->data.real : (num->type == EXPR_INTEGER) ? (double)num->data.integer : numeric_mpz_to_double(num->data.bigint);
        double vden = (den->type == EXPR_REAL) ? den->data.real : (den->type == EXPR_INTEGER) ? (double)den->data.integer : numeric_mpz_to_double(den->data.bigint);
        if (vden == 0.0) {
            arith_warn("Power::infy: Infinite expression 1/0 encountered.\n");
            /* 0/0 is 0 ComplexInfinity, Indeterminate -- as the rewrite
             * Times[num, Power[den, -1]] answers it (audit NDArray#2: the
             * delisted lanes of NDArray[{0.}] / NDArray[{0.}] land here). */
            if (vnum == 0.0) {
                arith_warn("Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        return expr_new_real(vnum / vden);
    }

    /* A machine complex quotient: one division (mth_cdiv_machine, the
     * operation Compile's DIV_C and the buffer quotient use), not the
     * rewrite's rounded reciprocal times the numerator (K1). An exact operand
     * meets the machine one only where the contagion conversion is exact. */
    {
        double ar, ai, br, bi;
        bool ac, bc;
        if (divide_machine_complex_operand(num, &ar, &ai, &ac) &&
            divide_machine_complex_operand(den, &br, &bi, &bc) && (ac || bc)) {
            if (br == 0.0 && bi == 0.0) {
                arith_warn("Power::infy: Infinite expression 1/0 encountered.\n");
                if (ar == 0.0 && ai == 0.0) {
                    arith_warn("Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                    return expr_new_symbol(SYM_Indeterminate);
                }
                return expr_new_symbol(SYM_ComplexInfinity);
            }
            double qr, qi;
            mth_cdiv_machine(ar, ai, br, bi, &qr, &qi);
            Expr* parts[2] = { expr_new_real(qr), expr_new_real(qi) };
            return expr_new_function(expr_new_symbol(SYM_Complex), parts, 2);
        }
    }

#ifdef USE_MPFR
    /* A real arbitrary-precision quotient: ONE mpfr_div at the operands'
     * combined precision (K1, extended to MPFR: numeric-tower review, found by
     * PARANOIA). The rewrite below rounds 1/den and then the product -- twice,
     * so N[3,50]/N[9,50] missed N[1/3,50] by an ulp -- and its reciprocal
     * overflows when den is in MPFR's bottom binades: z/z was Overflow[] at
     * z = 2.^-2^30, where the quotient is 1. An MPFR operand with an exact
     * integer meets it at that precision; a machine operand goes by contagion
     * (machine), and a Rational / Complex keeps the rewrite. */
    if ((num->type == EXPR_MPFR || den->type == EXPR_MPFR) &&
        (num->type == EXPR_MPFR || num->type == EXPR_INTEGER || num->type == EXPR_BIGINT) &&
        (den->type == EXPR_MPFR || den->type == EXPR_INTEGER || den->type == EXPR_BIGINT) &&
        !(den->type == EXPR_MPFR ? mpfr_zero_p(den->data.mpfr) || !mpfr_number_p(den->data.mpfr)
                                 : expr_numeric_sign(den) == 0) &&
        !(num->type == EXPR_MPFR && !mpfr_number_p(num->data.mpfr))) {
        /* The MPFR operand's precision is the result's (contagion); an EXACT
         * integer operand is converted exactly -- at its own bit length --
         * so the division is the only rounding (converting 10^400 at 53 bits
         * first made Divide[1.`53, 10^400] an ulp off). */
        mpfr_prec_t prec = num->type == EXPR_MPFR ? mpfr_get_prec(num->data.mpfr)
                                                  : mpfr_get_prec(den->data.mpfr);
        if (num->type == EXPR_MPFR && den->type == EXPR_MPFR &&
            mpfr_get_prec(den->data.mpfr) < prec) prec = mpfr_get_prec(den->data.mpfr);
        mpfr_t a, b, out;
        bool own_a = num->type != EXPR_MPFR, own_b = den->type != EXPR_MPFR;
        if (own_a || own_b) {
            const Expr* ex = own_a ? num : den;
            mpz_t z;
            expr_to_mpz((Expr*)ex, z);
            size_t zb = mpz_sizeinbase(z, 2);
            mpfr_prec_t ep = (mpfr_prec_t)(zb < (size_t)MPFR_PREC_MIN ? MPFR_PREC_MIN : zb);
            if (own_a) { mpfr_init2(a, ep); mpfr_set_z(a, z, MPFR_RNDN); }
            else       { mpfr_init2(b, ep); mpfr_set_z(b, z, MPFR_RNDN); }
            mpz_clear(z);
        }
        mpfr_init2(out, prec);
        mpfr_div(out, own_a ? a : num->data.mpfr, own_b ? b : den->data.mpfr, MPFR_RNDN);
        if (own_a) mpfr_clear(a);
        if (own_b) mpfr_clear(b);
        return expr_new_mpfr_move(out);
    }
#endif

    int64_t n1, d1, n2, d2;
    if (is_rational(num, &n1, &d1) && is_rational(den, &n2, &d2)) {
        if (n2 == 0) {
            /* x / 0 with rational/integer x: 0/0 -> Indeterminate (handled in
             * Times when 0 multiplies ComplexInfinity); otherwise emit the
             * Power::infy message and yield ComplexInfinity. */
            arith_warn("Power::infy: Infinite expression 1/0 encountered.\n");
            if (n1 == 0) {
                arith_warn("Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        /* Cross products overflow int64 for wide operands (Divide[1/3^39, 3^39]
         * gave -1/7710839681495117703); on overflow take the exact
         * Times/Power rewrite below, which promotes to GMP. */
        int64_t rn, rd;
        if (!ci_mul_i64(n1, d2, &rn) && !ci_mul_i64(d1, n2, &rd)) {
            Expr* r = make_rational(rn, rd);
            if (r) return r;
        }
    }

rewrite:
    {
        Expr* minus_one = expr_new_integer(-1);
        Expr* p_args[2] = { expr_copy(den), minus_one };
        Expr* power = expr_new_function(expr_new_symbol(SYM_Power), p_args, 2);
        /* Evaluated here, so the Times below carries no Power[den, -1] for the
         * evaluator's one-rounding fusion (K1) to send back to Divide: a
         * buffer pair ndarray_divide declined (shape, dtype) would loop. */
        Expr* rec = evaluate(power);
        expr_free(power);

        Expr* t_args[2] = { expr_copy(num), rec };
        return expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
    }
}

bool is_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

bool is_complex_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_ComplexInfinity;
}

bool is_indeterminate_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Indeterminate;
}

int expr_numeric_sign(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    /* Rational[num, den].  Recurse into the components so bignum-component
     * rationals (e.g. Rational[1, 10^25], whose denominator overflows int64)
     * get a correct sign — the old is_rational() int64 path silently reported
     * 0 for them, which propagated into Times/Plus canonical ordering, the
     * superficially-negative test, and Sign itself.  The denominator is
     * conventionally positive, but reading both signs keeps a hand-built
     * Rational with a negative denominator correct. */
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        int sn = expr_numeric_sign(e->data.function.args[0]);
        int sd = expr_numeric_sign(e->data.function.args[1]);
        return sn * sd;
    }
    return 0;
}

bool expr_is_superficially_negative(Expr* e) {
    if (!e) return false;
    int s = expr_numeric_sign(e);
    if (s < 0) return true;
    if (s > 0) return false;
    /* s == 0: either a zero numeric (not negative) or a non-numeric --
     * fall through to Complex and Times shape checks. */
    Expr* re; Expr* im;
    if (is_complex(e, &re, &im)) {
        int rs = expr_numeric_sign(re);
        if (rs < 0) return true;
        if (rs > 0) return false;
        /* Real part is zero: use the imaginary part's sign. Catches the
         * "pure imaginary with negative coefficient" case, e.g. -2 I =
         * Complex[0, -2]. */
        return expr_numeric_sign(im) < 0;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times &&
        e->data.function.arg_count > 0) {
        /* Leading factor carries the syntactic sign, per Times canonical
         * ordering (numerics sort first). */
        return expr_is_superficially_negative(e->data.function.args[0]);
    }
    return false;
}

bool is_neg_infinity_form(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol.name != SYM_Times) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (!is_infinity_sym(b)) return false;
    return expr_numeric_sign(a) < 0;
}
