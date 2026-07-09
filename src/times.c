#include "times.h"
#include "arithmetic.h"
#include "complex.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "trig_canon.h"
#include "series.h"
#include "ndarray.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

static bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol == SYM_Overflow;
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

static Expr* multiply_numbers(Expr* a, Expr* b) {
    if (is_overflow(a) || is_overflow(b)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
#ifdef USE_MPFR
    if (numeric_any_mpfr(a, b) && !is_complex(a, NULL, NULL) && !is_complex(b, NULL, NULL)) {
        Expr* r = numeric_mpfr_mul(a, b, 0);
        if (r) return r;
    }
#endif
    if (a->type == EXPR_REAL || b->type == EXPR_REAL) {
        double va = 0.0, vb = 0.0;
        int64_t n, d;
        
        if (a->type == EXPR_REAL) va = a->data.real;
        else if (a->type == EXPR_INTEGER) va = (double)a->data.integer;
        else if (a->type == EXPR_BIGINT) va = mpz_get_d(a->data.bigint);
        else if (is_rational(a, &n, &d)) va = (double)n / d;

        if (b->type == EXPR_REAL) vb = b->data.real;
        else if (b->type == EXPR_INTEGER) vb = (double)b->data.integer;
        else if (b->type == EXPR_BIGINT) vb = mpz_get_d(b->data.bigint);
        else if (is_rational(b, &n, &d)) vb = (double)n / d;
        
        return expr_new_real(va * vb);
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

        mpz_t av, bv, r;
        if (a_is_rat) mpz_init_set_si(av, n1);
        else expr_to_mpz(a, av);

        if (b_is_rat) mpz_init_set_si(bv, n2);
        else expr_to_mpz(b, bv);
        
        mpz_init(r);
        mpz_mul(r, av, bv);
        mpz_clear(av); mpz_clear(bv);
        
        // Now handle the denominators!
        int64_t den = d1 * d2;
        if (den == 1) {
            Expr* res = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clear(r);
            return res;
        }
        
        // Compute GCD of r and den
        mpz_t m_den, m_gcd;
        mpz_inits(m_den, m_gcd, NULL);
        mpz_set_si(m_den, den);
        mpz_gcd(m_gcd, r, m_den);
        
        mpz_divexact(r, r, m_gcd);
        mpz_divexact(m_den, m_den, m_gcd);
        
        int64_t final_den = mpz_get_si(m_den);
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
            mpz_init_set_si(mn1, n1); mpz_init_set_si(md1, d1);
            mpz_init_set_si(mn2, n2); mpz_init_set_si(md2, d2);
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
                   a->data.function.head->data.symbol == SYM_Rational &&
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
                       b->data.function.head->data.symbol == SYM_Rational &&
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

Expr* builtin_times(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_integer(1);
    if (n == 1) return expr_copy(res->data.function.args[0]);

    /* NDArray fast path: n same-shape NDArray operands multiply elementwise
     * over raw double buffers (mirrors Listable Times over Lists). If they are
     * all NDArrays but of disagreeing shape, warn (NDArray::shape) and leave
     * the product unevaluated, like Dot::dotsh. A NULL from a mixed
     * NDArray/scalar operand set falls through to the generic path, which
     * treats the NDArrays as opaque non-numeric factors. */
    {
        Expr* fast = ndarray_elementwise(res->data.function.args, n, false);
        if (fast) return fast;
        if (ndarray_warn_shape_mismatch(res->data.function.args, n, "multiplied"))
            return NULL;
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
        Expr** numed = numeric_contagion_args(res->data.function.args, n);
        if (numed) {
            for (size_t i = 0; i < n; i++) {
                expr_free(res->data.function.args[i]);
                res->data.function.args[i] = numed[i];
            }
            free(numed);
        }
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
                    if (!arith_warnings_muted())
                        fprintf(stderr,
                            "Infinity::indet: Indeterminate expression 0 %s encountered.\n", what);
                    return expr_new_symbol(SYM_Indeterminate);
                }
                if (cinf_count > 0) return expr_new_symbol(SYM_ComplexInfinity);
                /* inf_count > 0, coeff != 0 */
                if (sign > 0) return expr_new_symbol(SYM_Infinity);
                return expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
            }
            expr_free(coeff);
            /* fall through: keep symbolic, let the rest run */
        }
    }

    Expr* num_prod = expr_new_integer(1);
    Expr* complex_val = NULL;

    BasePower* groups = malloc(sizeof(BasePower) * n);
    size_t group_count = 0;

    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        if (is_overflow(arg)) {
            expr_free(num_prod); if (complex_val) expr_free(complex_val);
            for(size_t j=0; j<group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
            free(groups); return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        }

        if (expr_is_numeric_like(arg) && !is_complex(arg, NULL, NULL)) {
            Expr* next = multiply_numbers(num_prod, arg);
            if (next == NULL) {
                /* Could not fold this numeric factor into num_prod;
                 * stash it as its own group so we don't lose it (and so
                 * we don't propagate NULL to the type checks below). */
                groups[group_count].base = expr_copy(arg);
                groups[group_count].exponent = expr_new_integer(1);
                group_count++;
            } else {
                expr_free(num_prod);
                num_prod = next;
            }
        } else if (is_complex(arg, NULL, NULL) || (arg->type == EXPR_SYMBOL && arg->data.symbol == SYM_I)) {
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
        } else {
            Expr* base = arg; Expr* exponent;
            if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL && arg->data.function.head->data.symbol == SYM_Power && arg->data.function.arg_count == 2) {
                base = arg->data.function.args[0]; exponent = expr_copy(arg->data.function.args[1]);
            } else { exponent = expr_new_integer(1); }
            
            int found = -1;
            for (size_t j = 0; j < group_count; j++) { if (expr_eq(groups[j].base, base)) { found = (int)j; break; } }
            if (found != -1) {
                Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){groups[found].exponent, exponent}, 2));
                groups[found].exponent = new_exp;
            } else {
                groups[group_count].base = expr_copy(base);
                groups[group_count].exponent = (arg == base) ? exponent : exponent;
                group_count++;
            }
        }
    }

    if (num_prod->type == EXPR_INTEGER && num_prod->data.integer == 0) {
        if (complex_val) expr_free(complex_val);
        for(size_t j=0; j<group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
        free(groups); return num_prod;
    }

    /* Radical canonicalization: for each Power[b, q] group with b a positive
     * integer (>= 2) and q rational, fold factors of b appearing in the
     * accumulated rational coefficient num_prod into the exponent. This is
     * what turns
     *      Sqrt[2] * 1/2   ->  Power[2, -1/2]            (Sqrt[2]/2)
     *      2^(1/3) * 1/2   ->  Power[2, -2/3]            (2^(1/3)/2)
     *      8/Sqrt[2]       ->  4 Sqrt[2]                  (Power[2, 1/2] * 8)
     * Pull from den(num_prod) unconditionally — the canonical form has
     * den(num_prod) coprime to b. Pull from num(num_prod) only enough to
     * raise a still-negative exponent up toward 0; pulling more would push
     * the exponent >= 1 and trigger Power's integer-part extraction, which
     * would re-emit the factor and loop. The cap k_n_target is the minimum
     * number of pulls needed to make the exponent non-negative.
     *
     * We work directly on GMP num/den rather than constructing intermediate
     * Expr nodes per pull — divisibility tests and exact divisions on small
     * mpz values are essentially free, and there is no evaluator round-trip.
     * The whole pass is O(sum of factors-of-b in num_prod) per group with
     * a positive integer base.
     *
     * Runs unconditionally even with complex_val present, so that
     * Plus terms with mixed real/complex coefficients normalise their
     * Power exponents consistently and collect like-radical terms. */
    if (group_count > 0 &&
        (expr_is_integer_like(num_prod) || is_rational(num_prod, NULL, NULL))) {
        bool any_radical_group = false;
        for (size_t gi = 0; gi < group_count; gi++) {
            Expr* gb = groups[gi].base;
            if (gb->type == EXPR_INTEGER && gb->data.integer >= 2 &&
                is_rational(groups[gi].exponent, NULL, NULL)) {
                any_radical_group = true;
                break;
            }
        }
        if (any_radical_group) {
            mpz_t num_z, den_z;
            bool num_extracted = false;
            if (expr_is_integer_like(num_prod)) {
                expr_to_mpz(num_prod, num_z);
                mpz_init_set_ui(den_z, 1);
                num_extracted = true;
            } else if (num_prod->type == EXPR_FUNCTION &&
                       num_prod->data.function.head->type == EXPR_SYMBOL &&
                       num_prod->data.function.head->data.symbol == SYM_Rational &&
                       num_prod->data.function.arg_count == 2 &&
                       expr_is_integer_like(num_prod->data.function.args[0]) &&
                       expr_is_integer_like(num_prod->data.function.args[1])) {
                expr_to_mpz(num_prod->data.function.args[0], num_z);
                expr_to_mpz(num_prod->data.function.args[1], den_z);
                num_extracted = true;
            }

            if (num_extracted) {
                bool num_changed = false;
                for (size_t gi = 0; gi < group_count; gi++) {
                    Expr* gb = groups[gi].base;
                    Expr* ge = groups[gi].exponent;
                    if (gb->type != EXPR_INTEGER || gb->data.integer < 2) continue;
                    int64_t en, ed;
                    if (!is_rational(ge, &en, &ed)) continue;  /* also accepts Integer with ed=1 */

                    int64_t b = gb->data.integer;
                    mpz_t b_z;
                    mpz_init_set_si(b_z, b);

                    /* Pull factors of b from den into the exponent, but
                     * stop when the resulting exponent magnitude would
                     * exceed 1: pulling further would push |new_en/ed| > 1
                     * and Power's integer-part extraction would re-emit
                     * the factor (oscillation). The cap is
                     *    k_d_max = ceil(en/ed)  for en > 0
                     *            = 0            for en <= 0
                     * which keeps new_en > -ed (i.e., new exponent > -1).
                     * Examples:
                     *    Sqrt[2]/4 (en=1, ed=2):  k_d_max=1 -> 2^(-1/2)/2
                     *    1/4 * 2^(2/3) (en=2,ed=3): k_d_max=1 -> 2^(-1/3)/2
                     *    Power[2, -2/3] * 1 (en=-2, ed=3): k_d_max=0 (skip)
                     */
                    int64_t k_d_max = (en > 0) ? ((en + ed - 1) / ed) : 0;
                    int64_t k_d = 0;
                    while (k_d < k_d_max && mpz_divisible_p(den_z, b_z)) {
                        mpz_divexact(den_z, den_z, b_z);
                        k_d++;
                    }

                    /* en/ed - k_d. Compute new numerator over the same ed. */
                    int64_t new_en;
                    if (k_d > 0 && ed > 0 && k_d > (INT64_MAX / ed)) {
                        /* Pathological denominator power; back out the den
                         * pull rather than overflow the int64 exponent. */
                        for (int64_t i = 0; i < k_d; i++) mpz_mul(den_z, den_z, b_z);
                        mpz_clear(b_z);
                        continue;
                    }
                    new_en = en - k_d * ed;

                    /* If still negative, pull from num just enough to make
                     * the new exponent >= 0. ceil((-new_en)/ed) bounds the
                     * pull; we may pull fewer if num runs out of factors. */
                    int64_t k_n = 0;
                    if (new_en < 0) {
                        int64_t target = (-new_en + ed - 1) / ed;
                        for (int64_t i = 0; i < target; i++) {
                            if (!mpz_divisible_p(num_z, b_z)) break;
                            mpz_divexact(num_z, num_z, b_z);
                            k_n++;
                        }
                        new_en += k_n * ed;
                    }

                    /* If still negative AND we have a complex coefficient
                     * whose Re/Im are integer-like and *both* divisible by
                     * b, pull factors of b out of the complex too -- this
                     * keeps Complex[2k, 2m] * Sqrt[2]/2 from leaving a
                     * stale factor of 2 in the complex (gives us the
                     * canonical Complex[k, m] * Power[2, -1/2] rather than
                     * Complex[2k, 2m] * Power[2, -3/2]). */
                    int64_t k_c = 0;
                    if (new_en < 0 && complex_val) {
                        Expr *re = NULL, *im = NULL;
                        is_complex(complex_val, &re, &im);
                        if (re && im &&
                            expr_is_integer_like(re) && expr_is_integer_like(im)) {
                            mpz_t re_z, im_z;
                            expr_to_mpz(re, re_z);
                            expr_to_mpz(im, im_z);
                            int64_t target = (-new_en + ed - 1) / ed;
                            for (int64_t i = 0; i < target; i++) {
                                /* Each component must be either zero or
                                 * divisible by b -- 0/b = 0 leaves it
                                 * unchanged, which is fine for purely-
                                 * imaginary or purely-real complexes. */
                                bool re_ok = (mpz_sgn(re_z) == 0) ||
                                             mpz_divisible_p(re_z, b_z);
                                bool im_ok = (mpz_sgn(im_z) == 0) ||
                                             mpz_divisible_p(im_z, b_z);
                                if (!re_ok || !im_ok) break;
                                /* Stop if we'd reduce the whole complex to
                                 * 0+0i (Re=0 AND Im=0). */
                                if (mpz_sgn(re_z) == 0 && mpz_sgn(im_z) == 0) break;
                                if (mpz_sgn(re_z) != 0) mpz_divexact(re_z, re_z, b_z);
                                if (mpz_sgn(im_z) != 0) mpz_divexact(im_z, im_z, b_z);
                                k_c++;
                            }
                            if (k_c > 0) {
                                expr_free(complex_val);
                                Expr* re_e = expr_bigint_normalize(expr_new_bigint_from_mpz(re_z));
                                Expr* im_e = expr_bigint_normalize(expr_new_bigint_from_mpz(im_z));
                                complex_val = make_complex(re_e, im_e);
                                new_en += k_c * ed;
                            }
                            mpz_clear(re_z); mpz_clear(im_z);
                        }
                    }

                    if (k_d == 0 && k_n == 0 && k_c == 0) { mpz_clear(b_z); continue; }
                    num_changed = true;

                    /* Reduce new_en/ed by their gcd. */
                    int64_t abs_en = new_en >= 0 ? new_en : -new_en;
                    int64_t g_e = gcd(abs_en, ed);
                    if (g_e > 1) { new_en /= g_e; ed /= g_e; }

                    expr_free(groups[gi].exponent);
                    if (ed == 1) {
                        groups[gi].exponent = expr_new_integer(new_en);
                    } else {
                        groups[gi].exponent = make_rational(new_en, ed);
                    }
                    mpz_clear(b_z);
                }

                if (num_changed) {
                    /* Reduce num/den by gcd and rebuild num_prod. */
                    if (mpz_sgn(den_z) < 0) { mpz_neg(num_z, num_z); mpz_neg(den_z, den_z); }
                    mpz_t g_z;
                    mpz_init(g_z);
                    mpz_gcd(g_z, num_z, den_z);
                    if (mpz_cmp_ui(g_z, 1) > 0) {
                        mpz_divexact(num_z, num_z, g_z);
                        mpz_divexact(den_z, den_z, g_z);
                    }
                    mpz_clear(g_z);
                    expr_free(num_prod);
                    if (mpz_cmp_ui(den_z, 1) == 0) {
                        num_prod = expr_bigint_normalize(expr_new_bigint_from_mpz(num_z));
                    } else {
                        Expr* num_e = expr_bigint_normalize(expr_new_bigint_from_mpz(num_z));
                        Expr* den_e = expr_bigint_normalize(expr_new_bigint_from_mpz(den_z));
                        num_prod = expr_new_function(expr_new_symbol(SYM_Rational),
                                                     (Expr*[]){num_e, den_e}, 2);
                    }
                }
                mpz_clear(num_z);
                mpz_clear(den_z);
            }
        }
    }

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
        if (!is_positive_numeric_expr(groups[i].base)) continue;
        int64_t pi_n, pi_d;
        if (!is_rational(groups[i].exponent, &pi_n, &pi_d)) continue;
        if (pi_d == 1) continue;
        if (pi_n == 0) continue;

        for (size_t j = i + 1; j < group_count; j++) {
            if (!is_positive_numeric_expr(groups[j].base)) continue;
            int64_t pj_n, pj_d;
            if (!is_rational(groups[j].exponent, &pj_n, &pj_d)) continue;
            if (pj_d != pi_d) continue;
            if (pj_n == 0) continue;

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

    /* Trig / hyperbolic ratio canonicalization. After base grouping, look for
     * factors of the form Sin/Cos/Tan/Cot/Sec/Csc[arg] (or hyperbolic
     * counterparts) sharing an arg and re-emit the shortest form:
     *   Sin[x]/Cos[x] -> Tan[x],   1/Cos[x] -> Sec[x],   1/Tanh[x] -> Coth[x],
     *   Sin[x]*Csc[x] -> 1,        Cos[x]*Tan[x] -> Sin[x], etc. */
    trig_canon_groups(groups, &group_count);

    /* Sqrt-coefficient absorption: c * Power[r, +/-1/2] folds into
     * sign(c) * Sqrt[c^2 * r^eps] followed by perfect-square extraction.
     * The result is the canonical Mathematica form sign * (p_sq / q_sq) *
     * Sqrt[p_rest / q_rest]. Examples:
     *   14 / Sqrt[10]        -> 7 Sqrt[2/5]    (98/5 after merge,  49 extracted)
     *   78 / Sqrt[66]        -> 13 Sqrt[6/11]  (1014/11 after,    169 extracted)
     *   (3/5) / Sqrt[2/5]    -> 3 / Sqrt[10]   (9/10 after,         9 extracted)
     *   21/5 Sqrt[2/213]     -> 7/5 Sqrt[6/71] (294/1775 after, 7 and 5 extracted)
     *
     * The eps = -1 path covers the original "Sqrt in denominator" case
     * Mathematica rationalises. The eps = +1 path covers the symmetric
     * case `c * Sqrt[r]` for rational coefficients on Rational radicands
     * (e.g. 21/5 Sqrt[2/213] above) -- the tight `p_sq > 1 AND q_sq > 1`
     * gate keeps integer-coefficient siblings like 2 Sqrt[2/15] and
     * 30 Sqrt[2/15] (which dominate Gram-Schmidt residuals) unchanged.
     *
     * The unchanged-form pre-check keeps inputs already in canonical form
     * literally identical -- e.g. 2/Sqrt[3] stays as is because 4/3 has
     * no extractable square beyond the trivial 2/Sqrt[3]. */
    if (complex_val == NULL && group_count == 1 &&
        (expr_is_integer_like(num_prod) || is_rational(num_prod, NULL, NULL)) &&
        !(num_prod->type == EXPR_INTEGER &&
          (num_prod->data.integer == 0 ||
           num_prod->data.integer == 1 ||
           num_prod->data.integer == -1))) {
        int64_t en_sa, ed_sa;
        bool eps_ok = is_rational(groups[0].exponent, &en_sa, &ed_sa) &&
                       (en_sa == -1 || en_sa == 1) && ed_sa == 2;
        /* eps=+1 gating for c * Sqrt[Rational[a, b]]:
         *   - cd >= 2 (rational coefficient): integer coefficients can't
         *     produce q_sq > 1 from a squarefree rd, so the rewrite is
         *     a no-op for them.
         *   - |cn|, cd, rn, rd all small (<= 1024): Simplify's internal
         *     polynomial Euclidean iteration produces sums of products
         *     of Sqrt[Rational] coefficients whose numerators grow
         *     without bound; canonicalising those gives PolynomialGCD
         *     a moving target and the test_qr_rank_deficient_3x3
         *     Simplify call doesn't terminate.  User-visible inputs
         *     (21/5 Sqrt[2/213], -1/7 Sqrt[21/2]) sit well below the
         *     threshold. */
        int64_t cd_chk = 1, cn_chk = 0;
        bool num_prod_rational = is_rational(num_prod, &cn_chk, &cd_chk);
        int64_t cn_abs = cn_chk < 0 ? -cn_chk : cn_chk;
        const int64_t CANON_CAP = 1024;
        bool plus_eps_rational_base_and_coeff =
            (en_sa == 1) &&
            (groups[0].base->type == EXPR_FUNCTION) &&
            (groups[0].base->data.function.head->type == EXPR_SYMBOL) &&
            (groups[0].base->data.function.head->data.symbol == SYM_Rational) &&
            num_prod_rational && (cd_chk >= 2) &&
            (cn_abs <= CANON_CAP) && (cd_chk <= CANON_CAP);
        /* Also bound the radicand magnitudes. */
        if (plus_eps_rational_base_and_coeff) {
            Expr* base = groups[0].base;
            int64_t rn_chk, rd_chk;
            if (!(is_rational(base, &rn_chk, &rd_chk) &&
                  (rn_chk >= 0 ? rn_chk : -rn_chk) <= CANON_CAP &&
                  rd_chk <= CANON_CAP)) {
                plus_eps_rational_base_and_coeff = false;
            }
        }
        /* Disable the eps=+1 rewrite while inside cancel_recursive's
         * PolynomialGCD pass.  Together / PolynomialGCD iterates over
         * polynomials whose coefficients carry Sqrt[Rational] terms;
         * each multiplication produces fresh `c * Sqrt[a/b]` factors
         * that my rewrite canonicalises to `c' * Sqrt[a'/b']` --
         * different shape, so PolynomialGCD's coefficient explosion
         * cannot terminate.  Leaving the form alone inside that pass
         * keeps the polynomial Euclidean stable, while still applying
         * the canonical form at user-facing evaluation. */
        extern int cancel_recursive_inside_gcd;
        if (cancel_recursive_inside_gcd > 0 && en_sa == 1) {
            plus_eps_rational_base_and_coeff = false;
        }
        if (eps_ok && (en_sa == -1 || plus_eps_rational_base_and_coeff)) {
            int eps_sa = (int)en_sa;
            Expr* gb = groups[0].base;
            mpz_t rn, rd;
            bool base_ok = false;
            if (gb->type == EXPR_INTEGER && gb->data.integer >= 2) {
                mpz_init_set_si(rn, gb->data.integer);
                mpz_init_set_ui(rd, 1);
                base_ok = true;
            } else if (gb->type == EXPR_BIGINT &&
                       mpz_sgn(gb->data.bigint) > 0 &&
                       mpz_cmp_ui(gb->data.bigint, 1) > 0) {
                mpz_init_set(rn, gb->data.bigint);
                mpz_init_set_ui(rd, 1);
                base_ok = true;
            } else if (gb->type == EXPR_FUNCTION &&
                       gb->data.function.head->type == EXPR_SYMBOL &&
                       gb->data.function.head->data.symbol == SYM_Rational &&
                       gb->data.function.arg_count == 2 &&
                       expr_is_integer_like(gb->data.function.args[0]) &&
                       expr_is_integer_like(gb->data.function.args[1])) {
                expr_to_mpz(gb->data.function.args[0], rn);
                expr_to_mpz(gb->data.function.args[1], rd);
                if (mpz_sgn(rn) > 0 && mpz_sgn(rd) > 0 &&
                    !(mpz_cmp_ui(rn, 1) == 0 && mpz_cmp_ui(rd, 1) == 0)) {
                    base_ok = true;
                } else {
                    mpz_clear(rn); mpz_clear(rd);
                }
            }
            if (base_ok) {
                mpz_t cn, cd;
                int sign = 1;
                if (expr_is_integer_like(num_prod)) {
                    expr_to_mpz(num_prod, cn);
                    sign = mpz_sgn(cn);
                    mpz_abs(cn, cn);
                    mpz_init_set_ui(cd, 1);
                } else {
                    int64_t cn_i, cd_i;
                    is_rational(num_prod, &cn_i, &cd_i);
                    sign = (cn_i < 0) ? -1 : 1;
                    mpz_init_set_si(cn, cn_i < 0 ? -cn_i : cn_i);
                    mpz_init_set_si(cd, cd_i);
                }
                if (sign != 0) {
                    /* p/q = |c|^2 * r^eps:
                     *   eps = +1: (cn^2 * rn) / (cd^2 * rd)
                     *   eps = -1: (cn^2 * rd) / (cd^2 * rn) */
                    mpz_t p, q, g;
                    mpz_inits(p, q, g, NULL);
                    mpz_mul(p, cn, cn);
                    mpz_mul(p, p, eps_sa > 0 ? rn : rd);
                    mpz_mul(q, cd, cd);
                    mpz_mul(q, q, eps_sa > 0 ? rd : rn);
                    mpz_gcd(g, p, q);
                    if (mpz_cmp_ui(g, 1) > 0) {
                        mpz_divexact(p, p, g);
                        mpz_divexact(q, q, g);
                    }
                    /* Extract square parts in-place from p and q. After the
                     * call, p holds p_rest and p_sq holds the extracted root
                     * (similarly for q). Trial-division by primes is O(sqrt)
                     * which is acceptable here; large bases would be unusual
                     * in symbolic input. */
                    mpz_t p_sq, q_sq;
                    mpz_init_set_ui(p_sq, 1);
                    mpz_init_set_ui(q_sq, 1);
                    mpz_t pr_p, pr_psq;
                    mpz_inits(pr_p, pr_psq, NULL);
                    /* p */
                    if (mpz_cmp_ui(p, 1) > 0) {
                        mpz_set_ui(pr_p, 2);
                        mpz_set_ui(pr_psq, 4);
                        while (mpz_divisible_p(p, pr_psq)) {
                            mpz_divexact(p, p, pr_psq);
                            mpz_mul(p_sq, p_sq, pr_p);
                        }
                        mpz_set_ui(pr_p, 3);
                        mpz_set_ui(pr_psq, 9);
                        while (mpz_cmp(pr_psq, p) <= 0) {
                            if (mpz_divisible_p(p, pr_psq)) {
                                mpz_divexact(p, p, pr_psq);
                                mpz_mul(p_sq, p_sq, pr_p);
                            } else {
                                mpz_add_ui(pr_p, pr_p, 2);
                                mpz_mul(pr_psq, pr_p, pr_p);
                            }
                        }
                    }
                    /* q */
                    if (mpz_cmp_ui(q, 1) > 0) {
                        mpz_set_ui(pr_p, 2);
                        mpz_set_ui(pr_psq, 4);
                        while (mpz_divisible_p(q, pr_psq)) {
                            mpz_divexact(q, q, pr_psq);
                            mpz_mul(q_sq, q_sq, pr_p);
                        }
                        mpz_set_ui(pr_p, 3);
                        mpz_set_ui(pr_psq, 9);
                        while (mpz_cmp(pr_psq, q) <= 0) {
                            if (mpz_divisible_p(q, pr_psq)) {
                                mpz_divexact(q, q, pr_psq);
                                mpz_mul(q_sq, q_sq, pr_p);
                            } else {
                                mpz_add_ui(pr_p, pr_p, 2);
                                mpz_mul(pr_psq, pr_p, pr_p);
                            }
                        }
                    }
                    mpz_clears(pr_p, pr_psq, NULL);

                    /* Usefulness criteria differ by eps:
                     *   eps=-1: original criterion -- extract a perfect
                     *           square (p_sq or q_sq > 1) or collapse
                     *           to integer base (p_rest=1 or q_rest=1).
                     *           A pure shuffle (e.g. 2/Sqrt[30] -> Sqrt[2/15])
                     *           costs Plus's like-term collection in QR.
                     *   eps=+1: tight gate -- BOTH p_sq > 1 AND q_sq > 1.
                     *           This is the user-visible Mathematica
                     *           canonicalisation 21/5 Sqrt[2/213] ->
                     *           7/5 Sqrt[6/71].  Looser criteria here
                     *           rewrite QR intermediates like
                     *           30 Sqrt[2/15] -> Sqrt[60] in ways that
                     *           PolynomialGCD / Together can't reduce. */
                    bool integer_radical = (mpz_cmp_ui(p, 1) == 0) ||
                                           (mpz_cmp_ui(q, 1) == 0);
                    bool useful_extraction;
                    if (eps_sa < 0) {
                        useful_extraction = (mpz_cmp_ui(p_sq, 1) > 0) ||
                                            (mpz_cmp_ui(q_sq, 1) > 0) ||
                                            integer_radical;
                    } else {
                        useful_extraction = (mpz_cmp_ui(p_sq, 1) > 0) &&
                                            (mpz_cmp_ui(q_sq, 1) > 0);
                    }

                    /* Unchanged-form check: the rewrite is a no-op exactly
                     * when the extracted coefficient equals |c| and the
                     * residual radical matches the input radical:
                     *   eps=-1, integer base: output 1/Sqrt[q_rest]
                     *                          need p=1, rd=1, q=rn.
                     *   eps=+1, Rational base: output Sqrt[Rational[p,q]]
                     *                          need p=rn, q=rd, both >=2. */
                    bool unchanged =
                        (mpz_cmp(p_sq, cn) == 0) &&
                        (mpz_cmp(q_sq, cd) == 0) &&
                        ((eps_sa < 0 && mpz_cmp_ui(p, 1) == 0 &&
                                         mpz_cmp_ui(rd, 1) == 0 &&
                                         mpz_cmp(q, rn) == 0) ||
                         (eps_sa > 0 && mpz_cmp_ui(rd, 1) > 0 &&
                                         mpz_cmp(p, rn) == 0 &&
                                         mpz_cmp(q, rd) == 0 &&
                                         mpz_cmp_ui(p, 1) > 0 &&
                                         mpz_cmp_ui(q, 1) > 0));
                    if (!useful_extraction) unchanged = true;

                    if (!unchanged) {
                        /* New coefficient (sign * p_sq / q_sq). */
                        mpz_t coef_num;
                        mpz_init_set(coef_num, p_sq);
                        if (sign < 0) mpz_neg(coef_num, coef_num);
                        expr_free(num_prod);
                        if (mpz_cmp_ui(q_sq, 1) == 0) {
                            num_prod = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(coef_num));
                        } else {
                            Expr* num_e = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(coef_num));
                            Expr* den_e = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(q_sq));
                            num_prod = expr_new_function(
                                expr_new_symbol(SYM_Rational),
                                (Expr*[]){ num_e, den_e }, 2);
                        }
                        mpz_clear(coef_num);

                        /* Residual radical. Three shapes by (p_rest, q_rest):
                         *   (1, 1)    radical is 1 -> drop the group (exp=0)
                         *   (>=2, 1)  Sqrt[p_rest]      (exp = +1/2)
                         *   (1, >=2)  1/Sqrt[q_rest]    (exp = -1/2)
                         *   else      Sqrt[p_rest/q_rest] (exp = +1/2) */
                        bool p_is_one = (mpz_cmp_ui(p, 1) == 0);
                        bool q_is_one = (mpz_cmp_ui(q, 1) == 0);
                        expr_free(groups[0].base);
                        expr_free(groups[0].exponent);
                        if (p_is_one && q_is_one) {
                            groups[0].base = expr_new_integer(1);
                            groups[0].exponent = expr_new_integer(0);
                        } else if (p_is_one) {
                            groups[0].base = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(q));
                            groups[0].exponent = make_rational(-1, 2);
                        } else if (q_is_one) {
                            groups[0].base = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(p));
                            groups[0].exponent = make_rational(1, 2);
                        } else {
                            Expr* num_e = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(p));
                            Expr* den_e = expr_bigint_normalize(
                                expr_new_bigint_from_mpz(q));
                            groups[0].base = expr_new_function(
                                expr_new_symbol(SYM_Rational),
                                (Expr*[]){ num_e, den_e }, 2);
                            groups[0].exponent = make_rational(1, 2);
                        }
                    }

                    mpz_clears(p, q, g, p_sq, q_sq, NULL);
                }
                mpz_clear(cn); mpz_clear(cd);
                mpz_clear(rn); mpz_clear(rd);
            }
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
    for (size_t i = 0; i < group_count; i++) {
        if (!(groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 0)) final_count++;
    }

    if (final_count == 0) {
        expr_free(num_prod); if (complex_val) expr_free(complex_val);
        for(size_t j=0; j<group_count; j++) { expr_free(groups[j].base); expr_free(groups[j].exponent); }
        free(groups); return expr_new_integer(1);
    }

    Expr** final_args = malloc(sizeof(Expr*) * final_count); size_t idx = 0;
    if (!(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) final_args[idx++] = num_prod;
    else expr_free(num_prod);
    if (complex_val) final_args[idx++] = complex_val;
    for (size_t i = 0; i < group_count; i++) {
        if (groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 0) {
            expr_free(groups[i].base); expr_free(groups[i].exponent); continue;
        }
        if (groups[i].exponent->type == EXPR_INTEGER && groups[i].exponent->data.integer == 1) {
            final_args[idx++] = groups[i].base; expr_free(groups[i].exponent);
        } else {
            final_args[idx++] = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){groups[i].base, groups[i].exponent}, 2));
        }
    }
    free(groups);
    if (idx == 1) { Expr* res_final = final_args[0]; free(final_args); return res_final; }
    Expr* result = expr_new_function(expr_new_symbol(SYM_Times), final_args, idx);
    free(final_args); return result;
}
