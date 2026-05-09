#include "times.h"
#include "arithmetic.h"
#include "complex.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "trig_canon.h"
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
 * positive in PicoCAS). Used to guard the radical-fusion rewrite
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

/* Sum of two rational exponents equals zero, exactly? */
static bool exponents_sum_to_zero(Expr* p, Expr* q) {
    int64_t pn, pd, qn, qd;
    if (!is_rational(p, &pn, &pd) || !is_rational(q, &qn, &qd)) return false;
    /* pn/pd + qn/qd = 0  iff  pn*qd + qn*pd = 0 */
    __int128_t lhs = (__int128_t)pn * qd + (__int128_t)qn * pd;
    return lhs == 0;
}

/* Sign of a rational Expr -- returns -1, 0, or +1, or 0 if not rational. */
static int rational_sign(Expr* e) {
    int64_t n, d;
    if (!is_rational(e, &n, &d)) return 0;
    /* Denominators are conventionally positive in PicoCAS Rational[n, d]. */
    if (n > 0) return (d > 0) ? +1 : -1;
    if (n < 0) return (d > 0) ? -1 : +1;
    return 0;
}

static Expr* multiply_numbers(Expr* a, Expr* b) {
    if (is_overflow(a) || is_overflow(b)) return expr_new_function(expr_new_symbol("Overflow"), NULL, 0);
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
    if (a->type == EXPR_BIGINT || b->type == EXPR_BIGINT) {
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
        // Wait! In PicoCAS, Rational CAN hold BigInts?! Let's allow it:
        Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clear(r);
        Expr* r_den = expr_new_integer(final_den);
        Expr* r_args[2] = { r_num, r_den };
        return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
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
            return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
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
            return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
        }
        if (a_ok) { mpz_clear(an); mpz_clear(ad); }
        if (b_ok) { mpz_clear(bn); mpz_clear(bd); }
    }

    return NULL;
}

Expr* make_times(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol("Times"), args, 2);
}

/* BasePower (the (base, exponent) pair) is shared with trig_canon.h, which
 * mutates the array we build here. */

Expr* builtin_times(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_integer(1);
    if (n == 1) return expr_copy(res->data.function.args[0]);

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
        if (has_indet) return expr_new_symbol("Indeterminate");
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
                    return expr_new_symbol("Indeterminate");
                }
                if (cinf_count > 0) return expr_new_symbol("ComplexInfinity");
                /* inf_count > 0, coeff != 0 */
                if (sign > 0) return expr_new_symbol("Infinity");
                return expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_new_symbol("Infinity") }, 2);
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
            free(groups); return expr_new_function(expr_new_symbol("Overflow"), NULL, 0);
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
                Expr* re = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(re1), expr_copy(re2)}, 2),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(im1), expr_copy(im2)}, 2)}, 2)
                }, 2));
                Expr* im = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(re1), expr_copy(im2)}, 2),
                    expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(re2), expr_copy(im1)}, 2)
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
                Expr* new_exp = eval_and_free(expr_new_function(expr_new_symbol("Plus"), (Expr*[]){groups[found].exponent, exponent}, 2));
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

    /* Radical fusion: collapse Power[a, q] * Power[b, -q] into Power[a/b, q]
     * when a, b are both positive numeric (integer, bigint, rational, real)
     * -- a > 0 and b > 0 ensures we stay on the principal branch. Applied
     * after base-grouping so same-base factors have already merged; only
     * heterogeneous pairs reach here. Examples:
     *   Sqrt[6]/Sqrt[2]    -> Sqrt[3]         (ratio is a positive integer)
     *   2^(1/3)/3^(1/3)    -> (2/3)^(1/3)     (ratio is a positive rational)
     * The ratio is computed by delegating to the evaluator so integer,
     * bigint, rational, and real bases all compose correctly. */
    for (size_t i = 0; i < group_count; i++) {
        if (!is_positive_numeric_expr(groups[i].base)) continue;
        if (!is_rational(groups[i].exponent, NULL, NULL)) continue;
        for (size_t j = i + 1; j < group_count; j++) {
            if (!is_positive_numeric_expr(groups[j].base)) continue;
            if (!is_rational(groups[j].exponent, NULL, NULL)) continue;
            if (!exponents_sum_to_zero(groups[i].exponent, groups[j].exponent)) continue;

            /* Pick the positive-exponent factor as numerator so the fused
             * power has the positive exponent (otherwise we'd print the
             * result as 1/(b/a)^q rather than (a/b)^q). */
            size_t ni = i, nj = j;
            if (rational_sign(groups[i].exponent) < 0) { ni = j; nj = i; }

            Expr* ratio = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
                expr_copy(groups[ni].base),
                expr_new_function(expr_new_symbol("Power"), (Expr*[]){
                    expr_copy(groups[nj].base), expr_new_integer(-1)
                }, 2)
            }, 2));

            /* Guard: only keep the fusion if the ratio is a positive numeric
             * we can wrap in a single Power. If evaluation produced anything
             * else (it shouldn't with positive numeric bases, but be safe)
             * discard and leave the original pair intact. */
            if (!is_positive_numeric_expr(ratio)) { expr_free(ratio); continue; }

            Expr* new_exp = expr_copy(groups[ni].exponent);

            size_t keep = (ni == i) ? i : j;
            size_t drop = (ni == i) ? j : i;
            expr_free(groups[keep].base);
            expr_free(groups[keep].exponent);
            expr_free(groups[drop].base);
            expr_free(groups[drop].exponent);
            groups[keep].base = ratio;
            groups[keep].exponent = new_exp;
            for (size_t k = drop; k + 1 < group_count; k++) groups[k] = groups[k + 1];
            group_count--;
            /* Restart from the beginning -- the fused base may pair with an
             * earlier group. Each fusion strictly decreases group_count so
             * the overall loop still terminates in O(group_count) restarts. */
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

    if (complex_val && !(num_prod->type == EXPR_INTEGER && num_prod->data.integer == 1)) {
        Expr *re, *im; is_complex(complex_val, &re, &im);
        Expr* nr = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(num_prod), expr_copy(re)}, 2));
        Expr* ni = eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_copy(num_prod), expr_copy(im)}, 2));
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
            final_args[idx++] = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){groups[i].base, groups[i].exponent}, 2));
        }
    }
    free(groups);
    if (idx == 1) { Expr* res_final = final_args[0]; free(final_args); return res_final; }
    Expr* result = expr_new_function(expr_new_symbol("Times"), final_args, idx);
    free(final_args); return result;
}
