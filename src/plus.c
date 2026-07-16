
#include "plus.h"
#include "arithmetic.h"
#include "complex.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include "series.h"
#include "ndarray.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <gmp.h>

/* Classify a single Plus term for Infinity/Indeterminate handling.
 * Returns:
 *    0 = ordinary finite term
 *    1 = +Infinity (Infinity itself, or Times[c, Infinity, ...] with c > 0)
 *   -1 = -Infinity (Times[c, Infinity, ...] with c < 0)
 *    2 = ComplexInfinity (or any Times factor that contains it)
 *    3 = Indeterminate (or any Times factor that contains it)
 */
static int classify_plus_term(Expr* e) {
    if (is_indeterminate_sym(e)) return 3;
    if (is_complex_infinity_sym(e)) return 2;
    if (is_infinity_sym(e)) return 1;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times) {
        size_t ac = e->data.function.arg_count;
        bool has_inf = false, has_cinf = false, has_indet = false;
        for (size_t i = 0; i < ac; i++) {
            Expr* f = e->data.function.args[i];
            if (is_indeterminate_sym(f)) has_indet = true;
            else if (is_complex_infinity_sym(f)) has_cinf = true;
            else if (is_infinity_sym(f)) has_inf = true;
        }
        if (has_indet) return 3;
        if (has_cinf) return 2;
        if (has_inf) {
            /* Canonical Times has the numeric coefficient first. */
            Expr* f0 = e->data.function.args[0];
            if (expr_numeric_sign(f0) < 0) return -1;
            return 1;
        }
    }
    return 0;
}

static bool is_overflow(Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Overflow;
}

/* Helper: extract numeric coefficient and base expression from a term. 
   Always returns an OWNED coefficient. Base is NOT owned by default, unless *allocated_base is set to true. */
static void get_coeff_base(Expr* e, Expr** coeff, Expr** base, bool* allocated_base) {
    *allocated_base = false;
    if (is_overflow(e)) {
        *coeff = expr_copy(e);
        *base = NULL;
        return;
    }

    if (expr_is_numeric_like(e) || is_complex(e, NULL, NULL)) {
        *coeff = expr_copy(e);
        *base = NULL;
        return;
    }

    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol.name == SYM_Times) {
        if (e->data.function.arg_count >= 2) {
            Expr* first = e->data.function.args[0];
            if (expr_is_numeric_like(first) || is_complex(first, NULL, NULL)) {
                *coeff = expr_copy(first);
                if (e->data.function.arg_count == 2) {
                    *base = e->data.function.args[1];
                } else {
                    Expr** rest_args = malloc(sizeof(Expr*) * (e->data.function.arg_count - 1));
                    for (size_t i = 1; i < e->data.function.arg_count; i++) {
                        rest_args[i-1] = expr_copy(e->data.function.args[i]);
                    }
                    *base = expr_new_function(expr_new_symbol(SYM_Times), rest_args, e->data.function.arg_count - 1);
                    *allocated_base = true;
                    free(rest_args);
                }
                return;
            }
        }
    }
    
    *coeff = expr_new_integer(1);
    *base = e;
}

static Expr* add_numbers(Expr* a, Expr* b) {
    if (is_overflow(a) || is_overflow(b)) return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);

#ifdef USE_MPFR
    /* MPFR path: if either operand carries arbitrary precision, fold
     * through MPFR so precision is preserved (and, for Complex-of-MPFR,
     * the Complex case below handles real/imag separately). */
    if (numeric_any_mpfr(a, b) && !is_complex(a, NULL, NULL) && !is_complex(b, NULL, NULL)) {
        Expr* r = numeric_mpfr_add(a, b, 0);
        if (r) return r;
        /* Fall through if the operand structure wasn't purely real. */
    }
#endif

    Expr *re1 = NULL, *im1 = NULL, *re2 = NULL, *im2 = NULL;
    bool a_comp = is_complex(a, &re1, &im1);
    bool b_comp = is_complex(b, &re2, &im2);
    
    if (a_comp || b_comp) {
        if (!a_comp) re1 = a;
        if (!b_comp) re2 = b;
        Expr* p1 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){expr_copy(re1), expr_copy(re2)}, 2));
        
        Expr* zero = expr_new_integer(0);
        Expr* i1 = a_comp ? expr_copy(im1) : expr_copy(zero);
        Expr* i2 = b_comp ? expr_copy(im2) : expr_copy(zero);
        Expr* p2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){i1, i2}, 2));
        expr_free(zero);
        
        if (is_overflow(p1) || is_overflow(p2)) {
            expr_free(p1); expr_free(p2);
            return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        }

        Expr* res = make_complex(p1, p2);
        return res;
    }

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
        
        return expr_new_real(va + vb);
    }

    /* Mixed BigInt fast path: both sides must be plain integer-like
     * (Integer or BigInt) — Rational[BigInt, BigInt] would fail
     * is_rational() above and then crash expr_to_mpz here, so it has
     * to fall through to the generic GMP block at the bottom. */
    if ((a->type == EXPR_BIGINT || b->type == EXPR_BIGINT) &&
        expr_is_integer_like(a) && expr_is_integer_like(b)) {
        mpz_t av, bv, r;
        expr_to_mpz(a, av);
        expr_to_mpz(b, bv);
        mpz_init(r);
        mpz_add(r, av, bv);
        mpz_clear(av); mpz_clear(bv);
        Expr* res = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clear(r);
        return res;
    }

    if (a->type == EXPR_INTEGER && b->type == EXPR_INTEGER) {
        __int128_t res = (__int128_t)a->data.integer + b->data.integer;
        if (res > INT64_MAX || res < INT64_MIN) {
            mpz_t av, bv, r;
            expr_to_mpz(a, av);
            expr_to_mpz(b, bv);
            mpz_init(r);
            mpz_add(r, av, bv);
            mpz_clear(av); mpz_clear(bv);
            Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clear(r);
            return result;
        }
        return expr_new_integer((int64_t)res);
    }

    int64_t n1, d1, n2, d2;
    if (is_rational(a, &n1, &d1) && is_rational(b, &n2, &d2)) {
        __int128_t num = (__int128_t)n1 * d2 + (__int128_t)n2 * d1;
        __int128_t den = (__int128_t)d1 * d2;
        if (num > INT64_MAX || num < INT64_MIN || den > INT64_MAX || den < INT64_MIN) {
            /* Promote to a BigInt-backed Rational rather than flagging
             * overflow: this keeps symbolic series/polynomial computations
             * exact when intermediate numerators or denominators exceed
             * 64 bits. */
            mpz_t mn1, md1, mn2, md2, t1, t2, mnum, mden, g;
            mpz_init_set_si(mn1, n1); mpz_init_set_si(md1, d1);
            mpz_init_set_si(mn2, n2); mpz_init_set_si(md2, d2);
            mpz_inits(t1, t2, mnum, mden, g, NULL);
            mpz_mul(t1, mn1, md2);
            mpz_mul(t2, mn2, md1);
            mpz_add(mnum, t1, t2);
            mpz_mul(mden, md1, md2);
            if (mpz_sgn(mden) < 0) { mpz_neg(mden, mden); mpz_neg(mnum, mnum); }
            mpz_gcd(g, mnum, mden);
            if (mpz_sgn(g) != 0) { mpz_divexact(mnum, mnum, g); mpz_divexact(mden, mden, g); }
            Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(mnum));
            Expr* r_den = expr_bigint_normalize(expr_new_bigint_from_mpz(mden));
            mpz_clears(mn1, md1, mn2, md2, t1, t2, mnum, mden, g, NULL);
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
     * Without this, a Rational with BigInt components (produced by the
     * overflow branch above, or by intermediate polynomial arithmetic)
     * cannot be added to an Integer or another such Rational, and we
     * would return NULL — which the callers in builtin_plus dereference
     * via is_overflow(). */
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
            mpz_t t1, t2, mnum, mden, g;
            mpz_inits(t1, t2, mnum, mden, g, NULL);
            mpz_mul(t1, an, bd);
            mpz_mul(t2, bn, ad);
            mpz_add(mnum, t1, t2);
            mpz_mul(mden, ad, bd);
            if (mpz_sgn(mden) < 0) { mpz_neg(mden, mden); mpz_neg(mnum, mnum); }
            mpz_gcd(g, mnum, mden);
            if (mpz_sgn(g) != 0) {
                mpz_divexact(mnum, mnum, g);
                mpz_divexact(mden, mden, g);
            }
            Expr* r_num = expr_bigint_normalize(expr_new_bigint_from_mpz(mnum));
            Expr* r_den = expr_bigint_normalize(expr_new_bigint_from_mpz(mden));
            mpz_clears(an, ad, bn, bd, t1, t2, mnum, mden, g, NULL);
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

Expr* make_plus(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
}

/* Detect Times[-1, Plus[...]] (canonical form for `-(a+b+...)`).
 * Any other Times shape (more than 2 args, non-Plus second factor,
 * or a coefficient other than -1) returns false. Limiting to -1
 * matches Mathematica's Plus normalisation: `a - (b+c)` distributes
 * to `a - b - c`, while `a + 2(b+c)` stays unexpanded. */
static bool is_neg_of_plus(Expr* e, Expr** inner_plus_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (!e->data.function.head ||
        e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Times)
        return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* c = e->data.function.args[0];
    Expr* p = e->data.function.args[1];
    if (!(c->type == EXPR_INTEGER && c->data.integer == -1)) return false;
    if (p->type != EXPR_FUNCTION ||
        !p->data.function.head ||
        p->data.function.head->type != EXPR_SYMBOL ||
        p->data.function.head->data.symbol.name != SYM_Plus)
        return false;
    *inner_plus_out = p;
    return true;
}

Expr* builtin_plus(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t n = res->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);
    if (n == 1) return expr_copy(res->data.function.args[0]);

    /* NDArray fast path: same-shape NDArray operands add elementwise over raw
     * buffers, with numpy-style broadcasting of numeric scalars (1 + NDArray).
     * If the array operands disagree in shape, warn (NDArray::shape) and leave
     * the sum unevaluated, like Dot::dotsh. A NULL (a symbolic operand, or no
     * NDArray at all) falls through to the generic grouping, which treats any
     * NDArrays as opaque non-numeric terms. */
    {
        Expr* fast = ndarray_elementwise(res->data.function.args, n, true);
        if (fast) return fast;
        if (ndarray_warn_shape_mismatch(res->data.function.args, n, "added"))
            return NULL;
        /* NDArray combined with a symbolic term (NDArray + a): purely numeric,
         * so it can't be added elementwise. Warn, then fall through to leave
         * the sum unevaluated. */
        ndarray_warn_symbolic(res->data.function.args, n, "added");
    }

    /* Distribute Times[-1, Plus[...]] over the outer Plus. This is
     * the cancellation-enabling step that lets `Plus[A, -A]` evaluate
     * to 0 even when A is a Plus. Without it,
     *
     *     Plus[x, x^2, Times[-1, Plus[x, x^2]]]
     *
     * has bases {x, x^2, Plus[x,x^2]} -- three distinct -- and never
     * collects, leaving the input unchanged. Distributing the -1
     * gives `Plus[x, x^2, Times[-1, x], Times[-1, x^2]]`, whose
     * grouping by (coeff, base) collapses to 0.
     *
     * Mathematica-equivalent: matches the user-visible behaviour of
     * `a + b - (a+b) -> 0` without distributing arbitrary `c·Plus[..]`
     * (those stay unexpanded, just like in Mathematica). */
    bool any_neg_plus = false;
    size_t expanded_count = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* inner = NULL;
        if (is_neg_of_plus(res->data.function.args[i], &inner)) {
            any_neg_plus = true;
            expanded_count += inner->data.function.arg_count;
        } else {
            expanded_count++;
        }
    }
    if (any_neg_plus) {
        Expr** new_args = malloc(sizeof(Expr*) * expanded_count);
        size_t k = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* inner = NULL;
            if (is_neg_of_plus(res->data.function.args[i], &inner)) {
                /* Replace the Times[-1, Plus[t_1, ..., t_m]] arg with
                 * m fresh Times[-1, t_i] args. */
                for (size_t j = 0; j < inner->data.function.arg_count; j++) {
                    Expr* tj = inner->data.function.args[j];
                    new_args[k++] = eval_and_free(expr_new_function(
                        expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_new_integer(-1), expr_copy(tj) }, 2));
                }
                expr_free(res->data.function.args[i]);
                res->data.function.args[i] = NULL;
            } else {
                new_args[k++] = res->data.function.args[i];
                res->data.function.args[i] = NULL;
            }
        }
        free(res->data.function.args);
        res->data.function.args = new_args;
        res->data.function.arg_count = expanded_count;
        n = expanded_count;
        if (n == 0) return expr_new_integer(0);
        if (n == 1) {
            Expr* sole = expr_copy(res->data.function.args[0]);
            return sole;
        }
    }

    /* SeriesData arithmetic: if any summand is a power series, fold the whole
     * Plus through the series algebra. Runs before inexact contagion, which
     * would otherwise apply N to the whole SeriesData (turning its integer
     * nmin/nmax/den into reals and defeating recognition); mixed-precision
     * coefficients still combine correctly because so_add evaluates each sum.
     * A NULL result means the operands are incompatible (e.g. different
     * expansion points); fall through and keep them as a symbolic Plus. */
    for (size_t i = 0; i < n; i++) {
        if (is_series_data(res->data.function.args[i])) {
            Expr* r = series_combine_plus(res->data.function.args, n);
            if (r) return r;
            break;
        }
    }

    /* Inexact contagion: if any summand is an inexact Real/MPFR,
     * numericalize exact numeric parts in-place so `1. + Pi` collapses to
     * `4.14159` instead of staying as a frozen `1. + Pi` Plus. */
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
     * Mathematica semantics for Plus:
     *   Indeterminate + anything                      -> Indeterminate
     *   Infinity + (-Infinity)                        -> Indeterminate (with msg)
     *   ComplexInfinity + ComplexInfinity             -> Indeterminate (with msg)
     *   ComplexInfinity + (any other infinity)        -> Indeterminate (with msg)
     *   ComplexInfinity + finite                      -> ComplexInfinity
     *   Infinity + finite (and no -Infinity)          -> Infinity
     *   -Infinity + finite (and no +Infinity)         -> -Infinity
     */
    {
        bool has_indet = false;
        int pos_inf = 0, neg_inf = 0, cinf = 0;
        for (size_t i = 0; i < n; i++) {
            int c = classify_plus_term(res->data.function.args[i]);
            if (c == 3) has_indet = true;
            else if (c == 1) pos_inf++;
            else if (c == -1) neg_inf++;
            else if (c == 2) cinf++;
        }
        if (has_indet) return expr_new_symbol(SYM_Indeterminate);
        if (pos_inf > 0 && neg_inf > 0) {
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression -Infinity + Infinity encountered.\n");
            return expr_new_symbol(SYM_Indeterminate);
        }
        if (cinf > 1 || (cinf > 0 && (pos_inf > 0 || neg_inf > 0))) {
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression involving ComplexInfinity encountered.\n");
            return expr_new_symbol(SYM_Indeterminate);
        }
        if (cinf == 1) return expr_new_symbol(SYM_ComplexInfinity);
        if (pos_inf > 0) return expr_new_symbol(SYM_Infinity);
        if (neg_inf > 0) {
            return expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), expr_new_symbol(SYM_Infinity) }, 2);
        }
    }

    Expr* num_sum = expr_new_integer(0);
    
    typedef struct {
        Expr* base;
        Expr* coeff;
        bool temp_base;
    } TermGroup;
    
    TermGroup* groups = malloc(sizeof(TermGroup) * n);
    size_t group_count = 0;

    /* Open-addressing hash table mapping base-hash -> group index, so that
     * collecting like terms is O(1) expected per term rather than a linear
     * scan of every existing group. Without it, a Plus of many distinct
     * monomials (e.g. Total of 10^5 terms over ~10^3 powers) is O(n*groups),
     * i.e. quadratic. Capacity is a power of two >= 2n for a low load factor;
     * slot_group holds the group index (-1 = empty), slot_hash the base hash
     * so most probes reject without touching expr_eq. */
    size_t ht_cap = 8;
    while (ht_cap < n * 2) ht_cap <<= 1;
    size_t ht_mask = ht_cap - 1;
    int64_t* slot_group = malloc(sizeof(int64_t) * ht_cap);
    uint64_t* slot_hash = malloc(sizeof(uint64_t) * ht_cap);
    for (size_t s = 0; s < ht_cap; s++) slot_group[s] = -1;

    for (size_t i = 0; i < n; i++) {
        Expr* arg = res->data.function.args[i];
        Expr* c = NULL;
        Expr* b = NULL;
        bool is_temp_base = false;
        get_coeff_base(arg, &c, &b, &is_temp_base);
        
        if (is_overflow(c)) {
            expr_free(num_sum);
            for (size_t j = 0; j < group_count; j++) {
                expr_free(groups[j].coeff);
                if (groups[j].temp_base) expr_free(groups[j].base);
            }
            free(groups);
            free(slot_group);
            free(slot_hash);
            expr_free(c);
            if (is_temp_base && b != arg) expr_free(b);
            return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
        }

        if (b == NULL) {
            Expr* next_sum = add_numbers(num_sum, c);
            if (next_sum == NULL) {
                /* Unrecognised coefficient form — leave the sum
                 * symbolic by treating `c` as its own group with a
                 * unit base, and keep `num_sum` unchanged. */
                groups[group_count].base = expr_copy(c);
                groups[group_count].coeff = expr_new_integer(1);
                groups[group_count].temp_base = true;
                group_count++;
                expr_free(c);
                if (is_temp_base && b != arg) expr_free(b);
                continue;
            }
            expr_free(num_sum);
            num_sum = next_sum;
            expr_free(c);
            if (is_temp_base && b != arg) expr_free(b);
            if (is_overflow(num_sum)) {
                for (size_t j = 0; j < group_count; j++) {
                    expr_free(groups[j].coeff);
                    if (groups[j].temp_base) expr_free(groups[j].base);
                }
                free(groups);
                free(slot_group);
                free(slot_hash);
                return num_sum;
            }
        } else {
            uint64_t hb = expr_hash(b);
            size_t slot = (size_t)hb & ht_mask;
            int found = -1;
            while (slot_group[slot] != -1) {
                if (slot_hash[slot] == hb &&
                    expr_eq(groups[slot_group[slot]].base, b)) {
                    found = (int)slot_group[slot];
                    break;
                }
                slot = (slot + 1) & ht_mask;
            }

            if (found >= 0) {
                Expr* next_c = add_numbers(groups[found].coeff, c);
                if (next_c == NULL) {
                    /* Unrecognised coefficient combination — keep the
                     * existing group and add a fresh group for c*b,
                     * preserving the input rather than crashing. */
                    Expr* base_copy = is_temp_base ? b : expr_copy(b);
                    groups[group_count].base = base_copy;
                    groups[group_count].coeff = c;  /* take ownership */
                    groups[group_count].temp_base = is_temp_base;
                    group_count++;
                    continue;
                }
                expr_free(groups[found].coeff);
                groups[found].coeff = next_c;
                expr_free(c);
                if (is_temp_base) expr_free(b);
                if (is_overflow(groups[found].coeff)) {
                    expr_free(num_sum);
                    for (size_t j = 0; j < group_count; j++) {
                        expr_free(groups[j].coeff);
                        if (groups[j].temp_base) expr_free(groups[j].base);
                    }
                    free(groups);
                    free(slot_group);
                    free(slot_hash);
                    return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
                }
            } else {
                /* `slot` is the empty slot the probe stopped on. */
                slot_group[slot] = (int64_t)group_count;
                slot_hash[slot] = hb;
                groups[group_count].base = b;
                groups[group_count].coeff = c; // Already owned
                groups[group_count].temp_base = is_temp_base;
                group_count++;
            }
        }
    }
    free(slot_group);
    free(slot_hash);
    
    size_t final_count = group_count;
    bool has_num = !(num_sum->type == EXPR_INTEGER && num_sum->data.integer == 0);
    if (has_num) final_count++;
    
    if (final_count == 0) {
        expr_free(num_sum);
        free(groups);
        return expr_new_integer(0);
    }
    
    Expr** final_args = malloc(sizeof(Expr*) * final_count);
    size_t idx = 0;
    if (has_num) {
        Expr* re = NULL, *im = NULL;
        if (is_complex(num_sum, &re, &im) && im->type == EXPR_INTEGER && im->data.integer == 0) {
            final_args[idx++] = expr_copy(re);
        } else {
            final_args[idx++] = expr_copy(num_sum);
        }
    }
    
    for (size_t j = 0; j < group_count; j++) {
        if (groups[j].coeff->type == EXPR_INTEGER && groups[j].coeff->data.integer == 0) {
            continue;
        }
        if (groups[j].coeff->type == EXPR_INTEGER && groups[j].coeff->data.integer == 1) {
            final_args[idx++] = expr_copy(groups[j].base);
        } else {
            Expr* t_args[] = {expr_copy(groups[j].coeff), expr_copy(groups[j].base)};
            final_args[idx++] = expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
        }
    }

    Expr* final_res = NULL;
    if (idx == 0) {
        final_res = expr_new_integer(0);
        free(final_args);
    } else if (idx == 1) {
        final_res = final_args[0];
        free(final_args);
    } else {
        final_res = expr_new_function(expr_new_symbol(SYM_Plus), final_args, idx);
        free(final_args);
    }
    
    expr_free(num_sum);
    for(size_t j=0; j<group_count; j++) {
        expr_free(groups[j].coeff);
        if (groups[j].temp_base) expr_free(groups[j].base);
    }
    free(groups);
    
    return final_res;
}
