
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

int g_arith_warnings_muted = 0;

int64_t gcd(int64_t a, int64_t b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        a %= b;
        int64_t tmp = a;
        a = b;
        b = tmp;
    }
    return a;
}

int64_t lcm(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    a = llabs(a);
    b = llabs(b);
    return (a / gcd(a, b)) * b;
}

/* Single-argument GCD/LCM: mathematically the GCD or LCM of one number
 * is its absolute value (GCD[-5] == 5, GCD[-1/2] == 1/2).  Symbolic
 * input is returned unchanged so user-supplied rules / OneIdentity
 * pattern matching still apply (GCD[x] == x). */
static Expr* single_arg_abs_or_copy(Expr* arg) {
    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        return make_rational(llabs(n), llabs(d));
    }
    if (arg->type == EXPR_BIGINT) {
        mpz_t a;
        mpz_init(a);
        mpz_abs(a, arg->data.bigint);
        Expr* r = expr_bigint_normalize(expr_new_bigint_from_mpz(a));
        mpz_clear(a);
        return r;
    }
    return expr_copy(arg);
}

Expr* builtin_gcd(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return expr_new_integer(0);
    if (count == 1) return single_arg_abs_or_copy(res->data.function.args[0]);

    // Single pass: detect any bigint while confirming all args are integer-like
    bool any_bigint = false, all_integer_like = true;
    for (size_t i = 0; i < count; i++) {
        Expr* arg = res->data.function.args[i];
        if (!expr_is_integer_like(arg)) { all_integer_like = false; break; }
        if (arg->type == EXPR_BIGINT) any_bigint = true;
    }

    if (all_integer_like && any_bigint) {
        mpz_t running, tmp;
        mpz_init_set_ui(running, 0);
        for (size_t i = 0; i < count; i++) {
            expr_to_mpz(res->data.function.args[i], tmp);
            mpz_abs(tmp, tmp);
            mpz_gcd(running, running, tmp);
            mpz_clear(tmp);
        }
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(running));
        mpz_clear(running);
        return result;
    }

    int64_t running_n = 0;
    int64_t running_d = 1;

    for (size_t i = 0; i < count; i++) {
        int64_t n, d;
        if (!is_rational(res->data.function.args[i], &n, &d)) {
            return NULL;
        }
        if (i == 0) {
            running_n = llabs(n);
            running_d = llabs(d);
        } else {
            running_n = gcd(running_n, n);
            running_d = lcm(running_d, d);
        }
    }

    return make_rational(running_n, running_d);
}

Expr* builtin_lcm(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t count = res->data.function.arg_count;
    if (count == 0) return expr_new_integer(1);
    if (count == 1) return single_arg_abs_or_copy(res->data.function.args[0]);

    int64_t running_n = 0;
    int64_t running_d = 0;

    for (size_t i = 0; i < count; i++) {
        int64_t n, d;
        if (!is_rational(res->data.function.args[i], &n, &d)) {
            return NULL;
        }
        if (i == 0) {
            running_n = llabs(n);
            running_d = llabs(d);
        } else {
            if (running_n == 0 || n == 0) {
                running_n = 0;
                running_d = 1;
            } else {
                running_n = lcm(running_n, n);
                running_d = gcd(running_d, d);
            }
        }
    }

    return make_rational(running_n, running_d);
}

Expr* make_rational(int64_t n, int64_t d) {
    if (d == 0) return NULL; // Error
    if (n == 0) return expr_new_integer(0);
    
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
    return expr_new_function(expr_new_symbol("Rational"), args, 2);
}

Expr* builtin_rational(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_expr = res->data.function.args[0];
    Expr* d_expr = res->data.function.args[1];
    
    if (n_expr->type == EXPR_INTEGER && d_expr->type == EXPR_INTEGER) {
        int64_t n = n_expr->data.integer;
        int64_t d = d_expr->data.integer;
        if (d == 0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol("Indeterminate");
            }
            return expr_new_symbol("ComplexInfinity");
        }
        
        Expr* r = make_rational(n, d);
        if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL && r->data.function.head->data.symbol == SYM_Rational) {
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
        e->data.function.head->data.symbol == SYM_Rational) {
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
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        return true;
    }
    return false;
}

bool is_complex(Expr* e, Expr** re, Expr** im) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex &&
        e->data.function.arg_count == 2) {
        if (re) *re = e->data.function.args[0];
        if (im) *im = e->data.function.args[1];
        return true;
    }
    return false;
}

Expr* make_complex(Expr* re, Expr* im) {
    Expr* args[2] = { re, im };
    return expr_new_function(expr_new_symbol("Complex"), args, 2);
}

Expr* builtin_subtract(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    Expr* minus_one = expr_new_integer(-1);
    Expr* mb_args[2] = { minus_one, expr_copy(b) };
    Expr* minus_b = expr_new_function(expr_new_symbol("Times"), mb_args, 2);

    Expr* p_args[2] = { expr_copy(a), minus_b };
    return expr_new_function(expr_new_symbol("Plus"), p_args, 2);
}

Expr* builtin_complex(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* re = res->data.function.args[0];
    Expr* im = res->data.function.args[1];

    if (im->type == EXPR_INTEGER && im->data.integer == 0) {
        return expr_copy(re);
    }
    if (im->type == EXPR_REAL && im->data.real == 0.0) {
        if (re->type == EXPR_INTEGER) {
            return expr_new_real((double)re->data.integer);
        }
        return expr_copy(re);
    }

    return NULL;
}

Expr* builtin_divide(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];

    if (num->type == EXPR_REAL || den->type == EXPR_REAL) {
        double vnum = (num->type == EXPR_REAL) ? num->data.real : (num->type == EXPR_INTEGER) ? (double)num->data.integer : (num->type == EXPR_BIGINT) ? mpz_get_d(num->data.bigint) : 0.0;
        double vden = (den->type == EXPR_REAL) ? den->data.real : (den->type == EXPR_INTEGER) ? (double)den->data.integer : (den->type == EXPR_BIGINT) ? mpz_get_d(den->data.bigint) : 0.0;
        if (vden == 0.0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            return expr_new_symbol("ComplexInfinity");
        }
        return expr_new_real(vnum / vden);
    }

    int64_t n1, d1, n2, d2;
    if (is_rational(num, &n1, &d1) && is_rational(den, &n2, &d2)) {
        if (n2 == 0) {
            /* x / 0 with rational/integer x: 0/0 -> Indeterminate (handled in
             * Times when 0 multiplies ComplexInfinity); otherwise emit the
             * Power::infy message and yield ComplexInfinity. */
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n1 == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol("Indeterminate");
            }
            return expr_new_symbol("ComplexInfinity");
        }
        Expr* r = make_rational(n1 * d2, d1 * n2);
        if (r) return r;
    }

    Expr* minus_one = expr_new_integer(-1);
    Expr* p_args[2] = { expr_copy(den), minus_one };
    Expr* power = expr_new_function(expr_new_symbol("Power"), p_args, 2);
    
    Expr* t_args[2] = { expr_copy(num), power };
    return expr_new_function(expr_new_symbol("Times"), t_args, 2);
}




Expr* builtin_powermod(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 3) return NULL;
    
    Expr* a_expr = res->data.function.args[0];
    Expr* b_expr = res->data.function.args[1];
    Expr* m_expr = res->data.function.args[2];
    
    if (!(a_expr->type == EXPR_INTEGER || a_expr->type == EXPR_BIGINT) || 
        !(m_expr->type == EXPR_INTEGER || m_expr->type == EXPR_BIGINT)) return NULL;

    mpz_t a, m;
    expr_to_mpz(a_expr, a);
    expr_to_mpz(m_expr, m);
    if (mpz_cmp_ui(m, 0) == 0) {
        mpz_clear(a); mpz_clear(m);
        return NULL;
    }
    
    if (b_expr->type == EXPR_INTEGER || b_expr->type == EXPR_BIGINT) {
        mpz_t b, r;
        expr_to_mpz(b_expr, b);
        mpz_init(r);
        if (mpz_cmp_ui(b, 0) < 0) {
            mpz_t inv, b_neg;
            mpz_inits(inv, b_neg, NULL);
            if (mpz_invert(inv, a, m) == 0) {
                mpz_clears(inv, b_neg, b, a, m, r, NULL);
                return expr_copy(res);
            }
            mpz_neg(b_neg, b);
            mpz_powm(r, inv, b_neg, m);
            mpz_clears(inv, b_neg, NULL);
        } else {
            mpz_powm(r, a, b, m);
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clears(b, a, m, r, NULL);
        return out;
    }
    mpz_clears(a, m, NULL);
    return NULL;
}

Expr* builtin_factorial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        if (d == 1) {
            if (n < 0) return expr_new_symbol("ComplexInfinity");
            if (n <= 20) {
                int64_t f = 1;
                for (int64_t i = 2; i <= n; i++) f *= i;
                return expr_new_integer(f);
            } else {
                mpz_t result;
                mpz_init(result);
                mpz_fac_ui(result, (unsigned long)n);
                Expr* r = expr_new_bigint_from_mpz(result);
                mpz_clear(result);
                return r;
            }
        } else if (d == 2 || d == -2) {
            if (d == -2) { n = -n; }
            int64_t num = 1;
            int64_t den = 1;
            
            if (n > 0) {
                for (int64_t i = n; i >= 1; i -= 2) num *= i;
                den = 1LL << ((n + 1) / 2);
            } else {
                for (int64_t i = n + 2; i <= -1; i += 2) {
                    num *= 2;
                    den *= i;
                }
            }
            
            Expr* coeff = make_rational(num, den);
            if (!coeff) coeff = expr_new_integer(0);
            
            Expr* pi_sym = expr_new_symbol("Pi");
            Expr* half = make_rational(1, 2);
            Expr* sqrt_pi = eval_and_free(expr_new_function(expr_new_symbol("Power"), (Expr*[]){pi_sym, half}, 2));
            
            if (coeff->type == EXPR_INTEGER && coeff->data.integer == 1) {
                expr_free(coeff);
                return sqrt_pi;
            } else {
                return eval_and_free(expr_new_function(expr_new_symbol("Times"), (Expr*[]){coeff, sqrt_pi}, 2));
            }
        }
    }

    return NULL;
}

/* FactorialPower[n, k] -- the falling factorial (Pochhammer-style).
 *
 *   FactorialPower[n, k] = n * (n - 1) * (n - 2) * ... * (n - k + 1)
 *                        = n! / (n - k)!     (for non-negative integer k)
 *                        = 1                  (k == 0)
 *
 * Concrete numeric paths:
 *   - Both arguments integer with k >= 0: compute the product (GMP for
 *     k > 30 to avoid int64 overflow).
 *   - Integer n with negative integer k: undefined / no closed-form
 *     elementary identity in MMA's symbolic core; leave unevaluated.
 *
 * Symbolic n with concrete non-negative integer k expands to a product
 * of (n - i) factors so downstream Expand / D can act on it -- e.g.
 *   FactorialPower[n, 3] -> n (n - 1) (n - 2)
 * Mathematically equal under Mathematica's definition and lets the
 * surrounding evaluator collapse it via Expand[] when the user wants.
 *
 * Anything else (k symbolic, or both arguments symbolic): NULL --
 * stays as `FactorialPower[n, k]` until more information is supplied.
 */
Expr* builtin_factorialpower(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_arg = res->data.function.args[0];
    Expr* k_arg = res->data.function.args[1];

    /* k must be a concrete non-negative integer for any reduction. */
    if (k_arg->type != EXPR_INTEGER) return NULL;
    int64_t k = k_arg->data.integer;
    if (k < 0) return NULL;
    if (k == 0) return expr_new_integer(1);

    /* Both n and k integer: compute exact product in GMP. */
    if (n_arg->type == EXPR_INTEGER || n_arg->type == EXPR_BIGINT) {
        mpz_t n_z, term, product;
        mpz_init(n_z);
        if (n_arg->type == EXPR_INTEGER) mpz_set_si(n_z, n_arg->data.integer);
        else                              mpz_set(n_z, n_arg->data.bigint);
        mpz_init_set_ui(product, 1);
        mpz_init(term);
        for (int64_t i = 0; i < k; i++) {
            mpz_set(term, n_z);
            mpz_sub_ui(term, term, (unsigned long)i);
            mpz_mul(product, product, term);
        }
        mpz_clear(n_z);
        mpz_clear(term);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(product));
        mpz_clear(product);
        return result;
    }

    /* Symbolic n, k a small concrete integer (<= 32 to keep tree size
     * bounded). Expand to Times[n, n-1, ..., n-k+1]. */
    if (k <= 32) {
        Expr** factors = malloc(sizeof(Expr*) * (size_t)k);
        for (int64_t i = 0; i < k; i++) {
            if (i == 0) {
                factors[i] = expr_copy(n_arg);
            } else {
                Expr* shift = expr_new_integer(-i);
                factors[i] = expr_new_function(expr_new_symbol("Plus"),
                    (Expr*[]){ expr_copy(n_arg), shift }, 2);
            }
        }
        Expr* product = expr_new_function(expr_new_symbol("Times"), factors, (size_t)k);
        free(factors);
        return eval_and_free(product);
    }
    return NULL;
}

/* Factorial2[n] -- the double factorial.
 *
 * For non-negative integer n:
 *   n!! = n * (n-2) * (n-4) * ...    (terminating at 1 if n odd, 2 if n even).
 *   0!! = 1, (-1)!! = 1.
 * For negative odd integers n <= -3:
 *   n!! = ComplexInfinity (poles, like Factorial of negative integers).
 *
 * Symbolic / non-integer arguments return NULL (the call stays as-is)
 * so the surrounding evaluator can keep them in expression form.
 *
 * Implemented in two paths: int64 fast-path for small n (<= 30 -- fits in
 * int64 for all valid inputs since 31!! ~= 7.9e9 fits), and a GMP path
 * for larger n. */
Expr* builtin_factorial2(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    /* Special-case the small literal -1 / 0 inputs. */
    if (arg->type == EXPR_INTEGER) {
        int64_t n = arg->data.integer;
        if (n == -1 || n == 0) return expr_new_integer(1);
        if (n < 0) {
            /* Even negative integers: ComplexInfinity (pole). Odd negative
             * integers below -1 also diverge. */
            return expr_new_symbol("ComplexInfinity");
        }
        if (n <= 30) {
            int64_t f = 1;
            for (int64_t i = n; i >= 1; i -= 2) f *= i;
            return expr_new_integer(f);
        }
        /* GMP path. */
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (int64_t i = n; i >= 1; i -= 2) {
            mpz_set_si(k, (long)i);
            mpz_mul(r, r, k);
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    if (arg->type == EXPR_BIGINT) {
        /* For arguments large enough to be BigInt, use GMP. We only
         * handle non-negative inputs; negative BigInt is ComplexInfinity. */
        if (mpz_sgn(arg->data.bigint) < 0) {
            return expr_new_symbol("ComplexInfinity");
        }
        if (mpz_sgn(arg->data.bigint) == 0) return expr_new_integer(1);
        /* Cap at a sane size: if mpz_fits_ulong_p fails the input is too
         * large for our explicit loop. Return NULL to leave it symbolic. */
        if (!mpz_fits_ulong_p(arg->data.bigint)) return NULL;
        unsigned long n = mpz_get_ui(arg->data.bigint);
        mpz_t r, k;
        mpz_init_set_ui(r, 1);
        mpz_init(k);
        for (unsigned long i = n; i >= 1; i -= 2) {
            mpz_set_ui(k, i);
            mpz_mul(r, r, k);
            if (i < 2) break; /* guard underflow when i==1 */
        }
        Expr* out = expr_new_bigint_from_mpz(r);
        mpz_clears(r, k, NULL);
        return out;
    }

    return NULL; /* symbolic -- leave Factorial2[expr] alone */
}

Expr* builtin_binomial(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* arg_n = res->data.function.args[0];
    Expr* arg_m = res->data.function.args[1];

    if (arg_n->type == EXPR_REAL || arg_m->type == EXPR_REAL) {
        double n_val = (arg_n->type == EXPR_REAL) ? arg_n->data.real : (double)arg_n->data.integer;
        double m_val = (arg_m->type == EXPR_REAL) ? arg_m->data.real : (double)arg_m->data.integer;
        double val = tgamma(n_val + 1.0) / (tgamma(m_val + 1.0) * tgamma(n_val - m_val + 1.0));
        return expr_new_real(val);
    }

    int64_t n_num = 0, n_den = 1, m_num = 0, m_den = 1;
    bool n_is_rat = false;
    if (arg_n->type == EXPR_INTEGER || arg_n->type == EXPR_BIGINT) { n_is_rat = true; }
    else if (is_rational(arg_n, &n_num, &n_den)) { n_is_rat = true; }

    bool m_is_rat = false;
    if (arg_m->type == EXPR_INTEGER || arg_m->type == EXPR_BIGINT) { m_is_rat = true; }
    else if (is_rational(arg_m, &m_num, &m_den)) { m_is_rat = true; }
    
    if (n_is_rat && m_is_rat) {
        if ((arg_n->type == EXPR_INTEGER || arg_n->type == EXPR_BIGINT) && 
            (arg_m->type == EXPR_INTEGER || arg_m->type == EXPR_BIGINT)) {
            
            mpz_t n, m;
            expr_to_mpz(arg_n, n);
            expr_to_mpz(arg_m, m);
            
            if (mpz_cmp_ui(m, 0) < 0) {
                mpz_clear(n); mpz_clear(m);
                return expr_new_integer(0);
            }
            if (mpz_cmp_ui(n, 0) >= 0 && mpz_cmp(m, n) > 0) {
                mpz_clear(n); mpz_clear(m);
                return expr_new_integer(0);
            }
            
            mpz_t r;
            mpz_init(r);
            
            if (mpz_cmp_ui(n, 0) < 0) {
                // Binomial[n, m] = (-1)^m * Binomial[-n+m-1, m] for n < 0
                mpz_t temp_n;
                mpz_init(temp_n);
                mpz_neg(temp_n, n);
                mpz_add(temp_n, temp_n, m);
                mpz_sub_ui(temp_n, temp_n, 1);
                
                mpz_bin_ui(r, temp_n, mpz_get_ui(m));
                
                if (mpz_odd_p(m)) {
                    mpz_neg(r, r);
                }
                mpz_clear(temp_n);
            } else {
                mpz_bin_ui(r, n, mpz_get_ui(m));
            }
            
            Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
            mpz_clears(n, m, r, NULL);
            return out;
        } else {
            // Handle fractional case...
            if (arg_n->type == EXPR_INTEGER) { n_num = arg_n->data.integer; n_den = 1; }
            if (arg_m->type == EXPR_INTEGER) { m_num = arg_m->data.integer; m_den = 1; }
            int64_t n = n_num; int64_t m = m_num;
            if (m_den == 2 && n_den == 2) {
                if (m > 0 && n > 0 && m <= n) {
                    int64_t num = 1;
                    for (int64_t i = 1; i <= m; i += 2) {
                        num *= (n - i + 2);
                    }
                    int64_t den = 1;
                    for (int64_t i = 1; i <= m; i += 2) {
                        den *= (i + 1);
                    }
                    return make_rational(num, den);
                }
            } else if (n_den == 1 && m_den == 1) {
                // Fallback for huge rationals?
            }
        }
    }
    return NULL;
}

bool is_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

bool is_complex_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_ComplexInfinity;
}

bool is_indeterminate_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Indeterminate;
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
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        // d is conventionally positive in Mathilda Rational[n, d].
        if (n > 0) return (d > 0) ? 1 : -1;
        if (n < 0) return (d > 0) ? -1 : 1;
        return 0;
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
        e->data.function.head->data.symbol == SYM_Times &&
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
    if (e->data.function.head->data.symbol != SYM_Times) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (!is_infinity_sym(b)) return false;
    return expr_numeric_sign(a) < 0;
}
