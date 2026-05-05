#include "eval.h"
#include "power.h"
#include "arithmetic.h"
#include "times.h"
#include "numeric.h"
#include "sym_names.h"
#include "trig_canon.h"
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

static int64_t ipow(int64_t base, int64_t exp, bool* overflow) {
    if (exp < 0) return 0;
    if (exp == 0) return 1;
    if (base == 0) return 0;
    if (base == 1) return 1;
    if (base == -1) return (exp % 2 == 0) ? 1 : -1;

    __int128_t res = 1;
    __int128_t b = base;
    while (exp > 0) {
        if (exp % 2 == 1) {
            res *= b;
            if (res > INT64_MAX || res < INT64_MIN) {
                *overflow = true;
                return 0;
            }
        }
        exp /= 2;
        if (exp > 0) {
            if (b > INT64_MAX || b < INT64_MIN) {
                *overflow = true;
                return 0;
            }
            b *= b;
        }
    }
    return (int64_t)res;
}

static void factor_out_kth_power(int64_t n, int64_t k, int64_t* out_m, int64_t* out_r) {
    int64_t m = 1;
    int64_t r = n;
    int64_t d = 2;
    while (d * d <= r) {
        bool ov = false;
        int64_t dk = ipow(d, k, &ov);
        if (ov || dk <= 0 || dk > r) {
            if (ov || dk > r) break;
            d++;
            continue;
        }
        while (r % dk == 0) {
            m *= d;
            r /= dk;
        }
        d++;
    }
    *out_m = m;
    *out_r = r;
}

static Expr* bigint_pow(const Expr* base, int64_t exp) {
    if (exp < 0) return NULL;
    mpz_t b, r;
    expr_to_mpz(base, b);
    mpz_init(r);
    mpz_pow_ui(r, b, (unsigned long)exp);
    mpz_clear(b);
    Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
    mpz_clear(r);
    return result;
}

Expr* make_power(Expr* base, Expr* exp) {
    Expr* args[2] = { base, exp };
    return expr_new_function(expr_new_symbol("Power"), args, 2);
}

static bool is_head_call(Expr* e, const char* name, size_t argc) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0 &&
           e->data.function.arg_count == argc;
}

/* Cancel Log inside a Power exponent:
 *   Exp[c Log[a]]                  -> a^c                   (base == E)
 *   Power[base, c Log[base, a]]    -> a^c                   (general)
 *
 * Internally Log[b, a] is represented as Log[a] * Log[b]^(-1), so case 2
 * reduces to matching Log[a] * Power[Log[base], -1] among the exponent's
 * Times-factors. We require exactly one Log[a] factor and, when base != E,
 * exactly one matching Power[Log[base], -1] factor; the remaining factors
 * become the coefficient c in the rewritten Power[a, c]. */
static Expr* simplify_exp_log(Expr* base, Expr* exp) {
    bool exp_is_times = exp->type == EXPR_FUNCTION &&
                        exp->data.function.head->type == EXPR_SYMBOL &&
                        exp->data.function.head->data.symbol == SYM_Times;
    size_t nf = exp_is_times ? exp->data.function.arg_count : 1;
    Expr** factors = exp_is_times ? exp->data.function.args : &exp;

    bool base_is_E = (base->type == EXPR_SYMBOL && base->data.symbol == SYM_E);

    int log_idx = -1;
    Expr* a = NULL;
    for (size_t i = 0; i < nf; i++) {
        if (!is_head_call(factors[i], "Log", 1)) continue;
        log_idx = (int)i;
        a = factors[i]->data.function.args[0];
        break;
    }
    if (log_idx < 0) return NULL;

    int inv_log_idx = -1;
    if (!base_is_E) {
        for (size_t i = 0; i < nf; i++) {
            if ((int)i == log_idx) continue;
            Expr* f = factors[i];
            if (!is_head_call(f, "Power", 2)) continue;
            Expr* ib = f->data.function.args[0];
            Expr* ie = f->data.function.args[1];
            if (!(ie->type == EXPR_INTEGER && ie->data.integer == -1)) continue;
            if (!is_head_call(ib, "Log", 1)) continue;
            if (!expr_eq(ib->data.function.args[0], base)) continue;
            inv_log_idx = (int)i;
            break;
        }
        if (inv_log_idx < 0) return NULL;
    }

    Expr* coeff;
    if (!exp_is_times) {
        coeff = expr_new_integer(1);
    } else {
        Expr** rest = malloc(sizeof(Expr*) * nf);
        size_t kept = 0;
        for (size_t i = 0; i < nf; i++) {
            if ((int)i == log_idx) continue;
            if ((int)i == inv_log_idx) continue;
            rest[kept++] = expr_copy(factors[i]);
        }
        if (kept == 0) {
            free(rest);
            coeff = expr_new_integer(1);
        } else if (kept == 1) {
            coeff = rest[0];
            free(rest);
        } else {
            coeff = expr_new_function(expr_new_symbol("Times"), rest, kept);
            free(rest);
        }
    }

    return eval_and_free(expr_new_function(expr_new_symbol("Power"),
        (Expr*[]){ expr_copy(a), coeff }, 2));
}

Expr* builtin_power(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    size_t n_args = res->data.function.arg_count;
    if (n_args == 0) return NULL;
    if (n_args == 1) return expr_copy(res->data.function.args[0]);

    if (n_args > 2) {
        Expr** sub_args = malloc(sizeof(Expr*) * (n_args - 1));
        for (size_t i = 0; i < n_args - 1; i++) {
            sub_args[i] = expr_copy(res->data.function.args[i+1]);
        }
        Expr* sub_power = expr_new_function(expr_new_symbol("Power"), sub_args, n_args - 1);
        free(sub_args);
        
        Expr* final_args[2] = { expr_copy(res->data.function.args[0]), sub_power };
        return expr_new_function(expr_new_symbol("Power"), final_args, 2);
    }

    if (n_args == 1) return NULL;

    Expr* base = res->data.function.args[0];
    Expr* exp = res->data.function.args[1];

    /* Infinity / Indeterminate preprocessing.
     *
     * Mathematica semantics for Power (selected):
     *   Indeterminate ^ x, x ^ Indeterminate          -> Indeterminate
     *   1 ^ Infinity, 1 ^ -Infinity, 1 ^ ComplexInfinity -> Indeterminate (msg)
     *   Infinity ^ 0, (-Infinity) ^ 0, ComplexInfinity ^ 0 -> Indeterminate (msg)
     *   0 ^ Infinity                                   -> 0
     *   0 ^ -Infinity                                  -> ComplexInfinity (msg)
     *   0 ^ ComplexInfinity                            -> Indeterminate (msg)
     *   Infinity ^ Infinity                            -> ComplexInfinity
     *   Infinity ^ -Infinity                           -> 0
     *   Infinity ^ n   (n numeric, n > 0)              -> Infinity
     *   Infinity ^ n   (n numeric, n < 0)              -> 0
     *   ComplexInfinity ^ n (n numeric, n > 0)         -> ComplexInfinity
     *   ComplexInfinity ^ n (n numeric, n < 0)         -> 0
     *
     * Note: 0^positive and 0^negative for ordinary exponents are still handled
     * by the existing logic below; only the infinity-flavoured cases are
     * intercepted here.
     */
    {
        bool b_indet = is_indeterminate_sym(base);
        bool e_indet = is_indeterminate_sym(exp);
        if (b_indet || e_indet) return expr_new_symbol("Indeterminate");

        bool b_inf  = is_infinity_sym(base);
        bool b_ninf = is_neg_infinity_form(base);
        bool b_cinf = is_complex_infinity_sym(base);
        bool e_inf  = is_infinity_sym(exp);
        bool e_ninf = is_neg_infinity_form(exp);
        bool e_cinf = is_complex_infinity_sym(exp);

        bool base_is_one  = (base->type == EXPR_INTEGER && base->data.integer == 1);
        bool base_is_zero_lit = (base->type == EXPR_INTEGER && base->data.integer == 0) ||
                                (base->type == EXPR_REAL && base->data.real == 0.0) ||
                                (base->type == EXPR_BIGINT && mpz_sgn(base->data.bigint) == 0);
        bool exp_is_zero_lit  = (exp->type == EXPR_INTEGER && exp->data.integer == 0) ||
                                (exp->type == EXPR_REAL && exp->data.real == 0.0) ||
                                (exp->type == EXPR_BIGINT && mpz_sgn(exp->data.bigint) == 0);

        bool base_is_inf = b_inf || b_ninf || b_cinf;
        bool exp_is_inf  = e_inf || e_ninf || e_cinf;

        if (base_is_one && exp_is_inf) {
            const char* what = e_inf ? "Infinity" : (e_ninf ? "-Infinity" : "ComplexInfinity");
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression 1^%s encountered.\n", what);
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_inf && exp_is_zero_lit) {
            const char* what = b_inf ? "Infinity" : (b_ninf ? "-Infinity" : "ComplexInfinity");
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression %s^0 encountered.\n", what);
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_zero_lit && e_cinf) {
            if (!arith_warnings_muted())
                fprintf(stderr,
                    "Infinity::indet: Indeterminate expression 0^ComplexInfinity encountered.\n");
            return expr_new_symbol("Indeterminate");
        }

        if (base_is_zero_lit && e_inf) {
            return expr_new_integer(0);
        }

        if (base_is_zero_lit && e_ninf) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            return expr_new_symbol("ComplexInfinity");
        }

        if (b_inf) {
            if (e_inf)  return expr_new_symbol("ComplexInfinity");
            if (e_ninf) return expr_new_integer(0);
            if (e_cinf) return expr_new_symbol("Indeterminate");
            int es = expr_numeric_sign(exp);
            if (es > 0) return expr_new_symbol("Infinity");
            if (es < 0) return expr_new_integer(0);
        }

        if (b_cinf) {
            if (e_inf)  return expr_new_symbol("ComplexInfinity");
            if (e_ninf) return expr_new_integer(0);
            int es = expr_numeric_sign(exp);
            if (es > 0) return expr_new_symbol("ComplexInfinity");
            if (es < 0) return expr_new_integer(0);
        }
    }

    if (exp->type == EXPR_INTEGER && exp->data.integer == 0) return expr_new_integer(1);
    if (exp->type == EXPR_INTEGER && exp->data.integer == 1) return expr_copy(base);
    if (base->type == EXPR_INTEGER && base->data.integer == 1) return expr_new_integer(1);

    bool base_is_zero = false;
    if (base->type == EXPR_INTEGER && base->data.integer == 0) base_is_zero = true;
    if (base->type == EXPR_REAL && base->data.real == 0.0) base_is_zero = true;

    bool exp_is_negative = false;
    bool exp_is_positive = false;
    if (exp->type == EXPR_INTEGER) {
        if (exp->data.integer < 0) exp_is_negative = true;
        else if (exp->data.integer > 0) exp_is_positive = true;
    }
    if (exp->type == EXPR_REAL) {
        if (exp->data.real < 0.0) exp_is_negative = true;
        else if (exp->data.real > 0.0) exp_is_positive = true;
    }
    int64_t rn, rd;
    if (is_rational(exp, &rn, &rd)) {
        if (rn < 0) exp_is_negative = true;
        else if (rn > 0) exp_is_positive = true;
    }

    if (base_is_zero && exp_is_negative) {
        if (!arith_warnings_muted())
            fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
        return expr_new_symbol("ComplexInfinity");
    }
    /* 0^positive -> 0 (includes 0^(1/2) = Sqrt[0] = 0, 0^0.5 = 0, etc.).
     * Preserves the base's numeric type: real 0.0 for real exponents or a
     * real-base zero, integer 0 otherwise. */
    if (base_is_zero && exp_is_positive) {
        if (base->type == EXPR_REAL || exp->type == EXPR_REAL) return expr_new_real(0.0);
        return expr_new_integer(0);
    }

    Expr *re_b = NULL, *im_b = NULL, *re_e = NULL, *im_e = NULL;
    bool base_comp = is_complex(base, &re_b, &im_b);
    bool exp_comp = is_complex(exp, &re_e, &im_e);

    if (base_comp || exp_comp ||
        base->type == EXPR_REAL || exp->type == EXPR_REAL ||
        base->type == EXPR_BIGINT || exp->type == EXPR_BIGINT
#ifdef USE_MPFR
        || base->type == EXPR_MPFR || exp->type == EXPR_MPFR
#endif
        ) {
        bool base_num = expr_is_numeric_like(base) || base_comp;
        bool exp_num  = expr_is_numeric_like(exp)  || exp_comp;

        if (base_num && exp_num) {
#ifdef USE_MPFR
            /* MPFR real-valued power: fast path when both sides are real
             * and at least one carries MPFR precision. The complex case
             * (negative base with fractional exponent, or Complex operands)
             * still falls through to the cpow path below — true MPC
             * complex arithmetic is out of scope for Phase 2. */
            if (!base_comp && !exp_comp && numeric_any_mpfr(base, exp)) {
                /* Only take the real MPFR path when the result is real:
                 * base > 0 OR exponent is an integer. */
                Expr *b_re_check = base;
                bool base_nonneg = (b_re_check->type == EXPR_INTEGER && b_re_check->data.integer >= 0)
                                   || (b_re_check->type == EXPR_REAL && b_re_check->data.real >= 0.0)
                                   || (b_re_check->type == EXPR_BIGINT && mpz_sgn(b_re_check->data.bigint) >= 0)
                                   || (b_re_check->type == EXPR_MPFR && mpfr_sgn(b_re_check->data.mpfr) >= 0);
                int64_t en, ed;
                bool exp_is_int = (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT
                                   || (is_rational(exp, &en, &ed) && ed == 1));
                if (base_nonneg || exp_is_int) {
                    Expr* r = numeric_mpfr_pow(base, exp, 0);
                    if (r) return r;
                }
            }
#endif
            bool has_real = (base->type == EXPR_REAL || exp->type == EXPR_REAL ||
                            (base_comp && (re_b->type == EXPR_REAL || im_b->type == EXPR_REAL)) ||
                            (exp_comp && (re_e->type == EXPR_REAL || im_e->type == EXPR_REAL)));
#ifdef USE_MPFR
            if (!has_real) {
                /* Treat an MPFR operand as "real" for the double-complex
                 * fallback path too. */
                has_real = (base->type == EXPR_MPFR || exp->type == EXPR_MPFR);
            }
#endif
            if (has_real) {
                double vbase_re = 0, vbase_im = 0, vexp_re = 0, vexp_im = 0;
                int64_t n, d;
                
                if (base_comp) {
                    vbase_re = (re_b->type == EXPR_REAL) ? re_b->data.real : (re_b->type == EXPR_INTEGER) ? (double)re_b->data.integer : (re_b->type == EXPR_BIGINT) ? mpz_get_d(re_b->data.bigint) : 0;
                    vbase_im = (im_b->type == EXPR_REAL) ? im_b->data.real : (im_b->type == EXPR_INTEGER) ? (double)im_b->data.integer : (im_b->type == EXPR_BIGINT) ? mpz_get_d(im_b->data.bigint) : 0;
                    if (is_rational(re_b, &n, &d)) vbase_re = (double)n / d;
                    if (is_rational(im_b, &n, &d)) vbase_im = (double)n / d;
                } else {
                    vbase_re = (base->type == EXPR_REAL) ? base->data.real : (base->type == EXPR_INTEGER) ? (double)base->data.integer : (base->type == EXPR_BIGINT) ? mpz_get_d(base->data.bigint) : 0;
                    if (is_rational(base, &n, &d)) vbase_re = (double)n / d;
                }
                
                if (exp_comp) {
                    vexp_re = (re_e->type == EXPR_REAL) ? re_e->data.real : (re_e->type == EXPR_INTEGER) ? (double)re_e->data.integer : (re_e->type == EXPR_BIGINT) ? mpz_get_d(re_e->data.bigint) : 0;
                    vexp_im = (im_e->type == EXPR_REAL) ? im_e->data.real : (im_e->type == EXPR_INTEGER) ? (double)im_e->data.integer : (im_e->type == EXPR_BIGINT) ? mpz_get_d(im_e->data.bigint) : 0;
                    if (is_rational(re_e, &n, &d)) vexp_re = (double)n / d;
                    if (is_rational(im_e, &n, &d)) vexp_im = (double)n / d;
                } else {
                    vexp_re = (exp->type == EXPR_REAL) ? exp->data.real : (exp->type == EXPR_INTEGER) ? (double)exp->data.integer : (exp->type == EXPR_BIGINT) ? mpz_get_d(exp->data.bigint) : 0;
                    if (is_rational(exp, &n, &d)) vexp_re = (double)n / d;
                }

                double complex z = vbase_re + vbase_im * I;
                double complex w = vexp_re + vexp_im * I;
                
                if (vbase_re < 0.0 && vbase_im == 0.0 && vexp_im == 0.0 && floor(vexp_re) != vexp_re) {
                    double r = pow(-vbase_re, vexp_re);
                    double theta = vexp_re * 3.14159265358979323846;
                    Expr* c_re = expr_new_real(r * cos(theta));
                    Expr* c_im = expr_new_real(r * sin(theta));
                    Expr* comp = make_complex(c_re, c_im);
                    return comp;
                }
                
                double complex result = cpow(z, w);
                double r_re = creal(result);
                double r_im = cimag(result);
                
                if (r_im == 0.0 || (fabs(r_im) < 1e-14 && fabs(r_re) > 1e-14)) {
                    return expr_new_real(r_re);
                }
                
                Expr* c_re = expr_new_real(r_re);
                Expr* c_im = expr_new_real(r_im);
                Expr* comp = make_complex(c_re, c_im);
                return comp;
            }
        }
    }
    
    if (base_comp && exp->type == EXPR_INTEGER) {
        int64_t e = exp->data.integer;
        if (e < 0 && e > -1000) {
            // e < 0, compute (base^-e)^-1
            Expr* pos_pow = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(base), expr_new_integer(-e)}, 2);
            Expr* denom = evaluate(pos_pow);
            expr_free(pos_pow);
            
            // Now we have denom = a + b I. We need 1 / denom.
            // We can emit: Divide[1, denom], which is: Times[1, Power[denom, -1]] -> WAIT that loops!
            // No, we need to construct the conjugate division manually!
            Expr* re; Expr* im;
            if (is_complex(denom, &re, &im)) {
                // conj = a - b I
                Expr* conj = make_complex(expr_copy(re), expr_new_function(expr_new_symbol("Times"), (Expr*[]){expr_new_integer(-1), expr_copy(im)}, 2));
                // mag_sq = a^2 + b^2
                Expr* a2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(re), expr_new_integer(2)}, 2);
                Expr* b2 = expr_new_function(expr_new_symbol("Power"), (Expr*[]){expr_copy(im), expr_new_integer(2)}, 2);
                Expr* mag_sq = expr_new_function(expr_new_symbol("Plus"), (Expr*[]){a2, b2}, 2);
                
                // Result = conj / mag_sq
                Expr* res_ast = expr_new_function(expr_new_symbol("Divide"), (Expr*[]){conj, mag_sq}, 2);
                Expr* result = evaluate(res_ast);
                expr_free(res_ast);
                expr_free(denom);
                return result;
            } else {
                // If it evaluated to Real or Integer, just do Power[denom, -1] natively
                Expr* res_ast = expr_new_function(expr_new_symbol("Divide"), (Expr*[]){expr_new_integer(1), expr_copy(denom)}, 2);
                Expr* result = evaluate(res_ast);
                expr_free(res_ast);
                expr_free(denom);
                return result;
            }
        }
        if (e >= 0 && e < 1000) { 
            Expr** prod_args = malloc(sizeof(Expr*) * (size_t)e);
            for(int64_t i=0; i<e; i++) prod_args[i] = expr_copy(base);
            Expr* prod = expr_new_function(expr_new_symbol("Times"), prod_args, (size_t)e);
            free(prod_args);
            return prod;
        }
    }

    if (base->type == EXPR_BIGINT && exp->type == EXPR_INTEGER) {
        int64_t e = exp->data.integer;
        if (e > 0) {
            return bigint_pow(base, e);
        }
        if (e < 0) {
            Expr* denom = bigint_pow(base, -e);
            if (!denom) return NULL;
            Expr* r_args[2] = { expr_new_integer(1), denom };
            return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
        }
    }

    if (base->type == EXPR_INTEGER && exp->type == EXPR_INTEGER) {
        int64_t b = base->data.integer;
        int64_t e = exp->data.integer;
        if (e > 0) {
            bool overflow = false;
            int64_t res_val = ipow(b, e, &overflow);
            if (overflow) return bigint_pow(base, e);
            return expr_new_integer(res_val);
        }
        if (e < 0) {
            bool overflow = false;
            int64_t res_val = ipow(b, -e, &overflow);
            if (overflow) {
                Expr* denom = bigint_pow(base, -e);
                if (!denom) return NULL;
                Expr* r_args[2] = { expr_new_integer(1), denom };
                return expr_new_function(expr_new_symbol("Rational"), r_args, 2);
            }
            return make_rational(1, res_val);
        }
    }

    int64_t bn, bd, e;
    if (is_rational(base, &bn, &bd) && exp->type == EXPR_INTEGER) {
        e = exp->data.integer;
        if (e == 0) return expr_new_integer(1);
        if (e > 0) {
            bool overflow = false;
            int64_t res_n = ipow(bn, e, &overflow);
            int64_t res_d = ipow(bd, e, &overflow);
            if (!overflow) return make_rational(res_n, res_d);
        } else {
            bool overflow = false;
            int64_t res_n = ipow(bn, -e, &overflow);
            int64_t res_d = ipow(bd, -e, &overflow);
            if (!overflow) return make_rational(res_d, res_n);
        }
    }

    if (exp->type == EXPR_INTEGER && base->type == EXPR_FUNCTION && base->data.function.head->data.symbol == SYM_Times) {
        size_t bc = base->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * bc);
        for (size_t i = 0; i < bc; i++) {
            Expr* p_args[2] = { expr_copy(base->data.function.args[i]), expr_copy(exp) };
            new_args[i] = expr_new_function(expr_new_symbol("Power"), p_args, 2);
        }
        Expr* result = expr_new_function(expr_new_symbol("Times"), new_args, bc);
        free(new_args);
        return result;
    }

    if (exp->type == EXPR_INTEGER && base->type == EXPR_FUNCTION && base->data.function.head->data.symbol == SYM_Power) {
        if (base->data.function.arg_count == 2) {
            Expr* inner_base = base->data.function.args[0];
            Expr* inner_exp = base->data.function.args[1];
            Expr* t_args[2] = { expr_copy(inner_exp), expr_copy(exp) };
            Expr* new_exp = expr_new_function(expr_new_symbol("Times"), t_args, 2);
            Expr* p_args[2] = { expr_copy(inner_base), new_exp };
            return expr_new_function(expr_new_symbol("Power"), p_args, 2);
        }
    }

    {
        Expr* simp = simplify_exp_log(base, exp);
        if (simp) return simp;
    }

    int64_t p, q;
    if (base->type == EXPR_INTEGER && is_rational(exp, &p, &q) && q > 1) {
        int64_t n = base->data.integer;
        if (n < 0) {
            if (q == 2) {
                Expr* i_val = expr_new_function(expr_new_symbol("Complex"), (Expr*[]){expr_new_integer(0), expr_new_integer(1)}, 2);
                Expr* pos_base = expr_new_integer(-n);
                Expr* tmp_p_args[2] = { pos_base, expr_copy(exp) };
                Expr* tmp_power = expr_new_function(expr_new_symbol("Power"), tmp_p_args, 2);
                Expr* rest = builtin_power(tmp_power);
                if (!rest) rest = tmp_power; else expr_free(tmp_power);
                if (p == 1) {
                    Expr* t_args[2] = { i_val, rest };
                    return expr_new_function(expr_new_symbol("Times"), t_args, 2);
                } else {
                    expr_free(i_val);
                    expr_free(rest);
                    return NULL;
                }
            }
            return NULL;
        }
        int64_t m, r;
        factor_out_kth_power(n, q, &m, &r);
        if (m > 1) {
            Expr* coeff = NULL;
            if (p == 1) {
                coeff = expr_new_integer(m);
            } else {
                Expr* m_p_args[2] = { expr_new_integer(m), expr_new_integer(p) };
                Expr* tmp_m_power = expr_new_function(expr_new_symbol("Power"), m_p_args, 2);
                coeff = builtin_power(tmp_m_power);
                if (!coeff) coeff = tmp_m_power; else expr_free(tmp_m_power);
            }
            Expr* residue = NULL;
            if (r == 1) {
                residue = expr_new_integer(1);
            } else {
                Expr* r_p_args[2] = { expr_new_integer(r), expr_copy(exp) };
                residue = expr_new_function(expr_new_symbol("Power"), r_p_args, 2);
            }
            Expr* t_args[2] = { coeff, residue };
            return expr_new_function(expr_new_symbol("Times"), t_args, 2);
        }
    }

    /* Trig / hyperbolic reciprocal naming: Power[Cos[x], -k] -> Sec[x]^k, etc.
     * The parser collapses 1/Cos[x] to Power[Cos[x], -1] without going through
     * Times, so the Times-level canonicalizer never sees these solo cases. */
    if (exp->type == EXPR_INTEGER) {
        Expr* tc = trig_canon_power(base, exp->data.integer);
        if (tc) return tc;
    }

    return NULL;
}

Expr* builtin_sqrt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr* half = make_rational(1, 2);
    Expr* p_args[2] = { expr_copy(arg), half };
    return expr_new_function(expr_new_symbol("Power"), p_args, 2);
}
