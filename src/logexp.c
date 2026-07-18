#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "logexp.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "complex.h"
#include "numeric.h"
#include "numeric_complex.h"
#include "sym_names.h"
#include "common.h"

/*
 * logexp_init:
 * Initializes the logarithmic and exponential built-in functions in Mathilda.
 * Adds the 'Log' and 'Exp' C functions to the symbol table and protects the constant 'E'.
 */
void logexp_init(void) {
    symtab_add_builtin("Log", builtin_log);
    symtab_add_builtin("Exp", builtin_exp);
    symtab_get_def("E")->attributes |= ATTR_PROTECTED;
}

/*
 * get_approx:
 * Attempts to extract a numeric approximation (double complex) from a given expression.
 * Handles EXPR_REAL, EXPR_INTEGER, Rational[n, d], and Complex[re, im] formats.
 * Returns true if a valid numeric approximation could be extracted, storing it in *out.
 */
static bool get_approx(Expr* e, double complex* out, bool* is_inexact) {
    if (is_inexact) *is_inexact = true;
    if (e->type == EXPR_REAL) {
        *out = e->data.real + 0.0 * I;
        return true;
    }
    
    Expr* re; Expr* im;
    // Check if the expression is a complex number representation
    if (is_complex(e, &re, &im)) {
        bool has_real = false;
        double r = 0.0, i = 0.0;
        int64_t n, d;
        
        // Extract real part
        if (re->type == EXPR_REAL) { r = re->data.real; has_real = true; }
        else if (re->type == EXPR_INTEGER) r = (double)re->data.integer;
        else if (is_rational(re, &n, &d)) r = (double)n / d;
        else return false;
        
        // Extract imaginary part
        if (im->type == EXPR_REAL) { i = im->data.real; has_real = true; }
        else if (im->type == EXPR_INTEGER) i = (double)im->data.integer;
        else if (is_rational(im, &n, &d)) i = (double)n / d;
        else return false;
        
        if (has_real) {
            *out = r + i * I;
            return true;
        }
    }
    return false;
}

/*
 * is_infinity:
 * Checks if the given expression represents the symbol 'Infinity'.
 */
static bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

/*
 * is_minus_infinity:
 * Checks if the given expression represents '-Infinity',
 * which is represented internally as Times[-1, Infinity].
 */
static bool is_minus_infinity(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 && 
        e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol.name == SYM_Times) {
        Expr* a1 = e->data.function.args[0];
        Expr* a2 = e->data.function.args[1];
        if (a1->type == EXPR_INTEGER && a1->data.integer == -1 && is_infinity(a2)) return true;
        if (a2->type == EXPR_INTEGER && a2->data.integer == -1 && is_infinity(a1)) return true;
    }
    return false;
}

/*
 * make_minus_infinity:
 * Constructs and returns a new expression representing '-Infinity'.
 */
static Expr* make_minus_infinity() {
    Expr* args[2] = { expr_new_integer(-1), expr_new_symbol(SYM_Infinity) };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* True when e is a concrete real numeric value (int, bigint, rational, real).
 * Used to guard Log[E^k] -> k and Log[b, b^k] -> k: on the principal branch
 * these identities hold for any real k, but for complex k the result can
 * differ from k by 2 pi i / Log[b]. Restricting k to real numerics keeps us
 * safe without a full assumption system. */
static bool is_real_numeric_expr(Expr* e) {
    if (e->type == EXPR_INTEGER) return true;
    if (e->type == EXPR_BIGINT)  return true;
    if (e->type == EXPR_REAL)    return true;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return true;
    return false;
}

/* True when e is known to be strictly positive: a positive numeric value or
 * a symbol whose value is intrinsically positive (E, Pi). Used to guard
 * Log[b, b^k] -> k so we don't cross a branch cut. */
static bool is_positive_known(Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) > 0;
    if (e->type == EXPR_REAL)    return e->data.real > 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return n > 0;
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        if (s == SYM_E)  return true;
        if (s == SYM_Pi) return true;
    }
    return false;
}

/* True when e has the shape Power[_, _] (a two-argument Power call). */
static bool is_power_call(Expr* e) {
    return e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Power &&
           e->data.function.arg_count == 2;
}

/*
 * builtin_log:
 * Implements the evaluation logic for the 'Log' function.
 * Supports Log[z] (natural logarithm) and Log[b, z] (base-b logarithm).
 */
Expr* builtin_log(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    // Wrong arity: emit Mathematica's `Log::argt` diagnostic and leave
    // the call unevaluated. Log accepts 1 (natural log) or 2 (base-b log).
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2)
        return builtin_arg_error("Log", argc, 1, 2);

    // Log[z] - Natural logarithm
    if (argc == 1) {
        Expr* z = res->data.function.args[0];

        // Exact evaluations for special constants
        if (z->type == EXPR_INTEGER && z->data.integer == 0) {
            Expr* ret = make_minus_infinity(); // Log[0] = -Infinity
            return ret;
        }
        // Inexact zero: direction is ambiguous in floating point, so the
        // limit is Indeterminate (matches Mathematica's Log[0.] behaviour).
        // Catches both +0.0 and -0.0 by IEEE 754 equality.
        if (z->type == EXPR_REAL && z->data.real == 0.0) {
            return expr_new_symbol(SYM_Indeterminate);
        }
        if (z->type == EXPR_INTEGER && z->data.integer == 1) {
            Expr* ret = expr_new_integer(0); // Log[1] = 0
            return ret;
        }

        // Negative integer: Log[n] = I*Pi + Log[-n] for n < 0
        if ((z->type == EXPR_INTEGER && z->data.integer < 0) ||
            (z->type == EXPR_BIGINT && mpz_sgn(z->data.bigint) < 0)) {
            Expr* neg_z;
            if (z->type == EXPR_INTEGER) {
                if (z->data.integer == INT64_MIN) {
                    mpz_t tmp;
                    mpz_init_set_si(tmp, z->data.integer);
                    mpz_neg(tmp, tmp);
                    neg_z = expr_new_bigint_from_mpz(tmp);
                    mpz_clear(tmp);
                } else {
                    neg_z = expr_new_integer(-z->data.integer);
                }
            } else {
                mpz_t tmp;
                mpz_init(tmp);
                mpz_neg(tmp, z->data.bigint);
                neg_z = expr_new_bigint_from_mpz(tmp);
                mpz_clear(tmp);
                neg_z = expr_bigint_normalize(neg_z);
            }
            Expr* log_args[1] = { neg_z };
            Expr* log_neg = expr_new_function(expr_new_symbol(SYM_Log), log_args, 1);
            Expr* times_args[2] = { expr_new_symbol(SYM_I), expr_new_symbol(SYM_Pi) };
            Expr* i_pi = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
            Expr* plus_args[2] = { i_pi, log_neg };
            return expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
        }

        // Pure imaginary exact argument: Log[Complex[0, b]] with b an exact
        // real numeric. On the principal branch arg(b I) = +Pi/2 for b > 0 and
        // -Pi/2 for b < 0, so
        //   Log[b I] = Log[|b|] +/- (I Pi)/2.
        // Restricting b to exact real numerics keeps us on the principal
        // branch (inexact imaginaries fall through to the clog path below) and
        // matches Mathematica: Log[I] = (I Pi)/2, Log[-3 I] = -((I Pi)/2)+Log[3].
        {
            Expr *cre, *cim;
            if (is_complex(z, &cre, &cim) &&
                cre->type == EXPR_INTEGER && cre->data.integer == 0 &&
                is_real_numeric_expr(cim)) {
                bool positive = is_positive_known(cim);
                // Magnitude Log[|b|]. For b < 0 wrap Times[-1, b]; the evaluator
                // folds it to a positive numeric before recursing into Log.
                Expr* mag_arg;
                if (positive) {
                    mag_arg = expr_copy(cim);
                } else {
                    Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(cim) };
                    mag_arg = expr_new_function(expr_new_symbol(SYM_Times), neg_args, 2);
                }
                Expr* log_arg[1] = { mag_arg };
                Expr* log_mag = expr_new_function(expr_new_symbol(SYM_Log), log_arg, 1);
                // +/- (I Pi)/2 = Times[Rational[+/-1, 2], I, Pi].
                Expr* half = make_rational(positive ? 1 : -1, 2);
                Expr* ipi_args[3] = { half, expr_new_symbol(SYM_I), expr_new_symbol(SYM_Pi) };
                Expr* ipi = expr_new_function(expr_new_symbol(SYM_Times), ipi_args, 3);
                Expr* plus_args[2] = { log_mag, ipi };
                return expr_new_function(expr_new_symbol(SYM_Plus), plus_args, 2);
            }
        }

        if (is_infinity(z)) {
            Expr* ret = expr_new_symbol(SYM_Infinity); // Log[Infinity] = Infinity
            return ret;
        }
        if (z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_E) {
            Expr* ret = expr_new_integer(1); // Log[E] = 1
            return ret;
        }

        // Log[E^k] -> k for real numeric k (int, bigint, rational, real).
        // E > 0 puts us on the principal branch; restricting k to real
        // numerics avoids the branch-cut mismatch that can occur for
        // complex k.
        if (is_power_call(z)) {
            Expr* pbase = z->data.function.args[0];
            Expr* pexp  = z->data.function.args[1];
            if (pbase->type == EXPR_SYMBOL &&
                pbase->data.symbol.name == SYM_E &&
                is_real_numeric_expr(pexp)) {
                return expr_copy(pexp);
            }
        }

#ifdef USE_MPFR
        if (numeric_expr_is_mpfr(z)) {
            /* Real MPFR path: positive real input goes through mpfr_log
             * directly for an EXPR_MPFR result. */
            long bits = numeric_combined_bits(z, NULL, 0);
            mpfr_t rr, ii;
            mpfr_init2(rr, bits); mpfr_init2(ii, bits);
            bool ok = get_approx_mpfr(z, rr, ii, NULL);
            if (ok && mpfr_zero_p(ii) && mpfr_sgn(rr) > 0) {
                mpfr_t out; mpfr_init2(out, bits);
                mpfr_log(out, rr, MPFR_RNDN);
                mpfr_clear(rr); mpfr_clear(ii);
                return expr_new_mpfr_move(out);
            }
            mpfr_clear(rr); mpfr_clear(ii);
            /* Complex MPFR path: negative real MPFR, or Complex with an
             * MPFR component. Goes through log(|z|) + I arg(z) at the
             * working precision via the helper. The helper handles the
             * zero-input case by producing -Infinity-style MPFR output
             * — but the symbolic Log[0] = ComplexInfinity path above
             * has already filtered exact-zero inputs. */
            Expr* r = numeric_mpfr_apply_complex_unary(z, 0, mpfr_complex_log);
            if (r) return r;
        }
#endif
        // Approximate numerical evaluation
        double complex c;
        bool inexact = false;
    if (get_approx(z, &c, &inexact) && inexact) {
            double complex s = clog(c);
            // Return real result if output is purely real, otherwise return complex
            Expr* ret = NULL;
            if (cimag(c) == 0.0 && creal(c) > 0.0) ret = expr_new_real(creal(s));
            else ret = make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
            return ret;
        }
    } 
    // Log[b, z] - Logarithm to base b
    else if (argc == 2) {
        Expr* b = res->data.function.args[0];
        Expr* z = res->data.function.args[1];

        // Inexact zero argument: Log[b, 0.] = Indeterminate, matching the
        // 1-arg float-zero case. Direction is ambiguous in floating point.
        if (z->type == EXPR_REAL && z->data.real == 0.0) {
            return expr_new_symbol(SYM_Indeterminate);
        }
        // Exact integer zero argument with a real-positive base b != 1:
        //   Log[b, 0] = Log[0] / Log[b] = -Infinity / Log[b].
        // For b > 1 (Log[b] > 0) the directed limit is -Infinity;
        // for 0 < b < 1 (Log[b] < 0) it is +Infinity. Other bases
        // (b <= 0, b == 1, symbolic) fall through to the default rewrite.
        if (z->type == EXPR_INTEGER && z->data.integer == 0) {
            int sgn = 0;   /* -1: 0<b<1, 0: unknown/skip, +1: b>1 */
            if (b->type == EXPR_INTEGER) {
                if (b->data.integer > 1) sgn = +1;
            } else if (b->type == EXPR_BIGINT) {
                /* Bigints are always |n| > INT64_MAX, so positive means > 1. */
                if (mpz_sgn(b->data.bigint) > 0) sgn = +1;
            } else if (b->type == EXPR_REAL) {
                if (b->data.real > 1.0) sgn = +1;
                else if (b->data.real > 0.0 && b->data.real < 1.0) sgn = -1;
            } else if (b->type == EXPR_SYMBOL) {
                /* E (~2.718) and Pi (~3.14) are both > 1. */
                if (b->data.symbol.name == SYM_E || b->data.symbol.name == SYM_Pi) sgn = +1;
            } else {
                int64_t n, d;
                if (is_rational(b, &n, &d)) {
                    /* Denominator d is always positive in canonical form. */
                    if (n > d) sgn = +1;
                    else if (n > 0 && n < d) sgn = -1;
                }
            }
            if (sgn == +1) return make_minus_infinity();
            if (sgn == -1) return expr_new_symbol(SYM_Infinity);
            /* else fall through to the default rewrite */
        }

        // Log[b, b] = 1
        if (expr_eq(b, z)) {
            Expr* ret = expr_new_integer(1);
            return ret;
        }

        // Log[b, b^k] -> k when b is known positive (positive numeric or
        // a symbol like E, Pi) and k is a real numeric. Same branch-cut
        // reasoning as Log[E^k]: the identity holds for all real k but
        // can fail on the principal branch for complex k.
        if (is_power_call(z)) {
            Expr* pbase = z->data.function.args[0];
            Expr* pexp  = z->data.function.args[1];
            if (expr_eq(pbase, b) &&
                is_positive_known(b) &&
                is_real_numeric_expr(pexp)) {
                return expr_copy(pexp);
            }
        }

        // Attempt to return exact rational results for integer bases and arguments (e.g. Log[2, 8] = 3)
        if (b->type == EXPR_INTEGER && z->type == EXPR_INTEGER) {
            int64_t bv = b->data.integer;
            int64_t zv = z->data.integer;
            if (bv > 1 && zv > 0) {
                int64_t temp = zv;
                int64_t p = 0;
                while (temp > 1 && temp % bv == 0) {
                    temp /= bv;
                    p++;
                }
                if (temp == 1) {
                    Expr* ret = expr_new_integer(p);
                    return ret;
                }
            }
        }

        // Default rewrite: Log[b, z] -> Log[z] / Log[b]
        Expr* num_args[1] = { expr_copy(z) };
        Expr* den_args[1] = { expr_copy(b) };
        Expr* num = expr_new_function(expr_new_symbol(SYM_Log), num_args, 1);
        Expr* den = expr_new_function(expr_new_symbol(SYM_Log), den_args, 1);

        Expr* pow_args[2] = { den, expr_new_integer(-1) };
        Expr* inv_den = expr_new_function(expr_new_symbol(SYM_Power), pow_args, 2);

        Expr* times_args[2] = { num, inv_den };
        Expr* ret = expr_new_function(expr_new_symbol(SYM_Times), times_args, 2);
        return ret;
    }

    // Remains unevaluated if it doesn't match above rules
    return NULL;
}

/*
 * builtin_exp:
 * Implements the evaluation logic for the 'Exp' function.
 */
Expr* builtin_exp(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Exp", res->data.function.arg_count, 1, 1);
    Expr* z = res->data.function.args[0];

    // Exact evaluations for special constants
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        Expr* ret = expr_new_integer(1); // Exp[0] = 1
        return ret;
    }
    if (is_minus_infinity(z)) {
        Expr* ret = expr_new_integer(0); // Exp[-Infinity] = 0
        return ret;
    }
    if (is_infinity(z)) {
        Expr* ret = expr_new_symbol(SYM_Infinity); // Exp[Infinity] = Infinity
        return ret;
    }

    // Exact evaluation of Exp[I * q * Pi] using Euler's formula (e^{i x} = Cos[x] + I Sin[x])
    // The argument z may be internally structured as Times[Complex[0, q], Pi]
    if (z->type == EXPR_FUNCTION && z->data.function.head->data.symbol.name == SYM_Times) {
        bool has_pi = false;
        Expr* im_coeff = NULL;

        // Scan the arguments of Times to identify a pure imaginary coefficient and Pi
        for (size_t i = 0; i < z->data.function.arg_count; i++) {
            Expr* arg = z->data.function.args[i];
            if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_Pi) {
                has_pi = true;
            } else {
                Expr *re, *im;
                if (is_complex(arg, &re, &im)) {
                    // We only rewrite if the complex number is purely imaginary (real part == 0)
                    if (re->type == EXPR_INTEGER && re->data.integer == 0) {
                        im_coeff = im;
                    }
                }
            }
        }

        // If we found Pi and exactly one pure imaginary coefficient (total args = 2)
        if (has_pi && im_coeff && z->data.function.arg_count == 2) {
            int64_t n, d;
            // Only perform expansion if the coefficient is rational or integer (e.g. q * Pi)
            if (im_coeff->type == EXPR_INTEGER || is_rational(im_coeff, &n, &d)) {
                // Canonicalize Exp[I q Pi] -> (-1)^q (matching Mathematica).
                // Power[-1, q] already reduces the special cases itself:
                // (-1)^1 = -1, (-1)^(1/2) = I, (-1)^2 = 1, etc., while
                // leaving non-reducible roots like (-1)^(1/5) intact rather
                // than over-eagerly expanding into trig radicals.
                Expr* pow_args[2] = { expr_new_integer(-1), expr_copy(im_coeff) };
                Expr* ret = expr_new_function(expr_new_symbol(SYM_Power), pow_args, 2);
                return ret;
            }
        }
    }

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(z)) {
        /* Real MPFR path first: pure real input takes mpfr_exp directly
         * and returns EXPR_MPFR. */
        Expr* r = numeric_mpfr_apply_unary(z, 0, mpfr_exp);
        if (r) return r;
        /* Otherwise complex MPFR: Complex[MPFR, MPFR] (or a real with
         * inexact imag component) goes through the complex helper,
         * which uses exp(a)(cos(b) + i sin(b)) at the working
         * precision and collapses on zero-imag-out. */
        r = numeric_mpfr_apply_complex_unary(z, 0, mpfr_complex_exp);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
    if (get_approx(z, &c, &inexact) && inexact) {
        double complex s = cexp(c);
        // Return real result if output is purely real, otherwise return complex
        Expr* ret = NULL;
        if (cimag(c) == 0.0) ret = expr_new_real(creal(s));
        else ret = make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
        return ret;
    }

    // Remains unevaluated if it doesn't match above rules
    return expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){expr_new_symbol(SYM_E), expr_copy(z)}, 2);
}
