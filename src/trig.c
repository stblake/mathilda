#include <complex.h>
#include "trig.h"
#include "arithmetic.h"
#include "plus.h"
#include "times.h"
#include "power.h"
#include "complex.h"
#include "symtab.h"
#include "eval.h"
#include "common.h"
#include "numeric.h"
#include "numeric_complex.h"
#include "sym_intern.h"
#include "sym_names.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#include <string.h>
#include <stdint.h>

#ifdef USE_MPFR
/* MPFR helpers for the reciprocal-input inverses ArcCot, ArcSec, ArcCsc.
 * Each computes f(1/x) at the input's precision. Exposed as
 * MpfrUnaryOp-compatible function pointers so we can reuse
 * numeric_mpfr_apply_unary below. */
static int mpfr_atan_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_atan(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
static int mpfr_acos_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_acos(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
static int mpfr_asin_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_asin(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
#endif

/* Helper function to construct Sqrt[n] as Power[n, 1/2] */
static Expr* make_sqrt(int64_t n) {
    Expr* base = expr_new_integer(n);
    Expr* exponent = make_rational(1, 2);
    Expr* res = make_power(base, exponent);
    return res;
}

/* Helper function to construct Sqrt[expr] as Power[expr, 1/2] */
static Expr* make_sqrt_expr(Expr* e) {
    Expr* exponent = make_rational(1, 2);
    Expr* res = make_power(e, exponent);
    return res;
}

/* 
 * extract_pi_multiplier:
 * Checks if the expression is of the form n/d * Pi or simply Pi (1/1 * Pi).
 * Returns true if match found, and sets n and d.
 */
static bool extract_pi_multiplier(Expr* e, int64_t* n, int64_t* d) {
    // Case 1: Pi
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Pi) {
        *n = 1; *d = 1;
        return true;
    }
    
    // Case 2: n/d * Pi (Times[Rational[n, d], Pi])
    if (e->type == EXPR_FUNCTION && e->data.function.head->data.symbol.name == SYM_Times && e->data.function.arg_count == 2) {
        Expr* first = e->data.function.args[0];
        Expr* second = e->data.function.args[1];
        
        if (second->type == EXPR_SYMBOL && second->data.symbol.name == SYM_Pi) {
            if (is_rational(first, n, d)) return true;
        }
    }
    return false;
}

/*
 * flip_sign:
 * Returns a newly allocated, evaluated expression equal to -e. Caller
 * owns the result; `e` is untouched. Uses the evaluator so that e.g.
 * Times[-1, Complex[0, -2]] canonicalises back to Complex[0, 2] and
 * Times[-1, Times[-k, x]] collapses to Times[k, x].
 */
static Expr* flip_sign(Expr* e) {
    Expr* args[2] = { expr_new_integer(-1), expr_copy(e) };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), args, 2));
}

/*
 * even_fold / odd_fold:
 * If `arg` is superficially negative (see expr_is_superficially_negative
 * in arithmetic.h), rewrite f[arg] as the simplified form -- f[|arg|]
 * for even f, -f[|arg|] for odd f -- and return the new expression.
 * Return NULL otherwise so the caller can continue with other
 * simplifications. The caller's `res` is freed by evaluate_step once
 * this helper returns a non-NULL replacement, per the builtin calling
 * convention -- do NOT free it here.
 */
static Expr* even_fold(Expr* arg, const char* head_name) {
    if (!expr_is_superficially_negative(arg)) return NULL;
    Expr* neg = flip_sign(arg);
    Expr* inner_args[1] = { neg };
    return expr_new_function(expr_new_symbol(head_name), inner_args, 1);
}

static Expr* odd_fold(Expr* arg, const char* head_name) {
    if (!expr_is_superficially_negative(arg)) return NULL;
    Expr* neg = flip_sign(arg);
    Expr* inner_args[1] = { neg };
    Expr* inner = expr_new_function(expr_new_symbol(head_name), inner_args, 1);
    return make_times(expr_new_integer(-1), inner);
}

/*
 * peel_imaginary_unit:
 * If `arg` is a pure-imaginary-unit-scaled expression -- Complex[0, k]
 * (k > 0) or Times[Complex[0, k], rest...] (k > 0) -- return a newly
 * allocated, evaluated expression equal to `arg / I` (the "real"
 * coefficient). Otherwise NULL. Assumes sign-extraction has already
 * canonicalised the leading imaginary coefficient to positive.
 */
static Expr* peel_imaginary_unit(Expr* arg) {
    Expr* re; Expr* im;
    if (is_complex(arg, &re, &im)) {
        if (expr_numeric_sign(re) != 0) return NULL;
        if (expr_numeric_sign(im) <= 0) return NULL;
        return expr_copy(im);
    }
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_Times &&
        arg->data.function.arg_count > 0) {
        Expr* first = arg->data.function.args[0];
        Expr* fre; Expr* fim;
        if (!is_complex(first, &fre, &fim)) return NULL;
        if (expr_numeric_sign(fre) != 0) return NULL;
        if (expr_numeric_sign(fim) <= 0) return NULL;
        size_t n = arg->data.function.arg_count;
        if (n == 1) return expr_copy(fim);
        Expr** new_args = (Expr**)malloc(sizeof(Expr*) * n);
        new_args[0] = expr_copy(fim);
        for (size_t i = 1; i < n; i++) new_args[i] = expr_copy(arg->data.function.args[i]);
        Expr* t = expr_new_function(expr_new_symbol(SYM_Times), new_args, n);
        free(new_args);
        return eval_and_free(t);
    }
    return NULL;
}

/*
 * trig_i_fold:
 * Rewrites f[I*y] using the Cos/Cosh -- Sin/Sinh -- Tan/Tanh bridge:
 *   Cos[I y]  -> Cosh[y]          Sec[I y]  -> Sech[y]
 *   Sin[I y]  -> I Sinh[y]        Csc[I y]  -> -I Csch[y]
 *   Tan[I y]  -> I Tanh[y]        Cot[I y]  -> -I Coth[y]
 * `counterpart` is the hyperbolic head to emit. `coeff_kind` selects
 * the multiplier prepended to the result: 0 = none, +1 = I, -1 = -I.
 * Returns NULL when `arg` is not of the form I*y. The caller's `res`
 * is freed by evaluate_step after a non-NULL return.
 */
static Expr* trig_i_fold(Expr* arg, const char* counterpart, int coeff_kind) {
    Expr* y = peel_imaginary_unit(arg);
    if (!y) return NULL;
    Expr* a[1] = { y };
    Expr* inner = expr_new_function(expr_new_symbol(counterpart), a, 1);
    if (coeff_kind == 0) return inner;
    /* +I coefficient: Complex[0, 1]; -I coefficient: Complex[0, -1]. */
    Expr* coeff = make_complex(expr_new_integer(0), expr_new_integer(coeff_kind));
    return make_times(coeff, inner);
}

/*
 * arc_pi_minus_fold:
 * Negation reflection for inverse functions whose negation introduces a Pi
 * shift rather than a sign flip:  ArcCos[-x] -> Pi - ArcCos[x],
 * ArcSec[-x] -> Pi - ArcSec[x]. Returns NULL when arg is not superficially
 * negative; otherwise returns Plus[Pi, Times[-1, head_name[|arg|]]] for
 * the canonicalising evaluator to flatten.
 */
static Expr* arc_pi_minus_fold(Expr* arg, const char* head_name) {
    if (!expr_is_superficially_negative(arg)) return NULL;
    Expr* neg = flip_sign(arg);
    Expr* inner_args[1] = { neg };
    Expr* inner = expr_new_function(expr_new_symbol(head_name), inner_args, 1);
    Expr* neg_inner = make_times(expr_new_integer(-1), inner);
    Expr* sum_args[2] = { expr_new_symbol(SYM_Pi), neg_inner };
    return expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
}

/*
 * arccos_i_fold:
 * Imaginary-axis bridge for ArcCos: ArcCos[I y] -> Pi/2 - I*ArcSinh[y].
 * Mirrors trig_i_fold's I-extraction but emits a constant + scaled hyperbolic
 * inverse instead of a single scaled call. Returns NULL when arg is not of
 * the form I*y with positive imaginary coefficient.
 */
static Expr* arccos_i_fold(Expr* arg) {
    Expr* y = peel_imaginary_unit(arg);
    if (!y) return NULL;
    Expr* a[1] = { y };
    Expr* asinh_call = expr_new_function(expr_new_symbol(SYM_ArcSinh), a, 1);
    Expr* neg_i = make_complex(expr_new_integer(0), expr_new_integer(-1));
    Expr* neg_i_asinh = make_times(neg_i, asinh_call);
    Expr* half = make_rational(1, 2);
    Expr* pi = expr_new_symbol(SYM_Pi);
    Expr* half_pi = make_times(half, pi);
    Expr* sum_args[2] = { half_pi, neg_i_asinh };
    return expr_new_function(expr_new_symbol(SYM_Plus), sum_args, 2);
}

/*
 * try_simp_forward_of_inverse:
 * Universal forward-of-inverse identities for trig functions. `outer` is the
 * forward head ("Sin","Cos","Tan",...), `arg` is the body of that call. If
 * `arg` is an inverse-trig call whose pairing with `outer` yields a closed
 * algebraic form, return the rewritten Expr*; otherwise NULL. Identities
 * follow the principal-branch conventions Mathematica exposes:
 *   Sin[ArcCos[x]] -> Sqrt[1-x^2]      Cos[ArcSin[x]] -> Sqrt[1-x^2]
 *   Sin[ArcTan[x]] -> x/Sqrt[1+x^2]    Cos[ArcTan[x]] -> 1/Sqrt[1+x^2]
 *   Tan[ArcSin[x]] -> x/Sqrt[1-x^2]    Tan[ArcCos[x]] -> Sqrt[1-x^2]/x
 *   Tan[ArcCot[x]] -> 1/x              Cot[ArcTan[x]] -> 1/x
 * The inverse-of-forward direction (e.g. ArcSin[Sin[x]]) is intentionally
 * NOT folded -- those reduce to x only on the principal domain of each f.
 */
static Expr* try_simp_forward_of_inverse(const char* outer, Expr* arg) {
    if (arg->type != EXPR_FUNCTION || arg->data.function.arg_count != 1) return NULL;
    if (!arg->data.function.head ||
        arg->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* inner = arg->data.function.head->data.symbol.name;
    Expr* x = arg->data.function.args[0];

    /* Helpers to build common subexpressions: Sqrt[1 - x^2], Sqrt[1 + x^2],
     * x / Sqrt[1 - x^2], etc. We construct the unevaluated tree and let the
     * outer evaluator canonicalise (it folds Plus, Power, Times). */
    #define SQRT_OF(arg_e)                                                     \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),              \
            (Expr*[]){ (arg_e), make_rational(1, 2) }, 2))
    #define X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),              \
            (Expr*[]){ expr_copy(x), expr_new_integer(2) }, 2))
    #define ONE_MINUS_X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),               \
            (Expr*[]){ expr_new_integer(1),                                    \
                       eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), \
                           (Expr*[]){ expr_new_integer(-1), X_SQ() }, 2)) }, 2))
    #define ONE_PLUS_X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),               \
            (Expr*[]){ expr_new_integer(1), X_SQ() }, 2))

    if (strcmp(outer, "Sin") == 0 && strcmp(inner, "ArcCos") == 0) {
        return SQRT_OF(ONE_MINUS_X_SQ());
    }
    if (strcmp(outer, "Cos") == 0 && strcmp(inner, "ArcSin") == 0) {
        return SQRT_OF(ONE_MINUS_X_SQ());
    }
    if (strcmp(outer, "Sin") == 0 && strcmp(inner, "ArcTan") == 0) {
        Expr* den = SQRT_OF(ONE_PLUS_X_SQ());
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
        return make_times(expr_copy(x), inv);
    }
    if (strcmp(outer, "Cos") == 0 && strcmp(inner, "ArcTan") == 0) {
        Expr* den = SQRT_OF(ONE_PLUS_X_SQ());
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
    }
    if (strcmp(outer, "Tan") == 0 && strcmp(inner, "ArcSin") == 0) {
        Expr* den = SQRT_OF(ONE_MINUS_X_SQ());
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
        return make_times(expr_copy(x), inv);
    }
    if (strcmp(outer, "Tan") == 0 && strcmp(inner, "ArcCos") == 0) {
        Expr* num = SQRT_OF(ONE_MINUS_X_SQ());
        Expr* inv_x = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(x), expr_new_integer(-1) }, 2));
        return make_times(num, inv_x);
    }
    if ((strcmp(outer, "Tan") == 0 && strcmp(inner, "ArcCot") == 0) ||
        (strcmp(outer, "Cot") == 0 && strcmp(inner, "ArcTan") == 0)) {
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(x), expr_new_integer(-1) }, 2));
    }

    #undef SQRT_OF
    #undef X_SQ
    #undef ONE_MINUS_X_SQ
    #undef ONE_PLUS_X_SQ
    return NULL;
}

/*
 * get_approx:
 * Tries to get a numeric complex approximation of the expression.
 */
static bool get_approx(Expr* e, double complex* out, bool* is_inexact) {
    if (e->type == EXPR_INTEGER) {
        *out = (double)e->data.integer + 0.0 * I;
        if (is_inexact) *is_inexact = false;
        return true;
    }
    if (e->type == EXPR_REAL) {
        *out = e->data.real + 0.0 * I;
        if (is_inexact) *is_inexact = true;
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        *out = (double)n / d + 0.0 * I;
        if (is_inexact) *is_inexact = false;
        return true;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        double complex r, i;
        bool rex = false, imx = false;
        if (get_approx(re, &r, &rex) && get_approx(im, &i, &imx)) {
            *out = creal(r) + creal(i) * I;
            if (is_inexact) *is_inexact = (rex || imx);
            return true;
        }
    }
    return false;
}

// --- Exact Symbolic Trigonometric Values ---

/*
 * exact_sin:
 * Computes exact symbolic evaluations for Sin[n/d * Pi].
 * We reduce the argument using 2Pi periodicity and symmetry properties.
 * Currently supports denominators: 1, 2, 3, 4, 5, 6, 10, 12.
 */
static Expr* exact_sin(int64_t n, int64_t d) {
    // Reduce to [0, 2Pi)
    n = n % (2 * d);
    if (n < 0) n += 2 * d;
    
    int64_t sign = 1;
    // Map to [0, Pi]
    if (n >= d) {
        sign = -1;
        n -= d;
    }
    // Map to [0, Pi/2]
    if (n > d / 2.0) {
        n = d - n;
    }
    
    // Reduce fraction
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_integer(0);
    else if (d == 2) res = expr_new_integer(1);
    else if (d == 3) res = make_times(make_rational(1, 2), make_sqrt(3));
    else if (d == 4) res = make_power(expr_new_integer(2), make_rational(-1, 2));
    else if (d == 5) {
        if (n == 1) res = make_times(make_rational(1, 4), make_sqrt_expr(make_plus(expr_new_integer(10), make_times(expr_new_integer(-2), make_sqrt(5)))));
        else if (n == 2) res = make_times(make_rational(1, 4), make_sqrt_expr(make_plus(expr_new_integer(10), make_times(expr_new_integer(2), make_sqrt(5)))));
    }
    else if (d == 6) res = make_rational(1, 2);
    else if (d == 10) {
        if (n == 1) res = make_times(make_rational(1, 4), make_plus(make_sqrt(5), expr_new_integer(-1)));
        else if (n == 3) res = make_times(make_rational(1, 4), make_plus(make_sqrt(5), expr_new_integer(1)));
    }
    else if (d == 12) {
        if (n == 1) res = make_times(make_rational(1, 4), make_plus(make_sqrt(6), make_times(expr_new_integer(-1), make_sqrt(2))));
        else if (n == 5) res = make_times(make_rational(1, 4), make_plus(make_sqrt(6), make_sqrt(2)));
    }
    
    // Apply resolved sign
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}

/*
 * exact_cos:
 * Computes exact symbolic evaluations for Cos[n/d * Pi].
 * Uses cosine's even symmetry and 2Pi periodicity.
 */
static Expr* exact_cos(int64_t n, int64_t d) {
    // Reduce to [0, 2Pi)
    n = n % (2 * d);
    if (n < 0) n += 2 * d;
    // Exploit symmetry Cos[x] = Cos[2Pi - x] to map to [0, Pi]
    if (n > d) n = 2 * d - n;
    
    int64_t sign = 1;
    // Map to [0, Pi/2]
    if (n > d / 2.0) {
        n = d - n;
        sign = -1;
    }
    
    // Reduce fraction
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_integer(1);
    else if (d == 2) res = expr_new_integer(0);
    else if (d == 3) res = make_rational(1, 2);
    else if (d == 4) res = make_power(expr_new_integer(2), make_rational(-1, 2));
    else if (d == 5) {
        if (n == 1) res = make_times(make_rational(1, 4), make_plus(make_sqrt(5), expr_new_integer(1)));
        else if (n == 2) res = make_times(make_rational(1, 4), make_plus(make_sqrt(5), expr_new_integer(-1)));
    }
    else if (d == 6) res = make_times(make_rational(1, 2), make_sqrt(3));
    else if (d == 10) {
        if (n == 1) res = make_times(make_rational(1, 4), make_sqrt_expr(make_plus(expr_new_integer(10), make_times(expr_new_integer(2), make_sqrt(5)))));
        else if (n == 3) res = make_times(make_rational(1, 4), make_sqrt_expr(make_plus(expr_new_integer(10), make_times(expr_new_integer(-2), make_sqrt(5)))));
    }
    else if (d == 12) {
        if (n == 1) res = make_times(make_rational(1, 4), make_plus(make_sqrt(6), make_sqrt(2)));
        else if (n == 5) res = make_times(make_rational(1, 4), make_plus(make_sqrt(6), make_times(expr_new_integer(-1), make_sqrt(2))));
    }
    
    // Apply resolved sign
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}

/*
 * exact_tan:
 * Computes exact symbolic evaluations for Tan[n/d * Pi].
 * Tangent is Pi-periodic. We reduce to [0, Pi/2] and resolve exact symbolic forms.
 */
static Expr* exact_tan(int64_t n, int64_t d) {
    n = n % d;
    if (n < 0) n += d;
    
    int64_t sign = 1;
    // Map to [0, Pi/2]
    if (n > d / 2.0) {
        n = d - n;
        sign = -1;
    }
    
    // Reduce fraction
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_integer(0);
    else if (d == 2) res = expr_new_symbol(SYM_ComplexInfinity);
    else if (d == 3) res = make_sqrt(3);
    else if (d == 4) res = expr_new_integer(1);
    else if (d == 5) {
        if (n == 1) res = make_sqrt_expr(make_plus(expr_new_integer(5), make_times(expr_new_integer(-2), make_sqrt(5))));
        else if (n == 2) res = make_sqrt_expr(make_plus(expr_new_integer(5), make_times(expr_new_integer(2), make_sqrt(5))));
    }
    else if (d == 6) res = make_power(expr_new_integer(3), make_rational(-1, 2));
    else if (d == 10) {
        if (n == 1) res = make_sqrt_expr(make_plus(expr_new_integer(1), make_times(make_rational(-2, 5), make_sqrt(5))));
        else if (n == 3) res = make_sqrt_expr(make_plus(expr_new_integer(1), make_times(make_rational(2, 5), make_sqrt(5))));
    }
    else if (d == 12) {
        if (n == 1) res = make_plus(expr_new_integer(2), make_times(expr_new_integer(-1), make_sqrt(3)));
        else if (n == 5) res = make_plus(expr_new_integer(2), make_sqrt(3));
    }
    
    // Apply resolved sign
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}

/*
 * exact_cot:
 * Computes exact symbolic evaluations for Cot[n/d * Pi].
 * Cotangent is Pi-periodic. We reduce to [0, Pi/2] and resolve exact symbolic forms.
 */
static Expr* exact_cot(int64_t n, int64_t d) {
    n = n % d;
    if (n < 0) n += d;
    
    int64_t sign = 1;
    // Map to [0, Pi/2]
    if (n > d / 2.0) {
        n = d - n;
        sign = -1;
    }
    
    // Reduce fraction
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_symbol(SYM_ComplexInfinity);
    else if (d == 2) res = expr_new_integer(0);
    else if (d == 3) res = make_power(expr_new_integer(3), make_rational(-1, 2));
    else if (d == 4) res = expr_new_integer(1);
    else if (d == 5) {
        if (n == 1) res = make_sqrt_expr(make_plus(expr_new_integer(1), make_times(make_rational(2, 1), make_power(expr_new_integer(5), make_rational(-1, 2)))));
        else if (n == 2) res = make_sqrt_expr(make_plus(expr_new_integer(1), make_times(make_rational(-2, 1), make_power(expr_new_integer(5), make_rational(-1, 2)))));
    }
    else if (d == 6) res = make_sqrt(3);
    else if (d == 10) {
        if (n == 1) res = make_sqrt_expr(make_plus(expr_new_integer(5), make_times(expr_new_integer(2), make_sqrt(5))));
        else if (n == 3) res = make_sqrt_expr(make_plus(expr_new_integer(5), make_times(expr_new_integer(-2), make_sqrt(5))));
    }
    else if (d == 12) {
        if (n == 1) res = make_plus(expr_new_integer(2), make_sqrt(3));
        else if (n == 5) res = make_plus(expr_new_integer(2), make_times(expr_new_integer(-1), make_sqrt(3)));
    }
    
    // Apply resolved sign
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}

/*
 * exact_sec:
 * Computes exact symbolic evaluations for Sec[n/d * Pi].
 * Uses cosine's 2Pi periodicity and even symmetry.
 */
static Expr* exact_sec(int64_t n, int64_t d) {
    n = n % (2 * d);
    if (n < 0) n += 2 * d;
    if (n > d) n = 2 * d - n;
    
    int64_t sign = 1;
    // Map to [0, Pi/2]
    if (n > d / 2.0) {
        n = d - n;
        sign = -1;
    }
    
    // Reduce fraction
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_integer(1);
    else if (d == 2) res = expr_new_symbol(SYM_ComplexInfinity);
    else if (d == 3) res = expr_new_integer(2);
    else if (d == 4) res = make_sqrt(2);
    else if (d == 5) {
        if (n == 1) res = make_plus(make_sqrt(5), expr_new_integer(-1));
        else if (n == 2) res = make_plus(make_sqrt(5), expr_new_integer(1));
    }
    else if (d == 6) res = make_times(expr_new_integer(2), make_power(expr_new_integer(3), make_rational(-1, 2)));
    else if (d == 10) {
        if (n == 1) res = make_power(make_sqrt_expr(make_plus(make_rational(5, 8), make_times(make_rational(1, 8), make_sqrt(5)))), expr_new_integer(-1));
        else if (n == 3) res = make_power(make_sqrt_expr(make_plus(make_rational(5, 8), make_times(make_rational(-1, 8), make_sqrt(5)))), expr_new_integer(-1));
    }
    else if (d == 12) {
        if (n == 1) res = make_times(make_sqrt(2), make_plus(make_sqrt(3), expr_new_integer(-1)));
        else if (n == 5) res = make_times(make_sqrt(2), make_plus(make_sqrt(3), expr_new_integer(1)));
    }
    
    // Apply resolved sign
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}

/*
 * exact_csc:
 * Computes exact symbolic evaluations for Csc[n/d * Pi].
 * Uses sine's 2Pi periodicity and odd symmetry.
 */
static Expr* exact_csc(int64_t n, int64_t d) {
    n = n % (2 * d);
    if (n < 0) n += 2 * d;
    
    int64_t sign = 1;
    if (n >= d) {
        sign = -1;
        n -= d;
    }
    if (n > d / 2.0) {
        n = d - n;
    }
    
    int64_t g = gcd(n, d);
    n /= g; d /= g;
    
    Expr* res = NULL;
    if (d == 1) res = expr_new_symbol(SYM_ComplexInfinity);
    else if (d == 2) res = expr_new_integer(1);
    else if (d == 3) res = make_times(expr_new_integer(2), make_power(expr_new_integer(3), make_rational(-1, 2)));
    else if (d == 4) res = make_sqrt(2);
    else if (d == 5) {
        if (n == 1) res = make_power(make_sqrt_expr(make_plus(make_rational(5, 8), make_times(make_rational(1, 8), make_sqrt(5)))), expr_new_integer(-1));
        else if (n == 2) res = make_power(make_sqrt_expr(make_plus(make_rational(5, 8), make_times(make_rational(1, 8), make_sqrt(5)))), expr_new_integer(-1));
    }
    else if (d == 6) res = expr_new_integer(2);
    else if (d == 10) {
        if (n == 1) res = make_plus(make_sqrt(5), expr_new_integer(1));
        else if (n == 3) res = make_plus(make_sqrt(5), expr_new_integer(-1));
    }
    else if (d == 12) {
        if (n == 1) res = make_times(make_sqrt(2), make_plus(make_sqrt(3), expr_new_integer(1)));
        else if (n == 5) res = make_times(make_sqrt(2), make_plus(make_sqrt(3), expr_new_integer(-1)));
    }
    
    if (res && sign == -1) {
        return make_times(expr_new_integer(-1), res);
    }
    return res;
}


// --- Built-in Functions ---

/* If arg is a one-argument call whose head is `inverse_name`, return a deep
 * copy of its single argument; otherwise NULL. Used to fold the direct
 * forward/inverse identities Sin[ArcSin[x]] -> x, Cos[ArcCos[x]] -> x, ...
 * These hold identically over the complex numbers because each ArcX is a
 * right inverse of X by construction. The two-argument form ArcTan[x, y]
 * is excluded by the arg_count guard (Tan[ArcTan[x, y]] = y/x, not a single
 * variable). We deliberately do NOT fold the opposite direction
 * (ArcSin[Sin[x]], etc.) because those only reduce to x on each function's
 * principal domain. */
static Expr* strip_inverse_call(Expr* arg, const char* inverse_name) {
    if (head_is(arg, intern_symbol(inverse_name)) &&
        arg->data.function.arg_count == 1) {
        return expr_copy(arg->data.function.args[0]);
    }
    return NULL;
}

/*
 * builtin_sin:
 * Implements the standard evaluation logic for Sin.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (csin) if possible.
 */
Expr* builtin_sin(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Sin", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcSin"); if (inv) return inv; }

    // Sin[ArcCos[x]] -> Sqrt[1-x^2], Sin[ArcTan[x]] -> x/Sqrt[1+x^2]
    { Expr* f = try_simp_forward_of_inverse("Sin", arg); if (f) return f; }

    // Sin is odd: Sin[-x] -> -Sin[x]
    { Expr* f = odd_fold(arg, "Sin"); if (f) return f; }

    // Sin[I y] -> I Sinh[y]
    { Expr* f = trig_i_fold(arg, "Sinh", +1); if (f) return f; }

    // Sin[0] = 0
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_sin(n, d);
        if (exact) return exact;
    }
    
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_sin);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_sin);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = csin(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_cos:
 * Implements the standard evaluation logic for Cos.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (ccos) if possible.
 */
Expr* builtin_cos(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Cos", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCos"); if (inv) return inv; }

    // Cos[ArcSin[x]] -> Sqrt[1-x^2], Cos[ArcTan[x]] -> 1/Sqrt[1+x^2]
    { Expr* f = try_simp_forward_of_inverse("Cos", arg); if (f) return f; }

    // Cos is even: Cos[-x] -> Cos[x]
    { Expr* f = even_fold(arg, "Cos"); if (f) return f; }

    // Cos[I y] -> Cosh[y]
    { Expr* f = trig_i_fold(arg, "Cosh", 0); if (f) return f; }

    // Cos[0] = 1
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(1);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_cos(n, d);
        if (exact) return exact;
    }
    
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_cos);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_cos);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = ccos(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_tan:
 * Implements the standard evaluation logic for Tan.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (ctan) if possible.
 */
Expr* builtin_tan(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Tan", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcTan"); if (inv) return inv; }

    // Tan[ArcSin[x]] -> x/Sqrt[1-x^2], Tan[ArcCos[x]] -> Sqrt[1-x^2]/x,
    // Tan[ArcCot[x]] -> 1/x.
    { Expr* f = try_simp_forward_of_inverse("Tan", arg); if (f) return f; }

    // Tan is odd: Tan[-x] -> -Tan[x]
    { Expr* f = odd_fold(arg, "Tan"); if (f) return f; }

    // Tan[I y] -> I Tanh[y]
    { Expr* f = trig_i_fold(arg, "Tanh", +1); if (f) return f; }

    // Tan[0] = 0
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_tan(n, d);
        if (exact) return exact;
    }
    
    // Approximate numerical evaluation
    double complex cplx;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_tan);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_tan);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &cplx, &inexact) && inexact) {
        double complex s = ctan(cplx);
        if (cimag(cplx) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_cot:
 * Implements the standard evaluation logic for Cot.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (1/ctan) if possible.
 */
Expr* builtin_cot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Cot", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCot"); if (inv) return inv; }

    // Cot[ArcTan[x]] -> 1/x
    { Expr* f = try_simp_forward_of_inverse("Cot", arg); if (f) return f; }

    // Cot is odd: Cot[-x] -> -Cot[x]
    { Expr* f = odd_fold(arg, "Cot"); if (f) return f; }

    // Cot[I y] -> -I Coth[y]
    { Expr* f = trig_i_fold(arg, "Coth", -1); if (f) return f; }

    // Cot[0] = ComplexInfinity
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_symbol(SYM_ComplexInfinity);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_cot(n, d);
        if (exact) return exact;
    }
    
    // Approximate numerical evaluation
    double complex cplx;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_cot);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_cot);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &cplx, &inexact) && inexact) {
        double complex s = 1.0 / ctan(cplx);
        if (cimag(cplx) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_sec:
 * Implements the standard evaluation logic for Sec.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (1/ccos) if possible.
 */
Expr* builtin_sec(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Sec", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcSec"); if (inv) return inv; }

    // Sec is even: Sec[-x] -> Sec[x]
    { Expr* f = even_fold(arg, "Sec"); if (f) return f; }

    // Sec[I y] -> Sech[y]
    { Expr* f = trig_i_fold(arg, "Sech", 0); if (f) return f; }

    // Sec[0] = 1
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(1);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_sec(n, d);
        if (exact) return exact;
    }
    
    // Approximate numerical evaluation
    double complex cplx;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_sec);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_sec);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &cplx, &inexact) && inexact) {
        double complex s = 1.0 / ccos(cplx);
        if (cimag(cplx) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_csc:
 * Implements the standard evaluation logic for Csc.
 * Handles exact evaluations for integer zeros and specific Pi rational multiples.
 * Falls back to approximate numeric evaluation via C standard math (1/csin) if possible.
 */
Expr* builtin_csc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Csc", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCsc"); if (inv) return inv; }

    // Csc is odd: Csc[-x] -> -Csc[x]
    { Expr* f = odd_fold(arg, "Csc"); if (f) return f; }

    // Csc[I y] -> -I Csch[y]
    { Expr* f = trig_i_fold(arg, "Csch", -1); if (f) return f; }

    // Csc[0] = ComplexInfinity
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_symbol(SYM_ComplexInfinity);
    
    // Attempt exact evaluation if argument is a rational multiple of Pi
    int64_t n, d;
    if (extract_pi_multiplier(arg, &n, &d)) {
        Expr* exact = exact_csc(n, d);
        if (exact) return exact;
    }
    
    // Approximate numerical evaluation
    double complex cplx;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_csc);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_csc);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &cplx, &inexact) && inexact) {
        double complex s = 1.0 / csin(cplx);
        if (cimag(cplx) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

static Expr* exact_arcsin(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = -d; n <= d; n++) {
            Expr* val = exact_sin(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}

static Expr* exact_arccos(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n <= d; n++) {
            Expr* val = exact_cos(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}

static Expr* exact_arctan(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = -d/2; n <= d/2; n++) {
            Expr* val = exact_tan(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}

static Expr* exact_arccot(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n < d; n++) {
            Expr* val = exact_cot(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}

static Expr* exact_arcsec(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n <= d; n++) {
            Expr* val = exact_sec(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}

static Expr* exact_arccsc(Expr* arg) {
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = -d; n <= d; n++) {
            Expr* val = exact_csc(n, d);
            if (val) {
                if (expr_eq(arg, val)) {
                    expr_free(val);
                    return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
                }
                expr_free(val);
            }
        }
    }
    return NULL;
}


/*
 * builtin_arcsin:
 * Implements the evaluation logic for the inverse sine (ArcSin).
 * Reverses exact algebraic mappings for supported inputs.
 * Falls back to C standard math (casin) for approximate numerical inputs.
 */
Expr* builtin_arcsin(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcSin", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcSin is odd: ArcSin[-x] -> -ArcSin[x]
    { Expr* f = odd_fold(arg, "ArcSin"); if (f) return f; }

    // ArcSin[I y] -> I ArcSinh[y]  (principal-branch identity)
    { Expr* f = trig_i_fold(arg, "ArcSinh", +1); if (f) return f; }

    // Attempt exact inverse evaluation
    Expr* exact = exact_arcsin(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_asin);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_asin);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation - only if input is already inexact
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = casin(c);
        if (cimag(c) == 0.0 && creal(c) >= -1.0 && creal(c) <= 1.0) return expr_new_real(creal(s));
        double im = cimag(s);
        /* On the (1,inf) branch cut, C99 casin(x+0i) lands on the upper side
         * while Mathematica uses the lower side; flip imag for real x>1.
         * The (-inf,-1) cut already agrees between the two conventions. */
        if (cimag(c) == 0.0 && creal(c) > 1.0) im = -im;
        return make_complex(expr_new_real(creal(s)), expr_new_real(im));
    }
    return NULL;
}

/*
 * builtin_arccos:
 * Implements the evaluation logic for the inverse cosine (ArcCos).
 * Reverses exact algebraic mappings for supported inputs.
 * Falls back to C standard math (cacos) for approximate numerical inputs.
 */
Expr* builtin_arccos(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCos", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcCos[-x] -> Pi - ArcCos[x]
    { Expr* f = arc_pi_minus_fold(arg, "ArcCos"); if (f) return f; }

    // ArcCos[I y] -> Pi/2 - I ArcSinh[y]
    { Expr* f = arccos_i_fold(arg); if (f) return f; }

    // Attempt exact inverse evaluation
    Expr* exact = exact_arccos(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_acos);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acos);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation - only if input is already inexact
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = cacos(c);
        if (cimag(c) == 0.0 && creal(c) >= -1.0 && creal(c) <= 1.0) return expr_new_real(creal(s));
        double im = cimag(s);
        /* C99 cacos(x+0i) for x>1 returns -i*acosh(x); Mathematica returns
         * +i*acosh(x). The x<-1 cut already agrees. */
        if (cimag(c) == 0.0 && creal(c) > 1.0) im = -im;
        return make_complex(expr_new_real(creal(s)), expr_new_real(im));
    }
    return NULL;
}

/*
 * builtin_arctan:
 * Implements the evaluation logic for the inverse tangent (ArcTan).
 * Supports both ArcTan[z] (1 argument) and ArcTan[x, y] (2 arguments).
 * Reverses exact algebraic mappings for supported inputs.
 * Uses C standard math (catan and atan2) for numeric evaluations.
 */
Expr* builtin_arctan(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    // ArcTan accepts 1 (arctangent) or 2 (two-argument / quadrant) arguments.
    if (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)
        return builtin_arg_error("ArcTan", res->data.function.arg_count, 1, 2);

    // Single argument ArcTan[z]
    if (res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];

        // ArcTan is odd: ArcTan[-x] -> -ArcTan[x]
        { Expr* f = odd_fold(arg, "ArcTan"); if (f) return f; }

        // ArcTan[I y] -> I ArcTanh[y]
        { Expr* f = trig_i_fold(arg, "ArcTanh", +1); if (f) return f; }

        // ArcTan[Infinity] = Pi/2  (ArcTan[-Infinity] -> -Pi/2 via the odd fold above)
        if (is_infinity_sym(arg))
            return make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi));

        // Attempt exact inverse evaluation
        Expr* exact = exact_arctan(arg);
        if (exact) return exact;

#ifdef USE_MPFR
        if (numeric_expr_is_mpfr(arg)) {
            Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_atan);
            if (r) return r;
            r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_atan);
            if (r) return r;
        }
#endif
        // Approximate numerical evaluation - only if input is already inexact
        double complex c;
        bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
            double complex s = catan(c);
            if (cimag(c) == 0.0) return expr_new_real(creal(s));
            return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
        }
    } 
    // Two argument ArcTan[x, y] (computes argument of x + I*y)
    else if (res->data.function.arg_count == 2) {
        Expr* x = res->data.function.args[0];
        Expr* y = res->data.function.args[1];
        
        // Exact discrete evaluations based on sign quadrant mapping
        if (x->type == EXPR_INTEGER && y->type == EXPR_INTEGER) {
            int64_t xv = x->data.integer;
            int64_t yv = y->data.integer;
            if (xv == 0 && yv == 0) return NULL; // Indeterminate form
            if (yv == 0) {
                if (xv > 0) return expr_new_integer(0);
                if (xv < 0) return expr_new_symbol(SYM_Pi);
            }
            if (xv == 0) {
                if (yv > 0) return make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi));
                if (yv < 0) return make_times(make_rational(-1, 2), expr_new_symbol(SYM_Pi));
            }
            if (xv == yv) {
                if (xv > 0) return make_times(make_rational(1, 4), expr_new_symbol(SYM_Pi));
                if (xv < 0) return make_times(make_rational(-3, 4), expr_new_symbol(SYM_Pi));
            }
            if (xv == -yv) {
                if (xv > 0) return make_times(make_rational(-1, 4), expr_new_symbol(SYM_Pi));
                if (xv < 0) return make_times(make_rational(3, 4), expr_new_symbol(SYM_Pi));
            }
        }
        
#ifdef USE_MPFR
        /* MPFR two-argument ArcTan[x, y] for real inputs (either carrying
         * MPFR precision). Uses mpfr_atan2. */
        if (numeric_expr_is_mpfr(x) || numeric_expr_is_mpfr(y)) {
            long bits = numeric_combined_bits(x, y, 0);
            mpfr_t rx, ix, ry, iy, out;
            mpfr_init2(rx, bits); mpfr_init2(ix, bits);
            mpfr_init2(ry, bits); mpfr_init2(iy, bits);
            bool ok_x = get_approx_mpfr(x, rx, ix, NULL);
            bool ok_y = get_approx_mpfr(y, ry, iy, NULL);
            if (ok_x && ok_y && mpfr_zero_p(ix) && mpfr_zero_p(iy)) {
                mpfr_init2(out, bits);
                mpfr_atan2(out, ry, rx, MPFR_RNDN);
                mpfr_clear(rx); mpfr_clear(ix);
                mpfr_clear(ry); mpfr_clear(iy);
                return expr_new_mpfr_move(out);
            }
            mpfr_clear(rx); mpfr_clear(ix);
            mpfr_clear(ry); mpfr_clear(iy);
        }
#endif
        // Approximate numerical evaluation using atan2 for strictly real inputs
        double complex cx, cy;
        if ((x->type == EXPR_REAL || y->type == EXPR_REAL) && get_approx(x, &cx, NULL) && get_approx(y, &cy, NULL)) {
            if (cimag(cx) == 0.0 && cimag(cy) == 0.0) {
                return expr_new_real(atan2(creal(cy), creal(cx)));
            }
        }
    }
    return NULL;
}

/*
 * builtin_arccot:
 * Implements the evaluation logic for the inverse cotangent (ArcCot).
 * Reverses exact algebraic mappings for supported inputs.
 * Falls back to C standard math (catan) for approximate numerical inputs.
 */
Expr* builtin_arccot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCot", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcCot is odd (Mathematica convention): ArcCot[-x] -> -ArcCot[x]
    { Expr* f = odd_fold(arg, "ArcCot"); if (f) return f; }

    // ArcCot[I y] -> -I ArcCoth[y]
    { Expr* f = trig_i_fold(arg, "ArcCoth", -1); if (f) return f; }

    // ArcCot[Infinity] = 0  (ArcCot[-Infinity] -> 0 via the odd fold above)
    if (is_infinity_sym(arg)) return expr_new_integer(0);

    // Attempt exact inverse evaluation
    Expr* exact = exact_arccot(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_atan_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acot);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        if (c == 0.0) return expr_new_real(M_PI / 2.0); // ArcCot[0] = Pi/2
        double complex s = catan(1.0 / c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_arcsec:
 * Implements the evaluation logic for the inverse secant (ArcSec).
 * Reverses exact algebraic mappings for supported inputs.
 * Falls back to C standard math (cacos) for approximate numerical inputs.
 */
Expr* builtin_arcsec(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcSec", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcSec[-x] -> Pi - ArcSec[x]
    { Expr* f = arc_pi_minus_fold(arg, "ArcSec"); if (f) return f; }

    // ArcSec[+-Infinity] = Pi/2  (ArcSec[-Infinity] -> Pi - Pi/2 via the fold above)
    if (is_infinity_sym(arg))
        return make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi));

    // Attempt exact inverse evaluation
    Expr* exact = exact_arcsec(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_acos_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_asec);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        if (c == 0.0) return expr_new_symbol(SYM_ComplexInfinity); // ArcSec[0] = ComplexInfinity
        double complex s = cacos(1.0 / c);
        if (cimag(c) == 0.0 && (creal(c) <= -1.0 || creal(c) >= 1.0)) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * builtin_arccsc:
 * Implements the evaluation logic for the inverse cosecant (ArcCsc).
 * Reverses exact algebraic mappings for supported inputs.
 * Falls back to C standard math (casin) for approximate numerical inputs.
 */
Expr* builtin_arccsc(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCsc", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcCsc is odd: ArcCsc[-x] -> -ArcCsc[x]
    { Expr* f = odd_fold(arg, "ArcCsc"); if (f) return f; }

    // ArcCsc[I y] -> -I ArcCsch[y]
    { Expr* f = trig_i_fold(arg, "ArcCsch", -1); if (f) return f; }

    // ArcCsc[Infinity] = 0  (ArcCsc[-Infinity] -> 0 via the odd fold above)
    if (is_infinity_sym(arg)) return expr_new_integer(0);

    // Attempt exact inverse evaluation
    Expr* exact = exact_arccsc(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_asin_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acsc);
        if (r) return r;
    }
#endif
    // Approximate numerical evaluation
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        if (c == 0.0) return expr_new_symbol(SYM_ComplexInfinity); // ArcCsc[0] = ComplexInfinity
        double complex s = casin(1.0 / c);
        if (cimag(c) == 0.0 && (creal(c) <= -1.0 || creal(c) >= 1.0)) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

void trig_init(void) {
    symtab_add_builtin("Sin", builtin_sin);
    symtab_add_builtin("Cos", builtin_cos);
    symtab_add_builtin("Tan", builtin_tan);
    symtab_add_builtin("Cot", builtin_cot);
    symtab_add_builtin("Sec", builtin_sec);
    symtab_add_builtin("Csc", builtin_csc);
    
    symtab_add_builtin("ArcSin", builtin_arcsin);
    symtab_add_builtin("ArcCos", builtin_arccos);
    symtab_add_builtin("ArcTan", builtin_arctan);
    symtab_add_builtin("ArcCot", builtin_arccot);
    symtab_add_builtin("ArcSec", builtin_arcsec);
    symtab_add_builtin("ArcCsc", builtin_arccsc);
}
