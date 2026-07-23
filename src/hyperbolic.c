#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
#include "hyperbolic.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "complex.h"
#include "times.h"
#include "common.h"
#include "numeric.h"
#include "numeric_complex.h"
#include "sym_intern.h"
#include "sym_names.h"

#ifdef USE_MPFR
/* Reciprocal-input helpers for ArcCoth, ArcSech, ArcCsch. Each computes
 * f(1/x) in MPFR. */
static int mpfr_atanh_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_atanh(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
static int mpfr_acosh_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_acosh(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
static int mpfr_asinh_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    mpfr_t tmp; mpfr_init2(tmp, mpfr_get_prec(in));
    mpfr_ui_div(tmp, 1, in, rnd);
    int r = mpfr_asinh(out, tmp, rnd);
    mpfr_clear(tmp);
    return r;
}
#endif

void hyperbolic_init(void) {
    symtab_add_builtin("Sinh", builtin_sinh);
    symtab_add_builtin("Cosh", builtin_cosh);
    symtab_add_builtin("Tanh", builtin_tanh);
    symtab_add_builtin("Coth", builtin_coth);
    symtab_add_builtin("Sech", builtin_sech);
    symtab_add_builtin("Csch", builtin_csch);
    
    symtab_add_builtin("ArcSinh", builtin_arcsinh);
    symtab_add_builtin("ArcCosh", builtin_arccosh);
    symtab_add_builtin("ArcTanh", builtin_arctanh);
    symtab_add_builtin("ArcCoth", builtin_arccoth);
    symtab_add_builtin("ArcSech", builtin_arcsech);
    symtab_add_builtin("ArcCsch", builtin_arccsch);

    const char* hypers[] = {"Sinh", "Cosh", "Tanh", "Coth", "Sech", "Csch", 
                            "ArcSinh", "ArcCosh", "ArcTanh", "ArcCoth", "ArcSech", "ArcCsch", NULL};
    for (int i = 0; hypers[i] != NULL; i++) {
        symtab_get_def(hypers[i])->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    }
}

static bool get_approx(Expr* e, double complex* out, bool* is_inexact) {
    if (is_inexact) *is_inexact = true;
    if (e->type == EXPR_REAL) {
        *out = e->data.real + 0.0 * I;
        return true;
    }
    Expr* re; Expr* im;
    if (is_complex(e, &re, &im)) {
        bool has_real = false;
        double r = 0.0, i = 0.0;
        int64_t n, d;
        if (re->type == EXPR_REAL) { r = re->data.real; has_real = true; }
        else if (re->type == EXPR_INTEGER) r = (double)re->data.integer;
        else if (is_rational(re, &n, &d)) r = (double)n / d;
        else return false;
        
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

static bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Infinity;
}

/*
 * flip_sign:
 * Returns a newly allocated, evaluated expression equal to -e. Uses
 * the evaluator so nested Complex / Times structures canonicalise
 * (e.g. -1 * Complex[0, -2] -> Complex[0, 2]; -1 * Times[-k, x] -> k x).
 */
static Expr* flip_sign(Expr* e) {
    Expr* args[2] = { expr_new_integer(-1), expr_copy(e) };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), args, 2));
}

/*
 * even_fold / odd_fold:
 * Rewrite f[arg] as f[|arg|] (even) or -f[|arg|] (odd) when `arg` is
 * superficially negative (see expr_is_superficially_negative). Return
 * NULL otherwise. The caller's `res` is freed by evaluate_step after a
 * non-NULL return -- do NOT free it here.
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
 * If `arg` is of the form I*y (Complex[0, k] or Times[Complex[0, k], ...])
 * with k > 0, return a new evaluated expression equal to arg / I. Otherwise
 * NULL. Sign-extraction runs first, so k is always positive here.
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
 * hyp_i_fold:
 * Rewrites hyperbolic-of-I*y using the bridge:
 *   Cosh[I y]  -> Cos[y]          Sech[I y]  -> Sec[y]
 *   Sinh[I y]  -> I Sin[y]        Csch[I y]  -> -I Csc[y]
 *   Tanh[I y]  -> I Tan[y]        Coth[I y]  -> -I Cot[y]
 * `counterpart` is the trig head to emit. `coeff_kind`: 0 = no prefix,
 * +1 = I, -1 = -I.
 */
static Expr* hyp_i_fold(Expr* arg, const char* counterpart, int coeff_kind) {
    Expr* y = peel_imaginary_unit(arg);
    if (!y) return NULL;
    Expr* a[1] = { y };
    Expr* inner = expr_new_function(expr_new_symbol(counterpart), a, 1);
    if (coeff_kind == 0) return inner;
    Expr* coeff = make_complex(expr_new_integer(0), expr_new_integer(coeff_kind));
    return make_times(coeff, inner);
}

/*
 * try_simp_forward_of_inverse_hyp:
 * Universal forward-of-inverse identities for hyperbolic functions, mirroring
 * try_simp_forward_of_inverse in trig.c. `outer` is the forward head
 * ("Sinh","Cosh","Tanh",...), `arg` the body of that call. Returns the
 * rewritten Expr* on a hit, NULL otherwise. Identities follow the
 * principal-branch conventions:
 *   Sinh[ArcCosh[x]] -> Sqrt[x-1]*Sqrt[x+1]    (NB: not Sqrt[x^2-1])
 *   Cosh[ArcSinh[x]] -> Sqrt[1+x^2]
 *   Sinh[ArcTanh[x]] -> x/Sqrt[1-x^2]
 *   Cosh[ArcTanh[x]] -> 1/Sqrt[1-x^2]
 *   Tanh[ArcSinh[x]] -> x/Sqrt[1+x^2]
 *   Tanh[ArcCosh[x]] -> Sqrt[x-1]*Sqrt[x+1]/x
 *   Tanh[ArcCoth[x]] -> 1/x                    Coth[ArcTanh[x]] -> 1/x
 * The Sqrt[x-1]*Sqrt[x+1] form (rather than Sqrt[x^2-1]) is the
 * principal-branch identity Mathematica uses; the two forms differ on
 * the branch cut x in (-1,1) where Sqrt[x-1]*Sqrt[x+1] = i*Sqrt[1-x^2]
 * agrees with the Sinh[ArcCosh[x]] continuation but Sqrt[x^2-1] does not.
 */
static Expr* try_simp_forward_of_inverse_hyp(const char* outer, Expr* arg) {
    if (arg->type != EXPR_FUNCTION || arg->data.function.arg_count != 1) return NULL;
    if (!arg->data.function.head ||
        arg->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* inner = arg->data.function.head->data.symbol.name;
    Expr* x = arg->data.function.args[0];

    #define SQRT_OF(arg_e)                                                     \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),              \
            (Expr*[]){ (arg_e), make_rational(1, 2) }, 2))
    #define X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),              \
            (Expr*[]){ expr_copy(x), expr_new_integer(2) }, 2))
    #define ONE_PLUS_X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),               \
            (Expr*[]){ expr_new_integer(1), X_SQ() }, 2))
    #define ONE_MINUS_X_SQ() \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),               \
            (Expr*[]){ expr_new_integer(1),                                    \
                       eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), \
                           (Expr*[]){ expr_new_integer(-1), X_SQ() }, 2)) }, 2))
    #define X_PLUS_K(k) \
        eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus),               \
            (Expr*[]){ expr_new_integer(k), expr_copy(x) }, 2))
    /* Sqrt[(x-1)(x+1)] expressed as Sqrt[x-1]*Sqrt[x+1] -- the principal-
     * branch form Mathematica returns for Sinh[ArcCosh[x]]. Storing it as
     * a product of two single-arg sqrts (rather than Sqrt[x^2-1]) preserves
     * the right branch cut and lets the assumption-aware rules in simp.c
     * collapse it cleanly when x > 1. */
    #define SQRT_X_MINUS_1_TIMES_SQRT_X_PLUS_1() \
        make_times(SQRT_OF(X_PLUS_K(-1)), SQRT_OF(X_PLUS_K(1)))

    if (strcmp(outer, "Sinh") == 0 && strcmp(inner, "ArcCosh") == 0) {
        return SQRT_X_MINUS_1_TIMES_SQRT_X_PLUS_1();
    }
    if (strcmp(outer, "Cosh") == 0 && strcmp(inner, "ArcSinh") == 0) {
        return SQRT_OF(ONE_PLUS_X_SQ());
    }
    if (strcmp(outer, "Sinh") == 0 && strcmp(inner, "ArcTanh") == 0) {
        Expr* den = SQRT_OF(ONE_MINUS_X_SQ());
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
        return make_times(expr_copy(x), inv);
    }
    if (strcmp(outer, "Cosh") == 0 && strcmp(inner, "ArcTanh") == 0) {
        Expr* den = SQRT_OF(ONE_MINUS_X_SQ());
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
    }
    if (strcmp(outer, "Tanh") == 0 && strcmp(inner, "ArcSinh") == 0) {
        Expr* den = SQRT_OF(ONE_PLUS_X_SQ());
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ den, expr_new_integer(-1) }, 2));
        return make_times(expr_copy(x), inv);
    }
    if (strcmp(outer, "Tanh") == 0 && strcmp(inner, "ArcCosh") == 0) {
        Expr* num = SQRT_X_MINUS_1_TIMES_SQRT_X_PLUS_1();
        Expr* inv_x = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(x), expr_new_integer(-1) }, 2));
        return make_times(num, inv_x);
    }
    if ((strcmp(outer, "Tanh") == 0 && strcmp(inner, "ArcCoth") == 0) ||
        (strcmp(outer, "Coth") == 0 && strcmp(inner, "ArcTanh") == 0)) {
        return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power),
            (Expr*[]){ expr_copy(x), expr_new_integer(-1) }, 2));
    }

    #undef SQRT_OF
    #undef X_SQ
    #undef ONE_PLUS_X_SQ
    #undef ONE_MINUS_X_SQ
    #undef X_PLUS_K
    #undef SQRT_X_MINUS_1_TIMES_SQRT_X_PLUS_1
    return NULL;
}

/* If arg is a one-argument call whose head is `inverse_name`, return a deep
 * copy of its single argument; otherwise NULL. Folds the direct
 * hyperbolic forward/inverse identities (Sinh[ArcSinh[x]] -> x, etc.):
 * each ArcXh is a right inverse of Xh by construction, so the identity
 * holds over the complex numbers. The opposite direction is NOT folded
 * because those only reduce on each function's principal domain. */
static Expr* strip_inverse_call(Expr* arg, const char* inverse_name) {
    if (head_is(arg, intern_symbol(inverse_name)) &&
        arg->data.function.arg_count == 1) {
        return expr_copy(arg->data.function.args[0]);
    }
    return NULL;
}

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

/* True for any flavour of infinite argument: Infinity, -Infinity,
 * ComplexInfinity, or DirectedInfinity[_]. The reciprocal inverse-hyperbolic
 * functions (ArcCoth, ArcSech, ArcCsch) all satisfy f[z] = g[1/z], and
 * 1/z -> 0 for every one of these, so f attains its value at 0 regardless of
 * the direction of approach. */
static bool is_any_infinity(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL &&
        (e->data.symbol.name == SYM_Infinity ||
         e->data.symbol.name == SYM_ComplexInfinity)) return true;
    if (is_minus_infinity(e)) return true;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity) return true;
    return false;
}

/* Builds the constant  s * I * Pi/2  (s = +1 or -1). Used for the
 * inverse-hyperbolic values at infinity that land on the imaginary axis:
 * ArcSech[oo] = I Pi/2, ArcTanh[Infinity] = -I Pi/2. */
static Expr* half_i_pi(int s) {
    return make_times(make_complex(expr_new_integer(0), expr_new_integer(s)),
                      make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi)));
}

Expr* builtin_sinh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Sinh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcSinh"); if (inv) return inv; }

    // Sinh[ArcCosh[x]] -> Sqrt[x-1]*Sqrt[x+1], Sinh[ArcTanh[x]] -> x/Sqrt[1-x^2]
    { Expr* f = try_simp_forward_of_inverse_hyp("Sinh", arg); if (f) return f; }

    // Sinh is odd: Sinh[-x] -> -Sinh[x]
    { Expr* f = odd_fold(arg, "Sinh"); if (f) return f; }

    // Sinh[I y] -> I Sin[y]
    { Expr* f = hyp_i_fold(arg, "Sin", +1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    if (is_infinity(arg)) return expr_new_symbol(SYM_Infinity);
    if (is_minus_infinity(arg)) {
        Expr* args[2] = { expr_new_integer(-1), expr_new_symbol(SYM_Infinity) };
        return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
    }

#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_sinh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_sinh);
        if (r) return r;
    }
#endif
    double complex c;
    bool inexact = false;
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = csinh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_cosh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Cosh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCosh"); if (inv) return inv; }

    // Cosh[ArcSinh[x]] -> Sqrt[1+x^2], Cosh[ArcTanh[x]] -> 1/Sqrt[1-x^2]
    { Expr* f = try_simp_forward_of_inverse_hyp("Cosh", arg); if (f) return f; }

    // Cosh is even: Cosh[-x] -> Cosh[x]
    { Expr* f = even_fold(arg, "Cosh"); if (f) return f; }

    // Cosh[I y] -> Cos[y]
    { Expr* f = hyp_i_fold(arg, "Cos", 0); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(1);
    if (is_infinity(arg) || is_minus_infinity(arg)) return expr_new_symbol(SYM_Infinity);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_cosh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_cosh);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = ccosh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_tanh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Tanh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcTanh"); if (inv) return inv; }

    // Tanh[ArcSinh[x]] -> x/Sqrt[1+x^2], Tanh[ArcCosh[x]] -> Sqrt[x-1]*Sqrt[x+1]/x,
    // Tanh[ArcCoth[x]] -> 1/x.
    { Expr* f = try_simp_forward_of_inverse_hyp("Tanh", arg); if (f) return f; }

    // Tanh is odd: Tanh[-x] -> -Tanh[x]
    { Expr* f = odd_fold(arg, "Tanh"); if (f) return f; }

    // Tanh[I y] -> I Tan[y]
    { Expr* f = hyp_i_fold(arg, "Tan", +1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    if (is_infinity(arg)) return expr_new_integer(1);
    if (is_minus_infinity(arg)) return expr_new_integer(-1);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_tanh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_tanh);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = ctanh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_coth(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Coth", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCoth"); if (inv) return inv; }

    // Coth[ArcTanh[x]] -> 1/x
    { Expr* f = try_simp_forward_of_inverse_hyp("Coth", arg); if (f) return f; }

    // Coth is odd: Coth[-x] -> -Coth[x]
    { Expr* f = odd_fold(arg, "Coth"); if (f) return f; }

    // Coth[I y] -> -I Cot[y]
    { Expr* f = hyp_i_fold(arg, "Cot", -1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_symbol(SYM_ComplexInfinity);
    if (is_infinity(arg)) return expr_new_integer(1);
    if (is_minus_infinity(arg)) return expr_new_integer(-1);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_coth);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_coth);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = 1.0 / ctanh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_sech(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Sech", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcSech"); if (inv) return inv; }

    // Sech is even: Sech[-x] -> Sech[x]
    { Expr* f = even_fold(arg, "Sech"); if (f) return f; }

    // Sech[I y] -> Sec[y]
    { Expr* f = hyp_i_fold(arg, "Sec", 0); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(1);
    if (is_infinity(arg) || is_minus_infinity(arg)) return expr_new_integer(0);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_sech);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_sech);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = 1.0 / ccosh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_csch(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("Csch", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    { Expr* inv = strip_inverse_call(arg, "ArcCsch"); if (inv) return inv; }

    // Csch is odd: Csch[-x] -> -Csch[x]
    { Expr* f = odd_fold(arg, "Csch"); if (f) return f; }

    // Csch[I y] -> -I Csc[y]
    { Expr* f = hyp_i_fold(arg, "Csc", -1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_symbol(SYM_ComplexInfinity);
    if (is_infinity(arg) || is_minus_infinity(arg)) return expr_new_integer(0);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_csch);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_csch);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = 1.0 / csinh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_arcsinh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcSinh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcSinh is odd: ArcSinh[-x] -> -ArcSinh[x]
    { Expr* f = odd_fold(arg, "ArcSinh"); if (f) return f; }

    // ArcSinh[I y] -> I ArcSin[y]  (principal-branch identity)
    { Expr* f = hyp_i_fold(arg, "ArcSin", +1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    if (is_infinity(arg)) return expr_new_symbol(SYM_Infinity);
    // ArcSinh grows without bound; a directionless infinite argument maps to
    // ComplexInfinity (magnitude infinite, direction undefined).
    if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_ComplexInfinity)
        return expr_new_symbol(SYM_ComplexInfinity);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_asinh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_asinh);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = casinh(c);
        if (cimag(c) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_arccosh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCosh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];
    
    if (arg->type == EXPR_INTEGER && arg->data.integer == 1) return expr_new_integer(0);
    if (is_infinity(arg)) return expr_new_symbol(SYM_Infinity);
    // ArcCosh grows without bound; a directionless infinite argument maps to
    // ComplexInfinity.
    if (arg->type == EXPR_SYMBOL && arg->data.symbol.name == SYM_ComplexInfinity)
        return expr_new_symbol(SYM_ComplexInfinity);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_acosh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acosh);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = cacosh(c);
        if (cimag(c) == 0.0 && creal(c) >= 1.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_arctanh(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcTanh", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcTanh is odd: ArcTanh[-x] -> -ArcTanh[x]
    { Expr* f = odd_fold(arg, "ArcTanh"); if (f) return f; }

    // ArcTanh[I y] -> I ArcTan[y]
    { Expr* f = hyp_i_fold(arg, "ArcTan", +1); if (f) return f; }

    if (arg->type == EXPR_INTEGER && arg->data.integer == 0) return expr_new_integer(0);
    /* ArcTanh[1] = Infinity (a genuine pole).  ArcTanh[-1] is reached as
     * -ArcTanh[1] via the odd fold above, giving -Infinity. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 1)
        return expr_new_symbol(SYM_Infinity);
    /* ArcTanh[Infinity] = -I Pi/2 (real part vanishes, imaginary -> -Pi/2).
     * ArcTanh[-Infinity] = +I Pi/2 is reached via the odd fold above.
     * ComplexInfinity / DirectedInfinity are left unevaluated: the value
     * (+/- I Pi/2) depends on the direction of approach. */
    if (is_infinity(arg)) return half_i_pi(-1);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_atanh);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_atanh);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = catanh(c);
        if (cimag(c) == 0.0 && creal(c) > -1.0 && creal(c) < 1.0) return expr_new_real(creal(s));
        double im = cimag(s);
        /* C99 catanh(x+0i) for x>1 lands on the +iPi/2 side of the branch cut;
         * Mathematica uses -iPi/2. The x<-1 cut already agrees. */
        if (cimag(c) == 0.0 && creal(c) > 1.0) im = -im;
        return make_complex(expr_new_real(creal(s)), expr_new_real(im));
    }
    return NULL;
}

Expr* builtin_arccoth(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCoth", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcCoth is odd: ArcCoth[-x] -> -ArcCoth[x]
    { Expr* f = odd_fold(arg, "ArcCoth"); if (f) return f; }

    // ArcCoth[I y] -> -I ArcCot[y]
    { Expr* f = hyp_i_fold(arg, "ArcCot", -1); if (f) return f; }

    /* ArcCoth[z] = ArcTanh[1/z] -> ArcTanh[0] = 0 for every infinite z. */
    if (is_any_infinity(arg)) return expr_new_integer(0);
    /* ArcCoth[1] = Infinity (a genuine pole).  ArcCoth[-1] is reached as
     * -ArcCoth[1] via the odd fold above, giving -Infinity. */
    if (arg->type == EXPR_INTEGER && arg->data.integer == 1)
        return expr_new_symbol(SYM_Infinity);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_atanh_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acoth);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = catanh(1.0 / c);
        if (cimag(c) == 0.0 && (creal(c) > 1.0 || creal(c) < -1.0)) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_arcsech(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcSech", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER && arg->data.integer == 1) return expr_new_integer(0);
    /* ArcSech[z] = ArcCosh[1/z] -> ArcCosh[0] = I Pi/2 for every infinite z. */
    if (is_any_infinity(arg)) return half_i_pi(+1);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_acosh_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_asech);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = cacosh(1.0 / c);
        if (cimag(c) == 0.0 && creal(c) > 0.0 && creal(c) <= 1.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

Expr* builtin_arccsch(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcCsch", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    // ArcCsch is odd: ArcCsch[-x] -> -ArcCsch[x]
    { Expr* f = odd_fold(arg, "ArcCsch"); if (f) return f; }

    // ArcCsch[I y] -> -I ArcCsc[y]
    { Expr* f = hyp_i_fold(arg, "ArcCsc", -1); if (f) return f; }

    /* ArcCsch[z] = ArcSinh[1/z] -> ArcSinh[0] = 0 for every infinite z. */
    if (is_any_infinity(arg)) return expr_new_integer(0);

    double complex c;
    bool inexact = false;
#ifdef USE_MPFR
    if (numeric_expr_is_mpfr(arg)) {
        Expr* r = numeric_mpfr_apply_unary(arg, 0, mpfr_asinh_recip_op);
        if (r) return r;
        r = numeric_mpfr_apply_complex_unary(arg, 0, mpfr_complex_acsch);
        if (r) return r;
    }
#endif
    if (get_approx(arg, &c, &inexact) && inexact) {
        double complex s = casinh(1.0 / c);
        if (cimag(c) == 0.0 && creal(c) != 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}
