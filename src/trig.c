#include <complex.h>
#include "trig.h"
#include "mth_machine_elem.h"   /* the one machine elementary function (scalar/ND/VM) */
#include "arithmetic.h"
#include "interval.h"
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
#include "message.h"   /* mth_message: ArcTan::indet */
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
/* ArcSec / ArcCsc of a real x, |x| >= 1, without forming 1/x: acos(1/x) and
 * asin(1/x) are infinitely sensitive at x = +-1, so the rounding of 1/x cost
 * ~half the digits next to them (ArcSec[1 + 10^-25] at 60 digits had 41).
 * With r = sqrt((x - 1)(x + 1)), x - 1 exact there:
 *   acos(1/x) = atan2(r, sign x),   asin(1/x) = atan2(sign x, r).
 * |x| < 1 makes r NaN, the domain failure the complex fallback expects. */
static int mpfr_recip_acos_asin(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd, bool asin_) {
    mpfr_prec_t wp = mpfr_get_prec(out) + 32;
    mpfr_t r, t, sg;
    mpfr_inits2(wp, r, t, sg, (mpfr_ptr)NULL);
    mpfr_sub_ui(r, in, 1, MPFR_RNDN);
    mpfr_add_ui(t, in, 1, MPFR_RNDN);
    mpfr_mul(r, r, t, MPFR_RNDN);
    mpfr_sqrt(r, r, MPFR_RNDN);
    mpfr_abs(r, r, MPFR_RNDN);           /* sqrt(-0) = -0 at x = -1: atan2 -> -Pi */
    mpfr_set_si(sg, mpfr_sgn(in) < 0 ? -1 : 1, MPFR_RNDN);
    int ret = asin_ ? mpfr_atan2(out, sg, r, rnd) : mpfr_atan2(out, r, sg, rnd);
    mpfr_clears(r, t, sg, (mpfr_ptr)NULL);
    return ret;
}
static int mpfr_acos_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    return mpfr_recip_acos_asin(out, in, rnd, false);
}
static int mpfr_asin_recip_op(mpfr_t out, const mpfr_t in, mpfr_rnd_t rnd) {
    return mpfr_recip_acos_asin(out, in, rnd, true);
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
 * Recognises a rational multiple of Pi and returns n, d with the argument
 * equal to (n/d) * Pi. Two spellings fold to the same normal form:
 *   - directly against Pi:      Pi -> 1/1 * Pi,   c * Pi -> c * Pi
 *   - against Degree (= Pi/180): Degree -> Pi/180, c * Degree -> c/180 * Pi
 *
 * The Degree case is what makes Sin[30 Degree] fold to 1/2. Every exact-value
 * trig head (Sin/Cos/Tan/Cot/Sec/Csc) funnels its argument through this one
 * detector, so Degree is interpreted in a single place and costs only a symbol
 * pointer comparison on the common Pi path -- no evaluator-wide substitution.
 * Degree remains an opaque numeric constant everywhere else (30 Degree does
 * not auto-expand to Pi/6), matching Mathematica.
 *
 * The (n, d) returned for a Degree argument need not be in lowest terms; the
 * exact_* callees reduce via gcd. Denominators are bounded before scaling by
 * 180 so the 2*d periodicity reduction downstream cannot overflow int64 -- a
 * pathological coefficient falls through to the numeric path rather than
 * folding a wrapped value into a wrong exact answer.
 */
/*
 * reduce_big_pi_coefficient:
 * c * Pi with c a big exact rational -- a BigInt, or a Rational whose
 * numerator is a BigInt -- is reduced by the 2 Pi period before the exact
 * tables see it: the numerator is replaced by its residue mod 2 d (GMP), so
 * Sin[10^30 Pi] -> Sin[0 Pi] and Sin[(10^30 + 1/2) Pi] -> Sin[Pi/2].  Every
 * trig head here is 2 Pi- (or Pi-) periodic, so the residue is exact.  The
 * denominator must itself be a machine integer small enough that the
 * downstream 2 d arithmetic cannot overflow; otherwise decline.
 */
static bool reduce_big_pi_coefficient(Expr* c, int64_t* n, int64_t* d) {
    Expr* num = NULL;
    int64_t den = 1;
    if (c->type == EXPR_BIGINT) {
        num = c;
    } else if (c->type == EXPR_FUNCTION && head_is(c, SYM_Rational) &&
               c->data.function.arg_count == 2 &&
               c->data.function.args[0]->type == EXPR_BIGINT &&
               c->data.function.args[1]->type == EXPR_INTEGER) {
        num = c->data.function.args[0];
        den = c->data.function.args[1]->data.integer;
    } else {
        return false;
    }
    if (den <= 0 || den > INT64_MAX / 8) return false;
    mpz_t z, m;
    expr_to_mpz(num, z);                                 /* initialises z */
    mth_mpz_init_set_i64(m, den);
    mpz_mul_2exp(m, m, 1);                               /* m = 2 den */
    mpz_fdiv_r(z, z, m);                       /* 0 <= z < 2 den */
    int64_t r = 0;
    bool ok = expr_mpz_get_int64(z, &r);
    mpz_clear(z); mpz_clear(m);
    if (!ok) return false;
    *n = r; *d = den;
    return true;
}

static bool extract_pi_multiplier(Expr* e, int64_t* n, int64_t* d) {
    // Case 1: Pi
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Pi) {
        *n = 1; *d = 1;
        return true;
    }
    // Case 2: Degree = Pi/180
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_Degree) {
        *n = 1; *d = 180;
        return true;
    }

    // Case 3/4: c * Pi or c * Degree with c a rational literal.
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Times &&
        e->data.function.arg_count == 2) {
        Expr* first  = e->data.function.args[0];
        Expr* second = e->data.function.args[1];

        if (second->type == EXPR_SYMBOL) {
            const char* s = second->data.symbol.name;
            int64_t rn, rd;
            if (s == SYM_Pi && is_rational(first, &rn, &rd)) {
                *n = rn; *d = rd;
                return true;
            }
            if (s == SYM_Pi && reduce_big_pi_coefficient(first, &rn, &rd)) {
                *n = rn; *d = rd;
                return true;
            }
            if (s == SYM_Degree && is_rational(first, &rn, &rd)) {
                // c * Degree = c/180 * Pi. Bound rd so that neither rd*180 nor
                // the downstream 2*d (in exact_*) overflows int64.
                if (rd > INT64_MAX / 360 || rd < INT64_MIN / 360) return false;
                *n = rn; *d = rd * 180;
                return true;
            }
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
    if (!expr_is_superficially_negative(arg) && !expr_is_negated_sum(arg)) return NULL;
    Expr* neg = flip_sign(arg);
    Expr* inner_args[1] = { neg };
    return expr_new_function(expr_new_symbol(head_name), inner_args, 1);
}

static Expr* odd_fold(Expr* arg, const char* head_name) {
    if (!expr_is_superficially_negative(arg) && !expr_is_negated_sum(arg)) return NULL;
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
 * Forward-of-inverse identities for the six circular functions. `outer` is the
 * forward head ("Sin", "Cos", "Tan", "Cot", "Sec", "Csc"), `arg` the body of
 * that call. When `arg` is g[x] for an inverse head g, the pair
 * (s, c) = (Sin[g[x]], Cos[g[x]]) has a closed algebraic form on the principal
 * branch, and every outer head is one of s, c, s/c, c/s, 1/c, 1/s:
 *   g = ArcSin:  s = x,                    c = Sqrt[1 - x^2]
 *   g = ArcCos:  s = Sqrt[1 - x^2],        c = x
 *   g = ArcTan:  s = x/Sqrt[1 + x^2],      c = 1/Sqrt[1 + x^2]
 *   g = ArcCot:  s = 1/(x Sqrt[1 + 1/x^2]), c = 1/Sqrt[1 + 1/x^2]
 *   g = ArcSec:  s = Sqrt[1 - 1/x^2],      c = 1/x
 *   g = ArcCsc:  s = 1/x,                  c = Sqrt[1 - 1/x^2]
 * The last three are the first three at 1/x, by the definitions
 * ArcCot[x] = ArcTan[1/x], ArcSec[x] = ArcCos[1/x], ArcCsc[x] = ArcSin[1/x].
 * Each identity holds for every complex x off the branch points
 * (Sin[ArcCos[x]] = Sqrt[1 - x^2] is the principal square root throughout).
 * The pairing f[f^-1[x]] = x is handled by the callers (strip_inverse_call).
 * The inverse-of-forward direction (e.g. ArcSin[Sin[x]]) is intentionally
 * NOT folded -- those reduce to x only on the principal domain of each f.
 */
static Expr* try_simp_forward_of_inverse(const char* outer, Expr* arg) {
    if (arg->type != EXPR_FUNCTION || arg->data.function.arg_count != 1) return NULL;
    if (!arg->data.function.head ||
        arg->data.function.head->type != EXPR_SYMBOL) return NULL;
    const char* inner = arg->data.function.head->data.symbol.name;
    int g = inner == SYM_ArcSin ? 0 : inner == SYM_ArcCos ? 1
          : inner == SYM_ArcTan ? 2 : inner == SYM_ArcCot ? 3
          : inner == SYM_ArcSec ? 4 : inner == SYM_ArcCsc ? 5 : -1;
    if (g < 0) return NULL;
    int f = strcmp(outer, "Sin") == 0 ? 0 : strcmp(outer, "Cos") == 0 ? 1
          : strcmp(outer, "Tan") == 0 ? 2 : strcmp(outer, "Cot") == 0 ? 3
          : strcmp(outer, "Sec") == 0 ? 4 : strcmp(outer, "Csc") == 0 ? 5 : -1;
    if (f < 0) return NULL;

    /* t = x, or 1/x for the reciprocal inverses (ArcCot/ArcSec/ArcCsc) */
    Expr* x = arg->data.function.args[0];
    Expr* t = (g >= 3)
        ? expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(x), expr_new_integer(-1) }, 2)
        : expr_copy(x);
    /* r = Sqrt[1 - t^2] (ArcSin/ArcCos/ArcSec/ArcCsc) or Sqrt[1 + t^2] (ArcTan/ArcCot) */
    int sign = (g == 2 || g == 3) ? 1 : -1;
    Expr* t2 = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(t), expr_new_integer(2) }, 2);
    Expr* rad = expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_new_integer(1),
                    expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ expr_new_integer(sign), t2 }, 2) }, 2);
    Expr* r = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ rad, make_rational(1, 2) }, 2);
    Expr *sn, *cs;
    switch (g) {
        case 0: case 5: sn = t; cs = r; break;                   /* ArcSin, ArcCsc */
        case 1: case 4: sn = r; cs = t; break;                   /* ArcCos, ArcSec */
        default:                                                 /* ArcTan, ArcCot */
            sn = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ t,
                     expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ expr_copy(r), expr_new_integer(-1) }, 2) }, 2);
            cs = expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ r, expr_new_integer(-1) }, 2);
            break;
    }
    #define RECIP(e_) expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){ (e_), expr_new_integer(-1) }, 2)
    Expr* out;
    switch (f) {
        case 0:  out = sn; expr_free(cs); break;
        case 1:  out = cs; expr_free(sn); break;
        case 2:  out = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ sn, RECIP(cs) }, 2); break;
        case 3:  out = expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ cs, RECIP(sn) }, 2); break;
        case 4:  out = RECIP(cs); expr_free(sn); break;
        default: out = RECIP(sn); expr_free(cs); break;
    }
    #undef RECIP
    return eval_and_free(out);
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
        *out = mth_cmplx(e->data.real, 0.0);   /* keeps -0. (G7) */
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
            *out = mth_cmplx(creal(r), creal(i));
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
    else if (d == 8) {
        if (n == 1) res = make_plus(expr_new_integer(-1), make_sqrt(2));
        else if (n == 3) res = make_plus(expr_new_integer(1), make_sqrt(2));
    }
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
    else if (d == 8) {
        if (n == 1) res = make_plus(expr_new_integer(1), make_sqrt(2));
        else if (n == 3) res = make_plus(expr_new_integer(-1), make_sqrt(2));
    }
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
        /* Sin[Pi/5] = Sqrt[5/8 - Sqrt[5]/8], Sin[2 Pi/5] = Sqrt[5/8 + Sqrt[5]/8] */
        if (n == 1) res = make_power(make_sqrt_expr(make_plus(make_rational(5, 8), make_times(make_rational(-1, 8), make_sqrt(5)))), expr_new_integer(-1));
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
/*
 * trig_quarter_shift:
 * f[theta + k Pi/2] for an integer k (a Plus argument whose one exact
 * Pi-term is a multiple of Pi/2) is rewritten by the quarter-period
 * identities: Sin[x + Pi] -> -Sin[x], Sin[x + Pi/2] -> Cos[x],
 * Cos[Pi/2 - x] -> Sin[x], Tan[x + Pi/2] -> -Cot[x], Sin[x + 2 Pi] -> Sin[x].
 * A shift that is not a multiple of Pi/2 (x + Pi/3) is left alone, as is a
 * Degree term.  `head` is the interned name of the calling head; returns
 * NULL when no rewrite applies.  The result is unevaluated: the evaluator
 * finishes it (folding Sin[-x] etc.).
 */
static Expr* trig_quarter_shift(Expr* arg, const char* head) {
    if (!head_is(arg, SYM_Plus) || arg->data.function.arg_count < 2) return NULL;
    size_t nargs = arg->data.function.arg_count;
    size_t pi_at = nargs;
    for (size_t i = 0; i < nargs; i++) {
        Expr* t = arg->data.function.args[i];
        bool is_pi_term =
            (t->type == EXPR_SYMBOL && t->data.symbol.name == SYM_Pi) ||
            (head_is(t, SYM_Times) && t->data.function.arg_count == 2 &&
             t->data.function.args[1]->type == EXPR_SYMBOL &&
             t->data.function.args[1]->data.symbol.name == SYM_Pi);
        if (!is_pi_term) continue;
        if (pi_at != nargs) return NULL;        /* two Pi terms: not canonical */
        pi_at = i;
    }
    if (pi_at == nargs) return NULL;
    int64_t n, d;
    if (!extract_pi_multiplier(arg->data.function.args[pi_at], &n, &d)) return NULL;
    if (d != 1 && d != 2) return NULL;
    int64_t k = (d == 1) ? (n % 2) * 2 : (n % 4);   /* shift = k Pi/2 mod 2 Pi */
    k = ((k % 4) + 4) % 4;

    /* Per head: the result head and sign for k = 0, 1, 2, 3. */
    const char* to[4]; int sg[4];
    if (head == SYM_Sin) {
        to[0] = SYM_Sin; to[1] = SYM_Cos; to[2] = SYM_Sin; to[3] = SYM_Cos;
        sg[0] = 1; sg[1] = 1; sg[2] = -1; sg[3] = -1;
    } else if (head == SYM_Cos) {
        to[0] = SYM_Cos; to[1] = SYM_Sin; to[2] = SYM_Cos; to[3] = SYM_Sin;
        sg[0] = 1; sg[1] = -1; sg[2] = -1; sg[3] = 1;
    } else if (head == SYM_Csc) {
        to[0] = SYM_Csc; to[1] = SYM_Sec; to[2] = SYM_Csc; to[3] = SYM_Sec;
        sg[0] = 1; sg[1] = 1; sg[2] = -1; sg[3] = -1;
    } else if (head == SYM_Sec) {
        to[0] = SYM_Sec; to[1] = SYM_Csc; to[2] = SYM_Sec; to[3] = SYM_Csc;
        sg[0] = 1; sg[1] = -1; sg[2] = -1; sg[3] = 1;
    } else if (head == SYM_Tan) {
        to[0] = SYM_Tan; to[1] = SYM_Cot; to[2] = SYM_Tan; to[3] = SYM_Cot;
        sg[0] = 1; sg[1] = -1; sg[2] = 1; sg[3] = -1;
    } else if (head == SYM_Cot) {
        to[0] = SYM_Cot; to[1] = SYM_Tan; to[2] = SYM_Cot; to[3] = SYM_Tan;
        sg[0] = 1; sg[1] = -1; sg[2] = 1; sg[3] = -1;
    } else {
        return NULL;
    }

    /* theta = the Plus without its Pi term */
    Expr* theta;
    if (nargs == 2) {
        theta = expr_copy(arg->data.function.args[pi_at == 0 ? 1 : 0]);
    } else {
        Expr** rest = (Expr**)malloc((nargs - 1) * sizeof(Expr*));
        if (!rest) return NULL;
        size_t j = 0;
        for (size_t i = 0; i < nargs; i++)
            if (i != pi_at) rest[j++] = expr_copy(arg->data.function.args[i]);
        theta = expr_new_function(expr_new_symbol(SYM_Plus), rest, nargs - 1);
        free(rest);
    }
    Expr* inner_args[1] = { theta };
    Expr* out = expr_new_function(expr_new_symbol(to[k]), inner_args, 1);
    if (sg[k] < 0) out = make_times(expr_new_integer(-1), out);
    return out;
}

static Expr* strip_inverse_call(Expr* arg, const char* inverse_name) {
    /* inverse_name is an interned SYM_* pointer: no string hash per call */
    if (head_is(arg, inverse_name) &&
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Sin", arg); if (r) return r; }
#ifdef USE_MPFR
    /* Sin[I y] = I Sinh[y]: the zero real part stays exact */
    { Expr* ax = numeric_mpfr_imag_axis(arg, "Sinh"); if (ax) return ax; }
#endif

    { Expr* inv = strip_inverse_call(arg, SYM_ArcSin); if (inv) return inv; }

    // Sin[ArcCos[x]] -> Sqrt[1-x^2], Sin[ArcTan[x]] -> x/Sqrt[1+x^2]
    { Expr* f = try_simp_forward_of_inverse("Sin", arg); if (f) return f; }

    // Sin is odd: Sin[-x] -> -Sin[x]
    { Expr* f = odd_fold(arg, "Sin"); if (f) return f; }

    // Sin[I y] -> I Sinh[y]
    { Expr* f = trig_i_fold(arg, "Sinh", +1); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Sin); if (f) return f; }

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
        double complex s = mth_sin_machine_z(c);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Cos", arg); if (r) return r; }

    { Expr* inv = strip_inverse_call(arg, SYM_ArcCos); if (inv) return inv; }

    // Cos[ArcSin[x]] -> Sqrt[1-x^2], Cos[ArcTan[x]] -> 1/Sqrt[1+x^2]
    { Expr* f = try_simp_forward_of_inverse("Cos", arg); if (f) return f; }

    // Cos is even: Cos[-x] -> Cos[x]
    { Expr* f = even_fold(arg, "Cos"); if (f) return f; }

    // Cos[I y] -> Cosh[y]
    { Expr* f = trig_i_fold(arg, "Cosh", 0); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Cos); if (f) return f; }

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
        double complex s = mth_cos_machine_z(c);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Tan", arg); if (r) return r; }
#ifdef USE_MPFR
    /* Tan[I y] = I Tanh[y]: the zero real part stays exact */
    { Expr* ax = numeric_mpfr_imag_axis(arg, "Tanh"); if (ax) return ax; }
#endif

    { Expr* inv = strip_inverse_call(arg, SYM_ArcTan); if (inv) return inv; }

    // Tan[ArcSin[x]] -> x/Sqrt[1-x^2], Tan[ArcCos[x]] -> Sqrt[1-x^2]/x,
    // Tan[ArcCot[x]] -> 1/x.
    { Expr* f = try_simp_forward_of_inverse("Tan", arg); if (f) return f; }

    // Tan is odd: Tan[-x] -> -Tan[x]
    { Expr* f = odd_fold(arg, "Tan"); if (f) return f; }

    // Tan[I y] -> I Tanh[y]
    { Expr* f = trig_i_fold(arg, "Tanh", +1); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Tan); if (f) return f; }

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
        double complex s = mth_tan_machine_z(cplx);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Cot", arg); if (r) return r; }

    { Expr* inv = strip_inverse_call(arg, SYM_ArcCot); if (inv) return inv; }

    // Cot[ArcTan[x]] -> 1/x
    { Expr* f = try_simp_forward_of_inverse("Cot", arg); if (f) return f; }

    // Cot is odd: Cot[-x] -> -Cot[x]
    { Expr* f = odd_fold(arg, "Cot"); if (f) return f; }

    // Cot[I y] -> -I Coth[y]
    { Expr* f = trig_i_fold(arg, "Coth", -1); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Cot); if (f) return f; }

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
        double complex s = mth_cot_machine_z(cplx);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Sec", arg); if (r) return r; }

    { Expr* inv = strip_inverse_call(arg, SYM_ArcSec); if (inv) return inv; }
    { Expr* f = try_simp_forward_of_inverse("Sec", arg); if (f) return f; }

    // Sec is even: Sec[-x] -> Sec[x]
    { Expr* f = even_fold(arg, "Sec"); if (f) return f; }

    // Sec[I y] -> Sech[y]
    { Expr* f = trig_i_fold(arg, "Sech", 0); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Sec); if (f) return f; }

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
        double complex s = mth_sec_machine_z(cplx);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("Csc", arg); if (r) return r; }

    { Expr* inv = strip_inverse_call(arg, SYM_ArcCsc); if (inv) return inv; }
    { Expr* f = try_simp_forward_of_inverse("Csc", arg); if (f) return f; }

    // Csc is odd: Csc[-x] -> -Csc[x]
    { Expr* f = odd_fold(arg, "Csc"); if (f) return f; }

    // Csc[I y] -> -I Csch[y]
    { Expr* f = trig_i_fold(arg, "Csch", -1); if (f) return f; }

    // f[theta + k Pi/2]: quarter-period shift identities
    { Expr* f = trig_quarter_shift(arg, SYM_Csc); if (f) return f; }

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
        double complex s = mth_csc_machine_z(cplx);
        if (cimag(cplx) == 0.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/* True when `arg` could be one of the exact values Sin/Cos/Tan/... take at
 * rational multiples of Pi that the inverse tables below compare it with:
 * those are built from integers by Plus, Times, Power and Rational (radicals
 * and rationals), or are ComplexInfinity at a pole.  Anything else -- an
 * approximate number, a symbol, a call such as ArcTan[1/10^8] -- cannot
 * match, and building the ~70 table values to find that out cost ~8 us of
 * every machine ArcSin[0.1] (60x the libm call), of every MPFR call N makes,
 * and of every unevaluated ArcSin[ArcTan[x]] re-evaluated. */
static bool exact_inverse_table_candidate(const Expr* e, int depth) {
    if (depth > 8) return false;            /* deeper than any table value */
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return true;
        case EXPR_SYMBOL: return e->data.symbol.name == SYM_ComplexInfinity;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL) return false;
            const char* hn = h->data.symbol.name;
            if (hn != SYM_Plus && hn != SYM_Times && hn != SYM_Power &&
                hn != SYM_Rational && hn != SYM_Complex) return false;
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!exact_inverse_table_candidate(e->data.function.args[i], depth + 1))
                    return false;
            return true;
        }
        default: return false;              /* Real, MPFR, String, ... */
    }
}

static bool exact_inverse_table_cannot_match(const Expr* arg) {
    /* The only bare rationals among the values: 0, +-1, +-2 (Sec, Csc at
     * Pi/3, Pi/6), +-1/2. */
    if (arg->type == EXPR_INTEGER)
        return arg->data.integer < -2 || arg->data.integer > 2;
    if (arg->type == EXPR_BIGINT) return true;
    int64_t n, d;
    if (is_rational((Expr*)arg, &n, &d)) return !((n == 1 || n == -1) && d == 2);
    return !exact_inverse_table_candidate(arg, 0);
}

/* The table values exact_sin / exact_cos produce are not all in evaluated
 * form (Sin[Pi/12] is built as a difference that evaluates to
 * (-Sqrt[2] + Sqrt[6])/4), and the argument of ArcSin / ArcCos is: the
 * comparison is against the EVALUATED value, computed once per (n, d) and
 * kept for the process (ArcSin[Sin[Pi/12]] = Pi/12, ArcSin#3). */
#define EXACT_INV_D_MAX 12
static Expr* g_exact_inv_cache[2][EXACT_INV_D_MAX + 1][2 * EXACT_INV_D_MAX + 1];
static bool  g_exact_inv_done[2][EXACT_INV_D_MAX + 1][2 * EXACT_INV_D_MAX + 1];

static const Expr* exact_inverse_value(int which, int64_t n, int64_t d) {
    size_t i = (size_t)d, j = (size_t)(n + EXACT_INV_D_MAX);
    if (!g_exact_inv_done[which][i][j]) {
        Expr* v = which == 0 ? exact_sin(n, d) : exact_cos(n, d);
        g_exact_inv_cache[which][i][j] = v ? eval_and_free(v) : NULL;
        g_exact_inv_done[which][i][j] = true;
    }
    return g_exact_inv_cache[which][i][j];
}

static Expr* exact_arcsin(Expr* arg) {
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = -d / 2; n <= d / 2; n++) {         /* the range [-Pi/2, Pi/2] */
            const Expr* val = exact_inverse_value(0, n, d);
            if (val && expr_eq(arg, (Expr*)val))
                return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
        }
    }
    return NULL;
}

static Expr* exact_arccos(Expr* arg) {
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n <= d; n++) {
            const Expr* val = exact_inverse_value(1, n, d);
            if (val && expr_eq(arg, (Expr*)val))
                return make_times(make_rational(n, d), expr_new_symbol(SYM_Pi));
        }
    }
    return NULL;
}

static Expr* exact_arctan(Expr* arg) {
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 8, 10, 12};
    for (int i = 0; i < 9; i++) {
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
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 8, 10, 12};
    for (int i = 0; i < 9; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n < d; n++) {
            Expr* val = exact_cot(n, d);
            if (val) val = eval_and_free(val);   /* canonical form, as arg is */
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
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = 0; n <= d; n++) {
            Expr* val = exact_sec(n, d);
            if (val) val = eval_and_free(val);   /* canonical form, as arg is */
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
    if (exact_inverse_table_cannot_match(arg)) return NULL;
    int64_t ds[] = {1, 2, 3, 4, 5, 6, 10, 12};
    for (int i = 0; i < 8; i++) {
        int64_t d = ds[i];
        for (int64_t n = -d; n <= d; n++) {
            Expr* val = exact_csc(n, d);
            if (val) val = eval_and_free(val);   /* canonical form, as arg is */
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
/* True for any flavour of infinite argument: Infinity, -Infinity,
 * ComplexInfinity, or DirectedInfinity[_]. The reciprocal inverse-trig
 * functions (ArcCot, ArcSec, ArcCsc) satisfy f[z] = g[1/z] and 1/z -> 0 for
 * every one of these, so f attains its value at 0 for any approach. */
static bool is_any_infinity_arg(Expr* e) {
    if (!e) return false;
    if (is_infinity_sym(e) || is_complex_infinity_sym(e) || is_neg_infinity_form(e))
        return true;
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_DirectedInfinity;
}

Expr* builtin_arcsin(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1)
        return builtin_arg_error("ArcSin", res->data.function.arg_count, 1, 1);
    Expr* arg = res->data.function.args[0];

    if (is_interval(arg)) { Expr* r = interval_apply_function("ArcSin", arg); if (r) return r; }

    // ArcSin is odd: ArcSin[-x] -> -ArcSin[x]
    { Expr* f = odd_fold(arg, "ArcSin"); if (f) return f; }

    // ArcSin grows without bound; a directionless infinite argument maps to
    // ComplexInfinity (magnitude infinite, direction undefined).
    if (is_complex_infinity_sym(arg)) return expr_new_symbol(SYM_ComplexInfinity);

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
        /* the one machine ArcSin (mth_machine_elem.h): asin on [-1, 1], the
         * lower side of the (1, inf) cut, shared with the kernel and the VM */
        double complex s = mth_asin_machine_z(c);
        if (cimag(c) == 0.0 && creal(c) >= -1.0 && creal(c) <= 1.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("ArcCos", arg); if (r) return r; }

    // ArcCos[-x] -> Pi - ArcCos[x]
    { Expr* f = arc_pi_minus_fold(arg, "ArcCos"); if (f) return f; }

    // ArcCos grows without bound; a directionless infinite argument maps to
    // ComplexInfinity.
    if (is_complex_infinity_sym(arg)) return expr_new_symbol(SYM_ComplexInfinity);

    // ArcCos[I y] -> Pi/2 - I ArcSinh[y]
    { Expr* f = arccos_i_fold(arg); if (f) return f; }

    // Attempt exact inverse evaluation
    Expr* exact = exact_arccos(arg);
    if (exact) return exact;

#ifdef USE_MPFR
    /* ArcCos[x] = I ArcCosh[x] for a real x > 1: the real part is an exact 0
     * (the complex kernel left an inexact one) */
    if (arg->type == EXPR_MPFR && mpfr_cmp_ui(arg->data.mpfr, 1) > 0) {
        Expr* v = expr_new_mpfr_bits(mpfr_get_prec(arg->data.mpfr));
        mpfr_acosh(v->data.mpfr, arg->data.mpfr, MPFR_RNDN);
        return make_complex(expr_new_integer(0), v);
    }
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
        /* the one machine ArcCos (mth_machine_elem.h): acos on [-1, 1],
         * +i acosh(x) for real x > 1, shared with the kernel and the VM */
        double complex s = mth_acos_machine_z(c);
        if (cimag(c) == 0.0 && creal(c) >= -1.0 && creal(c) <= 1.0) return expr_new_real(creal(s));
        return make_complex(expr_new_real(creal(s)), expr_new_real(cimag(s)));
    }
    return NULL;
}

/*
 * arctan_direction_sign:
 * Quadrant sign for ArcTan at DirectedInfinity[dir]. ArcTan[z] -> Pi/2 for
 * Re[z] > 0 and -Pi/2 for Re[z] < 0; on the imaginary axis (Re == 0) the
 * on-cut value is +Pi/2 for the upper half (Im > 0) and -Pi/2 for the lower
 * half (matching the direct evaluation of ArcTan at large imaginary
 * arguments). Returns +1 / -1, or 0 when the quadrant cannot be determined.
 */
static int arctan_direction_sign(Expr* dir) {
    /* Numerically approximate the direction and read off the quadrant. Going
     * through N resolves symbolic-but-numeric directions (I Sqrt[2], 1 + I,
     * ...) that expr_numeric_sign's literal-only view would miss. */
    Expr* nd = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
                   (Expr*[]){ expr_copy(dir) }, 1));
    double complex c;
    bool inexact = false;
    bool ok = get_approx(nd, &c, &inexact);
    expr_free(nd);
    if (!ok) return 0;
    const double tol = 1e-9;
    double re = creal(c), im = cimag(c);
    if (re >  tol) return +1;
    if (re < -tol) return -1;
    if (im >  tol) return +1;
    if (im < -tol) return -1;
    return 0;
}

/*
 * builtin_arctan:
 * Implements the evaluation logic for the inverse tangent (ArcTan).
 * Supports both ArcTan[z] (1 argument) and ArcTan[x, y] (2 arguments).
 * Reverses exact algebraic mappings for supported inputs.
 * Uses C standard math (catan and atan2) for numeric evaluations.
 */
/* The sign of an exact real numeric quantity for the two-argument ArcTan:
 * -1, 0 (a literal zero), +1, or 2 when it is not an exact real number or
 * its sign is not resolved by its machine value. */
static int atan2_exact_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer > 0 ? 1 : e->data.integer < 0 ? -1 : 0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    if (arg_is_inexact(e) || !expr_is_numeric_quantity(e)) return 2;
    double v;
    if (!numeric_machine_value(e, &v)) {
        Expr* nv = numericalize_machine_unplanned(e);
        bool ok = nv && nv->type == EXPR_REAL;
        if (ok) v = nv->data.real;
        if (nv) expr_free(nv);
        if (!ok) return 2;
    }
    if (!isfinite(v) || v == 0.0) return 2;
    return v > 0 ? 1 : -1;
}

static bool atan2_mentions_arctan(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_ArcTan)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (atan2_mentions_arctan(e->data.function.args[i])) return true;
    return false;
}

Expr* builtin_arctan(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;

    // ArcTan accepts 1 (arctangent) or 2 (two-argument / quadrant) arguments.
    if (res->data.function.arg_count != 1 && res->data.function.arg_count != 2)
        return builtin_arg_error("ArcTan", res->data.function.arg_count, 1, 2);

    // Single argument ArcTan[z]
    if (res->data.function.arg_count == 1) {
        Expr* arg = res->data.function.args[0];

        if (is_interval(arg)) { Expr* r = interval_apply_function("ArcTan", arg); if (r) return r; }

        // Infinite arguments. Handle these BEFORE the odd / exact folds:
        // exact_arctan otherwise spuriously matches Tan[-Pi/2] ==
        // ComplexInfinity (n = -1 is tried before n = +1) and returns -Pi/2
        // for ArcTan[ComplexInfinity]. The genuine values are:
        //   ArcTan[ComplexInfinity]      -> Indeterminate (direction unknown)
        //   ArcTan[DirectedInfinity[d]]  -> +/-Pi/2 by the quadrant of d.
        if (arg->type == EXPR_SYMBOL &&
            arg->data.symbol.name == SYM_ComplexInfinity) {
            return expr_new_symbol(SYM_Indeterminate);
        }
        if (arg->type == EXPR_FUNCTION &&
            arg->data.function.head->type == EXPR_SYMBOL &&
            arg->data.function.head->data.symbol.name == SYM_DirectedInfinity &&
            arg->data.function.arg_count == 1) {
            int s = arctan_direction_sign(arg->data.function.args[0]);
            if (s > 0) return make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi));
            if (s < 0) return make_times(make_rational(-1, 2), expr_new_symbol(SYM_Pi));
            // s == 0: quadrant undetermined; leave ArcTan[...] unevaluated.
        }

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
            double complex s = mth_atan_machine_z(c);   /* atan for real c (G4) */
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
            /* The angle of the origin is undefined (audit ArcTan#3). */
            if (xv == 0 && yv == 0) {
                mth_message("ArcTan", "indet", "Indeterminate expression ArcTan[0, 0] encountered.");
                return expr_new_symbol(SYM_Indeterminate);
            }
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
        
        /* Exact real x, y beyond machine integers -- rationals, bignums,
         * surds, Pi (audit ArcTan#1): the quadrant of (x, y) and the exact
         * ArcTan[y/x] when it has a closed form.  ArcTan[1/2, 1/2] = Pi/4,
         * ArcTan[-1/2, 0] = Pi, ArcTan[Pi, 0] = 0, ArcTan[1, -Sqrt[3]] =
         * -Pi/3.  Without a closed form for ArcTan[y/x] the call is left. */
        {
            int sx = atan2_exact_sign(x), sy = atan2_exact_sign(y);
            if (sx != 2 && sy != 2 && !(sx == 0 && sy == 0)) {
                if (sy == 0) return sx > 0 ? expr_new_integer(0) : expr_new_symbol(SYM_Pi);
                if (sx == 0)
                    return make_times(make_rational(sy > 0 ? 1 : -1, 2), expr_new_symbol(SYM_Pi));
                Expr* q = eval_and_free(make_times(expr_copy(y),
                              make_power(expr_copy(x), expr_new_integer(-1))));
                Expr* t = eval_and_free(expr_new_function(expr_new_symbol(SYM_ArcTan), &q, 1));
                if (t && !atan2_mentions_arctan(t)) {
                    if (sx > 0) return t;
                    return eval_and_free(make_plus(t, make_times(expr_new_integer(sy > 0 ? 1 : -1),
                                                                 expr_new_symbol(SYM_Pi))));
                }
                if (t) expr_free(t);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("ArcCot", arg); if (r) return r; }

    // ArcCot is odd (Mathematica convention): ArcCot[-x] -> -ArcCot[x]
    { Expr* f = odd_fold(arg, "ArcCot"); if (f) return f; }

    // ArcCot[I y] -> -I ArcCoth[y]
    { Expr* f = trig_i_fold(arg, "ArcCoth", -1); if (f) return f; }

    // ArcCot[z] = ArcTan[1/z] -> ArcTan[0] = 0 for every infinite z.
    if (is_any_infinity_arg(arg)) return expr_new_integer(0);

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
        /* the one machine ArcCot (mth_machine_elem.h): atan(1/x) for real x,
         * so ArcCot[0.] is Pi/2 and ArcCot[-0.] is -Pi/2 (IEEE signed zero on
         * every surface, numeric-tower review G7) */
        double complex s = mth_acot_machine_z(c);
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("ArcSec", arg); if (r) return r; }

    // ArcSec[-x] -> Pi - ArcSec[x]
    { Expr* f = arc_pi_minus_fold(arg, "ArcSec"); if (f) return f; }

    // ArcSec[z] = ArcCos[1/z] -> ArcCos[0] = Pi/2 for every infinite z.
    if (is_any_infinity_arg(arg))
        return make_times(make_rational(1, 2), expr_new_symbol(SYM_Pi));
    // ArcSec[0]: 1/z runs off to infinity in every direction -- a pole of
    // undetermined direction, as the machine ArcSec[0.] says (audit ArcSec#1).
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_symbol(SYM_ComplexInfinity);

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
        double im = cimag(s);
        /* The same branch-cut correction ArcCos carries: for real 0 < x < 1
         * the reciprocal lands on the (1, inf) cut, where C99 cacos returns
         * -i*acosh and Mathematica returns +i*acosh. Without this the machine
         * path contradicted both the MPFR path and ArcCos[1/x]. */
        if (cimag(c) == 0.0 && creal(c) != 0.0 && fabs(creal(c)) < 1.0) im = -im;
        return make_complex(expr_new_real(creal(s)), expr_new_real(im));
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

    if (is_interval(arg)) { Expr* r = interval_apply_function("ArcCsc", arg); if (r) return r; }

    // ArcCsc is odd: ArcCsc[-x] -> -ArcCsc[x]
    { Expr* f = odd_fold(arg, "ArcCsc"); if (f) return f; }

    // ArcCsc[I y] -> -I ArcCsch[y]
    { Expr* f = trig_i_fold(arg, "ArcCsch", -1); if (f) return f; }

    // ArcCsc[z] = ArcSin[1/z] -> ArcSin[0] = 0 for every infinite z.
    if (is_any_infinity_arg(arg)) return expr_new_integer(0);
    // ArcCsc[0] = ComplexInfinity, as ArcCsc[0.] (audit ArcCsc#1).
    if (arg->type == EXPR_INTEGER && arg->data.integer == 0)
        return expr_new_symbol(SYM_ComplexInfinity);

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
        double im = cimag(s);
        /* See ArcSec above: for real 0 < x < 1 the reciprocal sits on the
         * (1, inf) cut, where C99 and Mathematica pick opposite sides. */
        if (cimag(c) == 0.0 && creal(c) != 0.0 && fabs(creal(c)) < 1.0) im = -im;
        return make_complex(expr_new_real(creal(s)), expr_new_real(im));
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
