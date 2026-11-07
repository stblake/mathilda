/* Mathilda -- Surd[x, n], the real-valued n-th root of a real x.
 *
 *   Surd[x, n] = Sign[x] Abs[x]^(1/n)   for odd n,
 *              = x^(1/n)                for even n and x >= 0,
 *
 * unlike x^(1/n), which is the principal (complex) root for x < 0. Each kind of
 * argument takes the cheapest route:
 *
 *   machine x           -> sf_surd_r (sqrt / cbrt / pow + one Newton step), the
 *                          kernel the NDArray, packed and Compile surfaces share
 *   arbitrary precision -> Power on |x| at the input precision, sign restored
 *   exact x >= 0        -> Power[x, 1/n]          (Surd[8, 3] = 2)
 *   exact x < 0, odd n  -> -Power[-x, 1/n]        (Surd[-2, 3] = -2^(1/3))
 *   symbolic x          -> stays unevaluated
 *
 * Messages, as in the reference:
 *   Surd::int    n is not an integer          (call stays unevaluated)
 *   Surd::noneg  even root of a negative x    (Indeterminate)
 *   Surd::indet  n == 0                       (Indeterminate)
 *   Surd::preal  x is a non-real number       (call stays unevaluated)
 *
 * Deliberate difference: a machine real n with an integer value (Surd[8., 3.])
 * is accepted as that integer. The reference rejects it with Surd::int, but the
 * packed / NDArray kernel receives every operand as a double and cannot tell 3
 * from 3., and the scalar path must agree with it.
 *
 * CubeRoot[x] is Surd[x, 3] under its own name (the reference rewrites a
 * symbolic CubeRoot[x] to Surd[x, 3]): machine x takes the same sf_surd_r
 * kernel, registered as CubeRoot's own unary NDArray / packed / Compile
 * kernel; everything else is Surd[x, 3]; a non-real number is CubeRoot::preal.
 *
 * Attributes: Listable, NumericFunction, Protected (both heads).
 */
#include "surd.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "numeric.h"
#include "print.h"
#include "symtab.h"
#include "sf_scalar.h"

/* ------------------------------------------------------------------ */
/* Machine kernel                                                     */
/* ------------------------------------------------------------------ */

bool sf_surd_r(double x, double n, double* out) {
    if (!isfinite(x) || !isfinite(n) || n != floor(n) || n == 0.0 || fabs(n) > 1e15)
        return false;
    int64_t k = (int64_t)fabs(n);
    bool odd = (k & 1) != 0;
    if (x < 0.0 && !odd) return false;          /* Surd::noneg: the List path */
    double a = fabs(x), r;
    if (a == 0.0) {
        if (n < 0.0) return false;               /* ComplexInfinity: the List path */
        *out = x;                                /* keeps the sign of -0. */
        return true;
    }
    if (k == 1)      r = a;
    else if (k == 2) r = sqrt(a);
    else if (k == 3) r = cbrt(a);
    else {
        r = pow(a, 1.0 / (double)k);
        /* One Newton step on r^k = a recovers the ulp pow loses to the
         * rounded exponent 1/k. */
        double rk1 = pow(r, (double)(k - 1));
        double corr = (rk1 * r - a) / ((double)k * rk1);
        if (isfinite(corr)) r -= corr;
    }
    if (n < 0.0) r = 1.0 / r;
    if (!isfinite(r)) return false;
    *out = (x < 0.0) ? -r : r;
    return true;
}

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

static Expr* surd_indeterminate(const char* tag, const char* fmt, Expr* call) {
    if (call) {
        char* s = expr_to_string(call);
        mth_message("Surd", tag, fmt, s ? s : "");
        free(s);
    } else {
        mth_message("Surd", tag, "%s", fmt);
    }
    return expr_new_symbol(SYM_Indeterminate);
}

static Expr* surd_power(Expr* x, Expr* expo) {
    Expr* a[2] = { expr_copy(x), expo };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), a, 2));
}

static Expr* surd_negate(Expr* x) {
    Expr* a[2] = { expr_new_integer(-1), x };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), a, 2));
}

/* -1, 0 or +1 for a real numeric x; 2 for a non-real number; 3 undecided. */
static int surd_sign(Expr* x) {
    if (x->type == EXPR_INTEGER) return (x->data.integer > 0) - (x->data.integer < 0);
    if (x->type == EXPR_BIGINT) return mpz_sgn(x->data.bigint);
    if (x->type == EXPR_REAL) return (x->data.real > 0) - (x->data.real < 0);
    if (!expr_is_numeric_quantity(x)) {
        if (x->type == EXPR_SYMBOL && x->data.symbol.name == SYM_Infinity) return 1;
        return 3;
    }
    Expr *re, *im;
    if (is_complex(x, &re, &im)) return 2;
    Expr* v = numericalize(x, numeric_machine_spec());
    int s = 3;
    if (v) {
        if (v->type == EXPR_REAL) {
            if (v->data.real > 0) s = 1;
            else if (v->data.real < 0) s = -1;
        } else if (is_complex(v, &re, &im)) {
            s = 2;
        }
        expr_free(v);
    }
    if (s == 3) {
        /* Too close to zero for a machine double: ask the comparisons. */
        Expr* a[1] = { expr_copy(x) };
        Expr* p = eval_and_free(expr_new_function(expr_new_symbol(SYM_Positive), a, 1));
        if (p->type == EXPR_SYMBOL && p->data.symbol.name == SYM_True) s = 1;
        expr_free(p);
        if (s == 3) {
            Expr* b[1] = { expr_copy(x) };
            Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_Negative), b, 1));
            if (q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True) s = -1;
            expr_free(q);
        }
    }
    return s;
}

static bool surd_is_neg_infinity(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_Times &&
           e->data.function.args[0]->type == EXPR_INTEGER &&
           e->data.function.args[0]->data.integer == -1 &&
           e->data.function.args[1]->type == EXPR_SYMBOL &&
           e->data.function.args[1]->data.symbol.name == SYM_Infinity;
}

/* ------------------------------------------------------------------ */
/* Builtin                                                            */
/* ------------------------------------------------------------------ */

Expr* builtin_surd(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("Surd", argc, 2, 2);
    Expr* x = res->data.function.args[0];
    Expr* n = res->data.function.args[1];

    /* The index: an exact integer, or a machine real with an integer value. */
    int64_t k;
    bool n_inexact = false;
    if (n->type == EXPR_INTEGER) {
        k = n->data.integer;
    } else if (n->type == EXPR_REAL && isfinite(n->data.real) &&
               n->data.real == floor(n->data.real) && fabs(n->data.real) <= 1e15) {
        k = (int64_t)n->data.real;
        n_inexact = true;
    } else if (n->type == EXPR_BIGINT) {
        k = mpz_sgn(n->data.bigint) > 0 ? INT64_MAX : INT64_MIN;  /* parity unused below */
    } else {
        if (n->type == EXPR_REAL ||
#ifdef USE_MPFR
            n->type == EXPR_MPFR ||
#endif
            (n->type == EXPR_FUNCTION && expr_is_numeric_quantity(n))) {
            char* s = expr_to_string(res);
            mth_message("Surd", "int", "Integer expected at position 2 in %s.", s ? s : "");
            free(s);
        }
        return NULL;
    }
    if (k == 0)
        return surd_indeterminate("indet", "Indeterminate expression %s encountered.", res);

    /* Machine x (or an inexact index): the shared kernel. */
    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && !ma.any_complex) {
        double v;
        if (sf_surd_r(ma.re[0], ma.re[1], &v)) return expr_new_real(v);
    }

    if (x->type == EXPR_SYMBOL && x->data.symbol.name == SYM_Infinity)
        return expr_copy(x);
    if (surd_is_neg_infinity(x)) {
        if (n->type == EXPR_BIGINT ? mpz_odd_p(n->data.bigint) : (k & 1))
            return expr_copy(x);
        return surd_indeterminate("noneg", "Surd is not defined for even roots of negative values.", NULL);
    }

    int s = surd_sign(x);
    if (s == 2) {
        char* t = expr_to_string(x);
        mth_message("Surd", "preal", "The parameter %s should be real-valued.", t ? t : "");
        free(t);
        return NULL;
    }
    if (s == 3) return NULL;
    bool odd = n->type == EXPR_BIGINT ? mpz_odd_p(n->data.bigint) != 0 : (k & 1) != 0;
    if (s < 0 && !odd)
        return surd_indeterminate("noneg", "Surd is not defined for even roots of negative values.", NULL);

    /* 1/n, exact; an inexact index makes the result inexact. */
    Expr* expo;
    if (n->type == EXPR_BIGINT) {
        Expr* a[2] = { expr_new_integer(1), expr_copy(n) };
        expo = eval_and_free(expr_new_function(expr_new_symbol(SYM_Rational), a, 2));
    } else if (n_inexact) {
        expo = expr_new_real(1.0 / (double)k);
    } else {
        expo = make_rational(1, k);
    }
    Expr* r;
    if (s >= 0) {
        r = surd_power(x, expo);
    } else {
        Expr* mx = surd_negate(expr_copy(x));
        r = surd_negate(surd_power(mx, expo));
        expr_free(mx);
    }
    if (n_inexact && !arg_is_inexact(r)) {
        Expr* v = numericalize(r, numeric_machine_spec());
        expr_free(r);
        r = v;
    }
    return r;
}

/* CubeRoot[x]: the real cube root. Machine reals answer through the shared
 * kernel (bit-identical to the packed / NDArray / Compile surfaces); a
 * non-real number is CubeRoot::preal and stays unevaluated; anything else
 * becomes Surd[x, 3] (exact roots, -2^(1/3), symbolic Surd[x, 3]). */
bool sf_cuberoot_r(double x, double* out) { return sf_surd_r(x, 3.0, out); }

Expr* builtin_cuberoot(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("CubeRoot", argc, 1, 1);
    Expr* x = res->data.function.args[0];

    SFMachineArgs ma;
    if (sf_machine_scalar_args(res, &ma) && !ma.any_complex) {
        double v;
        if (sf_cuberoot_r(ma.re[0], &v)) return expr_new_real(v);
    }
    Expr *re, *im;
    if (is_complex(x, &re, &im) || (expr_is_numeric_quantity(x) && surd_sign(x) == 2)) {
        char* t = expr_to_string(x);
        mth_message("CubeRoot", "preal", "The parameter %s should be real-valued.", t ? t : "");
        free(t);
        return NULL;
    }
    Expr* a[2] = { expr_copy(x), expr_new_integer(3) };
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Surd), a, 2));
}

void surd_init(void) {
    symtab_add_builtin("CubeRoot", builtin_cuberoot);
    symtab_get_def("CubeRoot")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_add_builtin("Surd", builtin_surd);
    symtab_get_def("Surd")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
