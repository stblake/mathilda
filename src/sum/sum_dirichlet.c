/*
 * sum_dirichlet.c -- Sum`Dirichlet: Dirichlet series of arithmetic functions.
 *
 * Closes  Sum[g(k)/k^s, {k, 1, Infinity}]  when g is a recognised multiplicative
 * arithmetic function, mapping to the standard Dirichlet generating function in
 * terms of the Riemann zeta function:
 *
 *     MoebiusMu[k]            -> 1/Zeta[s]
 *     MoebiusMu[k]^2 / |mu|   -> Zeta[s]/Zeta[2 s]
 *     LiouvilleLambda[k]      -> Zeta[2 s]/Zeta[s]
 *     EulerPhi[k]             -> Zeta[s-1]/Zeta[s]
 *     DivisorSigma[0,k] (d)   -> Zeta[s]^2
 *     DivisorSigma[0,k]^2     -> Zeta[s]^4/Zeta[2 s]
 *     DivisorSigma[a,k]       -> Zeta[s] Zeta[s-a]
 *
 * The summand is normalised, its k^{-s} factor and arithmetic-function core are
 * split out, the table entry is built as a Zeta expression, and a final
 * evaluate() folds even zetas (Zeta[6] -> Pi^6/945) while odd zetas stay
 * symbolic.  s may be a symbol or a number; it must be free of k.
 *
 * Convergence.  Each identity holds only where the series converges.  For a
 * numeric s the series is required to converge absolutely, Re s > sigma_a:
 * sigma_a = 1 for mu, |mu|, mu^2, lambda, d and d^2; 2 for phi; 1 + max(0, Re a)
 * for sigma_a.  Outside that half-plane the zeta ratio is the analytic
 * continuation, not the sum (EulerPhi[k]/k^(3/2) diverges, while
 * Zeta[1/2]/Zeta[3/2] is a finite negative number), so the stage declines.  The
 * one boundary point kept is s = 1 for mu and lambda, where the series
 * converges conditionally to 0 (the prime number theorem), which 1/Zeta[1] and
 * Zeta[2]/Zeta[1] evaluate to.  A symbolic s keeps the closed form, valid on
 * the convergence region.
 *
 * Placed before Sum`Rational in the cascade: these summands carry a
 * non-polynomial factor, so Sum`Rational declines them anyway, but Dirichlet
 * recognises them directly and avoids a wasted Apart.
 *
 * Memory contract: takes ownership of res but must not free it; returns an
 * owned closed form or NULL to fall through.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include "numeric.h"      /* numericalize: the convergence test on a numeric s */
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* --- small builders (all adopt their Expr* arguments) --- */

static Expr* d_zeta(Expr* arg) {
    return expr_new_function(expr_new_symbol(SYM_Zeta), (Expr*[]){ arg }, 1);
}
static Expr* d_pow(Expr* base, int e) {
    return expr_new_function(expr_new_symbol(SYM_Power),
               (Expr*[]){ base, sum_int(e) }, 2);
}
static Expr* d_times2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ a, b }, 2);
}
/* 2*s */
static Expr* d_two_s(Expr* s) {
    return d_times2(sum_int(2), expr_copy(s));
}
/* s - n */
static Expr* d_s_minus(Expr* s, Expr* n) {
    Expr* neg = d_times2(sum_int(-1), expr_copy(n));
    return expr_new_function(expr_new_symbol(SYM_Plus),
               (Expr*[]){ expr_copy(s), neg }, 2);
}

/* Is e the call HEAD[var]? (single argument equal to var) */
static bool d_is_call1(Expr* e, const char* head_sym, Expr* var) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != head_sym) return false;
    if (e->data.function.arg_count != 1) return false;
    return expr_eq(e->data.function.args[0], var);
}
/* Is e EulerPhi[var]? (EulerPhi has no SYM_ constant -> match by name) */
static bool d_is_eulerphi(Expr* e, Expr* var) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_EulerPhi) return false;
    if (e->data.function.arg_count != 1) return false;
    return expr_eq(e->data.function.args[0], var);
}
/* Is e DivisorSigma[a, var]?  On success *a_out aliases the order argument. */
static bool d_is_divisorsigma(Expr* e, Expr* var, Expr** a_out) {
    if (e->type != EXPR_FUNCTION) return false;
    Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_DivisorSigma) return false;
    if (e->data.function.arg_count != 2) return false;
    if (!expr_eq(e->data.function.args[1], var)) return false;
    *a_out = e->data.function.args[0];
    return true;
}

/* Build the zeta-ratio closed form for the recognised (g_base, g_pow, abs).
 * Returns an owned, still-unevaluated tree, or NULL if unrecognised.
 * s and a are copied, never consumed. */
static Expr* d_closed_form(Expr* g_base, int g_pow, bool abs_flag,
                           Expr* var, Expr* s) {
    /* MoebiusMu */
    if (d_is_call1(g_base, SYM_MoebiusMu, var)) {
        if (abs_flag || g_pow == 2) {
            /* |mu| or mu^2:  Zeta[s]/Zeta[2 s]. */
            return d_times2(d_zeta(expr_copy(s)), d_pow(d_zeta(d_two_s(s)), -1));
        }
        if (g_pow == 1) {
            /* 1/Zeta[s]. */
            return d_pow(d_zeta(expr_copy(s)), -1);
        }
        return NULL;
    }
    /* LiouvilleLambda:  Zeta[2 s]/Zeta[s]. */
    if (g_pow == 1 && !abs_flag && d_is_call1(g_base, SYM_LiouvilleLambda, var)) {
        return d_times2(d_zeta(d_two_s(s)), d_pow(d_zeta(expr_copy(s)), -1));
    }
    /* EulerPhi:  Zeta[s-1]/Zeta[s]. */
    if (g_pow == 1 && !abs_flag && d_is_eulerphi(g_base, var)) {
        Expr* one = sum_int(1);
        Expr* zsm1 = d_zeta(d_s_minus(s, one));
        expr_free(one);
        return d_times2(zsm1, d_pow(d_zeta(expr_copy(s)), -1));
    }
    /* DivisorSigma */
    Expr* a = NULL;
    if (!abs_flag && d_is_divisorsigma(g_base, var, &a)) {
        bool a_is_zero = (a->type == EXPR_INTEGER && a->data.integer == 0);
        if (a_is_zero && g_pow == 1) {
            /* d(k):  Zeta[s]^2. */
            return d_pow(d_zeta(expr_copy(s)), 2);
        }
        if (a_is_zero && g_pow == 2) {
            /* d(k)^2:  Zeta[s]^4/Zeta[2 s]. */
            return d_times2(d_pow(d_zeta(expr_copy(s)), 4),
                            d_pow(d_zeta(d_two_s(s)), -1));
        }
        if (g_pow == 1) {
            /* sigma_a(k):  Zeta[s] Zeta[s-a]. */
            if (!sum_free_of(a, var)) return NULL;
            return d_times2(d_zeta(expr_copy(s)), d_zeta(d_s_minus(s, a)));
        }
        return NULL;
    }
    return NULL;
}

/* Real part of a numeric expression as a double; false when not numeric. */
static bool d_real_part(const Expr* e, double* out) {
    Expr* re = sum_eval("Re", (Expr*[]){ expr_copy((Expr*)e) }, 1);
    Expr* num = re ? numericalize(re, numeric_machine_spec()) : NULL;
    expr_free(re);
    bool ok = false;
    if (num && num->type == EXPR_REAL) { *out = num->data.real; ok = isfinite(*out); }
    else if (num && num->type == EXPR_INTEGER) { *out = (double)num->data.integer; ok = true; }
    expr_free(num);
    return ok;
}

/* Does the recognised series converge at the exponent s?  True for a symbolic
 * s (the closed form then states the value on the convergence region); for a
 * numeric s, absolute convergence Re s > sigma_a, plus the conditionally
 * convergent boundary s = 1 of mu and lambda (see the header). */
static bool d_converges(Expr* g_base, int g_pow, bool abs_flag, Expr* var, Expr* s) {
    double rs;
    if (!d_real_part(s, &rs)) return true;           /* symbolic s */
    double sigma = 1.0;
    bool boundary_ok = false;
    Expr* a = NULL;
    if (d_is_call1(g_base, SYM_MoebiusMu, var)) {
        boundary_ok = (g_pow == 1 && !abs_flag);
    } else if (d_is_call1(g_base, SYM_LiouvilleLambda, var)) {
        boundary_ok = (g_pow == 1);
    } else if (d_is_eulerphi(g_base, var)) {
        sigma = 2.0;
    } else if (d_is_divisorsigma(g_base, var, &a)) {
        double ra;
        if (!d_real_part(a, &ra)) return false;      /* order symbolic, s numeric */
        if (g_pow == 1 && ra > 0.0) sigma = 1.0 + ra;
    }
    if (rs > sigma) return true;
    if (boundary_ok && rs == 1.0) {
        Expr* im = sum_eval("Im", (Expr*[]){ expr_copy(s) }, 1);
        bool zero_im = im && im->type == EXPR_INTEGER && im->data.integer == 0;
        expr_free(im);
        return zero_im;
    }
    return false;
}

Expr* builtin_sum_dirichlet(Expr* res);

Expr* builtin_sum_dirichlet(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    if (!(imin->type == EXPR_INTEGER && imin->data.integer == 1)) return NULL;

    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    bool is_times = (fn->type == EXPR_FUNCTION
                     && fn->data.function.head->type == EXPR_SYMBOL
                     && fn->data.function.head->data.symbol.name == SYM_Times);
    size_t n = is_times ? fn->data.function.arg_count : 1;

    Expr* s = NULL;          /* the exponent s in k^{-s} (owned) */
    Expr* g_base = NULL;     /* arithmetic-function core (alias into fn) */
    int g_pow = 0;
    bool abs_flag = false;
    Expr* coeff = NULL;      /* accumulated k-free constant (owned) or NULL=1 */
    bool ok = true;

    for (size_t i = 0; i < n && ok; i++) {
        Expr* x = is_times ? fn->data.function.args[i] : fn;

        /* Power[base, exp] */
        if (x->type == EXPR_FUNCTION
            && x->data.function.head->type == EXPR_SYMBOL
            && x->data.function.head->data.symbol.name == SYM_Power
            && x->data.function.arg_count == 2) {
            Expr* base = x->data.function.args[0];
            Expr* exp  = x->data.function.args[1];
            if (expr_eq(base, var)) {
                /* k^exp  ->  s = -exp. */
                if (s) { ok = false; break; }
                s = sum_eval("Times", (Expr*[]){ sum_int(-1), expr_copy(exp) }, 2);
                continue;
            }
            /* g(k)^m with m a positive integer. */
            if (exp->type == EXPR_INTEGER && exp->data.integer > 0
                && !sum_free_of(base, var)) {
                if (g_base) { ok = false; break; }
                g_base = base;
                g_pow = (int)exp->data.integer;
                continue;
            }
            ok = false; break;
        }

        /* Abs[MoebiusMu[var]] */
        if (x->type == EXPR_FUNCTION
            && x->data.function.head->type == EXPR_SYMBOL
            && x->data.function.head->data.symbol.name == SYM_Abs
            && x->data.function.arg_count == 1
            && d_is_call1(x->data.function.args[0], SYM_MoebiusMu, var)) {
            if (g_base) { ok = false; break; }
            g_base = x->data.function.args[0];
            g_pow = 1;
            abs_flag = true;
            continue;
        }

        /* Bare arithmetic-function factor g(k). */
        if (!sum_free_of(x, var)) {
            if (g_base) { ok = false; break; }
            g_base = x;
            g_pow = 1;
            continue;
        }

        /* k-free constant coefficient. */
        if (coeff) {
            Expr* nc = sum_eval("Times", (Expr*[]){ coeff, expr_copy(x) }, 2);
            coeff = nc;
        } else {
            coeff = expr_copy(x);
        }
    }

    if (!ok || !g_base || !s || !sum_free_of(s, var)) {
        expr_free(fn); if (s) expr_free(s); if (coeff) expr_free(coeff);
        return NULL;
    }

    Expr* form = d_converges(g_base, g_pow, abs_flag, var, s)
                 ? d_closed_form(g_base, g_pow, abs_flag, var, s) : NULL;
    expr_free(s);
    expr_free(fn);          /* g_base aliased into fn; safe now that form is built */
    if (!form) { if (coeff) expr_free(coeff); return NULL; }

    if (coeff) form = d_times2(coeff, form);   /* adopts both */

    Expr* out = evaluate(form);
    expr_free(form);
    return out;
}

void sum_dirichlet_init(void) {
    symtab_add_builtin("Sum`Dirichlet", builtin_sum_dirichlet);
    symtab_get_def("Sum`Dirichlet")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Sum`Dirichlet",
        "Sum`Dirichlet[f, i, 1, Infinity] evaluates Dirichlet series Sum[g[i]/i^s] "
        "for arithmetic functions g (MoebiusMu, LiouvilleLambda, EulerPhi, "
        "DivisorSigma) as ratios of Riemann zeta values.");
}
