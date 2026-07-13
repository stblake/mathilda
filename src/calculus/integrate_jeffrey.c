/* integrate_jeffrey.c
 *
 * Jeffrey-Rich continuous Weierstrass-substitution integrator.  See
 * integrate_jeffrey.h for the high-level description and the reference.
 *
 * ---------------------------------------------------------------------------
 * Algorithm (paper section 5)
 * ---------------------------------------------------------------------------
 * Given a rational function f of the trig kernels {Sin,Cos,Tan,Cot,Sec,Csc}[x]
 * (trig mode) or the hyperbolic kernels {Sinh,Cosh,Tanh,Coth,Sech,Csch}[x]
 * (hyperbolic mode), with every occurrence of x inside such a kernel:
 *
 *   1. Substitute (Table I, choice (a) for trig, u = tan(x/2); the analogous
 *      u = tanh(x/2) for hyperbolic): each kernel becomes a rational function
 *      of a fresh symbol u; multiply by the Jacobian dx.  Cancel[Together[.]]
 *      yields a rational function g(u) free of x.
 *   2. Integrate g in u (recurses, closing through Integrate`BronsteinRational).
 *      Call the result ghat_u(u).  For trig, if choice (a) does not close, retry
 *      with choice (c), u = cot(x/2), b = 0.
 *   3. (Trig only.)  Compute the jump introduced by the pole of tan(x/2):
 *          K = lim_{u->+inf} ghat_u - lim_{u->-inf} ghat_u
 *      (negated for choice (c), whose half-angle runs the other way).  The core
 *      Limit engine does not distribute over Plus/Times at infinity, so the jump
 *      is summed term by term here: constant factors are pulled out and the
 *      limit of each ArcTan/Log core is resolved through the limit of its inner
 *      argument.  If any term diverges (a genuine singularity of the integrand),
 *      K is dropped and no correction is applied -- paper section 4.
 *   4. ghat(x) = ghat_u with u -> Tan[x/2] (a) / Cot[x/2] (c) / Tanh[x/2] (hyp).
 *      The continuous antiderivative is
 *          g(x) = ghat(x) + K * Floor[(x - b)/p],   b = Pi (a) / 0 (c), p = 2 Pi.
 *      For hyperbolic integrands no correction is needed (tanh(x/2) is a monotone
 *      bijection with no poles), so g(x) = ghat(x).
 *
 * The substitution is an exact identity and the rational sub-integral is closed
 * by a verified integrator, and Floor' = 0 almost everywhere, so the result is
 * correct by construction; no differentiate-back verification is performed (it
 * is unnecessary, and the Floor term defeats symbolic differentiation anyway --
 * cf. integrate_linrad.c).
 *
 * Memory: builtins take ownership of `res`; helpers own every Expr they build
 * and free intermediates explicitly.  `eval_take` consumes its argument.  We
 * never expr_free(res).  The small builder/eval helpers mirror the private copy
 * kept by the other integration-method modules.
 */

#include "integrate_jeffrey.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "poly.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers                                    */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(name), args, 2);
}

/* Unevaluated Plus[a, b] / Times[a, b] / Power[a, b]. */
static Expr* mk_plus2(Expr* a, Expr* b)  { return mk_fn2(SYM_Plus,  a, b); }
static Expr* mk_times2(Expr* a, Expr* b) { return mk_fn2(SYM_Times, a, b); }
static Expr* mk_pow(Expr* a, Expr* b)    { return mk_fn2(SYM_Power, a, b); }
/* a / b  ==  Times[a, Power[b, -1]]. */
static Expr* mk_div(Expr* a, Expr* b)    { return mk_times2(a, mk_pow(b, mk_int(-1))); }

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* True if `f` contains no subexpression structurally equal to `x`. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq((Expr*)f, (Expr*)x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* True if any symbol with interned name `s` occurs anywhere in `e`. */
static bool contains_sym(const Expr* e, const char* s) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == s;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_sym(e->data.function.head, s)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_sym(e->data.function.args[i], s)) return true;
    return false;
}

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (contains_unintegrated(e->data.function.args[i])) return true;
    return false;
}

/* Cancel[Together[e]], consuming `e`. */
static Expr* cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* ReplaceAll[expr, rules]; consumes `expr` and `rules`; returns owned. */
static Expr* replace_all_take(Expr* expr, Expr* rules) {
    return eval_take(internal_replace_all((Expr*[]){ expr, rules }, 2));
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2(SYM_Integrate, expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define WJ_MAX_DEPTH 6
static int wj_depth = 0;
static unsigned long wj_sym_counter = 0;

/* Mode / substitution-choice selectors. */
enum { WMODE_TRIG, WMODE_HYP };
enum { WCHOICE_A, WCHOICE_C, WCHOICE_TANH };

/* ---------------------------------------------------------------------- */
/* Kernel detection                                                       */
/* ---------------------------------------------------------------------- */

static bool is_trig_head(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Tan ||
           h == SYM_Cot || h == SYM_Sec || h == SYM_Csc;
}
static bool is_hyp_head(const char* h) {
    return h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh ||
           h == SYM_Coth || h == SYM_Sech || h == SYM_Csch;
}

/* Walk `e` classifying how the integration variable x appears.  A trig/hyp
 * kernel head applied to exactly x shields its argument (not recursed into).
 * Any bare x reached outside such a kernel sets *bad. */
static void wj_scan(const Expr* e, const char* xsym,
                    bool* has_trig, bool* has_hyp, bool* bad) {
    if (*bad || !e) return;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == xsym) *bad = true;  /* x outside a kernel */
        return;
    }
    if (e->type != EXPR_FUNCTION) return;

    const Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* h = head->data.symbol.name;
        const Expr* arg = e->data.function.args[0];
        bool arg_is_x = (arg->type == EXPR_SYMBOL && arg->data.symbol.name == xsym);
        if (arg_is_x && is_trig_head(h)) { *has_trig = true; return; }
        if (arg_is_x && is_hyp_head(h))  { *has_hyp  = true; return; }
    }
    wj_scan(head, xsym, has_trig, has_hyp, bad);
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        wj_scan(e->data.function.args[i], xsym, has_trig, has_hyp, bad);
}

/* Automatic-cascade gate: true if some Power[base, neg] with base depending on
 * x occurs in f, i.e. a kernel sits in a denominator.  Polynomial trig such as
 * Integrate[Sin[x], x] fails this test and is left to the cleaner table/Risch
 * forms; only genuine rational integrands (1/(a + cos x), ...) take this route. */
static bool wj_has_kernel_in_denominator(const Expr* e, const char* xsym) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        bool neg = false;
        if (exp->type == EXPR_INTEGER) neg = exp->data.integer < 0;
        else if (is_rational(exp, &p, &q)) neg = (p < 0) ^ (q < 0);
        if (neg && contains_sym(base, xsym)) return true;
    }
    if (wj_has_kernel_in_denominator(e->data.function.head, xsym)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (wj_has_kernel_in_denominator(e->data.function.args[i], xsym)) return true;
    return false;
}

/* ---------------------------------------------------------------------- */
/* Substitution table                                                     */
/* ---------------------------------------------------------------------- */

/* Reusable rational building blocks in u. */
static Expr* B_2u(Expr* u)    { return mk_times2(mk_int(2), expr_copy(u)); }            /* 2u     */
static Expr* B_1pu2(Expr* u)  { return mk_plus2(mk_int(1), mk_pow(expr_copy(u), mk_int(2))); }                  /* 1+u^2  */
static Expr* B_1mu2(Expr* u)  { return mk_plus2(mk_int(1), mk_times2(mk_int(-1), mk_pow(expr_copy(u), mk_int(2)))); } /* 1-u^2 */
static Expr* B_u2m1(Expr* u)  { return mk_plus2(mk_pow(expr_copy(u), mk_int(2)), mk_int(-1)); }                 /* u^2-1  */

/* Rule  head[x] -> rhs.  Consumes rhs. */
static Expr* wj_rule(const char* head, Expr* x, Expr* rhs) {
    return mk_fn2(SYM_Rule, mk_fn1(head, expr_copy(x)), rhs);
}

/* Build the six kernel substitution rules and the Jacobian for (mode, choice).
 * Fills rules[0..5]; returns the Jacobian dx expression (owned). */
static Expr* wj_build_rules(Expr* x, Expr* u, int mode, int choice, Expr* rules[6]) {
    if (mode == WMODE_TRIG && choice == WCHOICE_A) {
        rules[0] = wj_rule(SYM_Sin, x, mk_div(B_2u(u),   B_1pu2(u)));
        rules[1] = wj_rule(SYM_Cos, x, mk_div(B_1mu2(u), B_1pu2(u)));
        rules[2] = wj_rule(SYM_Tan, x, mk_div(B_2u(u),   B_1mu2(u)));
        rules[3] = wj_rule(SYM_Cot, x, mk_div(B_1mu2(u), B_2u(u)));
        rules[4] = wj_rule(SYM_Sec, x, mk_div(B_1pu2(u), B_1mu2(u)));
        rules[5] = wj_rule(SYM_Csc, x, mk_div(B_1pu2(u), B_2u(u)));
        return mk_div(mk_int(2), B_1pu2(u));                 /* dx = 2/(1+u^2) */
    }
    if (mode == WMODE_TRIG && choice == WCHOICE_C) {
        rules[0] = wj_rule(SYM_Sin, x, mk_div(B_2u(u),   B_1pu2(u)));
        rules[1] = wj_rule(SYM_Cos, x, mk_div(B_u2m1(u), B_1pu2(u)));
        rules[2] = wj_rule(SYM_Tan, x, mk_div(B_2u(u),   B_u2m1(u)));
        rules[3] = wj_rule(SYM_Cot, x, mk_div(B_u2m1(u), B_2u(u)));
        rules[4] = wj_rule(SYM_Sec, x, mk_div(B_1pu2(u), B_u2m1(u)));
        rules[5] = wj_rule(SYM_Csc, x, mk_div(B_1pu2(u), B_2u(u)));
        return mk_div(mk_int(-2), B_1pu2(u));                /* dx = -2/(1+u^2) */
    }
    /* WMODE_HYP, u = tanh(x/2). */
    rules[0] = wj_rule(SYM_Sinh, x, mk_div(B_2u(u),   B_1mu2(u)));
    rules[1] = wj_rule(SYM_Cosh, x, mk_div(B_1pu2(u), B_1mu2(u)));
    rules[2] = wj_rule(SYM_Tanh, x, mk_div(B_2u(u),   B_1pu2(u)));
    rules[3] = wj_rule(SYM_Coth, x, mk_div(B_1pu2(u), B_2u(u)));
    rules[4] = wj_rule(SYM_Sech, x, mk_div(B_1mu2(u), B_1pu2(u)));
    rules[5] = wj_rule(SYM_Csch, x, mk_div(B_1mu2(u), B_2u(u)));
    return mk_div(mk_int(2), B_1mu2(u));                     /* dx = 2/(1-u^2) */
}

/* The half-angle kernel u = phi(x) used for back-substitution. */
static const char* wj_half_head(int mode, int choice) {
    if (mode == WMODE_HYP) return SYM_Tanh;
    return (choice == WCHOICE_A) ? SYM_Tan : SYM_Cot;
}

/* Apply the substitution: f(sin x, ...) dx -> g(u), Cancel[Together] reduced.
 * Returns owned g(u), or NULL.  Borrows f, x, u. */
static Expr* wj_substitute(Expr* f, Expr* x, Expr* u, int mode, int choice) {
    Expr* rules[6];
    Expr* jac = wj_build_rules(x, u, mode, choice, rules);
    Expr* rulelist = expr_new_function(expr_new_symbol(SYM_List), rules, 6);
    Expr* sub = replace_all_take(expr_copy(f), rulelist);
    if (!sub) { expr_free(jac); return NULL; }
    Expr* integrand = mk_times2(sub, jac);
    return cancel_together(integrand);
}

/* True if e is a rational function of u (Numerator and Denominator both
 * polynomial in u). */
static bool wj_is_rational_in(Expr* e, Expr* u) {
    Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(e)));
    Expr* den = eval_take(mk_fn1("Denominator", expr_copy(e)));
    Expr* vars[1] = { u };
    bool ok = num && den && is_polynomial(num, vars, 1) && is_polynomial(den, vars, 1);
    if (num) expr_free(num);
    if (den) expr_free(den);
    return ok;
}

/* ---------------------------------------------------------------------- */
/* Jump (K) computation via limits at +/- infinity                        */
/* ---------------------------------------------------------------------- */

typedef enum { LIM_FINITE, LIM_POSINF, LIM_NEGINF, LIM_BAIL } LimKind;

/* Numeric sign of a real constant c: -1, 0, +1, or 2 (cannot decide). */
static int const_sign(Expr* c) {
    Expr* n = eval_take(mk_fn1("N", expr_copy(c)));
    int r = 2;
    if (n) {
        if (n->type == EXPR_INTEGER)
            r = n->data.integer > 0 ? 1 : (n->data.integer < 0 ? -1 : 0);
        else if (n->type == EXPR_REAL)
            r = n->data.real > 0.0 ? 1 : (n->data.real < 0.0 ? -1 : 0);
        expr_free(n);
    }
    return r;
}

/* Classify an evaluated limit result.  On LIM_FINITE, *fin receives an owned
 * copy of the finite value; otherwise *fin is left NULL. */
static LimKind classify_inf(Expr* e, Expr** fin) {
    *fin = NULL;
    if (!e) return LIM_BAIL;

    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        if (s == SYM_Infinity) return LIM_POSINF;
        if (s == SYM_ComplexInfinity || s == SYM_Indeterminate) return LIM_BAIL;
        *fin = expr_copy(e);
        return LIM_FINITE;
    }
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT) {
        *fin = expr_copy(e);
        return LIM_FINITE;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if (h == SYM_DirectedInfinity && e->data.function.arg_count == 1) {
            int s = const_sign(e->data.function.args[0]);
            return s > 0 ? LIM_POSINF : (s < 0 ? LIM_NEGINF : LIM_BAIL);
        }
        if (h == SYM_Times) {
            /* k * Infinity  (e.g. -Infinity = Times[-1, Infinity]). */
            int infc = 0;
            bool dirty = false;
            Expr* dir = mk_int(1);
            for (size_t i = 0; i < e->data.function.arg_count; i++) {
                Expr* a = e->data.function.args[i];
                if (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_Infinity) {
                    infc++;
                } else if (contains_sym(a, SYM_Infinity) ||
                           contains_sym(a, SYM_DirectedInfinity) ||
                           contains_sym(a, SYM_ComplexInfinity)) {
                    dirty = true;
                } else {
                    dir = mk_times2(dir, expr_copy(a));
                }
            }
            if (infc == 1 && !dirty) {
                int s = const_sign(dir);
                expr_free(dir);
                return s > 0 ? LIM_POSINF : (s < 0 ? LIM_NEGINF : LIM_BAIL);
            }
            expr_free(dir);
            return LIM_BAIL;
        }
        if (head_is(e, "Limit")) return LIM_BAIL;  /* unresolved */
    }
    /* Any residual infinity head means the limit did not resolve cleanly. */
    if (contains_sym(e, SYM_Infinity) || contains_sym(e, SYM_DirectedInfinity) ||
        contains_sym(e, SYM_ComplexInfinity) || contains_sym(e, SYM_Indeterminate))
        return LIM_BAIL;

    *fin = expr_copy(e);
    return LIM_FINITE;
}

/* Limit[e, u -> sgn*Infinity].  Returns owned evaluated result. */
static Expr* limit_at(Expr* e, Expr* u, int sgn) {
    Expr* pt = (sgn > 0) ? expr_new_symbol(SYM_Infinity)
                         : mk_times2(mk_int(-1), expr_new_symbol(SYM_Infinity));
    Expr* rule = mk_fn2(SYM_Rule, expr_copy(u), pt);
    return eval_take(mk_fn2("Limit", expr_copy(e), rule));
}

/* Limit of a single core (no u-free factors) at sgn*infinity.  On success sets
 * *out to the owned finite limit and returns true; false means the term
 * diverges or could not be resolved. */
static bool core_limit(Expr* core, Expr* u, int sgn, Expr** out) {
    *out = NULL;

    /* ArcTan[arg]: arctan of the inner-argument limit. */
    if (head_is(core, SYM_ArcTan) && core->data.function.arg_count == 1) {
        Expr* il = limit_at(core->data.function.args[0], u, sgn);
        Expr* fin = NULL;
        LimKind k = classify_inf(il, &fin);
        if (il) expr_free(il);
        if (k == LIM_POSINF) {
            *out = mk_div(expr_new_symbol(SYM_Pi), mk_int(2));        /* +Pi/2 */
            return true;
        }
        if (k == LIM_NEGINF) {
            *out = mk_times2(mk_int(-1), mk_div(expr_new_symbol(SYM_Pi), mk_int(2)));
            return true;
        }
        if (k == LIM_FINITE) {
            *out = eval_take(mk_fn1(SYM_ArcTan, fin));               /* ArcTan[L] */
            return true;
        }
        return false;
    }

    /* Log[arg]: log of the inner-argument limit, provided it is finite. */
    if (head_is(core, SYM_Log) && core->data.function.arg_count == 1) {
        Expr* il = limit_at(core->data.function.args[0], u, sgn);
        Expr* fin = NULL;
        LimKind k = classify_inf(il, &fin);
        if (il) expr_free(il);
        if (k != LIM_FINITE) { if (fin) expr_free(fin); return false; }
        Expr* ll = eval_take(mk_fn1(SYM_Log, fin));                  /* Log[L] */
        Expr* fin2 = NULL;
        LimKind k2 = classify_inf(ll, &fin2);
        if (ll) expr_free(ll);
        if (k2 == LIM_FINITE) { *out = fin2; return true; }
        if (fin2) expr_free(fin2);
        return false;                                               /* Log diverged */
    }

    /* Otherwise: rational / algebraic core -- must converge to a finite value. */
    Expr* l = limit_at(core, u, sgn);
    Expr* fin = NULL;
    LimKind k = classify_inf(l, &fin);
    if (l) expr_free(l);
    if (k == LIM_FINITE) { *out = fin; return true; }
    if (fin) expr_free(fin);
    return false;
}

/* Limit of one additive term at sgn*infinity.  Splits off u-free factors as a
 * constant multiplier, then takes the limit of the u-core.  On success sets
 * *out (owned) and returns true. */
static bool term_limit(Expr* term, Expr* u, int sgn, Expr** out) {
    *out = NULL;

    if (expr_free_of(term, u)) {           /* constant additive term */
        *out = expr_copy(term);
        return true;
    }

    /* Partition factors into a u-free constant product C and a u-core. */
    Expr* C = mk_int(1);
    Expr* core = NULL;
    bool is_times = head_is(term, SYM_Times);
    size_t nf = is_times ? term->data.function.arg_count : 1;
    for (size_t i = 0; i < nf; i++) {
        Expr* fac = is_times ? term->data.function.args[i] : term;
        if (expr_free_of(fac, u)) C = mk_times2(C, expr_copy(fac));
        else core = core ? mk_times2(core, expr_copy(fac)) : expr_copy(fac);
    }
    if (!core) { expr_free(C); return false; }  /* shouldn't happen */

    Expr* cl = NULL;
    bool ok = core_limit(core, u, sgn, &cl);
    expr_free(core);
    if (!ok) { expr_free(C); return false; }

    *out = eval_take(mk_times2(C, cl));
    return true;
}

/* Total limit of g at sgn*infinity, summed over additive terms.  On success
 * sets *out (owned) and returns true; false if any term diverges. */
static bool total_limit(Expr* g, Expr* u, int sgn, Expr** out) {
    *out = NULL;
    Expr* acc = NULL;
    bool g_is_plus = head_is(g, SYM_Plus);
    size_t n = g_is_plus ? g->data.function.arg_count : 1;
    for (size_t i = 0; i < n; i++) {
        Expr* term = g_is_plus ? g->data.function.args[i] : g;
        Expr* c = NULL;
        if (!term_limit(term, u, sgn, &c)) {
            if (acc) expr_free(acc);
            return false;
        }
        acc = acc ? eval_take(mk_plus2(acc, c)) : c;
    }
    *out = acc ? acc : mk_int(0);
    return true;
}

/* The jump K = lim_{+inf} g - lim_{-inf} g (negated for choice (c)), Simplified.
 * Returns NULL if either one-sided limit fails to resolve to a finite value. */
static Expr* wj_jump(Expr* g, Expr* u, int choice) {
    Expr* lp = NULL;
    Expr* ln = NULL;
    if (!total_limit(g, u, +1, &lp)) return NULL;
    if (!total_limit(g, u, -1, &ln)) { expr_free(lp); return NULL; }
    Expr* k = eval_take(mk_plus2(lp, mk_times2(mk_int(-1), ln)));
    if (choice == WCHOICE_C) k = eval_take(mk_times2(mk_int(-1), k));
    return eval_take(mk_fn1("Simplify", k));
}

static bool is_zero_expr(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer == 0) ||
           (e->type == EXPR_REAL && e->data.real == 0.0);
}

/* ---------------------------------------------------------------------- */
/* One substitution attempt                                               */
/* ---------------------------------------------------------------------- */

static Expr* wj_attempt(Expr* f, Expr* x, int mode, int choice) {
    char uname[96];
    snprintf(uname, sizeof(uname),
             "Integrate`Weierstrass`u$%lu", wj_sym_counter++);
    Expr* u = expr_new_symbol(uname);

    Expr* result = NULL;

    /* 1. Substitute to a rational function of u. */
    Expr* g = wj_substitute(f, x, u, mode, choice);
    if (!g) goto done;
    if (!expr_free_of(g, x) || !wj_is_rational_in(g, u)) { expr_free(g); goto done; }

    /* 2. Integrate in u. */
    wj_depth++;
    Expr* G = integrate_in(g, u);
    wj_depth--;
    expr_free(g);
    if (!G) goto done;

    /* 3. Back-substitute u -> phi(x). */
    Expr* half = mk_fn1(wj_half_head(mode, choice), mk_div(expr_copy(x), mk_int(2)));
    Expr* rule = mk_fn2(SYM_Rule, expr_copy(u), half);
    Expr* ghat_x = replace_all_take(expr_copy(G), rule);

    if (mode == WMODE_HYP) {
        /* tanh(x/2) is a monotone bijection: no spurious jump, no correction. */
        result = ghat_x;
        expr_free(G);
        goto done;
    }

    /* 4. Continuity correction K * Floor[(x - b)/p]. */
    Expr* K = wj_jump(G, u, choice);
    expr_free(G);
    if (!K || is_zero_expr(K)) {
        if (K) expr_free(K);
        result = ghat_x;                       /* already continuous */
        goto done;
    }
    {
        Expr* b = (choice == WCHOICE_A) ? expr_new_symbol(SYM_Pi) : mk_int(0);
        Expr* p = mk_times2(mk_int(2), expr_new_symbol(SYM_Pi));
        Expr* arg = mk_div(mk_plus2(expr_copy(x), mk_times2(mk_int(-1), b)), p);
        Expr* corr = mk_times2(K, mk_fn1(SYM_Floor, arg));
        result = eval_take(mk_plus2(ghat_x, corr));
    }

done:
    expr_free(u);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

/* Detect the mode, apply the Automatic-cascade gate, and run the substitution
 * attempt(s) on the (already-normalised) integrand `f`. */
static Expr* wj_run(Expr* f, Expr* x, const char* xsym, bool explicit_call) {
    bool has_trig = false, has_hyp = false, bad = false;
    wj_scan(f, xsym, &has_trig, &has_hyp, &bad);
    if (bad) return NULL;                 /* x appears outside a kernel */
    if (has_trig == has_hyp) return NULL; /* none, or mixed trig + hyperbolic */

    /* Automatic cascade: only intercept genuine rational integrands (kernel in a
     * denominator); leave polynomial trig to the cleaner table / Risch forms. */
    if (!explicit_call && !wj_has_kernel_in_denominator(f, xsym)) return NULL;

    if (has_hyp)
        return wj_attempt(f, x, WMODE_HYP, WCHOICE_TANH);

    Expr* r = wj_attempt(f, x, WMODE_TRIG, WCHOICE_A);
    if (!r) r = wj_attempt(f, x, WMODE_TRIG, WCHOICE_C);
    return r;
}

static Expr* wj_core(Expr* f, Expr* x, bool explicit_call) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (wj_depth >= WJ_MAX_DEPTH) return NULL;
    if (expr_free_of(f, x)) return NULL;

    const char* xsym = x->data.symbol.name;

    /* First try the integrand verbatim (the common, fast path).  If that does
     * not match -- typically because a kernel argument is a multiple or sum of x
     * (Cos[2 x], Sin[x + 1], Cosh[x] Cosh[2 x], ...) -- retry on TrigExpand[f],
     * which rewrites those into kernels of the bare variable x.  TrigExpand
     * leaves single-angle kernels untouched, so the verbatim path is unaffected. */
    Expr* result = wj_run(f, x, xsym, explicit_call);
    if (result) return result;

    Expr* fe = eval_take(mk_fn1("TrigExpand", expr_copy(f)));
    if (fe) {
        if (!expr_eq(fe, f)) result = wj_run(fe, x, xsym, explicit_call);
        expr_free(fe);
    }
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_jeffrey_try(Expr* f, Expr* x) {
    return wj_core(f, x, false);
}

Expr* integrate_jeffrey_full(Expr* f, Expr* x) {
    return wj_core(f, x, true);
}

Expr* builtin_integrate_jeffrey(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return wj_core(f, x, true);
}

void integrate_jeffrey_init(void) {
    symtab_add_builtin("Integrate`Weierstrass", builtin_integrate_jeffrey);
    symtab_get_def("Integrate`Weierstrass")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`Weierstrass",
        "Integrate`Weierstrass[f, x] integrates a rational function of the\n"
        "trigonometric kernels Sin/Cos/Tan/Cot/Sec/Csc[x] (or the hyperbolic\n"
        "kernels Sinh/Cosh/Tanh/Coth/Sech/Csch[x]) by the continuous Weierstrass\n"
        "substitution u = Tan[x/2] (Tanh[x/2] for hyperbolic). It integrates the\n"
        "resulting rational function of u, back-substitutes, and -- for the\n"
        "trigonometric case -- adds the Jeffrey-Rich secular correction\n"
        "K Floor[(x - b)/p] that removes the spurious discontinuities introduced\n"
        "at the poles of Tan[x/2]. Strict: returns unevaluated when f is not such\n"
        "a rational function or the reduced rational integral does not close.");
}
