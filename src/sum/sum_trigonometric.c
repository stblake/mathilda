/*
 * sum_trigonometric.c -- Sum`Trigonometric: infinite Fourier-type sums.
 *
 * Handles  Sum_{k=imin}^Infinity  T(k) / k^s , where T(k) is a trigonometric
 * polynomial in k (products/powers of Sin/Cos of arguments linear in k) and s
 * is a positive integer.  The general identity is, for integer s >= 1,
 *
 *     Sum_{k>=1} Sin[a k]/k^s = Im PolyLog[s, e^{i a}],
 *     Sum_{k>=1} Cos[a k]/k^s = Re PolyLog[s, e^{i a}].
 *
 * The trigonometric part is first linearised with TrigReduce, turning any
 * product/power of sines and cosines into a sum of single-angle terms
 * c_j {Sin|Cos}[a_j k + phi_j].  Each term maps to the (Im|Re) of a PolyLog,
 * and constant ("DC") terms c_j contribute c_j Zeta[s] (for s >= 2).
 *
 * Elementary collapse (Wolfram-Language parity): at s == 1 with no phase and a
 * numeric coefficient 0 < a < 2*Pi,
 *
 *     Sum_{k>=1} Sin[a k]/k = (Pi - a)/2,
 *     Sum_{k>=1} Cos[a k]/k = -Log[2 Sin[a/2]].
 *
 * so e.g. Sum[Sin[k]/k, {k,1,Infinity}] -> (Pi - 1)/2.  For s >= 2 (or a phase,
 * or a symbolic / out-of-range coefficient) the result is left as the
 * corresponding Im/Re of a PolyLog -- exactly as Wolfram leaves it -- never a
 * fabricated value.  Divergent inputs (a DC term at s == 1) stay unevaluated.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* --- tiny expression builders (all adopt their Expr* arguments) --- */

static Expr* tr_fn2(const char* sym, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(sym), (Expr*[]){ a, b }, 2);
}
static Expr* tr_fn1(const char* sym, Expr* a) {
    return expr_new_function(expr_new_symbol(sym), (Expr*[]){ a }, 1);
}

/* The imaginary unit I as a Complex literal. */
static Expr* tr_I(void) {
    Expr* c = tr_fn2(SYM_Complex, sum_int(0), sum_int(1));
    Expr* r = evaluate(c);
    expr_free(c);
    return r;
}

/* e^{i x}, adopting x. */
static Expr* tr_expi(Expr* x) {
    return tr_fn2(SYM_Power, expr_new_symbol(SYM_E), tr_fn2(SYM_Times, tr_I(), x));
}

/* Is g a var-dependent Sin[.]/Cos[.]?  Sets *head (1=Sin,2=Cos) and *arg (alias). */
static bool tr_is_trig(Expr* g, Expr* var, int* head, Expr** arg) {
    if (g->type != EXPR_FUNCTION || g->data.function.head->type != EXPR_SYMBOL
        || g->data.function.arg_count != 1) return false;
    const char* h = g->data.function.head->data.symbol;
    int hh = (h == SYM_Sin) ? 1 : (h == SYM_Cos) ? 2 : 0;
    if (!hh) return false;
    if (sum_free_of(g->data.function.args[0], var)) return false;  /* Sin[const] is a constant */
    *head = hh;
    *arg = g->data.function.args[0];
    return true;
}

/* True iff e contains a var-dependent Sin/Cos anywhere (cheap pre-gate). */
static bool tr_has_trig(Expr* e, Expr* var) {
    int hh; Expr* a;
    if (tr_is_trig(e, var, &hh, &a)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (tr_has_trig(e->data.function.head, var)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (tr_has_trig(e->data.function.args[i], var)) return true;
    }
    return false;
}

/* Classify one additive TrigReduce term as c * {1 | Sin|Cos}[a var + phi].
 * On success returns 1 and sets *c (owned coeff), *th (0 DC, 1 Sin, 2 Cos),
 * and (when th != 0) *a and *phi (owned).  Returns 0 if the term is not of the
 * required const * single-trig-of-linear-arg shape. */
static int tr_classify(Expr* term, Expr* var, Expr** c_out, int* th,
                       Expr** a_out, Expr** phi_out) {
    bool is_times = (term->type == EXPR_FUNCTION
                     && term->data.function.head->type == EXPR_SYMBOL
                     && term->data.function.head->data.symbol == SYM_Times);
    size_t n = is_times ? term->data.function.arg_count : 1;

    Expr* c = sum_int(1);
    Expr* trigarg = NULL;
    int head = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* g = is_times ? term->data.function.args[i] : term;
        int hh; Expr* ta;
        if (tr_is_trig(g, var, &hh, &ta)) {
            if (head) { expr_free(c); return 0; }   /* two trig factors: not linearised */
            head = hh; trigarg = ta;
        } else if (sum_free_of(g, var)) {
            Expr* t = tr_fn2(SYM_Times, c, expr_copy(g));
            c = evaluate(t); expr_free(t);
        } else {
            expr_free(c); return 0;                 /* var-dependent non-trig factor */
        }
    }

    if (!head) { *th = 0; *c_out = c; return 1; }    /* DC term */

    Expr* a = sum_eval("D", (Expr*[]){ expr_copy(trigarg), expr_copy(var) }, 2);
    if (!sum_free_of(a, var)) { expr_free(c); expr_free(a); return 0; }  /* nonlinear arg */
    Expr* zero = sum_int(0);
    Expr* phi = sum_subst(trigarg, var, zero);       /* arg(0) */
    expr_free(zero);
    *th = head; *c_out = c; *a_out = a; *phi_out = phi;
    return 1;
}

/* c * Sum_{k>=1} {Sin|Cos}[a k + phi]/k^s, owned.  Consumes c; copies a, phi. */
static Expr* tr_contribution(Expr* c, int th, Expr* a, Expr* phi, int s) {
    bool phi_zero = (phi->type == EXPR_INTEGER && phi->data.integer == 0);

    /* Elementary closed form: s == 1, no phase, numeric 0 < a < 2 Pi. */
    if (s == 1 && phi_zero) {
        Expr* na = sum_eval("N", (Expr*[]){ expr_copy(a) }, 1);
        bool numeric = (na->type == EXPR_REAL);
        double av = numeric ? na->data.real : 0.0;
        expr_free(na);
        if (numeric && av > 0.0 && av < 2.0 * M_PI) {
            Expr* core;
            if (th == 1) {   /* (Pi - a)/2 */
                core = tr_fn2(SYM_Times, make_rational(1, 2),
                           tr_fn2(SYM_Plus, expr_new_symbol(SYM_Pi),
                               tr_fn2(SYM_Times, sum_int(-1), expr_copy(a))));
            } else {          /* -Log[2 Sin[a/2]] */
                Expr* half_a = tr_fn2(SYM_Times, make_rational(1, 2), expr_copy(a));
                core = tr_fn2(SYM_Times, sum_int(-1),
                           tr_fn1(SYM_Log,
                               tr_fn2(SYM_Times, sum_int(2), tr_fn1(SYM_Sin, half_a))));
            }
            Expr* t = tr_fn2(SYM_Times, c, core);
            Expr* r = evaluate(t); expr_free(t);
            return r;
        }
    }

    /* General form: c * (Im|Re)[ e^{i phi} PolyLog[s, e^{i a}] ]. */
    Expr* pl = expr_new_function(expr_new_symbol(SYM_PolyLog),
                   (Expr*[]){ sum_int(s), tr_expi(expr_copy(a)) }, 2);
    Expr* inner = phi_zero ? pl : tr_fn2(SYM_Times, tr_expi(expr_copy(phi)), pl);
    Expr* reim = tr_fn1(th == 1 ? SYM_Im : SYM_Re, inner);
    Expr* t = tr_fn2(SYM_Times, c, reim);
    Expr* r = evaluate(t); expr_free(t);
    return r;
}

Expr* builtin_sum_trigonometric(Expr* res);

Expr* builtin_sum_trigonometric(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    /* k^{-s} is singular at k <= 0; require the natural lower bound k >= 1. */
    if (imin->type != EXPR_INTEGER || imin->data.integer < 1) return NULL;
    int64_t k0 = imin->data.integer;

    Expr* fn = evaluate(expr_copy(f));
    if (!tr_has_trig(fn, var)) { expr_free(fn); return NULL; }

    /* Split fn = (product of k-powers) * trigpart. */
    bool is_times = (fn->type == EXPR_FUNCTION
                     && fn->data.function.head->type == EXPR_SYMBOL
                     && fn->data.function.head->data.symbol == SYM_Times);
    size_t n = is_times ? fn->data.function.arg_count : 1;

    int64_t powsum = 0;              /* sum of integer exponents of var */
    Expr* trigpart = sum_int(1);
    for (size_t i = 0; i < n; i++) {
        Expr* g = is_times ? fn->data.function.args[i] : fn;
        if (g->type == EXPR_SYMBOL && g == var) {
            powsum += 1;
        } else if (g->type == EXPR_FUNCTION && g->data.function.head->type == EXPR_SYMBOL
                   && g->data.function.head->data.symbol == SYM_Power
                   && g->data.function.arg_count == 2
                   && expr_eq(g->data.function.args[0], var)
                   && g->data.function.args[1]->type == EXPR_INTEGER) {
            powsum += g->data.function.args[1]->data.integer;
        } else {
            Expr* t = tr_fn2(SYM_Times, trigpart, expr_copy(g));
            trigpart = evaluate(t); expr_free(t);
        }
    }
    expr_free(fn);

    int64_t s = -powsum;             /* summand ~ trigpart / k^s */
    if (s < 1 || s > 1000000) { expr_free(trigpart); return NULL; }

    /* Linearise the trigonometric part into single-angle terms.  TrigReduce
     * returns e.g. 1/2 (1 - Cos[2 k]) as Times[1/2, Plus[...]]; Expand pushes it
     * to a flat sum  c_1 + c_2 Cos[2 k] + ...  so the term splitter can see each. */
    Expr* trr = sum_eval("TrigReduce", (Expr*[]){ trigpart }, 1);
    Expr* tr = sum_eval("Expand", (Expr*[]){ trr }, 1);
    size_t nt = (tr->type == EXPR_FUNCTION && tr->data.function.head->type == EXPR_SYMBOL
                 && tr->data.function.head->data.symbol == SYM_Plus)
                ? tr->data.function.arg_count : 1;

    Expr** contribs = malloc(sizeof(Expr*) * nt);
    size_t nc = 0;
    bool ok = true, saw_trig = false;
    for (size_t i = 0; i < nt && ok; i++) {
        Expr* term = (nt > 1) ? tr->data.function.args[i] : tr;
        Expr *c = NULL, *a = NULL, *phi = NULL; int th = 0;
        if (!tr_classify(term, var, &c, &th, &a, &phi)) { ok = false; break; }
        if (th == 0) {
            /* DC term c/k^s: c Zeta[s] for s >= 2, divergent for s == 1. */
            if (s == 1) { expr_free(c); ok = false; break; }
            Expr* zt = tr_fn2(SYM_Times, c,
                           expr_new_function(expr_new_symbol(SYM_Zeta),
                               (Expr*[]){ sum_int((int64_t)s) }, 1));
            Expr* ev = evaluate(zt); expr_free(zt);
            contribs[nc++] = ev;
        } else {
            saw_trig = true;
            Expr* con = tr_contribution(c, th, a, phi, (int)s);
            expr_free(a); expr_free(phi);
            contribs[nc++] = con;
        }
    }
    expr_free(tr);

    if (!ok || !saw_trig) {          /* not a trig sum (or divergent): fall through */
        for (size_t i = 0; i < nc; i++) expr_free(contribs[i]);
        free(contribs);
        return NULL;
    }

    Expr* summed = (nc == 1) ? contribs[0]
                             : expr_new_function(expr_new_symbol(SYM_Plus), contribs, nc);
    free(contribs);

    /* Shift the lower limit: the closed form is for k >= 1; subtract the head. */
    if (k0 > 1) {
        Expr* head_spec = expr_new_function(expr_new_symbol(SYM_List),
                              (Expr*[]){ expr_copy(var), sum_int(1), sum_int(k0 - 1) }, 3);
        Expr* head_sum = sum_eval("Sum", (Expr*[]){ expr_copy(f), head_spec }, 2);
        Expr* diff = tr_fn2(SYM_Plus, summed,
                         tr_fn2(SYM_Times, sum_int(-1), head_sum));
        summed = evaluate(diff); expr_free(diff);
    }

    Expr* ev = evaluate(summed); expr_free(summed);
    return sum_eval("Simplify", (Expr*[]){ ev }, 1);
}

void sum_trigonometric_init(void) {
    symtab_add_builtin("Sum`Trigonometric", builtin_sum_trigonometric);
    symtab_get_def("Sum`Trigonometric")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Trigonometric",
        "Sum`Trigonometric[f, k, imin, Infinity] gives the closed form of an "
        "infinite Fourier-type sum T(k)/k^s with T a trigonometric polynomial in "
        "k, via TrigReduce linearisation and Im/Re of PolyLog[s, e^{i a}]. At s = 1 "
        "results collapse to elementary form ((Pi - a)/2, -Log[2 Sin[a/2]]). "
        "Returns unevaluated for non-trigonometric or divergent inputs.");
}
