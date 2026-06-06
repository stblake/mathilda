/* integrate_linrad.c
 *
 * Integration of rational functions of x and radicals of a common linear
 * argument -- the "linear radicals" substitution.  See integrate_linrad.h for
 * the high-level description.
 *
 * ---------------------------------------------------------------------
 * Algorithm
 * ---------------------------------------------------------------------
 * Identify the integrand as
 *
 *     R[x, (a x + b)^(m1/n1), (a x + b)^(m2/n2), ...]
 *
 * with R rational, a, b free of x, and a single shared linear base (a x + b):
 *
 *   1. Scan f for every Power[base, p/q] with q > 1.  Each base must be a
 *      degree-1 polynomial in x; all such bases must be structurally equal.
 *      n = LCM of the q's.  Reject otherwise.
 *   2. Substitute the radicals to a fresh symbol u via
 *      poly_subst_radical_to_gen (shared with the algebraic-factoring path):
 *      (a x + b)^(p/q) -> u^(p n / q), and the bare base (a x + b) -> u^n.
 *      The result exy = R(x, u) must be a rational function of {x, u}.
 *   3. Form the rationalised integrand
 *          (n/a) u^(n-1) exy  /.  x -> (u^n - b)/a,
 *      Cancel[Together[.]] it, and recurse: Integrate[., u].
 *   4. Back-substitute u -> (a x + b)^(1/n) and accept iff the result
 *      differentiates back to f (unconditional verification gate).
 *
 * Memory: builtins take ownership of `res`; helpers below own every Expr they
 * construct and free intermediates explicitly.  `eval_take` consumes its
 * argument.  We never expr_free(res).  Several small helpers mirror
 * integrate_derivdivides.c (the project keeps a private copy of this helper set
 * per integration method file).
 */

#include "integrate_linrad.h"

#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "internal.h"
#include "sym_intern.h"
#include "sym_names.h"
#include "arithmetic.h"
#include "poly.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------- */
/* Small builders / evaluation helpers (mirrors integrate_derivdivides.c) */
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

/* Evaluate `call` to a fixed point, freeing `call`. */
static Expr* eval_take(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static bool is_lit_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
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

/* True if `e` contains any unevaluated Integrate[...] call. */
static bool contains_unintegrated(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (head_is(e, SYM_Integrate)) return true;
    if (contains_unintegrated(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_unintegrated(e->data.function.args[i])) return true;
    }
    return false;
}

/* D[expr, x]; borrows, returns owned (evaluated). */
static Expr* deriv_dx(const Expr* expr, const Expr* x) {
    return eval_take(mk_fn2("D", expr_copy((Expr*)expr), expr_copy((Expr*)x)));
}

/* Cancel[Together[e]], consuming `e`. */
static Expr* cancel_together(Expr* e) {
    if (!e) return NULL;
    Expr* t = eval_take(internal_together((Expr*[]){ e }, 1));
    if (!t) return NULL;
    return eval_take(internal_cancel((Expr*[]){ t }, 1));
}

/* ReplaceAll[expr, from -> to]; consumes `expr`; returns owned (evaluated). */
static Expr* replace_one(Expr* expr, const Expr* from, const Expr* to) {
    Expr* rule = mk_fn2("Rule", expr_copy((Expr*)from), expr_copy((Expr*)to));
    Expr* call = internal_replace_all((Expr*[]){ expr, rule }, 2);
    return eval_take(call);
}

/* The differentiation-verification gate: true iff D[r, x] - f is zero.
 * PossibleZeroQ pre-screen, then Simplify, then Cancel[Together[Expand[.]]].
 * Borrows r and f. */
static bool differentiates_back(const Expr* r, const Expr* f, const Expr* x) {
    Expr* dr = deriv_dx(r, x);
    if (!dr) return false;
    Expr* diff = mk_fn2("Plus", dr,
                        mk_fn2("Times", mk_int(-1), expr_copy((Expr*)f)));

    Expr* pz = eval_take(mk_fn1("PossibleZeroQ", expr_copy(diff)));
    bool possible = pz && pz->type == EXPR_SYMBOL && pz->data.symbol == SYM_True;
    if (pz) expr_free(pz);
    if (!possible) { expr_free(diff); return false; }

    Expr* s = eval_take(mk_fn1("Simplify", expr_copy(diff)));
    bool zero = is_lit_zero(s);
    if (s) expr_free(s);
    if (!zero) {
        Expr* ex = eval_take(internal_expand((Expr*[]){ expr_copy(diff) }, 1));
        Expr* ct = cancel_together(ex);
        zero = is_lit_zero(ct);
        if (ct) expr_free(ct);
    }
    expr_free(diff);
    return zero;
}

/* Recurse into the full integrator: Integrate[g, u].  Returns the closed
 * antiderivative, or NULL if it does not close.  Borrows g and u. */
static Expr* integrate_in(const Expr* g, const Expr* u) {
    Expr* r = eval_take(mk_fn2("Integrate", expr_copy((Expr*)g),
                               expr_copy((Expr*)u)));
    if (!r) return NULL;
    if (contains_unintegrated(r)) { expr_free(r); return NULL; }
    return r;
}

/* ---------------------------------------------------------------------- */
/* Recursion guard + fresh-symbol counter                                 */
/* ---------------------------------------------------------------------- */

#define LR_MAX_DEPTH 8
static int lr_depth = 0;
static unsigned long lr_sym_counter = 0;

/* ---------------------------------------------------------------------- */
/* Identification: the shared linear base and its radical LCM             */
/* ---------------------------------------------------------------------- */

/* Walk `e` collecting every Power[base, p/q] (q > 1) whose base depends on x.
 * Each such base must be a degree-1 polynomial in x and all must be equal.
 *  *first_base : owned copy of the shared base (NULL until found).
 *  *n          : running LCM of the radical denominators.
 *  *found      : count flag (0 = none yet).
 *  *ok         : cleared on any rejecting condition. */
static void lr_scan(Expr* e, Expr* x, Expr** first_base,
                    int64_t* n, int* found, bool* ok) {
    if (!*ok || !e) return;
    if (e->type != EXPR_FUNCTION) return;

    if (head_is(e, SYM_Power) && e->data.function.arg_count == 2) {
        Expr* base = e->data.function.args[0];
        Expr* exp  = e->data.function.args[1];
        int64_t p, q;
        if (is_rational(exp, &p, &q) && (q > 1 || q < -1)) {
            if (!expr_free_of(base, x)) {
                /* A radical whose argument depends on x: it must be linear. */
                Expr* vars[1] = { x };
                if (is_polynomial(base, vars, 1) && get_degree_poly(base, x) == 1) {
                    int64_t qa = q < 0 ? -q : q;
                    if (*found == 0) {
                        *first_base = expr_copy(base);
                        *n = qa;
                        *found = 1;
                    } else if (!expr_eq(base, *first_base)) {
                        *ok = false;   /* distinct linear bases: not R[x,(ax+b)^.] */
                        return;
                    } else {
                        *n = lcm(*n, qa);
                    }
                } else {
                    *ok = false;       /* radical of a non-linear arg in x */
                    return;
                }
            }
            /* base free of x: a constant radical -- ignored. */
        }
    }

    lr_scan(e->data.function.head, x, first_base, n, found, ok);
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        lr_scan(e->data.function.args[i], x, first_base, n, found, ok);
    }
}

/* ---------------------------------------------------------------------- */
/* Core driver                                                            */
/* ---------------------------------------------------------------------- */

static Expr* lr_core(Expr* f, Expr* x) {
    if (x->type != EXPR_SYMBOL) return NULL;
    if (expr_free_of(f, x))      return NULL;   /* nothing to integrate in x */
    if (lr_depth >= LR_MAX_DEPTH) return NULL;

    /* --- 1. Identify the shared linear base and the radical LCM n. --- */
    Expr* base = NULL;
    int64_t n = 1;
    int found = 0;
    bool ok = true;
    lr_scan(f, x, &base, &n, &found, &ok);
    if (!ok || !found || n < 2) { if (base) expr_free(base); return NULL; }

    /* a x + b coefficients (a != 0 since deg == 1). */
    Expr* a = get_coeff(base, x, 1);
    Expr* b = get_coeff(base, x, 0);
    if (!a || !b) { expr_free(base); if (a) expr_free(a); if (b) expr_free(b); return NULL; }

    /* Fresh substitution variable u. */
    char uname[80];
    snprintf(uname, sizeof(uname), "Integrate`LinearRadicals`u$%lu", lr_sym_counter++);
    Expr* u = expr_new_symbol(uname);

    Expr* result = NULL;

    /* --- 2. Substitute radicals -> u and require a rational function. --- */
    Expr* one = mk_int(1);
    Expr* exy = poly_subst_radical_to_gen(f, base, one, n, uname);
    expr_free(one);
    if (exy) exy = cancel_together(exy);
    if (!exy) goto done;

    {
        Expr* num = eval_take(mk_fn1("Numerator",   expr_copy(exy)));
        Expr* den = eval_take(mk_fn1("Denominator", expr_copy(exy)));
        Expr* vars2[2] = { x, u };
        bool rat_ok = num && den &&
                      is_polynomial(num, vars2, 2) &&
                      is_polynomial(den, vars2, 2);
        if (num) expr_free(num);
        if (den) expr_free(den);
        if (!rat_ok) { expr_free(exy); goto done; }
    }

    /* --- 3. (n/a) u^(n-1) exy  /.  x -> (u^n - b)/a, then Cancel. --- */
    Expr* jac = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        mk_int(n),
        mk_fn2("Power", expr_copy(a), mk_int(-1)),
        mk_fn2("Power", expr_copy(u), mk_int(n - 1)),
        exy                                       /* consumes exy */
    }, 4);
    /* sub = (u^n - b)/a */
    Expr* sub = mk_fn2("Times",
        mk_fn2("Plus", mk_fn2("Power", expr_copy(u), mk_int(n)),
                       mk_fn2("Times", mk_int(-1), expr_copy(b))),
        mk_fn2("Power", expr_copy(a), mk_int(-1)));
    Expr* integrand = replace_one(jac, x, sub);   /* consumes jac */
    expr_free(sub);
    integrand = cancel_together(integrand);
    if (!integrand) goto done;

    /* --- 4. Recurse, back-substitute, verify. --- */
    lr_depth++;
    Expr* G = integrate_in(integrand, u);
    lr_depth--;
    expr_free(integrand);
    if (!G) goto done;

    Expr* back = mk_fn2("Power", expr_copy(base), make_rational(1, n));
    Expr* r = replace_one(G, u, back);            /* consumes G */
    expr_free(back);

    if (r && differentiates_back(r, f, x)) result = r;
    else if (r) expr_free(r);

done:
    expr_free(base);
    expr_free(a);
    expr_free(b);
    expr_free(u);
    return result;
}

/* ---------------------------------------------------------------------- */
/* Public entry points                                                    */
/* ---------------------------------------------------------------------- */

Expr* integrate_linrad_try(Expr* f, Expr* x) {
    return lr_core(f, x);
}

Expr* builtin_integrate_linrad(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 2) return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    return lr_core(f, x);
}

void integrate_linrad_init(void) {
    symtab_add_builtin("Integrate`LinearRadicals", builtin_integrate_linrad);
    symtab_get_def("Integrate`LinearRadicals")->attributes |=
        ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Integrate`LinearRadicals",
        "Integrate`LinearRadicals[f, x] integrates a rational function of x and\n"
        "radicals (a x + b)^(m/n) sharing one linear argument. It substitutes\n"
        "u = (a x + b)^(1/n) (n = LCM of the radical denominators), reducing f to\n"
        "a rational function of u, integrates that, back-substitutes, and verifies\n"
        "by differentiation. Strict: returns unevaluated when f is not of this\n"
        "form or the reduced integral does not close.");
}
